#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GenMap bedGraph -> (optional) smoothing -> valley(anchor) detection

This is a minimally invasive adaptation of your original "self-BLAST outfmt6 -> coverage"
script to work with GenMap bedGraph output (e.g. genout.bedgraph).

Key idea:
- Instead of building per-base coverage via BLAST HSP intervals (imos),
  we directly read the per-position GenMap signal from bedGraph and treat it as `cov`.

Expected input:
- GenMap command example:
    ./genmap-build/bin/genmap map -K 21 -E 2 -I index/ -O genout -t -w -bg -fl
  -> produces genout.bedgraph

Notes:
- bedGraph is usually sparse (runs of identical values). We expand to per-base array.
- We keep the rest of your pipeline intact: smoothing -> threshold -> valleys -> anchors/repeats -> plots/bigwig.
"""

import sys
import argparse
import gzip
import time
import os
import csv
import subprocess
import tempfile
import shutil
from bisect import bisect_left

# ========== 共通ログ/ユーティリティ ==========

def log(msg: str) -> None:
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()

def human(n):
    try:
        n = float(n)
    except Exception:
        return str(n)
    for u in ["", "K", "M", "G"]:
        if abs(n) < 1024:
            return f"{n:.1f}{u}"
        n /= 1024.0
    return f"{n:.1f}T"

# ========== I/O ==========

def read_fasta_len(fa_path, log_every_bp=50_000_000):
    """
    Read first record length only (as in your original).
    Works for .fa/.fasta and .gz.
    """
    t0 = time.perf_counter()
    name, L = None, 0
    op = gzip.open if fa_path.endswith(".gz") else open
    try:
        sizeB = os.path.getsize(fa_path)
        log(f"[INFO] Reading FASTA length: {fa_path} ({human(sizeB)}B)")
    except Exception:
        log(f"[INFO] Reading FASTA length: {fa_path}")

    with op(fa_path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    break
                name = line[1:].strip().split()[0]
            else:
                if name is not None:
                    L += len(line.strip())
                    if log_every_bp and (L % log_every_bp == 0):
                        log(f"[FASTA] counted {human(L)} bp so far…")
    if name is None or L == 0:
        raise RuntimeError("FASTA の長さが読めませんでした: " + fa_path)
    dt = time.perf_counter() - t0
    log(f"[INFO] FASTA: {name} length={L:,} (took {dt:.2f}s)")
    return name, L

def parse_bedgraph(path, log_every_lines=5_000_000):
    """
    bedGraph parser:
      chrom  start  end  value
    start/end are 0-based, half-open.

    Yields (chrom, start0, end0, value_float)
    Supports .gz
    """
    t0 = time.perf_counter()
    op = gzip.open if path.endswith(".gz") else open
    size = None
    try:
        size = os.path.getsize(path)
    except Exception:
        pass
    log(f"[INFO] Reading bedGraph: {path}" + (f" ({human(size)}B)" if size is not None else ""))

    n = 0
    with op(path, "rt") as f:
        for ln in f:
            n += 1
            if log_every_lines and (n % log_every_lines == 0):
                dt = time.perf_counter() - t0
                log(f"[bedGraph] parsed {human(n)} lines in {dt:.1f}s (avg {n/max(dt,1e-9):.0f} lines/s)")
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            if ln.lower().startswith(("track", "browser")):
                continue
            cols = ln.split()
            if len(cols) < 4:
                continue
            chrom = cols[0]
            try:
                s0 = int(cols[1])
                e0 = int(cols[2])
                v = float(cols[3])
            except ValueError:
                continue
            if e0 <= s0:
                continue
            yield chrom, s0, e0, v

    dt = time.perf_counter() - t0
    log(f"[INFO] bedGraph done: lines={n:,}, took {dt:.2f}s")

# ========== マスク区間の読み込み/判定 ==========

def _merge_intervals(iv, max_gap=0):
    if not iv:
        return []
    iv.sort()
    out = []
    cs, ce = iv[0]
    for s, e in iv[1:]:
        if s <= ce + max_gap:
            if e > ce:
                ce = e
        else:
            out.append((cs, ce))
            cs, ce = s, e
    out.append((cs, ce))
    return out

def _load_skip_intervals(bed_path, ref_name, label_col=4, labels=None, pad=0, seqlen=None, merge_gap=0, invert=False):
    if not bed_path:
        return []

    labset = set()
    if labels:
        labset = {s.strip().lower() for s in labels.split(",") if s.strip()}

    idx = label_col - 1
    ivs = []

    op = gzip.open if str(bed_path).endswith(".gz") else open

    with op(bed_path, "rt") as f:
        for ln in f:
            if not ln.strip():
                continue
            if ln[0] == "#":
                continue
            ls = ln.lower()
            if ls.startswith("track") or ls.startswith("browser"):
                continue
            cols = ln.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            if cols[0] != ref_name:
                continue
            try:
                s = int(cols[1])
                e = int(cols[2])  # 0-based 半開
            except (ValueError, IndexError):
                continue
            lbl = cols[idx] if idx < len(cols) else ""

            lab_ok = (not labset) or any(k in lbl.lower() for k in labset)
            take = (lab_ok and not invert) or ((not lab_ok) and invert)
            if not take:
                continue

            if e > s:
                ivs.append((s, e))

    if not ivs:
        return []

    if pad and pad > 0:
        tmp = []
        for s, e in ivs:
            ns = max(0, s - pad)
            ne = e + pad if seqlen is None else min(seqlen, e + pad)
            if ne > ns:
                tmp.append((ns, ne))
        ivs = tmp

    merged = _merge_intervals(ivs, merge_gap)
    return merged

def _overlaps_skip(start0, end0, merged_ivs):
    if not merged_ivs:
        return False
    starts = [s for s, _ in merged_ivs]
    i = bisect_left(starts, end0)
    for j in (i - 1, i):
        if 0 <= j < len(merged_ivs):
            s, e = merged_ivs[j]
            if max(s, start0) < min(e, end0):
                return True
    return False

# ========== 配列処理 ==========

def moving_average(arr, win):
    if win <= 1:
        return arr[:]
    n = len(arr)
    out = [0.0] * n
    csum = [0.0] * (n + 1)
    for i, v in enumerate(arr):
        csum[i + 1] = csum[i] + float(v)
        if i and (i % 10_000_000 == 0):
            log(f"[SMOOTH] prefix-sum progressed: {human(i)} / {human(n)}")
    half = win // 2
    for i in range(n):
        left = max(0, i - half)
        r = min(n - 1, i + half)
        s = csum[r + 1] - csum[left]
        out[i] = s / (r - left + 1)
        if i and (i % 10_000_000 == 0):
            log(f"[SMOOTH] moving-average progressed: {human(i)} / {human(n)}")
    return out

def find_valleys_maskaware(sm, radius, thr, masked):
    L = len(sm)
    valleys = []
    for i in range(L):
        if masked is not None and masked[i]:
            continue
        left = max(0, i - radius)
        r = min(L - 1, i + radius)
        v = sm[i]
        if v > thr:
            continue
        is_min = True
        for j in range(left, r + 1):
            if masked is not None and masked[j]:
                continue
            if sm[j] < v - 1e-12:
                is_min = False
                break
        if is_min:
            valleys.append((i, v, thr))
        if i and (i % 10_000_000 == 0):
            log(f"[VALLEY] scan progressed: {human(i)} / {human(L)} (current valleys={len(valleys):,})")
    return valleys

def suppress_close_valleys(valleys, merge_radius):
    kept = []
    for i, v, thr in sorted(valleys, key=lambda x: x[1]):
        if any(abs(i - ki) <= merge_radius for ki, _, _ in kept):
            continue
        kept.append((i, v, thr))
    log(f"[INFO] suppress_close_valleys: in={len(valleys):,}, out={len(kept):,}, merge_radius={merge_radius}")
    return kept

def detect_valley_span(sm, center_idx, thr_low, thr_high=None, max_expand=None):
    n = len(sm)
    th = thr_low if thr_high is None else thr_high
    L = center_idx
    steps = 0
    while L - 1 >= 0:
        if max_expand is not None and steps >= max_expand:
            break
        if sm[L - 1] <= th:
            L -= 1
            steps += 1
        else:
            break
    R = center_idx + 1
    steps = 0
    while R < n:
        if max_expand is not None and steps >= max_expand:
            break
        if sm[R] <= th:
            R += 1
            steps += 1
        else:
            break
    return L, R  # 0-based 半開

# ========== プロット ==========

def plot_valleys_png(sm, valleys, thr, out_png,
                     relax_factor=1.2,
                     plot_range=None,
                     figsize=(12, 4)):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    n = len(sm)
    if plot_range is None:
        a1, b1 = 1, n
    else:
        a1, b1 = plot_range
        a1 = max(1, min(n, a1))
        b1 = max(1, min(n, b1))
        if a1 > b1:
            a1, b1 = b1, a1

    a0 = a1 - 1
    b0 = b1
    xs = list(range(a1, b1 + 1))
    ys = sm[a0:b0]

    log(f"[PLOT] range={a1}-{b1} (len={len(xs)}), valleys={len(valleys):,}, thr={thr:.3f}")
    t0 = time.perf_counter()

    plt.figure(figsize=figsize)
    ax = plt.gca()
    ax.plot(xs, ys, label="Coverage (smoothed)")
    ax.axhline(y=thr, linestyle="--", label=f"threshold={thr:.2f}")

    if valleys:
        thr_low = thr
        thr_high = thr * relax_factor if relax_factor else thr
        for (i, v, _) in valleys:
            if not (a0 <= i < b0):
                continue
            start0, end0 = detect_valley_span(sm, i, thr_low, thr_high, max_expand=None)
            start1 = max(a1, start0 + 1)
            end1 = min(b1, end0)
            if start1 <= end1:
                ax.axvspan(start1, end1, alpha=0.25, linewidth=0)

        red_patch = mpatches.Patch(alpha=0.25, label=f"Valley spans (auto, x{relax_factor:.2f})")
        handles, labels = ax.get_legend_handles_labels()
        handles.append(red_patch)
        labels.append(f"Valley spans (auto, x{relax_factor:.2f})")
        ax.legend(handles, labels)
    else:
        ax.legend()

    ax.set_title("Coverage profile with auto-detected valley spans")
    ax.set_xlabel("Position (1-based)")
    ax.set_ylabel("Coverage")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    log(f"[PLOT] saved {out_png} (took {time.perf_counter()-t0:.2f}s)")

# ========== GenMap bedGraph -> cov array ==========

def bedgraph_to_cov_array(bedgraph_path, ref_name, L, default_value=0.0, log_every_lines=5_000_000):
    """
    Convert bedGraph to per-base array `cov` length L.
    - Only uses records matching chrom == ref_name.
    - Clips to [0, L).
    """
    log(f"[INFO] Building cov array from bedGraph (default={default_value})")
    t0 = time.perf_counter()
    cov = [float(default_value)] * L

    n_used = 0
    n_other_chrom = 0
    n_clipped = 0

    for chrom, s0, e0, v in parse_bedgraph(bedgraph_path, log_every_lines=log_every_lines):
        if chrom != ref_name:
            n_other_chrom += 1
            continue
        if s0 < 0:
            s0 = 0
            n_clipped += 1
        if e0 > L:
            e0 = L
            n_clipped += 1
        if e0 <= s0:
            continue
        # Expand run
        # NOTE: This is O(total_bases). bedGraph usually has long runs, so OK.
        for i in range(s0, e0):
            cov[i] = float(v)

        n_used += 1
        if log_every_lines and (n_used % log_every_lines == 0):
            dt = time.perf_counter() - t0
            log(f"[cov] filled {n_used:,} bedGraph intervals in {dt:.1f}s")

    dt = time.perf_counter() - t0
    log(f"[INFO] cov from bedGraph done: intervals_used={n_used:,}, other_chrom={n_other_chrom:,}, clipped={n_clipped:,}, took {dt:.2f}s")
    return cov

# ========== main ==========

def main():
    ap = argparse.ArgumentParser(description="GenMap bedGraph を被覆として読み込み、平滑化して谷をアンカー検出（skip-bed対応・詳細進捗ログ付き）")

    ap.add_argument("--fasta", required=True, help="FASTA (single record recommended; .gz ok)")
    ap.add_argument("--outprefix", required=True)

    # NEW: bedGraph input (GenMap)
    ap.add_argument("--bedgraph", required=True, help="GenMap bedGraph (chrom start end value). .gz ok")
    ap.add_argument("--bedgraph-default", type=float, default=0.0, help="Default value for positions not covered by bedGraph (default: 0.0)")
    ap.add_argument("--bedgraph-log-every", type=int, default=5_000_000, help="Progress interval for bedGraph parsing (lines)")

    # smoothing + valley
    ap.add_argument("--smooth", type=int, default=1000)
    ap.add_argument("--valley-radius", type=int, default=None)
    ap.add_argument("--valley-pct", type=float, default=30.0)

    ap.add_argument("--anchor-half", type=int, default=0, help="固定幅で出す場合の半幅(bp)")
    ap.add_argument("--auto-span", action="store_true", help="谷ごとに自動で幅を検出してBEDに出力")
    ap.add_argument("--auto-span-factor", type=float, default=1.2, help="thr_high = thr * factor")
    ap.add_argument("--auto-span-max", type=int, default=None, help="左右最大拡張幅（bp）")
    ap.add_argument("--score-scale", type=float, default=100.0, help="BED score のスケール係数")

    # plot
    ap.add_argument("--plot-png", type=str, default=None, help="PNG を保存するパス（例: out.png）。未指定なら出力しない")
    ap.add_argument("--plot-width", type=float, default=12.0)
    ap.add_argument("--plot-height", type=float, default=4.0)
    ap.add_argument("--plot-range", type=str, default=None, help="1-based 範囲（例: 100000-200000）。'auto' なら非マスク領域のうち最長区間を自動選択")
    ap.add_argument("--relax-factor", type=float, default=1.2)

    # skip
    ap.add_argument("--skip-bed", type=str, default=None, help="この BED の該当ラベル区間を coverage/閾値/探索/出力すべてから除外")
    ap.add_argument("--skip-labels", type=str, default=None, help="スキップ対象ラベル（カンマ区切り, 部分一致, 大文字小文字無視）例: 'ct'")
    ap.add_argument("--skip-col", type=int, default=4, help="ラベル列（1始まり; 既定=4）")
    ap.add_argument("--skip-merge-distance", type=int, default=0, help="除外区間マージ距離（bp, 既定:0）")
    ap.add_argument("--skip-pad", type=int, default=0, help="除外区間の前後パディング(bp)")
    ap.add_argument("--skip-invert", action="store_true", help="--skip-labels に一致しないラベルの区間を skip 対象にする（反転）")
    ap.add_argument("--focus-labels", action="store_true", help="指定ラベルの区間だけを解析対象にし、それ以外の全域をskip（=補集合をskip）")

    # filters + exports
    ap.add_argument("--min-anchor-len", type=int, default=20, help="Minimum anchor length (bp) to keep (default: 20)")
    ap.add_argument("--export-highcov-bed", action="store_true", help="Export BED of high-coverage spans (coverage_smooth > threshold) excluding skip regions.")
    ap.add_argument("--highcov-threshold", type=float, default=None, help="Override threshold for high-coverage BED. If unset, uses the auto 'thr' from valley_pct.")
    ap.add_argument("--highcov-min-len", type=int, default=5, help="Minimum span length (bp) to keep in high-coverage BED (default: 5)")
    ap.add_argument("--no-coverage-tsv", action="store_true", help="Do not write *_coverage.tsv.gz (default: write)")

    # BigWig
    ap.add_argument("--write-bigwig", action="store_true", help="Write BigWig files for IGV/JBrowse (smooth/raw)")
    ap.add_argument("--bigwig-prefix", type=str, default=None, help="Prefix for BigWig outputs (default: outprefix)")
    ap.add_argument("--bigwig-kind", type=str, choices=["smooth", "raw", "both"], default="smooth", help="Which track(s) to export: smoothed, raw, or both (default: smooth)")

    args = ap.parse_args()
    t_all = time.perf_counter()

    ref_name, L = read_fasta_len(args.fasta)

    # BigWig writer (same as yours)
    def _write_bigwig_from_array(values, kind_label, masked=None):
        out_prefix = args.bigwig_prefix if args.bigwig_prefix else args.outprefix
        out_bw = f"{out_prefix}_{kind_label}.bw"
        log(f"[INFO] BigWig export ({kind_label}): target={out_bw}")

        t0 = time.perf_counter()
        with tempfile.NamedTemporaryFile("w", delete=False, prefix=f"{kind_label}_", suffix=".bedGraph") as tmp:
            bedgraph_path = tmp.name
            for i in range(L):
                if masked is not None and masked[i]:
                    continue
                v = values[i]
                if isinstance(v, float):
                    v = f"{v:.6f}"
                tmp.write(f"{ref_name}\t{i}\t{i+1}\t{v}\n")
        log(f"[INFO] bedGraph temp written: {bedgraph_path} (took {time.perf_counter()-t0:.2f}s)")

        with tempfile.NamedTemporaryFile("w", delete=False, prefix="chrom_", suffix=".sizes") as tmp2:
            chrom_sizes = tmp2.name
            tmp2.write(f"{ref_name}\t{L}\n")

        tool = shutil.which("bedGraphToBigWig")
        if tool:
            try:
                t0 = time.perf_counter()
                subprocess.check_call([tool, bedgraph_path, chrom_sizes, out_bw])
                log(f"[INFO] BigWig written via bedGraphToBigWig: {out_bw} (took {time.perf_counter()-t0:.2f}s)")
                os.unlink(bedgraph_path)
                os.unlink(chrom_sizes)
                return
            except subprocess.CalledProcessError as e:
                log(f"[WARN] bedGraphToBigWig failed with code {e.returncode}; will try pyBigWig fallback")

        try:
            import pyBigWig
            t0 = time.perf_counter()
            bw = pyBigWig.open(out_bw, "w")
            bw.addHeader([(ref_name, L)])

            batch_starts, batch_ends, batch_vals = [], [], []
            BATCH = 1_000_000

            for i in range(L):
                if masked is not None and masked[i]:
                    continue
                batch_starts.append(i)
                batch_ends.append(i + 1)
                batch_vals.append(float(values[i]))
                if len(batch_starts) >= BATCH:
                    bw.addEntries([ref_name] * len(batch_starts), batch_starts, ends=batch_ends, values=batch_vals)
                    batch_starts.clear()
                    batch_ends.clear()
                    batch_vals.clear()

            if batch_starts:
                bw.addEntries([ref_name] * len(batch_starts), batch_starts, ends=batch_ends, values=batch_vals)

            bw.close()
            log(f"[INFO] BigWig written via pyBigWig: {out_bw} (took {time.perf_counter()-t0:.2f}s)")
            os.unlink(bedgraph_path)
            os.unlink(chrom_sizes)
            return
        except Exception as e:
            log(f"[ERROR] BigWig export failed (no bedGraphToBigWig and pyBigWig unavailable?): {e}")
            log(f"[HINT] You can manually convert: bedGraphToBigWig {bedgraph_path} {chrom_sizes} {out_bw}")
            return

    # --- skip intervals ---
    skip_ivs = []
    if args.skip_bed:
        log(f"[INFO] Loading skip intervals from {args.skip_bed} (labels={args.skip_labels}, col={args.skip_col}, pad={args.skip_pad}, merge_gap={args.skip_merge_distance})")

        if args.focus_labels:
            focus_ivs = _load_skip_intervals(
                args.skip_bed, ref_name,
                label_col=args.skip_col,
                labels=args.skip_labels,
                pad=args.skip_pad,
                seqlen=L,
                merge_gap=args.skip_merge_distance,
                invert=False,
            )
            skip_ivs = []
            prev = 0
            for s, e in focus_ivs:
                if s > prev:
                    skip_ivs.append((prev, s))
                prev = e
            if prev < L:
                skip_ivs.append((prev, L))
        else:
            skip_ivs = _load_skip_intervals(
                args.skip_bed, ref_name,
                label_col=args.skip_col,
                labels=args.skip_labels,
                pad=args.skip_pad,
                seqlen=L,
                merge_gap=args.skip_merge_distance,
                invert=args.skip_invert,
            )

        total_skip_bp = sum(e - s for s, e in skip_ivs)
        pct = 100.0 * total_skip_bp / max(1, L)
        log(f"[INFO] skip intervals merged: {len(skip_ivs)}, total_skip_bp={total_skip_bp:,} ({pct:.2f}% of {L:,})")

    # --- NEW: Build cov directly from GenMap bedGraph ---
    cov = bedgraph_to_cov_array(
        bedgraph_path=args.bedgraph,
        ref_name=ref_name,
        L=L,
        default_value=float(args.bedgraph_default),
        log_every_lines=int(args.bedgraph_log_every),
    )

    # smoothing
    win = max(1, args.smooth)
    log(f"[INFO] smoothing window={win}")
    t0 = time.perf_counter()
    sm = moving_average(cov, win)
    log(f"[INFO] smoothing done in {time.perf_counter()-t0:.2f}s")

    # masked array + make masked sm BIG (same logic)
    masked = None
    if skip_ivs:
        masked = [False] * L
        for s, e in skip_ivs:
            s = max(0, s)
            e = min(L, e)
            for i in range(s, e):
                masked[i] = True
        BIG = (max(sm) + 1.0) if sm else 1.0
        for s, e in skip_ivs:
            s = max(0, s)
            e = min(L, e)
            for i in range(s, e):
                sm[i] = BIG

    # radius
    if args.valley_radius is None:
        radius = int(max(1, 2 * win / 5))
    else:
        radius = int(args.valley_radius)
    min_radius = int(max(1, 2 * win / 10))
    if radius < min_radius:
        log(f"[WARNING] valley_radius={radius} < recommended minimum ({min_radius}); using {min_radius}")
        radius = min_radius

    # threshold on unmasked
    if masked:
        xs = [sm[i] for i in range(L) if not masked[i]]
    else:
        xs = sm[:]
    xs_sorted = sorted(xs)
    if xs_sorted:
        k = max(0, min(len(xs_sorted) - 1, int(len(xs_sorted) * args.valley_pct / 100.0)))
        thr = xs_sorted[k]
    else:
        thr = 0.0
    log(f"[INFO] threshold (percentile {args.valley_pct}% on unmasked) = {thr:.6f}")

    # valleys
    t0 = time.perf_counter()
    valleys = find_valleys_maskaware(sm, radius, thr, masked)
    log(f"[INFO] find_valleys: valleys={len(valleys):,}, scan_time={time.perf_counter()-t0:.2f}s")

    valleys = suppress_close_valleys(valleys, merge_radius=radius)
    log(f"[INFO] valleys after suppression: {len(valleys):,}")

    # span-based non-overlap filtering (same as yours)
    if args.auto_span:
        thr_low = thr
        thr_high = thr * args.auto_span_factor if args.auto_span_factor else None
        max_expand = args.auto_span_max
    else:
        half = max(0, args.anchor_half)

    spans = []
    for (i, v, _) in valleys:
        if args.auto_span:
            s0, e0 = detect_valley_span(sm, i, thr_low, thr_high, max_expand)
        else:
            s0, e0 = max(0, i - half), min(L, i + half + 1)
        spans.append((s0, e0, i, v))

    spans.sort(key=lambda x: x[3])
    kept = []
    for s0, e0, i, v in spans:
        overlap = any(not (e0 <= ks or s0 >= ke) for ks, ke, *_ in kept)
        if overlap:
            continue
        kept.append((s0, e0, i, v))

    log(f"[INFO] non-overlapping valleys kept: {len(kept):,} / {len(valleys):,} (span-based filtering)")

    min_anchor_len = int(args.min_anchor_len)
    filtered = [(s0, e0, i, v) for (s0, e0, i, v) in kept if (e0 - s0) >= min_anchor_len]
    removed = len(kept) - len(filtered)
    if removed > 0:
        log(f"[INFO] removed {removed} anchors shorter than {min_anchor_len} bp")
    kept = filtered

    valleys = [(i, v, thr) for s0, e0, i, v in kept]
    total_valleys = max(1, len(valleys))

    # coverage TSV
    cov_tsv = args.outprefix + "_coverage.tsv.gz"
    if not args.no_coverage_tsv:
        log(f"[INFO] writing coverage: {cov_tsv}")
        t0 = time.perf_counter()
        with gzip.open(cov_tsv, "wt") as w:
            cw = csv.writer(w, delimiter="\t")
            cw.writerow(["pos1", "coverage_raw", "coverage_smooth"])
            for i in range(L):
                if i and (i % 10_000_000 == 0):
                    pct = 100.0 * i / max(1, L)
                    log(f"[WRITE] coverage rows written: {human(i)} / {human(L)} ({pct:.1f}%)")
                cw.writerow([i + 1, f"{float(cov[i]):.6f}", f"{sm[i]:.6f}"])
        log(f"[INFO] coverage written in {time.perf_counter()-t0:.2f}s")
    else:
        log("[INFO] skipping *_coverage.tsv.gz as requested (--no-coverage-tsv)")

    # anchors BED
    anc_bed = args.outprefix + "_anchors.bed"
    skipped = 0
    kept_cnt = 0
    t0 = time.perf_counter()
    with open(anc_bed, "w") as w:
        if args.auto_span:
            thr_low = thr
            thr_high = thr * args.auto_span_factor if args.auto_span_factor else None
            max_expand = args.auto_span_max
            for idx, (i, v, _) in enumerate(valleys, 1):
                start0, end0 = detect_valley_span(sm, i, thr_low, thr_high, max_expand)
                start0 = max(0, start0)
                end0 = min(L, end0)
                if _overlaps_skip(start0, end0, skip_ivs):
                    skipped += 1
                else:
                    score = int(round(float(v) * args.score_scale))
                    w.write(f"{ref_name}\t{start0}\t{end0}\tanchor_{i+1}\t{score}\t+\n")
                    kept_cnt += 1
                if idx % 1_000_000 == 0:
                    pct = 100.0 * idx / total_valleys
                    log(f"[WRITE] anchors bed: processed {idx:,}/{total_valleys:,} ({pct:.1f}%), kept={kept_cnt:,}, skipped={skipped:,}")
        else:
            half = max(0, args.anchor_half)
            for idx, (i, v, _) in enumerate(valleys, 1):
                start0 = max(0, i - half)
                end0 = min(L, i + half + 1)
                if _overlaps_skip(start0, end0, skip_ivs):
                    skipped += 1
                else:
                    score = int(round(float(v) * args.score_scale))
                    w.write(f"{ref_name}\t{start0}\t{end0}\tanchor_{i+1}\t{score}\t+\n")
                    kept_cnt += 1
                if idx % 1_000_000 == 0:
                    pct = 100.0 * idx / total_valleys
                    log(f"[WRITE] anchors bed: processed {idx:,}/{total_valleys:,} ({pct:.1f}%), kept={kept_cnt:,}, skipped={skipped:,}")
    log(f"[INFO] anchors kept={kept_cnt:,}, skipped_by_labels={skipped:,}")
    log(f"[INFO] anchors BED written in {time.perf_counter()-t0:.2f}s")

    # high-coverage BED (repeat candidates)
    if args.export_highcov_bed:
        highcov_thr = args.highcov_threshold if (args.highcov_threshold is not None) else thr
        min_len = max(0, int(args.highcov_min_len))

        out_hc_bed = args.outprefix + "_repeats.bed"
        log(f"[INFO] writing high-coverage BED: {out_hc_bed} (thr={highcov_thr:.6f}, min_len={min_len})")

        t0 = time.perf_counter()
        with open(out_hc_bed, "w") as w:
            in_seg = False
            seg_s = 0
            seg_sum = 0.0
            seg_len = 0
            seg_idx = 0

            def _flush_segment(end_i):
                nonlocal in_seg, seg_s, seg_sum, seg_len, seg_idx
                seg_e = end_i
                if seg_e - seg_s >= min_len:
                    seg_idx += 1
                    avg_sm = (seg_sum / max(1, seg_len))
                    score = int(round(avg_sm * args.score_scale))
                    name = f"repeat_{seg_idx}"
                    w.write(f"{ref_name}\t{seg_s}\t{seg_e}\t{name}\t{score}\t+\n")
                in_seg = False
                seg_sum = 0.0
                seg_len = 0

            for i in range(L):
                if masked is not None and masked[i]:
                    if in_seg:
                        _flush_segment(i)
                    continue

                if sm[i] > highcov_thr:
                    if not in_seg:
                        in_seg = True
                        seg_s = i
                        seg_sum = 0.0
                        seg_len = 0
                    seg_sum += float(sm[i])
                    seg_len += 1
                else:
                    if in_seg:
                        _flush_segment(i)

            if in_seg:
                _flush_segment(L)

        log(f"[INFO] high-coverage BED written in {time.perf_counter()-t0:.2f}s")

    # BigWig export
    if args.write_bigwig:
        kinds = []
        if args.bigwig_kind in ("smooth", "both"):
            kinds.append(("smooth", sm))
        if args.bigwig_kind in ("raw", "both"):
            kinds.append(("raw", cov))
        for label, arr in kinds:
            _write_bigwig_from_array(arr, label, masked=masked)

    done_files = [anc_bed]
    if not args.no_coverage_tsv:
        done_files.insert(0, cov_tsv)
    log(f"[DONE] wrote: {', '.join(done_files)}")

    # PNG
    if args.plot_png:
        pr = None
        if args.plot_range:
            if isinstance(args.plot_range, str) and args.plot_range.lower() == "auto":
                if skip_ivs:
                    unmasked = []
                    prev = 0
                    for s, e in skip_ivs:
                        if s > prev:
                            unmasked.append((prev, s))
                        prev = e
                    if prev < L:
                        unmasked.append((prev, L))

                    if unmasked:
                        a0, b0 = max(unmasked, key=lambda x: x[1] - x[0])
                        view_len = 10_000
                        mid = (a0 + b0) // 2
                        a0 = max(0, mid - view_len // 2)
                        b0 = min(L, mid + view_len // 2)
                        pr = (a0 + 1, b0)
                        log(f"[INFO] auto plot-range (fixed ±{view_len//2}bp around center): {pr[0]}-{pr[1]}")
                    else:
                        log("[WARN] auto plot-range requested but all positions are masked; skip plotting range selection")
                        pr = None
                else:
                    a0, b0 = 0, L
                    view_len = 10_000
                    mid = (a0 + b0) // 2
                    a0 = max(0, mid - view_len // 2)
                    b0 = min(L, mid + view_len // 2)
                    pr = (a0 + 1, b0)
                    log(f"[INFO] auto plot-range (fixed ±{view_len//2}bp around center): {pr[0]}-{pr[1]}")
            else:
                try:
                    a1s, b1s = str(args.plot_range).split("-")
                    pr = (int(a1s), int(b1s))
                except Exception:
                    pr = None

        plot_valleys_png(
            sm=sm,
            valleys=valleys,
            thr=thr,
            out_png=args.plot_png,
            relax_factor=args.relax_factor,
            plot_range=pr,
            figsize=(args.plot_width, args.plot_height),
        )

    log(f"[TOTAL] elapsed {time.perf_counter()-t_all:.2f}s")


if __name__ == "__main__":
    main()