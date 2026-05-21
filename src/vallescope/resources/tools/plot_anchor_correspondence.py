#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# 09_plot_anchor_correspondence.py
# families.tsv（score_filter出力）から、サンプル間に対応付けられたアンカーの位置関係を可視化
# 拡張: --famcols / --worknames を追加し、同じ列名を複数レーンに割り当て可能に

import matplotlib
matplotlib.use("Agg")  # headless 環境でも保存可能

import argparse
import os
import glob
import sys
import traceback
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List

def log(msg: str):
    print(msg, flush=True)

def detect_samples_from_families(df: pd.DataFrame) -> List[str]:
    """families.tsv の anchor_id_* / center_* のサフィックスから候補を検出"""
    samples = []
    for c in df.columns:
        if c.startswith("anchor_id_"):
            s = c[len("anchor_id_"):]
            if s and s not in samples:
                samples.append(s)
        if c.startswith("center_"):
            s = c[len("center_"):]
            if s and s not in samples:
                samples.append(s)
    return samples

def load_bed_centers(work_root: str, workname: str, verbose=False) -> Dict[str, float]:
    """
    work/<workname>/**/anchor_scoring/*.anchors.bed
    に加えて、
    work/*/<workname>/**/anchor_scoring/*.anchors.bed も探索（実運用の階層構造に対応）
    """
    patterns = [
        os.path.join(work_root, workname, "**", "anchor_scoring", "*.anchors.bed"),
        os.path.join(work_root, "*", workname, "**", "anchor_scoring", "*.anchors.bed"),
    ]
    paths = []
    for pat in patterns:
        found = sorted(glob.glob(pat, recursive=True))
        if verbose:
            print(f"[scan] pattern={pat} -> {len(found)} files", flush=True)
        paths.extend(found)

    # 重複除去
    paths = sorted(set(paths))
    if verbose:
        print(f"[scan] {workname}: total {len(paths)} BEDs resolved", flush=True)

    centers: Dict[str, float] = {}
    for bed in paths:
        try:
            df = pd.read_csv(
                bed, sep="\t", header=None,
                names=["chrom","start","end","anchor_id","score","strand"],
                usecols=[0,1,2,3,4,5]
            )
        except Exception as e:
            print(f"[warn] failed to read BED: {bed} ({e})", flush=True)
            continue
        df["center"] = (df["start"].astype(float) + df["end"].astype(float)) / 2.0
        for _, row in df.iterrows():
            centers[str(row["anchor_id"])] = float(row["center"])
    return centers


def normalize_map(center_map: Dict[str, Dict[str, float]]) -> Dict[str, Dict[str, float]]:
    """サンプルごとに center を 0..1 に線形正規化"""
    out: Dict[str, Dict[str, float]] = {}
    for s, amap in center_map.items():
        if not amap:
            out[s] = {}
            continue
        vals = np.array(list(amap.values()), dtype=float)
        vmin, vmax = np.nanmin(vals), np.nanmax(vals)
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
            out[s] = dict(amap)
            continue
        scale = vmax - vmin
        out[s] = {aid: (x - vmin)/scale for aid, x in amap.items()}
    return out

def main():
    ap = argparse.ArgumentParser(description="Plot correspondence from families.tsv")
    ap.add_argument("--families", required=True, help="results_by_chr/<group>.families.tsv")
    ap.add_argument("--work-root", default="work", help="root dir that contains per-sample dirs (default: work)")
    ap.add_argument("--out", required=True, help="output PNG")

    # 重要: 3者分離
    #   --samples  : 図のレーン名（見出し）
    #   --famcols  : families.tsv 内の anchor_id 列名（例: anchor_id_HG00171）
    #   --worknames: work/<workname>/... でBEDを探すサブディレクトリ名
    ap.add_argument("--samples", default="", help="comma-separated lane labels (e.g., HG00171_hap1,HG00171_hap2,HG00096)")
    ap.add_argument("--famcols", default="", help="comma-separated families columns for anchor ids (e.g., anchor_id_HG00171,anchor_id_HG00171,anchor_id_HG00096)")
    ap.add_argument("--worknames", default="", help="comma-separated work subdirs to fetch BED centers (e.g., HG00171_haplotype1-...,HG00171_haplotype2-...,HG00096_haplotype1-...)")

    ap.add_argument("--normalize", dest="normalize", action="store_true", default=True)
    ap.add_argument("--no-normalize", dest="normalize", action="store_false")
    ap.add_argument("--min-members", type=int, default=2, help="minimum #lanes present in a family to draw")
    ap.add_argument("--width", type=float, default=14.0)
    ap.add_argument("--height", type=float, default=8.0)
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--alpha", type=float, default=0.5)
    ap.add_argument("--title", default="Anchor correspondence (families)")
    ap.add_argument("--max-rows", type=int, default=0)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    try:
        if args.verbose:
            log(f"[args] families={args.families}")
            log(f"[args] work-root={args.work_root}")
            log(f"[args] out={args.out}")

        if not os.path.exists(args.families):
            log(f"[ERROR] families not found: {args.families}")
            sys.exit(2)

        fam = pd.read_csv(args.families, sep="\t")
        if args.max_rows and len(fam) > args.max_rows:
            fam = fam.iloc[:args.max_rows].reset_index(drop=True)

        # ===== サンプル列の決定 =====
        auto_candidates = detect_samples_from_families(fam)  # ['HG00171','HG00096', ...] のようなサフィックス
        if args.samples.strip():
            lanes = [s.strip() for s in args.samples.split(",") if s.strip()]
        else:
            # 自動検出（列サフィックス）を使うが、重複があると 2レーンにしかならないので注意
            lanes = auto_candidates

        # famcols/worknames が指定されていればそれを使う（len は lanes と一致させる）
        famcols = [s.strip() for s in args.famcols.split(",")] if args.famcols.strip() else []
        worknames = [s.strip() for s in args.worknames.split(",")] if args.worknames.strip() else []

        if famcols and len(famcols) != len(lanes):
            log(f"[ERROR] --famcols length ({len(famcols)}) must match --samples length ({len(lanes)})")
            sys.exit(3)
        if worknames and len(worknames) != len(lanes):
            log(f"[ERROR] --worknames length ({len(worknames)}) must match --samples length ({len(lanes)})")
            sys.exit(4)

        # famcols 未指定 → 従来通り anchor_id_<lane> を参照
        if not famcols:
            famcols = [f"anchor_id_{s}" for s in lanes]

        # worknames 未指定 → 従来通り work/<lane>/... を探索
        if not worknames:
            worknames = lanes[:]

        if args.verbose:
            log(f"[lanes] {lanes}")
            log(f"[famcols] {famcols}")
            log(f"[worknames] {worknames}")
            log(f"[families] columns: {list(fam.columns)}")

        # famcols の存在チェック
        missing_cols = [c for c in famcols if c not in fam.columns]
        if missing_cols:
            log(f"[ERROR] columns not found in families.tsv: {missing_cols}")
            sys.exit(5)

        if len(lanes) < 2:
            log("[ERROR] need >=2 lanes to draw.")
            sys.exit(6)

        # ===== BED から座標収集 =====
        centers_raw: Dict[str, Dict[str, float]] = {}  # lane_label -> {anchor_id: center}
        for lane, wname in zip(lanes, worknames):
            cmap = load_bed_centers(args.work_root, wname, verbose=args.verbose)
            centers_raw[lane] = cmap
            if args.verbose:
                log(f"[centers] lane={lane} (work={wname}): {len(cmap)} IDs")

        centers_for_plot = normalize_map(centers_raw) if args.normalize else centers_raw
        x_label = "Normalized position (0..1)" if args.normalize else "Genomic position (bp)"
        subtitle = " (normalized)" if args.normalize else " (raw)"

        # ===== families → 折れ線データ =====
        y_positions = {lane: i for i, lane in enumerate(reversed(lanes))}
        lines = []
        present_counts = {lane: 0 for lane in lanes}

        for _, row in fam.iterrows():
            xs, ys = [], []
            present = 0
            for lane, col in zip(lanes, famcols):
                aid = row.get(col)
                if isinstance(aid, str) and aid in centers_for_plot.get(lane, {}):
                    xs.append(centers_for_plot[lane][aid])
                    ys.append(float(y_positions[lane]))
                    present += 1
                    present_counts[lane] += 1
            if present >= args.min_members and len(xs) >= 2:
                lines.append({"xs": xs, "ys": ys})

        if not lines:
            log("[ERROR] no drawable families.")
            for lane in lanes:
                log(f"  lane={lane}: centers={len(centers_for_plot.get(lane, {}))}")
            log(f"  famcols used: {famcols}")
            log(f"  (min-members={args.min_members})")
            sys.exit(7)

        # ===== 描画 =====
        fig, ax = plt.subplots(figsize=(args.width, args.height), dpi=args.dpi)
        for ln in lines:
            ax.plot(ln["xs"], ln["ys"], lw=0.8, alpha=args.alpha, color="gray", zorder=1)

        # ラインガイド
        for lane in lanes:
            ax.axhline(y=y_positions[lane], color="lightgray", lw=0.5, zorder=0)

        ax.set_yticks([y_positions[s] for s in lanes])
        ax.set_yticklabels(list(reversed(lanes)))
        ax.set_xlabel(x_label)
        ax.set_title(args.title + subtitle)
        ax.grid(axis="x", alpha=0.2)
        plt.tight_layout()

        outdir = os.path.dirname(args.out)
        if outdir and not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        plt.savefig(args.out, bbox_inches="tight")
        plt.close()

        if os.path.exists(args.out) and os.path.getsize(args.out) > 0:
            log(f"[OK] saved: {args.out}")
        else:
            log(f"[ERROR] failed to write: {args.out}")
            sys.exit(8)

        log(f"[info] lines drawn: {len(lines)}")
        for lane in lanes:
            log(f"[info] lane {lane}: points present = {present_counts.get(lane,0)}")

    except Exception:
        log("[FATAL] Uncaught exception:")
        log("".join(traceback.format_exc()))
        sys.exit(99)

if __name__ == "__main__":
    main()
