#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import sys
from collections import defaultdict

# ---------------- I/O helpers ----------------

def open_auto(path, mode="rt"):
    """Auto-handle .gz files."""
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

# ---------------- BED utilities ----------------

def parse_bed(path, label_col, labels):
    """
    BED: 0-based half-open.
    label_col: 1-based index of the column to match (e.g., 4).
    labels: iterable of lowercase substrings to match (case-insensitive, substring).
    Return: dict[seq] -> list[(start, end)]
    """
    label_idx = label_col - 1
    labset = {s.strip().lower() for s in labels if s.strip()}
    ivs = defaultdict(list)

    with open_auto(path, "rt") as f:
        for ln in f:
            if not ln.strip():
                continue
            low = ln.lower()
            if low.startswith("track") or low.startswith("browser") or ln.startswith("#"):
                continue

            cols = ln.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            try:
                seq = cols[0]
                s = int(cols[1])
                e = int(cols[2])
            except Exception:
                continue

            label_val = cols[label_idx] if label_idx < len(cols) else ""
            lab_ok = any(k in label_val.lower() for k in labset) if labset else True

            if lab_ok and e > s:
                ivs[seq].append((s, e))
    return ivs

def merge_intervals(ivals, max_gap=0):
    """
    Merge intervals if gap <= max_gap.
    ivals: list[(s,e)], s < e, 0-based half-open
    """
    if not ivals:
        return []
    ivals = sorted(ivals, key=lambda x: (x[0], x[1]))
    out = []
    cs, ce = ivals[0]
    for s, e in ivals[1:]:
        if s - ce <= max_gap:  # overlaps/adjacent/close enough
            if e > ce:
                ce = e
        else:
            out.append((cs, ce))
            cs, ce = s, e
    out.append((cs, ce))
    return out

def apply_padding(ivals, pad, seqlen):
    """Expand each interval by ±pad, clamped to [0, seqlen]."""
    if pad <= 0 or not ivals:
        return ivals
    padded = []
    for s, e in ivals:
        ns = 0 if s - pad < 0 else s - pad
        ne = seqlen if e + pad > seqlen else e + pad
        if ne > ns:
            padded.append((ns, ne))
    return padded

def complement_intervals(ivals, seqlen):
    """
    Given merged, sorted, non-overlapping intervals within [0,seqlen),
    return the complement intervals (also 0-based half-open).
    """
    if seqlen <= 0:
        return []
    if not ivals:
        return [(0, seqlen)]
    comp = []
    cur = 0
    for s, e in ivals:
        if cur < s:
            comp.append((cur, s))
        cur = max(cur, e)
    if cur < seqlen:
        comp.append((cur, seqlen))
    return comp

# ---------------- FASTA utilities ----------------

def fasta_iter(path):
    """
    Yield (name, seqstr) for each record.
    'name' = token right after '>' up to first whitespace.
    Sequences are uppercased for consistent softmasking behavior.
    """
    with open_auto(path, "rt") as f:
        name = None
        buf = []
        for ln in f:
            if ln.startswith(">"):
                if name is not None:
                    yield name, "".join(buf).upper()
                name = ln[1:].strip().split()[0]
                buf = []
            else:
                buf.append(ln.strip())
        if name is not None:
            yield name, "".join(buf).upper()

def write_fasta(path, records, width=60):
    with open_auto(path, "wt") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i+width] + "\n")

# ---------------- Masking ----------------

def mask_sequence(seq, intervals, softmask=False, mask_char="N"):
    """
    Apply masking on given intervals (0-based half-open).
    - softmask: lowercase A/T/G/C (leave N/others as-is)
    - hardmask: replace with mask_char[0] (default 'N')
    """
    if not intervals:
        return seq
    arr = list(seq)
    L = len(arr)

    if softmask:
        for s, e in intervals:
            if s < 0: 
                s = 0
            if e > L: 
                e = L
            for i in range(s, e):
                c = arr[i]
                if "A" <= c <= "Z":
                    arr[i] = c.lower()
    else:
        mchar = mask_char[0] if mask_char else "N"
        for s, e in intervals:
            if s < 0: 
                s = 0
            if e > L: 
                e = L
            for i in range(s, e):
                arr[i] = mchar
    return "".join(arr)

# ---------------- Main ----------------

def main():
    ap = argparse.ArgumentParser(
        description="Mask FASTA by BED labels (ct, HSat2, etc.). Supports .fa and .fa.gz."
    )
    ap.add_argument("--fasta", required=True, help="Input FASTA (.fa or .fa.gz)")
    ap.add_argument("--bed", required=True, help="BED with annotations to select from")
    ap.add_argument("--out", required=True, help="Output FASTA (.fa or .fa.gz)")
    ap.add_argument("--labels", default="", help="Comma-separated labels to match (case-insensitive substring). e.g. 'ct,HSat2'")
    ap.add_argument("--column", type=int, default=4, help="1-based BED column to search labels in (default: 4)")
    ap.add_argument("--invert", action="store_true", help="Mask everything EXCEPT matching labels (i.e., mask the complement of matched intervals)")
    ap.add_argument("--merge-distance", type=int, default=0, help="Merge gaps <= this bp (default: 0)")
    ap.add_argument("--padding", type=int, default=0, help="Expand each interval by +/- bp before merge (default: 0)")
    ap.add_argument("--softmask", action="store_true", help="Use soft-masking (lowercase) instead of 'N'")
    ap.add_argument("--mask-char", default="N", help="Hard mask character (default: 'N')")
    args = ap.parse_args()

    labels = [x.strip() for x in args.labels.split(",")] if args.labels else []

    # 1) Collect intervals from BED by label match (no invert here)
    bed_intervals = parse_bed(args.bed, args.column, labels)

    masked_records = []
    # 2) For each FASTA sequence: padding -> merge -> (invert? complement) -> mask
    for name, seq in fasta_iter(args.fasta):
        ivs = bed_intervals.get(name, [])
        if ivs:
            ivs = apply_padding(ivs, args.padding, len(seq))
            ivs = merge_intervals(ivs, max_gap=args.merge_distance)
        else:
            ivs = []

        if args.invert:
            # Keep matched intervals; therefore mask the complement.
            ivs_to_mask = complement_intervals(ivs, len(seq))
        else:
            # Mask only matched intervals.
            ivs_to_mask = ivs

        masked = mask_sequence(seq, ivs_to_mask, softmask=args.softmask, mask_char=args.mask_char)
        masked_records.append((name, masked))

    # 3) Write out FASTA
    write_fasta(args.out, masked_records, width=60)

if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        pass
    except Exception as e:
        sys.stderr.write(f"ERROR: {e}\n")
        sys.exit(1)
