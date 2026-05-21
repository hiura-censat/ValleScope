#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Filter anchors whose UniqueKmerFrac is NaN in <prefix>.anchors.limited.UniqueKmerFrac.tsv,
then emit a filtered TSV, BED, and FASTA that contain only the kept anchors.

Input TSV format (tab):
  anchor_id\tUniqueKmerFrac\tNumUnique\tNumKmers

Behavior:
  - Drop rows where UniqueKmerFrac is NaN (case-insensitive, supports 'NaN', 'nan', empty).
  - Optionally also apply a min-UKF threshold (--min-ukf).
  - Output:
      <out-prefix>.anchors.limited.UniqueKmerFrac.filtered.tsv
      <out-prefix>.filtered.anchors.bed
      <out-prefix>.filtered.anchors.fa
      <out-prefix>.filtered.ids.txt

No external dependencies required.
"""

import argparse
import math
from typing import Iterable, Tuple, Set

def log(msg: str) -> None:
    print(f"[filter-ukf] {msg}", flush=True)

# ---------- simple readers/writers ----------

def load_keep_ids_from_tsv(tsv_path: str, min_ukf: float = 0.0) -> Tuple[Set[str], int, int]:
    """
    Return (keep_ids, n_rows, n_dropped_nan).

    Keep if:
      - UniqueKmerFrac is a valid float (not NaN), and
      - (optional) >= min_ukf if provided
    """
    keep: Set[str] = set()
    n_rows = 0
    n_nan  = 0

    def _is_nan_token(s: str) -> bool:
        if s is None: 
            return True
        s = s.strip()
        if s == "": 
            return True
        # accept many NaN spellings
        return s.lower() in ("nan", "na", "null", "none")

    with open(tsv_path) as f:
        header = f.readline()
        if not header:
            raise RuntimeError(f"Empty TSV: {tsv_path}")
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            n_rows += 1
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            aid = parts[0].strip()
            frac_s = parts[1].strip()
            if _is_nan_token(frac_s):
                n_nan += 1
                continue
            try:
                frac = float(frac_s)
            except Exception:
                # malformed -> treat as NaN
                n_nan += 1
                continue
            if math.isnan(frac):
                n_nan += 1
                continue
            if (min_ukf is not None) and (frac < float(min_ukf)):
                continue
            if aid:
                keep.add(aid)
    return keep, n_rows, n_nan

def filter_bed_by_ids(in_bed: str, out_bed: str, keep_ids: Set[str]) -> int:
    kept = 0
    with open(in_bed) as fin, open(out_bed, "w") as fout:
        for line in fin:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4 and parts[3] in keep_ids:
                fout.write(line)
                kept += 1
    return kept

def stream_fasta_records(path: str) -> Iterable[Tuple[str, str]]:
    """
    Yield (name, seq) for each FASTA record.
    Name is the first token after '>'.
    """
    name = None
    buf: list[str] = [] 
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield (name, "".join(buf))
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        yield (name, "".join(buf))

def write_fasta(path: str, records: Iterable[Tuple[str, str]]) -> int:
    n = 0
    with open(path, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            # wrap at 60 for readability
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")
            n += 1
    return n

# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(description="Filter anchors with NaN UniqueKmerFrac from TSV, then emit filtered BED/FASTA.")
    ap.add_argument("--tsv", required=True, help="<prefix>.anchors.limited.UniqueKmerFrac.tsv")
    ap.add_argument("--bed", required=True, help="anchors BED6")
    ap.add_argument("--fa",  required=True, help="anchors FASTA (names must match BED col4)")
    ap.add_argument("--out-prefix", required=True, help="prefix of outputs")
    ap.add_argument("--min-ukf", type=float, default=None, help="optional threshold (keep if UKF >= this); NaN is always dropped")
    args = ap.parse_args()

    # 1) decide keep IDs
    keep_ids, n_rows, n_nan = load_keep_ids_from_tsv(args.tsv, args.min_ukf)
    log(f"TSV rows (excluding header): {n_rows}, NaN dropped: {n_nan}, keep IDs: {len(keep_ids)}")

    if not keep_ids:
        log("[WARN] No IDs to keep after filtering. Outputs will be empty files.")

    # 2) write filtered TSV (only kept IDs)
    tsv_out = f"{args.out_prefix}.tsv"
    kept_rows = 0
    with open(args.tsv) as fin, open(tsv_out, "w") as fout:
        header = fin.readline()
        fout.write(header if header else "anchor_id\tUniqueKmerFrac\tNumUnique\tNumKmers\n")
        for line in fin:
            if not line.strip(): 
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 1: 
                continue
            aid = parts[0].strip()
            if aid in keep_ids:
                fout.write(line)
                kept_rows += 1
    log(f"Wrote filtered TSV: {tsv_out} (rows: {kept_rows})")

    # 3) write filtered BED
    bed_out = f"{args.out_prefix}.bed"
    n_bed = filter_bed_by_ids(args.bed, bed_out, keep_ids)
    log(f"Wrote filtered BED: {bed_out} (rows: {n_bed})")

    # 4) write filtered FASTA
    fa_out = f"{args.out_prefix}.fa"
    def _iter_kept():
        for name, seq in stream_fasta_records(args.fa):
            if name in keep_ids:
                yield (name, seq)
    n_fa = write_fasta(fa_out, _iter_kept())
    log(f"Wrote filtered FASTA: {fa_out} (records: {n_fa})")

    # 5) write IDs list (便利)
    ids_out = f"{args.out_prefix}.ids.txt"
    with open(ids_out, "w") as w:
        for aid in sorted(keep_ids):
            w.write(aid + "\n")
    log(f"Wrote kept IDs: {ids_out}")

if __name__ == "__main__":
    main()
