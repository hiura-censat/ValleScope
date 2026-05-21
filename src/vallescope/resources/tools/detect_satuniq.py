#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
detect_satuniq.py (v1.2, mypy-friendly)

What it does
------------
  - Build k-mer DB with meryl
  - Build 'limited' DB with meryl (<= max-count; 1 means exactly unique)
  - (A) Compute UniqueKmerFrac per anchor (BED6)
  - (B) (optional) Genome-wide uniqueness track from `meryl-lookup -bed-runs`
        -> bedGraph (0/1) and bigWig

Key changes vs v1.1
-------------------
  - Stable file names (no "k31.le5" suffixes). Artifacts:
      <prefix>.meryl
      <prefix>.limited.meryl
      <prefix>.stats.txt
      <prefix>.limited.stats.txt
      <prefix>.anchors.fa
      <prefix>.anchors.limited.existence.txt
      <prefix>.anchors.limited.UniqueKmerFrac.tsv
      <prefix>.limited.runs.bed
      <prefix>.limited.runs.bedGraph
      <prefix>.limited.runs.bw
      <prefix>.limited.filtered.ids.txt                (if --min-ukf)
      <prefix>.limited.filtered.anchors.bed            (if --min-ukf)
  - Optional --tag lets you disambiguate runs without exposing k/le in names:
      <prefix>.<tag>.limited.*  etc.

Requirements
------------
  - Executables on PATH: meryl, meryl-lookup, samtools
  - Python: pysam  (pip install pysam)
  - For bigWig: bedGraphToBigWig (optional; bedGraphは常に出力可)
"""

from __future__ import annotations
import argparse
import os
import shutil
import subprocess
import tempfile
import math
from typing import Dict, List, Tuple, Optional, Set

# ------------------ optional deps ------------------
try:
    import pysam
except Exception:
    pysam = None


# ------------------ utils ------------------

def check_exec(name: str) -> bool:
    return shutil.which(name) is not None

def run_cmd(cmd: List[str], cwd: Optional[str] = None) -> None:
    p = subprocess.run(cmd, cwd=cwd)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed ({p.returncode}): {' '.join(cmd)}")

def log(msg: str) -> None:
    print(f"[satuniq] {msg}", flush=True)

def revcomp(seq: str) -> str:
    # 上流ツール(pysam)が大文字を返すことが多いので、基本upperに寄せる
    seq = seq.upper()
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(comp)[::-1]

def ensure_pysam() -> None:
    if pysam is None:
        raise RuntimeError("pysam is required. Install: pip install pysam")


# ---------- Step 1 & 2: meryl DB build ----------

def meryl_count_kmer(genome_fa: str, k: int, out_db: str, threads: int = 8, memory_gb: Optional[int] = None) -> None:
    cmd = ["meryl", "count", f"k={k}", f"threads={threads}"]
    if memory_gb is not None:
        # meryl "memory" は MB 指定
        cmd.append(f"memory={memory_gb*1024}")
    cmd += [genome_fa, "output", out_db]
    log(f"Running: {' '.join(cmd)}")
    run_cmd(cmd)

def meryl_at_most(in_db: str, max_count: int, out_db: str) -> None:
    """
    Extract k-mers with count <= max_count into 'out_db'.
    If max_count == 1, it's equivalent to unique-only (equal-to 1).
    """
    if max_count <= 1:
        cmd = ["meryl", "equal-to", in_db, "1", "output", out_db]
    else:
        cmd = ["meryl", "at-most", in_db, str(max_count), "output", out_db]
    log(f"Running: {' '.join(cmd)}")
    run_cmd(cmd)

def save_meryl_statistics(db_path: str, out_txt: str) -> None:
    """Dump 'meryl statistics' into text file."""
    log(f"[stats] meryl statistics {db_path} -> {out_txt}")
    p = subprocess.run(["meryl", "statistics", db_path], capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"meryl statistics failed on {db_path}\n{p.stderr}")
    with open(out_txt, "w") as w:
        w.write(p.stdout)


# ---------- Step 3A: UniqueKmerFrac per anchor ----------

def read_bed6(path: str) -> List[Tuple[str,int,int,str,str,str]]:
    rows: List[Tuple[str,int,int,str,str,str]] = []
    with open(path) as f:
        for ln, line in enumerate(f, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                raise ValueError(f"BED line {ln} has <4 columns: {line}")
            chrom = str(parts[0])
            start = int(parts[1])
            end   = int(parts[2])
            name  = str(parts[3])
            score = str(parts[4]) if len(parts) >= 5 else "0"
            strand= str(parts[5]) if len(parts) >= 6 else "+"
            rows.append((chrom, start, end, name, score, strand))
    return rows

def make_anchors_fasta(genome_fa: str, bed6_path: str, fasta_out: str) -> None:
    ensure_pysam()
    fa = pysam.FastaFile(genome_fa) 
    rows = read_bed6(bed6_path)
    with open(fasta_out, "w") as out:
        for chrom, start, end, name, score, strand in rows:
            seq = fa.fetch(chrom, start, end) 
            if strand == "-":
                seq = revcomp(seq)
            out.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")
    log(f"Wrote anchors FASTA: {fasta_out} ({len(rows)} sequences)")

def meryl_lookup_existence(unique_db: str, anchors_fa: str, out_txt: str) -> None:
    cmd = ["meryl-lookup", "-existence", "-sequence", anchors_fa, "-mers", unique_db]
    log(f"Running: {' '.join(cmd)} > {out_txt}")
    with open(out_txt, "w") as out:
        p = subprocess.run(cmd, stdout=out)
        if p.returncode != 0:
            raise RuntimeError(f"meryl-lookup -existence failed ({p.returncode})")

def parse_existence_to_unique_frac(existence_txt: str, k: int, tsv_out: str) -> None:
    """
    Convert 'meryl-lookup -existence' output into per-anchor UniqueKmerFrac TSV.
    Supported formats:
      A) 4-column: anchor, NumKmers, TotalUniqueInDB, NumUniqueInAnchor
      B) bit-string: >anchor header lines followed by 0/1 lines

    Output columns:
      anchor_id\tUniqueKmerFrac\tNumUnique\tNumKmers
    """
    import re

    def write_header(w) -> None:
        w.write("anchor_id\tUniqueKmerFrac\tNumUnique\tNumKmers\n")

    with open(existence_txt) as f:
        first_nonempty: Optional[str] = None
        pos = f.tell()
        for line in f:
            if line.strip():
                first_nonempty = line.strip()
                break
        f.seek(pos)

        with open(tsv_out, "w") as out:
            write_header(out)

            # ---- A) summarized 4-col format ----
            if first_nonempty and not first_nonempty.startswith(">"):
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    if len(parts) < 4:
                        continue
                    anc = parts[0]
                    try:
                        num_kmers   = int(parts[1])
                        num_unique  = int(parts[3])
                    except ValueError:
                        continue
                    frac = "NaN" if num_kmers <= 0 else f"{num_unique/num_kmers:.6f}"
                    out.write(f"{anc}\t{frac}\t{num_unique}\t{num_kmers}\n")

            # ---- B) bit-string format ----
            else:
                name: Optional[str] = None
                bits_count = 0
                ones_count = 0

                def flush() -> None:
                    if name is None:
                        return
                    frac = "NaN" if bits_count <= 0 else f"{ones_count/bits_count:.6f}"
                    out.write(f"{name}\t{frac}\t{ones_count}\t{bits_count}\n")

                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        flush()
                        name = line[1:].split()[0]
                        bits_count = 0
                        ones_count = 0
                    else:
                        for ch in re.findall(r"[01]", line):
                            bits_count += 1
                            if ch == "1":
                                ones_count += 1
                flush()

    log(f"Wrote UniqueKmerFrac per anchor: {tsv_out}")


# ---------- Step 3B: Genome-wide uniqueness track via -bed-runs ----------

def samtools_faidx_if_needed(genome_fa: str) -> str:
    fai = genome_fa + ".fai"
    if not os.path.exists(fai):
        if not check_exec("samtools"):
            raise RuntimeError("samtools not found; required to build FASTA index (.fai)")
        run_cmd(["samtools", "faidx", genome_fa])
    return fai

def load_chrom_sizes_from_fai(fai_path: str) -> Dict[str,int]:
    sizes: Dict[str,int] = {}
    with open(fai_path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.split("\t")
            sizes[str(parts[0])] = int(parts[1])
    return sizes

def meryl_lookup_bed_runs(unique_db: str, genome_fa: str, out_bed: str) -> None:
    # Outputs merged runs (overlapping unique k-mers combined) per contig as BED
    cmd = ["meryl-lookup", "-bed-runs", "-sequence", genome_fa, "-output", out_bed, "-mers", unique_db]
    log(f"Running: {' '.join(cmd)}")
    run_cmd(cmd)

def _merge_intervals(iv: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """mypy-friendly non-Optional interval merger"""
    if not iv:
        return []
    iv.sort(key=lambda x: (x[0], x[1]))
    merged: List[Tuple[int, int]] = []
    cs, ce = iv[0]
    for s, e in iv[1:]:
        if s <= ce:
            if e > ce:
                ce = e
        else:
            merged.append((cs, ce))
            cs, ce = s, e
    merged.append((cs, ce))
    return merged

def write_bedgraph_and_bigwig_from_runs_bed(
    runs_bed: str,
    chrom_sizes: Dict[str,int],
    bedgraph_out: str,
    bigwig_out: Optional[str] = None
) -> None:
    """
    Consume -bed-runs BED (unique segments), generate:
      - bedGraph (0/1) that spans the whole genome/contig (0 for non-unique, 1 for unique segments)
      - bigWig (optional) using bedGraphToBigWig
    """
    # Read and collect intervals
    per_chrom: Dict[str, List[Tuple[int,int]]] = {}
    with open(runs_bed) as f:
        for line in f:
            if not line.strip() or line.startswith(("#","track","browser")):
                continue
            chrom_s, s_s, e_s, *_rest = line.rstrip("\n").split("\t")
            s = int(s_s)
            e = int(e_s)
            if e <= s:
                continue
            chrom = str(chrom_s)
            per_chrom.setdefault(chrom, []).append((s, e))

    # sort & merge (safety)
    for chrom in list(per_chrom.keys()):
        per_chrom[chrom] = _merge_intervals(per_chrom[chrom])

    # emit bedGraph (0/1 across full contigs)
    with open(bedgraph_out, "w") as out:
        out.write('track type=bedGraph name="unique_kmer_runs" description="unique k-mer runs (1=unique,0=non-unique)"\n')
        for chrom, size in chrom_sizes.items():
            runs = per_chrom.get(chrom, [])
            cursor = 0
            for s, e in runs:
                if s > cursor:
                    out.write(f"{chrom}\t{cursor}\t{s}\t0\n")
                out.write(f"{chrom}\t{s}\t{e}\t1\n")
                cursor = e
            if cursor < size:
                out.write(f"{chrom}\t{cursor}\t{size}\t0\n")
    log(f"Wrote bedGraph: {bedgraph_out}")

    if bigwig_out:
        # chrom.sizes temp file
        with tempfile.NamedTemporaryFile("w", delete=False) as tmp:
            for c, sz in chrom_sizes.items():
                tmp.write(f"{c}\t{sz}\n")
            chrom_sizes_path = tmp.name
        try:
            if check_exec("bedGraphToBigWig"):
                run_cmd(["bedGraphToBigWig", bedgraph_out, chrom_sizes_path, bigwig_out])
                log(f"Wrote bigWig (bedGraphToBigWig): {bigwig_out}")
            else:
                log("[warn] bedGraphToBigWig not found; skipped bigWig creation (bedGraph is available).")
        finally:
            try:
                os.unlink(chrom_sizes_path)
            except Exception:
                pass

# ---------- UKF utilities ----------

def load_ukf_ids(unique_frac_tsv: str, min_ukf: float) -> Set[str]:
    """
    Read UniqueKmerFrac TSV ('anchor_id\tUniqueKmerFrac\tNumUnique\tNumKmers')
    and return anchor IDs with frac >= min_ukf. NaN/invalid are ignored.
    """
    keep: Set[str] = set()
    with open(unique_frac_tsv) as f:
        _ = f.readline()  # header
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            anc, frac = parts[0], parts[1]
            try:
                v = float(frac)
            except Exception:
                continue
            if not math.isnan(v) and v >= min_ukf:
                keep.add(anc)
    return keep

def write_ids(ids: Set[str], out_path: str) -> None:
    with open(out_path, "w") as w:
        for anc in sorted(ids):
            w.write(anc + "\n")

def filter_bed_by_anchor_ids(in_bed: str, out_bed: str, keep_ids: Set[str]) -> int:
    """
    Keep lines whose 4th column (name) is in keep_ids.
    Return number of kept rows.
    """
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

def run_ukf_filter(prefix: str,
                   unique_frac_tsv: str,
                   anchors_bed: str,
                   min_ukf: float,
                   tag: str = "") -> None:
    """
    Filter anchors by UniqueKmerFrac, write:
      - IDs:   <prefix><tag>.limited.filtered.ids.txt
      - BED:   <prefix><tag>.limited.filtered.anchors.bed
    """
    ids = load_ukf_ids(unique_frac_tsv, min_ukf)
    ids_out = f"{prefix}{tag}.limited.filtered.ids.txt"
    bed_out = f"{prefix}{tag}.limited.filtered.anchors.bed"
    write_ids(ids, ids_out)
    kept = filter_bed_by_anchor_ids(anchors_bed, bed_out, ids)
    log(f"[filter] UniqueKmerFrac >= {min_ukf:.2f}: kept {kept} anchors -> {bed_out} / IDs -> {ids_out}")


# ------------------ main ------------------

def main() -> None:
    ap = argparse.ArgumentParser(description="Build unique-ish k-mer DB with meryl; compute UniqueKmerFrac and optional bigWig (no k/le in filenames).")
    ap.add_argument("--genome-fa", required=True, help="Reference FASTA")
    ap.add_argument("--anchors-bed", required=True, help="Anchor regions (BED6)")
    ap.add_argument("--k", type=int, default=21, help="k-mer length for meryl count")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--memory-gb", type=int, default=None)
    ap.add_argument("--out-prefix", required=True, help="Prefix for all outputs (no parameter suffixes will be appended)")
    ap.add_argument("--make-bigwig", action="store_true", help="Also build genome-wide uniqueness track (bedGraph+bigWig)")
    ap.add_argument("--max-kmer-count", type=int, choices=[1, 3, 5], default=1,
                    help="Maximum allowed k-mer count in 'limited' DB (1=exactly unique; 3=<=3; 5=<=5)")
    ap.add_argument("--min-ukf", type=float, default=None,
                    help="If set, filter anchors by UniqueKmerFrac >= this value (e.g., 0.90)")
    ap.add_argument("--tag", type=str, default="",
                    help="Optional disambiguation tag; if set, outputs become <prefix>.<tag>.*")
    args = ap.parse_args()

    # executables
    for exe in ["meryl", "meryl-lookup", "samtools"]:
        if not check_exec(exe):
            raise RuntimeError(f"Required executable not found in PATH: {exe}")
    if pysam is None:
        raise RuntimeError("pysam is required. Install: pip install pysam")

    # convenience
    k = args.k
    prefix_root = args.out_prefix
    tag = f".{args.tag}" if args.tag else ""
    prefix = f"{prefix_root}{tag}"
    maxc = args.max_kmer_count

    # [1] meryl count
    meryl_db = f"{prefix}.meryl"
    if not os.path.exists(meryl_db):
        log(f"[1/6] meryl count -> {meryl_db}")
        meryl_count_kmer(args.genome_fa, k, meryl_db, threads=args.threads, memory_gb=args.memory_gb)
    else:
        log(f"[1/6] Skip meryl count (exists): {meryl_db}")

    # [2] build limited DB (<= maxc)
    limited_db = f"{prefix}.limited.meryl"
    if not os.path.exists(limited_db):
        if maxc == 1:
            log(f"[2/6] meryl equal-to 1 -> {limited_db}")
        else:
            log(f"[2/6] meryl at-most {maxc} -> {limited_db}")
        meryl_at_most(meryl_db, maxc, limited_db)
    else:
        log(f"[2/6] Skip limited DB (exists): {limited_db}")

    # stats
    save_meryl_statistics(meryl_db,   f"{prefix}.stats.txt")
    save_meryl_statistics(limited_db, f"{prefix}.limited.stats.txt")

    # [3] anchors FASTA
    anchors_fa = f"{prefix}.anchors.fa"
    log(f"[3/6] Extract anchors FASTA: {anchors_fa}")
    make_anchors_fasta(args.genome_fa, args.anchors_bed, anchors_fa)

    # [4] existence on limited DB
    existence_txt = f"{prefix}.anchors.limited.existence.txt"
    log(f"[4/6] meryl-lookup -existence -> {existence_txt}")
    meryl_lookup_existence(limited_db, anchors_fa, existence_txt)

    # [5] UniqueKmerFrac per anchor
    unique_frac_tsv = f"{prefix}.anchors.limited.UniqueKmerFrac.tsv"
    log(f"[5/6] Parse existence -> {unique_frac_tsv}")
    parse_existence_to_unique_frac(existence_txt, k, unique_frac_tsv)

    # optional filter
    if args.min_ukf is not None:
        run_ukf_filter(prefix, unique_frac_tsv, args.anchors_bed, args.min_ukf, tag="")

    # [6] genome-wide runs and bigWig
    if args.make_bigwig:
        log("[6/6] Genome-wide unique track via -bed-runs")
        fai_path = samtools_faidx_if_needed(args.genome_fa)
        sizes = load_chrom_sizes_from_fai(fai_path)
        runs_bed    = f"{prefix}.limited.runs.bed"
        if not os.path.exists(runs_bed):
            meryl_lookup_bed_runs(limited_db, args.genome_fa, runs_bed)
        bedgraph_out = f"{prefix}.limited.runs.bedGraph"
        bigwig_out   = f"{prefix}.limited.runs.bw"
        write_bedgraph_and_bigwig_from_runs_bed(runs_bed, sizes, bedgraph_out, bigwig_out)
    else:
        log("[6/6] Skip bigWig (use --make-bigwig to enable)")

    log("Done.")


if __name__ == "__main__":
    main()
