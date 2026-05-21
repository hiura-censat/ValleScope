#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
vallescope_bundles_to_cigar_paf.py  (AUTO-CONTIG SELECTION VERSION)

Goal
----
Take ValleScope bundles.tsv (block-level correspondences) and generate a "real"
base-level, CIGAR-tagged PAF suitable for feeding to SyRI (or similar),
by aligning each bundle interval pair with minimap2 and then lifting the
segment-local alignment coordinates back to chromosome-global coordinates.

Key points
----------
- PAF semantics: q = QUERY assembly, t = TARGET/REF assembly
- qname is query contig name, tname is ref contig name
- Handles bundle strand '-' by revcomping the extracted query interval before alignment,
  then lifting segment coords back to original query chromosome coordinates.
- Uses samtools faidx for random access (samtools required).
- Uses minimap2 -c --cs=long to emit cg:Z: and cs:Z: tags (minimap2 required).

NEW in this version
-------------------
- Auto contig selection / mapping:
    * If the FASTA has exactly one contig, it's used automatically.
    * Otherwise, a bundle chrom like "chr1" can match FASTA contigs like
        "chr1:133000000-134000000" or "chr1:...__mut_seed1"
      via a "base chrom" normalization:
        base = contig.split(':')[0].split('__')[0]
    * If multiple candidates share the same base, choose the longest contig.

Outputs
-------
- A merged PAF with cg/cs tags and global coordinates.
- Keeps original bundle id in tag: id:Z:<bundle_id>
- For failures (if --keep-failed), writes placeholder lines with na:i:<code>:
    na:i:1  contig name not found (after mapping)
    na:i:2  faidx returns empty (coords invalid/out-of-range etc.)
    na:i:3  minimap2 returned no alignments

Example
-------
python3 vallescope_bundles_to_cigar_paf.py \
  --bundles-tsv vs.bundles.tsv \
  --qry-fa chr1_HSat2_1M_1.masked.fasta \
  --ref-fa chr1_HSat2_1M_1.seed1.snv0p001.indel0p0000.mut.masked.fasta \
  --out-paf out.syri.paf \
  --threads 8 \
  --minimap2-preset asm5 \
  --flank 0 \
  --max-bundle-bp 20000000 \
  --keep-failed

Notes
-----
- This aligns each bundle interval pair as a single chunk. If bundles contain internal
  rearrangements, consider splitting bundles first (not implemented here).
"""

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


def die(msg: str, code: int = 2) -> None:
    print(f"[ERROR] {msg}", file=sys.stderr)
    raise SystemExit(code)


def check_exe(name: str) -> str:
    p = shutil.which(name)
    if not p:
        die(f"Executable not found in PATH: {name}")
    return p


def revcomp(seq: str) -> str:
    t = str.maketrans("ACGTacgtnN", "TGCAtgcaNn")
    return seq.translate(t)[::-1]


def ensure_fai(samtools: str, fasta: str) -> str:
    fai = fasta + ".fai"
    if os.path.exists(fai) and os.path.getsize(fai) > 0:
        return fai
    subprocess.run([samtools, "faidx", fasta], check=True)
    if not os.path.exists(fai) or os.path.getsize(fai) == 0:
        die(f"Failed to create FASTA index: {fai}")
    return fai


def load_fai_lengths(fai_path: str) -> Dict[str, int]:
    d: Dict[str, int] = {}
    with open(fai_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            d[cols[0]] = int(cols[1])
    return d


def faidx_fetch_sequence(
    samtools: str, fasta: str, contig: str, start0: int, end0: int
) -> str:
    """Fetch 0-based half-open [start0,end0) from fasta using samtools faidx.
    samtools region is 1-based inclusive. So request (start0+1)-(end0).
    """
    if end0 <= start0:
        return ""
    region = f"{contig}:{start0+1}-{end0}"
    p = subprocess.run(
        [samtools, "faidx", fasta, region],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if p.returncode != 0:
        return ""
    lines = []
    for ln in p.stdout.splitlines():
        if ln.startswith(">"):
            continue
        lines.append(ln.strip())
    return "".join(lines)


@dataclass
class BundleRow:
    bundle_id: str
    ref_name: str
    qry_name: str
    ref_chr: str
    qry_chr: str
    ref_start: int
    ref_end: int
    qry_start: int
    qry_end: int
    strand: str  # '+' or '-'


def _get(row: Dict[str, str], *keys: str, default: str = "") -> str:
    for k in keys:
        if k in row and row[k] != "":
            return row[k]
    return default


def read_bundles_tsv(path: str) -> List[BundleRow]:
    opener = gzip.open if path.endswith((".gz", ".bgz")) else open
    with opener(path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required_any = {"chrom", "ref_start", "ref_end", "qry_start", "qry_end", "strand"}
        for k in required_any:
            if k not in reader.fieldnames:
                die(f"bundles.tsv missing required column: {k}. Found: {reader.fieldnames}")

        out: List[BundleRow] = []
        for i, r in enumerate(reader):
            bundle_id = _get(r, "bundle_id", default=f"BUND{i+1:06d}")
            ref_name = _get(r, "ref_name", default="ref")
            qry_name = _get(r, "qry_name", default="qry")

            # ValleScope bundles.tsv usually has one "chrom" column.
            chrom = _get(r, "chrom", default="")
            ref_chr = _get(r, "ref_chrom", "ref_chr", default=chrom)
            qry_chr = _get(r, "qry_chrom", "qry_chr", default=chrom)

            try:
                rs = int(float(_get(r, "ref_start")))
                re = int(float(_get(r, "ref_end")))
                qs = int(float(_get(r, "qry_start")))
                qe = int(float(_get(r, "qry_end")))
            except Exception:
                die(f"Failed to parse coordinates at row {i+1}: {r}")

            strand = _get(r, "strand", default="+")
            if strand not in {"+", "-"}:
                strand = "+"

            out.append(
                BundleRow(
                    bundle_id=bundle_id,
                    ref_name=ref_name,
                    qry_name=qry_name,
                    ref_chr=str(ref_chr),
                    qry_chr=str(qry_chr),
                    ref_start=rs,
                    ref_end=re,
                    qry_start=qs,
                    qry_end=qe,
                    strand=strand,
                )
            )
    return out


def contig_base(name: str) -> str:
    """Normalize contig name to a base identifier for mapping bundles->FASTA contigs.
    Examples:
      chr1                         -> chr1
      chr1:133-134                 -> chr1
      chr1:133-134__mut_seed1      -> chr1
      chr1__mut_seed1              -> chr1
    """
    s = str(name)
    s = s.split(":")[0]
    s = s.split("__")[0]
    return s


def build_base_to_contigs(contig_len: Dict[str, int]) -> Dict[str, List[str]]:
    m: Dict[str, List[str]] = {}
    for c in contig_len.keys():
        b = contig_base(c)
        m.setdefault(b, []).append(c)
    # deterministic: sort by length desc, then name
    for b, cs in m.items():
        cs.sort(key=lambda x: (-int(contig_len.get(x, 0)), str(x)))
    return m


def pick_contig_for_bundle(
    bundle_chr: str,
    contig_len: Dict[str, int],
    base_to_contigs: Dict[str, List[str]],
    forced_contig: str = "",
    one_contig_fallback: Optional[str] = None,
) -> str:
    """Choose the actual FASTA contig name to use for a given bundle 'chrom'."""
    if forced_contig:
        return forced_contig

    # exact hit
    if bundle_chr in contig_len:
        return bundle_chr

    # single contig fallback
    if one_contig_fallback:
        return one_contig_fallback

    # base match
    b = contig_base(bundle_chr)
    cands = base_to_contigs.get(b, [])
    if cands:
        return cands[0]  # longest (pre-sorted)
    return ""


def parse_paf_line(line: str) -> Tuple[List[str], List[str]]:
    cols = line.rstrip("\n").split("\t")
    base = cols[:12]
    tags = cols[12:] if len(cols) > 12 else []
    return base, tags


def choose_best_paf(lines: List[str]) -> Optional[str]:
    """Pick a single best alignment from minimap2 output.
    Prefer:
      1) highest mapq (col 12)
      2) then largest alnlen (col 11)
    """
    best = None
    best_key = (-1, -1)
    for ln in lines:
        if not ln.strip():
            continue
        base, _ = parse_paf_line(ln)
        if len(base) < 12:
            continue
        try:
            alnlen = int(base[10])
            mapq = int(base[11])
        except Exception:
            continue
        key = (mapq, alnlen)
        if key > best_key:
            best_key = key
            best = ln
    return best


def lift_query_coords_from_segment(
    seg_qstart: int, seg_qend: int, bundle_qs: int, bundle_qe: int, strand: str
) -> Tuple[int, int]:
    """Convert segment-local qstart/qend to chromosome-global.
    If strand '+': global = bundle_qs + seg_*
    If strand '-': the query segment was reverse-complemented before alignment.
      segment coords [a,b) in reversed correspond to original:
        [bundle_qe - b, bundle_qe - a)
    """
    if strand == "+":
        return bundle_qs + seg_qstart, bundle_qs + seg_qend
    gstart = bundle_qe - seg_qend
    gend = bundle_qe - seg_qstart
    if gstart > gend:
        gstart, gend = gend, gstart
    return gstart, gend


def lift_ref_coords_from_segment(seg_tstart: int, seg_tend: int, bundle_rs: int) -> Tuple[int, int]:
    return bundle_rs + seg_tstart, bundle_rs + seg_tend


def run_minimap2_on_segments(
    minimap2: str,
    preset: str,
    threads: int,
    q_fa: str,
    t_fa: str,
) -> List[str]:
    cmd = [
        minimap2,
        "-x",
        preset,
        "-t",
        str(int(threads)),
        "--secondary=no",
        "--eqx",
        "-c",         # emit cg:Z:
        # "--cs=long",  # emit cs:Z:
        t_fa,
        q_fa,
    ]
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        return []
    outs = [ln for ln in p.stdout.splitlines() if ln and not ln.startswith("#")]
    return outs


def write_temp_fasta(path: str, name: str, seq: str) -> None:
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i : i + 80] + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Add real cg/cs to ValleScope bundles by aligning bundle intervals with minimap2 and lifting to global coords."
    )
    ap.add_argument("--bundles-tsv", required=True, help="ValleScope bundles.tsv (input)")
    ap.add_argument("--ref-fa", required=True, help="Reference assembly FASTA (target for PAF)")
    ap.add_argument("--qry-fa", required=True, help="Query assembly FASTA (query for PAF)")
    ap.add_argument("--out-paf", required=True, help="Output PAF (SyRI-ready)")
    ap.add_argument("--threads", type=int, default=4, help="Threads for minimap2")
    ap.add_argument("--minimap2-preset", default="asm5", help="minimap2 preset (asm5/asm10/asm20)")
    ap.add_argument("--flank", type=int, default=0, help="Add flank bp around bundle intervals before alignment")
    ap.add_argument("--max-bundle-bp", type=int, default=3000000, help="Skip bundles with min(span) > this")
    ap.add_argument("--keep-failed", action="store_true", help="Write placeholder PAF lines for failed bundles")
    ap.add_argument("--prefix-contig-names", action="store_true",
                    help="Use qname/tname as <sample>:<contig> instead of contig only")

    # NEW: optional forced contigs (override auto mapping)
    ap.add_argument("--ref-contig", default="", help="Force ref contig name (overrides auto mapping)")
    ap.add_argument("--qry-contig", default="", help="Force query contig name (overrides auto mapping)")

    args = ap.parse_args()

    samtools = check_exe("samtools")
    minimap2 = check_exe("minimap2")

    # Ensure indices and lengths
    ref_fai = ensure_fai(samtools, args.ref_fa)
    qry_fai = ensure_fai(samtools, args.qry_fa)
    ref_len = load_fai_lengths(ref_fai)
    qry_len = load_fai_lengths(qry_fai)

    # Build base mapping
    ref_base_to_contigs = build_base_to_contigs(ref_len)
    qry_base_to_contigs = build_base_to_contigs(qry_len)

    # One-contig fallback if available
    ref_single = next(iter(ref_len.keys())) if len(ref_len) == 1 else None
    qry_single = next(iter(qry_len.keys())) if len(qry_len) == 1 else None
    if ref_single:
        print(f"[info] ref FASTA has 1 contig -> will auto-use: {ref_single}", file=sys.stderr)
    if qry_single:
        print(f"[info] qry FASTA has 1 contig -> will auto-use: {qry_single}", file=sys.stderr)

    bundles = read_bundles_tsv(args.bundles_tsv)
    if not bundles:
        die("No bundles found.")

    # Output
    outdir = os.path.dirname(os.path.abspath(args.out_paf))
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    kept = 0
    skipped = 0
    failed = 0

    with open(args.out_paf, "w") as fout, tempfile.TemporaryDirectory(prefix="vallescope_cg_") as td:
        for b in bundles:
            rs, re = int(b.ref_start), int(b.ref_end)
            qs, qe = int(b.qry_start), int(b.qry_end)
            strand = b.strand

            # AUTO contig mapping
            ref_contig = pick_contig_for_bundle(
                bundle_chr=b.ref_chr,
                contig_len=ref_len,
                base_to_contigs=ref_base_to_contigs,
                forced_contig=str(args.ref_contig),
                one_contig_fallback=ref_single,
            )
            qry_contig = pick_contig_for_bundle(
                bundle_chr=b.qry_chr,
                contig_len=qry_len,
                base_to_contigs=qry_base_to_contigs,
                forced_contig=str(args.qry_contig),
                one_contig_fallback=qry_single,
            )

            if (not ref_contig) or (not qry_contig) or (ref_contig not in ref_len) or (qry_contig not in qry_len):
                failed += 1
                if args.keep_failed:
                    qname = f"{b.qry_name}:{(qry_contig or b.qry_chr)}" if args.prefix_contig_names else (qry_contig or b.qry_chr)
                    tname = f"{b.ref_name}:{(ref_contig or b.ref_chr)}" if args.prefix_contig_names else (ref_contig or b.ref_chr)
                    qL = int(qry_len.get(qry_contig or "", 0))
                    tL = int(ref_len.get(ref_contig or "", 0))
                    fout.write(
                        f"{qname}\t{qL}\t{qs}\t{qe}\t{strand}\t"
                        f"{tname}\t{tL}\t{rs}\t{re}\t"
                        f"0\t0\t0\tid:Z:{b.bundle_id}\tna:i:1\n"
                    )
                continue

            # Apply flank and clip
            flank = int(args.flank)
            rs2 = max(0, rs - flank)
            re2 = min(ref_len[ref_contig], re + flank)
            qs2 = max(0, qs - flank)
            qe2 = min(qry_len[qry_contig], qe + flank)

            ref_span = max(0, re2 - rs2)
            qry_span = max(0, qe2 - qs2)
            if min(ref_span, qry_span) <= 0:
                skipped += 1
                continue
            if min(ref_span, qry_span) > int(args.max_bundle_bp):
                skipped += 1
                continue

            # Fetch sequences
            ref_seq = faidx_fetch_sequence(samtools, args.ref_fa, ref_contig, rs2, re2)
            qry_seq = faidx_fetch_sequence(samtools, args.qry_fa, qry_contig, qs2, qe2)
            if not ref_seq or not qry_seq:
                failed += 1
                if args.keep_failed:
                    qname = f"{b.qry_name}:{qry_contig}" if args.prefix_contig_names else qry_contig
                    tname = f"{b.ref_name}:{ref_contig}" if args.prefix_contig_names else ref_contig
                    qL = int(qry_len.get(qry_contig, 0))
                    tL = int(ref_len.get(ref_contig, 0))
                    fout.write(
                        f"{qname}\t{qL}\t{qs2}\t{qe2}\t{strand}\t"
                        f"{tname}\t{tL}\t{rs2}\t{re2}\t"
                        f"0\t0\t0\tid:Z:{b.bundle_id}\tna:i:2\n"
                    )
                continue

            # If inversion bundle, revcomp query sequence before alignment
            qry_seq_aln = revcomp(qry_seq) if strand == "-" else qry_seq

            # Write temp FASTAs (segment-local coordinate system)
            q_fa = os.path.join(td, "q.fa")
            t_fa = os.path.join(td, "t.fa")
            write_temp_fasta(q_fa, "QSEG", qry_seq_aln)
            write_temp_fasta(t_fa, "TSEG", ref_seq)

            paf_lines = run_minimap2_on_segments(
                minimap2=minimap2,
                preset=str(args.minimap2_preset),
                threads=int(args.threads),
                q_fa=q_fa,
                t_fa=t_fa,
            )
            best = choose_best_paf(paf_lines)
            if best is None:
                failed += 1
                if args.keep_failed:
                    qname = f"{b.qry_name}:{qry_contig}" if args.prefix_contig_names else qry_contig
                    tname = f"{b.ref_name}:{ref_contig}" if args.prefix_contig_names else ref_contig
                    qL = int(qry_len.get(qry_contig, 0))
                    tL = int(ref_len.get(ref_contig, 0))
                    fout.write(
                        f"{qname}\t{qL}\t{qs2}\t{qe2}\t{strand}\t"
                        f"{tname}\t{tL}\t{rs2}\t{re2}\t"
                        f"0\t0\t0\tid:Z:{b.bundle_id}\tna:i:3\n"
                    )
                continue

            base, tags = parse_paf_line(best)
            if len(base) < 12:
                failed += 1
                continue

            # Segment-local coords from minimap2 output
            seg_qstart = int(base[2])
            seg_qend = int(base[3])
            seg_strand = base[4]  # usually '+', because we pre-revcomp for strand '-'
            seg_tstart = int(base[7])
            seg_tend = int(base[8])
            nmatch = base[9]
            alnlen = base[10]
            mapq = base[11]

            # Lift to original chromosome coordinates
            qg_s, qg_e = lift_query_coords_from_segment(seg_qstart, seg_qend, qs2, qe2, strand)
            tg_s, tg_e = lift_ref_coords_from_segment(seg_tstart, seg_tend, rs2)

            # Final names & lengths (use actual FASTA contig names)
            qname = f"{b.qry_name}:{qry_contig}" if args.prefix_contig_names else qry_contig
            tname = f"{b.ref_name}:{ref_contig}" if args.prefix_contig_names else ref_contig
            qL = int(qry_len[qry_contig])
            tL = int(ref_len[ref_contig])

            final_strand = strand

            # Ensure tags include id:Z:<bundle_id>. Remove existing id:Z to avoid duplicates.
            new_tags = []
            for tg in tags:
                if tg.startswith("id:Z:"):
                    continue
                new_tags.append(tg)
            new_tags.append(f"id:Z:{b.bundle_id}")
            new_tags.append(f"bs:i:{int(min(ref_span, qry_span))}")
            new_tags.append(f"fl:i:{int(flank)}")
            new_tags.append(f"ss:Z:{seg_strand}")
            # Extra debug: original bundles chrom strings (optional)
            new_tags.append(f"rb:Z:{ref_contig}")
            new_tags.append(f"qb:Z:{qry_contig}")


            fout.write(
                f"{qname}\t{qL}\t{qg_s}\t{qg_e}\t{final_strand}\t"
                f"{tname}\t{tL}\t{tg_s}\t{tg_e}\t"
                f"{nmatch}\t{alnlen}\t{mapq}\t" + "\t".join(new_tags) + "\n"
            )
            kept += 1

    print(f"[done] wrote: {args.out_paf}", file=sys.stderr)
    print(f"[stats] kept={kept} skipped={skipped} failed={failed}", file=sys.stderr)


if __name__ == "__main__":
    main()