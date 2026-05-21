#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
paf_gap_patch_unialigner_to_syri_paf.py

(A) 入力PAF（ValleScope/minimap2等のブロック）から隣接ブロック間ギャップ候補を抽出
(B) ギャップ両側にflankを付けて ref/qry から配列パッチを切り出し
(C) UniAlignerで局所アライン（パッチ同士）
(D) UniAlignerのCIGARを cg:Z として付与したPAF（SyRI用）を出力
    - 追加PAF行は tp:Z:gapua でタグ付け
    - ギャップ部分だけにCIGARをトリム（keep_overlapで少しアンカー側に食い込む）

raw / simplified の2層表現:
  - cg:Z:  = simplified cigar（SyRI/downstream 用）
  - cr:Z:  = raw trimmed cigar（UniAligner の生CIGAR）
  - sm:Z:  = simplify mode
      * raw
      * 1I
      * 1D
      * splitI

今回の重要変更:
  - split 判定は raw ではなく、**simplified が =I= なら split** する
  - つまり raw が複雑でも、simplify 後に pure insertion になれば split する
  - pure deletion の split はまだ実装しない（必要なら後で追加可能）

前提：
  - samtools がPATHにある
  - UniAlignerの実行コマンドは --unialigner-cmd にテンプレを渡す（{REF} {QRY} {OUT} を置換）
    UniAlignerは {OUT} に CIGAR を1行で出力すること
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


# ----------------------------- util -----------------------------

def die(msg: str, code: int = 2) -> None:
    print(f"[ERROR] {msg}", file=sys.stderr)
    raise SystemExit(code)

def check_exe(name: str) -> str:
    p = shutil.which(name)
    if not p:
        die(f"Executable not found in PATH: {name}")
    return p

def ensure_fai(samtools: str, fasta: str) -> str:
    fai = fasta + ".fai"
    subprocess.run([samtools, "faidx", fasta], check=True)
    if not os.path.exists(fai) or os.path.getsize(fai) == 0:
        die(f"Failed to create FASTA index: {fai}")
    return fai

def load_fai_lengths(fai_path: str) -> Dict[str, int]:
    out: Dict[str, int] = {}
    with open(fai_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 2:
                out[cols[0]] = int(cols[1])
    return out

def faidx_fetch_sequence(samtools: str, fasta: str, contig: str, start0: int, end0: int) -> str:
    """Fetch 0-based half-open [start0,end0) using samtools faidx (1-based inclusive region)."""
    if end0 <= start0:
        return ""
    region = f"{contig}:{start0+1}-{end0}"
    p = subprocess.run([samtools, "faidx", fasta, region],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        return ""
    seq = []
    for ln in p.stdout.splitlines():
        if ln.startswith(">"):
            continue
        seq.append(ln.strip())
    return "".join(seq)

def write_fasta(path: str, name: str, seq: str) -> None:
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")

def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


# ----------------------------- PAF model -----------------------------

@dataclass
class PafRec:
    qname: str
    qlen: int
    qstart: int
    qend: int
    strand: str
    tname: str
    tlen: int
    tstart: int
    tend: int
    nmatch: int
    alen: int
    mapq: int
    tags: List[str]
    raw: str

def parse_paf_line(line: str) -> Optional[PafRec]:
    if not line.strip() or line.startswith("#"):
        return None
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 12:
        return None
    try:
        qname = cols[0]
        qlen = int(cols[1])
        qstart = int(cols[2])
        qend = int(cols[3])
        strand = cols[4]
        tname = cols[5]
        tlen = int(cols[6])
        tstart = int(cols[7])
        tend = int(cols[8])
        nmatch = int(cols[9])
        alen = int(cols[10])
        mapq = int(cols[11])
        tags = cols[12:] if len(cols) > 12 else []
    except Exception:
        return None
    return PafRec(qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend,
                  nmatch, alen, mapq, tags, raw=line.rstrip("\n"))


# ----------------------------- strand helpers -----------------------------

def forwardize_query_interval(qs: int, qe: int, qlen: int, strand: str) -> Tuple[int, int]:
    """Map query interval into forward coordinate along the query contig."""
    if strand == "+":
        return qs, qe
    return qlen - qe, qlen - qs

def forward_to_original_interval(f0: int, f1: int, qlen: int, strand: str) -> Tuple[int, int]:
    """Convert forward-coordinate interval [f0,f1) to original coordinate interval for faidx."""
    if strand == "+":
        return f0, f1
    return qlen - f1, qlen - f0


# ----------------------------- group & gap candidate -----------------------------

def group_and_sort_blocks(recs: List[PafRec]) -> Dict[Tuple[str, str, str], List[PafRec]]:
    groups: Dict[Tuple[str, str, str], List[PafRec]] = {}
    for r in recs:
        if r.strand not in {"+", "-"}:
            continue
        groups.setdefault((r.tname, r.qname, r.strand), []).append(r)
    for k in groups:
        groups[k].sort(key=lambda x: (x.tstart, x.tend, x.qstart, x.qend))
    return groups


# ----------------------------- UniAligner IO -----------------------------

_cigar_re = re.compile(r"(\d+)([M=XID])")  # accept M,=,X,I,D

def parse_cigar(cg: str) -> Optional[List[Tuple[int, str]]]:
    if not cg:
        return None
    ops: List[Tuple[int, str]] = []
    pos = 0
    for m in _cigar_re.finditer(cg):
        if m.start() != pos:
            return None
        n = int(m.group(1))
        op = m.group(2)
        if n <= 0:
            return None
        ops.append((n, op))
        pos = m.end()
    if pos != len(cg):
        return None
    return ops

def cigar_to_string(ops: List[Tuple[int, str]]) -> str:
    # SyRI requires eqx-like cigar in cg:Z: use =/X/I/D (not M)
    out = []
    for n, op in ops:
        if op == "M":
            op = "="
        out.append(f"{n}{op}")
    return "".join(out)

def cigar_counts(ops: List[Tuple[int, str]]) -> Tuple[int, int, int, int]:
    """
    Return (ref_consume, qry_consume, nmatch, alen)
      '=' and 'M' consume both; treat as matches for nmatch
      'X' consumes both; not counted as match
      'D' consumes ref
      'I' consumes qry
    """
    ref_c = qry_c = nmatch = alen = 0
    for n, op in ops:
        if op in ("=", "M"):
            ref_c += n
            qry_c += n
            nmatch += n
            alen += n
        elif op == "X":
            ref_c += n
            qry_c += n
            alen += n
        elif op == "D":
            ref_c += n
            alen += n
        elif op == "I":
            qry_c += n
            alen += n
    return ref_c, qry_c, nmatch, alen

def unialigner_run(unialigner_cmd_tpl: str, ref_fa: str, qry_fa: str, out_cigar_path: str) -> Optional[str]:
    """
    Run UniAligner using a template command.
    Template can contain {REF} {QRY} {OUT}.
    UniAligner must write a single cigar line to OUT.
    """
    cmd = unialigner_cmd_tpl.format(REF=ref_fa, QRY=qry_fa, OUT=out_cigar_path)
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        return None
    if not os.path.exists(out_cigar_path) or os.path.getsize(out_cigar_path) == 0:
        return None
    with open(out_cigar_path, "r") as f:
        for ln in f:
            s = ln.strip()
            if s:
                return s
    return None


# ----------------------------- CIGAR trimming to window -----------------------------

def trim_ops_left(ops: List[Tuple[int, str]], need_ref: int, need_qry: int) -> Tuple[List[Tuple[int, str]], int, int]:
    out = ops[:]
    cons_r = cons_q = 0
    guard = 0
    while out and (need_ref > 0 or need_qry > 0):
        n, op = out[0]
        if op in ("=", "M", "X"):
            k = min(n, max(need_ref, need_qry))
            cons_r += k
            cons_q += k
            need_ref -= k
            need_qry -= k
            if k == n:
                out.pop(0)
            else:
                out[0] = (n - k, op)
        elif op == "D":
            k = min(n, max(need_ref, 0))
            if k == 0 and need_qry > 0:
                k = n
            cons_r += k
            need_ref -= k
            if k == n:
                out.pop(0)
            else:
                out[0] = (n - k, op)
        elif op == "I":
            k = min(n, max(need_qry, 0))
            if k == 0 and need_ref > 0:
                k = n
            cons_q += k
            need_qry -= k
            if k == n:
                out.pop(0)
            else:
                out[0] = (n - k, op)
        else:
            out.pop(0)
        guard += 1
        if guard > 10_000_000:
            break
    return out, cons_r, cons_q

def trim_ops_right(ops: List[Tuple[int, str]], need_ref: int, need_qry: int) -> Tuple[List[Tuple[int, str]], int, int]:
    out = ops[:]
    cons_r = cons_q = 0
    guard = 0
    while out and (need_ref > 0 or need_qry > 0):
        n, op = out[-1]
        if op in ("=", "M", "X"):
            k = min(n, max(need_ref, need_qry))
            cons_r += k
            cons_q += k
            need_ref -= k
            need_qry -= k
            if k == n:
                out.pop()
            else:
                out[-1] = (n - k, op)
        elif op == "D":
            k = min(n, max(need_ref, 0))
            cons_r += k
            need_ref -= k
            if k == n:
                out.pop()
            else:
                out[-1] = (n - k, op)
        elif op == "I":
            k = min(n, max(need_qry, 0))
            cons_q += k
            need_qry -= k
            if k == n:
                out.pop()
            else:
                out[-1] = (n - k, op)
        else:
            out.pop()
        guard += 1
        if guard > 10_000_000:
            break
    return out, cons_r, cons_q

def trim_cigar_to_window(
    ops: List[Tuple[int, str]],
    tstart: int, tend: int,
    qstart: int, qend: int,
    t_lo: int, t_hi: int,
    q_lo: int, q_hi: int,
) -> Optional[Tuple[List[Tuple[int, str]], int, int, int, int]]:
    """
    Given cigar ops and implied local coords [tstart,tend), [qstart,qend),
    trim from both ends so that resulting coords fall within [t_lo,t_hi] and [q_lo,q_hi].
    Returns (ops_trimmed, new_tstart, new_tend, new_qstart, new_qend) or None.
    """
    need_left_t = max(0, t_lo - tstart)
    need_left_q = max(0, q_lo - qstart)
    need_right_t = max(0, tend - t_hi)
    need_right_q = max(0, qend - q_hi)

    ops2, consL_t, consL_q = trim_ops_left(ops, need_left_t, need_left_q)
    ops3, consR_t, consR_q = trim_ops_right(ops2, need_right_t, need_right_q)

    if not ops3:
        return None

    new_tstart = tstart + consL_t
    new_qstart = qstart + consL_q
    new_tend = tend - consR_t
    new_qend = qend - consR_q

    if new_tend <= new_tstart or new_qend <= new_qstart:
        return None
    if not (new_tstart >= t_lo and new_tend <= t_hi and new_qstart >= q_lo and new_qend <= q_hi):
        return None

    return ops3, new_tstart, new_tend, new_qstart, new_qend


# ----------------------------- simplification helpers -----------------------------

def _collapse_eq(ops: List[Tuple[int, str]]) -> List[Tuple[int, str]]:
    """Collapse consecutive '=' into one; convert 'M'->'='."""
    out: List[Tuple[int, str]] = []
    for n, op in ops:
        if op == "M":
            op = "="
        if not out:
            out.append((n, op))
            continue
        pn, pop = out[-1]
        if op == pop and op == "=":
            out[-1] = (pn + n, pop)
        else:
            out.append((n, op))
    return out

def normalize_ops(ops: List[Tuple[int, str]]) -> List[Tuple[int, str]]:
    """
    Convert M->= and merge adjacent identical ops.
    """
    out: List[Tuple[int, str]] = []
    for n, op in ops:
        if n <= 0:
            continue
        if op == "M":
            op = "="
        if out and out[-1][1] == op:
            out[-1] = (out[-1][0] + n, op)
        else:
            out.append((n, op))
    return out

def detect_pure_insertion(ops: List[Tuple[int, str]]) -> Optional[Tuple[int, int, int]]:
    """
    Detect pattern: (=)+ I (=)+ only
    Return (L_eq, I_len, R_eq)
    """
    o = _collapse_eq(ops)
    o = [(n, op) for n, op in o if n > 0]
    if len(o) != 3:
        return None
    if o[0][1] != "=" or o[1][1] != "I" or o[2][1] != "=":
        return None
    return (o[0][0], o[1][0], o[2][0])

def conservative_simplify_ops(
    ops: List[Tuple[int, str]],
    min_net_gap: int = 100,
    max_opposite_gap: int = 50,
    max_total_x: int = 200,
) -> Tuple[List[Tuple[int, str]], str]:
    """
    Very conservative simplifier for gap-level SV representation.

    Returns:
      (simplified_ops, simplify_mode)

    simplify_mode:
      - "raw"
      - "1I"
      - "1D"
    """
    o = normalize_ops(ops)
    if not o:
        return o, "raw"

    # already pure insertion -> leave as-is here; split decision is done later on simplified ops
    ref_cons, qry_cons, _, _ = cigar_counts(o)
    delta = qry_cons - ref_cons  # >0 insertion-like, <0 deletion-like

    if abs(delta) < min_net_gap:
        return o, "raw"

    dom = "I" if delta > 0 else "D"
    opp = "D" if dom == "I" else "I"

    total_opp = sum(n for n, op in o if op == opp)
    total_x = sum(n for n, op in o if op == "X")

    if total_opp > max_opposite_gap:
        return o, "raw"
    if total_x > max_total_x:
        return o, "raw"

    dom_runs = [(idx, n) for idx, (n, op) in enumerate(o) if op == dom]
    if not dom_runs:
        return o, "raw"

    dom_idx, _ = max(dom_runs, key=lambda x: x[1])

    co_total = min(ref_cons, qry_cons)
    left_eq = sum(n for n, op in o[:dom_idx] if op in ("=", "X", "M"))
    left_eq = max(0, min(left_eq, co_total))
    right_eq = co_total - left_eq

    simp: List[Tuple[int, str]] = []
    if left_eq > 0:
        simp.append((left_eq, "="))
    simp.append((abs(delta), dom))
    if right_eq > 0:
        simp.append((right_eq, "="))

    return normalize_ops(simp), f"1{dom}"


# ----------------------------- PAF lifting helpers -----------------------------

def map_rc_interval_to_original(q0: int, seg_len: int, rc_s: int, rc_e: int) -> Tuple[int, int]:
    """
    Query patch was reverse-complemented before alignment.
    Local interval [rc_s, rc_e) in RC coords maps to original:
      [q0 + (seg_len - rc_e), q0 + (seg_len - rc_s))
    """
    return q0 + (seg_len - rc_e), q0 + (seg_len - rc_s)

def emit_paf_line(
    qname: str, qlen: int, qg_s: int, qg_e: int, strand: str,
    tname: str, tlen: int, tg_s: int, tg_e: int,
    nmatch: int, alen: int, mapq: int,
    cg: str,
    extra_tags: List[str],
    raw_cg: Optional[str] = None,
) -> str:
    # PAF optional tags must be 2-letter tags.
    # cg = simplified representation for downstream / SyRI
    # cr = raw trimmed cigar from UniAligner
    tags = [f"cg:Z:{cg}"]
    if raw_cg is not None:
        tags.append(f"cr:Z:{raw_cg}")
    tags.extend(extra_tags)

    return (
        f"{qname}\t{qlen}\t{qg_s}\t{qg_e}\t{strand}\t"
        f"{tname}\t{tlen}\t{tg_s}\t{tg_e}\t"
        f"{nmatch}\t{alen}\t{mapq}\t" +
        "\t".join(tags)
    )


# ----------------------------- dedup key + best-keeper -----------------------------

def _bin(x: int, step: int) -> int:
    if step <= 1:
        return x
    return x // step

def make_gapua_dedup_key(
    tname: str, qname: str, strand: str,
    tg_s: int, tg_e: int,
    qg_s: int, qg_e: int,
    step: int,
) -> Tuple[str, str, str, int, int, int, int]:
    """
    "同じギャップ（近傍）" を粗く同一視するためのキー。
    """
    return (
        tname, qname, strand,
        _bin(tg_s, step), _bin(tg_e, step),
        _bin(qg_s, step), _bin(qg_e, step),
    )


# ----------------------------- main -----------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description="Gap patch local realignment using UniAligner; output SyRI-ready PAF with cg:Z.")
    ap.add_argument("--in-paf", required=True, help="Input PAF blocks (ValleScope/minimap2 etc.)")
    ap.add_argument("--ref-fa", required=True, help="Reference FASTA")
    ap.add_argument("--qry-fa", required=True, help="Query FASTA")
    ap.add_argument("--out-paf", required=True, help="Output PAF = original + gapua additions")

    ap.add_argument("--flank", type=int, default=10000, help="bp margin on both sides of the gap")
    ap.add_argument("--keep-overlap", type=int, default=200, help="bp to keep into anchors when trimming to gap-only")
    ap.add_argument("--max-gapfill", type=int, default=70000, help="only attempt gapfill when max(ref_gap,qry_gap) <= this")

    ap.add_argument("--bridge-tol", type=int, default=500,
                    help="bridge tolerance: require alignment to reach within tol bp of both patch ends on BOTH ref and query")

    ap.add_argument("--unialigner-cmd", required=True,
                    help="Shell command template to run UniAligner. Use {REF} {QRY} {OUT}. Must write a single CIGAR line to OUT.")
    ap.add_argument("--mapq", type=int, default=60, help="MAPQ to write for accepted gapua PAF (UniAligner may not output MAPQ)")

    ap.add_argument("--split-on-pure-insertion", action="store_true",
                    help="If simplified cigar is pure insertion (=I= only) and flanks long enough, emit two '=' lines instead of cg with I.")
    ap.add_argument("--min-split-flank", type=int, default=500,
                    help="Minimum '=' length required on BOTH sides to do split-on-pure-insertion.")
    ap.add_argument("--keep-failed", action="store_true", help="Append placeholder PAF for failed gaps with na:i:<code> (debug)")

    args = ap.parse_args()

    samtools = check_exe("samtools")
    ensure_fai(samtools, args.ref_fa)
    ensure_fai(samtools, args.qry_fa)
    ref_len = load_fai_lengths(args.ref_fa + ".fai")
    qry_len = load_fai_lengths(args.qry_fa + ".fai")

    recs: List[PafRec] = []
    raw_lines: List[str] = []
    with open(args.in_paf, "r") as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"):
                continue
            raw_lines.append(ln.rstrip("\n"))
            r = parse_paf_line(ln)
            if r:
                if r.tname in ref_len:
                    r.tlen = ref_len[r.tname]
                if r.qname in qry_len:
                    r.qlen = qry_len[r.qname]
                recs.append(r)

    if not recs:
        die("No valid PAF records parsed from --in-paf")

    for r in recs:
        if r.tname not in ref_len:
            die(f"tname not in ref.fa.fai: {r.tname}")
        if r.qname not in qry_len:
            die(f"qname not in qry.fa.fai: {r.qname}")
        if r.strand not in {"+", "-"}:
            die(f"Unsupported strand: {r.strand}")

    groups = group_and_sort_blocks(recs)

    flank = int(args.flank)
    keep_ov = int(args.keep_overlap)
    max_gapfill = int(args.max_gapfill)
    bridge_tol = int(args.bridge_tol)
    min_split_flank = int(args.min_split_flank)

    outdir = os.path.dirname(os.path.abspath(args.out_paf))
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    kept_gapua = 0
    kept_split = 0
    failed = 0
    gap_id = 0

    best_gapua: Dict[Tuple[str, str, str, int, int, int, int], Tuple[float, int, List[str]]] = {}
    dedup_step = max(1, keep_ov)

    with open(args.out_paf, "w") as out, tempfile.TemporaryDirectory(prefix="gapua_") as td:
        for ln in raw_lines:
            out.write(ln + "\n")

        ref_patch_fa = os.path.join(td, "REF_PATCH.fa")
        qry_patch_fa = os.path.join(td, "QRY_PATCH.fa")
        cigar_path = os.path.join(td, "ua.cigar.txt")

        for (tname, qname, strand), lst in groups.items():
            for i in range(len(lst) - 1):
                a = lst[i]
                b = lst[i + 1]

                gap_t = max(0, b.tstart - a.tend)

                a_qf_s, a_qf_e = forwardize_query_interval(a.qstart, a.qend, a.qlen, strand)
                b_qf_s, b_qf_e = forwardize_query_interval(b.qstart, b.qend, b.qlen, strand)
                gap_q = max(0, b_qf_s - a_qf_e)

                if max(gap_t, gap_q) > max_gapfill:
                    continue

                gap_id += 1
                gid = f"GAPUA{gap_id:06d}"

                tL = ref_len[tname]
                qL = qry_len[qname]

                t0 = max(0, a.tend - flank)
                t1 = min(tL, b.tstart + flank)
                if t1 <= t0:
                    continue

                qf0 = max(0, a_qf_e - flank)
                qf1 = min(qL, b_qf_s + flank)
                if qf1 <= qf0:
                    continue

                q0, q1 = forward_to_original_interval(qf0, qf1, qL, strand)
                q0 = max(0, min(qL, q0))
                q1 = max(0, min(qL, q1))
                if q1 < q0:
                    q0, q1 = q1, q0
                if q1 <= q0:
                    continue

                ref_seq = faidx_fetch_sequence(samtools, args.ref_fa, tname, t0, t1)
                qry_seq = faidx_fetch_sequence(samtools, args.qry_fa, qname, q0, q1)
                if not ref_seq or not qry_seq:
                    failed += 1
                    if args.keep_failed:
                        out.write(
                            emit_paf_line(
                                qname, qL, q0, q1, strand,
                                tname, tL, t0, t1,
                                0, 0, 0,
                                "0=",
                                ["tp:Z:gapua", f"id:Z:{gid}", "na:i:2", f"gt:i:{gap_t}", f"gq:i:{gap_q}"],
                                raw_cg="0=",
                            ) + "\n"
                        )
                    continue

                qry_seq_for_align = qry_seq
                qry_is_rc = False
                if strand == "-":
                    qry_seq_for_align = revcomp(qry_seq)
                    qry_is_rc = True

                write_fasta(ref_patch_fa, "REF_PATCH", ref_seq)
                write_fasta(qry_patch_fa, "QRY_PATCH", qry_seq_for_align)

                if os.path.exists(cigar_path):
                    try:
                        os.remove(cigar_path)
                    except Exception:
                        pass

                cg = unialigner_run(args.unialigner_cmd, ref_patch_fa, qry_patch_fa, cigar_path)
                if cg is None:
                    failed += 1
                    if args.keep_failed:
                        out.write(
                            emit_paf_line(
                                qname, qL, q0, q1, strand,
                                tname, tL, t0, t1,
                                0, 0, 0,
                                "0=",
                                ["tp:Z:gapua", f"id:Z:{gid}", "na:i:3", f"gt:i:{gap_t}", f"gq:i:{gap_q}"],
                                raw_cg="0=",
                            ) + "\n"
                        )
                    continue

                ops = parse_cigar(cg)
                if ops is None:
                    failed += 1
                    if args.keep_failed:
                        out.write(
                            emit_paf_line(
                                qname, qL, q0, q1, strand,
                                tname, tL, t0, t1,
                                0, 0, 0,
                                "0=",
                                ["tp:Z:gapua", f"id:Z:{gid}", "na:i:7", f"gt:i:{gap_t}", f"gq:i:{gap_q}"],
                                raw_cg=cg,
                            ) + "\n"
                        )
                    continue

                ref_cons, qry_cons, nmatch, alen = cigar_counts(ops)

                tstart = 0
                qstart = 0
                if ops and ops[0][1] == "D":
                    tstart = ops[0][0]
                if ops and ops[0][1] == "I":
                    qstart = ops[0][0]
                tend = tstart + ref_cons
                qend = qstart + qry_cons

                ref_seg_len = len(ref_seq)
                qry_seg_len = len(qry_seq_for_align)

                left_ok_t = (tstart <= bridge_tol)
                right_ok_t = (tend >= ref_seg_len - bridge_tol)
                left_ok_q = (qstart <= bridge_tol)
                right_ok_q = (qend >= qry_seg_len - bridge_tol)
                if not (left_ok_t and right_ok_t and left_ok_q and right_ok_q):
                    failed += 1
                    if args.keep_failed:
                        out.write(
                            emit_paf_line(
                                qname, qL, q0, q1, strand,
                                tname, tL, t0, t1,
                                0, 0, 0,
                                cigar_to_string(ops),
                                ["tp:Z:gapua", f"id:Z:{gid}", "na:i:4", f"gt:i:{gap_t}", f"gq:i:{gap_q}"],
                                raw_cg=cigar_to_string(ops),
                            ) + "\n"
                        )
                    continue

                t_left = a.tend - t0
                t_right = b.tstart - t0
                t_lo = max(0, t_left - keep_ov)
                t_hi = min(ref_seg_len, t_right + keep_ov)

                q_left_f = a_qf_e - qf0
                q_right_f = b_qf_s - qf0
                q_lo = max(0, q_left_f - keep_ov)
                q_hi = min(qry_seg_len, q_right_f + keep_ov)

                trimmed = trim_cigar_to_window(ops, tstart, tend, qstart, qend, t_lo, t_hi, q_lo, q_hi)
                if trimmed is None:
                    failed += 1
                    if args.keep_failed:
                        out.write(
                            emit_paf_line(
                                qname, qL, q0, q1, strand,
                                tname, tL, t0, t1,
                                0, 0, 0,
                                cigar_to_string(ops),
                                ["tp:Z:gapua", f"id:Z:{gid}", "na:i:5", f"gt:i:{gap_t}", f"gq:i:{gap_q}"],
                                raw_cg=cigar_to_string(ops),
                            ) + "\n"
                        )
                    continue

                ops_t, nt0, nt1, nq0, nq1 = trimmed
                cg_raw = cigar_to_string(ops_t)

                tg_s = t0 + nt0
                tg_e = t0 + nt1

                def lift_query_interval(local_s: int, local_e: int) -> Tuple[int, int]:
                    if not qry_is_rc:
                        return q0 + local_s, q0 + local_e
                    seg_len = len(qry_seq_for_align)
                    return map_rc_interval_to_original(q0, seg_len, local_s, local_e)

                out_strand = strand
                extra_base_tags = [
                    "tp:Z:gapua",
                    f"id:Z:{gid}",
                    f"gt:i:{gap_t}",
                    f"gq:i:{gap_q}",
                    f"fl:i:{flank}",
                    f"ko:i:{keep_ov}",
                    f"bt:i:{bridge_tol}",
                ]

                candidate_lines: List[str] = []
                ref_cons2, qry_cons2, nmatch2, alen2 = cigar_counts(ops_t)
                ident = (float(nmatch2) / float(alen2)) if alen2 > 0 else 0.0

                # simplify first
                simp_ops, simp_mode = conservative_simplify_ops(ops_t)
                cg_simp = cigar_to_string(simp_ops)

                # split decision is now based on simplified ops
                split_info = None
                if args.split_on_pure_insertion:
                    pi = detect_pure_insertion(simp_ops)
                    if pi is not None:
                        L_eq, I_len, R_eq = pi
                        if L_eq >= min_split_flank and R_eq >= min_split_flank:
                            split_info = (L_eq, I_len, R_eq)

                if split_info is not None:
                    L_eq, I_len, R_eq = split_info

                    # left/right split on simplified =I=
                    # global intervals are derived from overall trimmed query/ref window
                    # left segment
                    ltg_s = tg_s
                    ltg_e = tg_s + L_eq
                    lqg_s = lift_query_interval(nq0, nq0 + L_eq)[0]
                    lqg_e = lift_query_interval(nq0, nq0 + L_eq)[1]

                    # right segment:
                    # on ref, right starts after left_eq
                    rtg_s = tg_s + L_eq
                    rtg_e = tg_s + L_eq + R_eq

                    # on query, right starts after left_eq + insertion_len
                    rqg_s = lift_query_interval(nq0 + L_eq + I_len, nq0 + L_eq + I_len + R_eq)[0]
                    rqg_e = lift_query_interval(nq0 + L_eq + I_len, nq0 + L_eq + I_len + R_eq)[1]

                    if ltg_e > ltg_s and rtg_e > rtg_s and lqg_e > lqg_s and rqg_e > rqg_s:
                        candidate_lines.append(
                            emit_paf_line(
                                qname, qL, lqg_s, lqg_e, out_strand,
                                tname, tL, ltg_s, ltg_e,
                                L_eq, L_eq, int(args.mapq),
                                f"{L_eq}=",
                                extra_base_tags + ["tr:Z:splitI_left", f"hi:i:{I_len}", "sm:Z:splitI"],
                                raw_cg=cg_raw,
                            )
                        )
                        candidate_lines.append(
                            emit_paf_line(
                                qname, qL, rqg_s, rqg_e, out_strand,
                                tname, tL, rtg_s, rtg_e,
                                R_eq, R_eq, int(args.mapq),
                                f"{R_eq}=",
                                extra_base_tags + ["tr:Z:splitI_right", f"hi:i:{I_len}", "sm:Z:splitI"],
                                raw_cg=cg_raw,
                            )
                        )

                        qg_s_all, qg_e_all = lift_query_interval(nq0, nq1)
                        key = make_gapua_dedup_key(
                            tname, qname, out_strand,
                            tg_s, tg_e,
                            qg_s_all, qg_e_all,
                            dedup_step
                        )

                        prev = best_gapua.get(key)
                        if (prev is None) or (ident > prev[0]) or (ident == prev[0] and alen2 > prev[1]):
                            best_gapua[key] = (ident, alen2, candidate_lines)
                    else:
                        failed += 1
                    continue

                # normal one-line candidate
                qg_s, qg_e = lift_query_interval(nq0, nq1)
                if tg_e <= tg_s or qg_e <= qg_s:
                    failed += 1
                    continue

                _, _, nmatch_s, alen_s = cigar_counts(simp_ops)

                candidate_lines.append(
                    emit_paf_line(
                        qname, qL, qg_s, qg_e, out_strand,
                        tname, tL, tg_s, tg_e,
                        nmatch_s, alen_s, int(args.mapq),
                        cg_simp,
                        extra_base_tags + ["tr:Z:gap_only", f"sm:Z:{simp_mode}"],
                        raw_cg=cg_raw,
                    )
                )

                key = make_gapua_dedup_key(tname, qname, out_strand, tg_s, tg_e, qg_s, qg_e, dedup_step)
                prev = best_gapua.get(key)
                if (prev is None) or (ident > prev[0]) or (ident == prev[0] and alen2 > prev[1]):
                    best_gapua[key] = (ident, alen2, candidate_lines)

        for _, (_, _, lines) in best_gapua.items():
            for ln in lines:
                out.write(ln + "\n")

        for _, (_, _, lines) in best_gapua.items():
            if len(lines) == 1:
                kept_gapua += 1
            else:
                kept_split += len(lines)

    print(f"[done] wrote: {args.out_paf}", file=sys.stderr)
    print(f"[stats] kept_gapua={kept_gapua} kept_split_lines={kept_split} failed={failed} dedup_keys={len(best_gapua)}", file=sys.stderr)
    print("[hint] Run SyRI like: syri -c <out.paf> -r <ref.fa> -q <qry.fa> -F P --cigar --prefix out", file=sys.stderr)
    print("[na:i codes] 2=faidx empty, 3=UniAligner fail/no cigar, 4=bridge fail, 5=trim fail, 7=cigar parse fail", file=sys.stderr)
    print(f"[dedup] key bin step = keep_overlap = {dedup_step}", file=sys.stderr)


if __name__ == "__main__":
    main()