#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import os
import re
import pandas as pd
from collections import defaultdict
from typing import Dict, List, Tuple, Iterable, Optional

def parse_args():
    ap = argparse.ArgumentParser(
        description="families.tsv だけから JBrowse2 synteny 用の簡易PAFを自動生成（列名を自動判定）。"
    )
    ap.add_argument("--families", required=True, help="multi.families.tsv など（タブ区切り）")
    ap.add_argument("--half", type=int, default=50, help="アンカー中心 ±half を線分長にする（既定: 50）")
    ap.add_argument("--outdir", default=None, help="出力先ディレクトリ（既定: families.tsv と同じ）")
    ap.add_argument(
        "--pairs",
        default=None,
        help="生成するペアを限定（例: REF1:TGT1,REF1:TGT2）。指定しなければ最左のサンプルを参照として他すべてに対して出力。"
    )
    return ap.parse_args()

# ───────────────────────────────────────────────────────────────
# 1) families.tsv の列からサンプル接尾辞（suffix）を抽出
#    例: chrom_<SUF>, center_<SUF> が共に存在するものをサンプルとみなす
# ───────────────────────────────────────────────────────────────
def detect_samples(columns: Iterable[str]) -> List[str]:
    chrom_sufs = []
    for c in columns:
        m = re.match(r"^chrom_(.+)$", c)
        if m:
            suf = m.group(1)
            # center_ もあるものだけ採用
            if f"center_{suf}" in columns:
                chrom_sufs.append(suf)
    # 列順に近い順序を保つため columns の順で抽出済み
    # 重複は生じないはずだが念のため
    seen = set()
    out = []
    for s in chrom_sufs:
        if s not in seen:
            out.append(s)
            seen.add(s)
    return out

# ───────────────────────────────────────────────────────────────
# 2) 染色体長の推定：各サンプル・染色体ごとに center の最大値 + half
# ───────────────────────────────────────────────────────────────
def infer_chrom_sizes(df: pd.DataFrame, suf: str, half: int) -> Dict[str, int]:
    chrom_col = f"chrom_{suf}"
    ctr_col   = f"center_{suf}"
    sizes: Dict[str, int] = defaultdict(int)
    if chrom_col not in df.columns or ctr_col not in df.columns:
        return {}
    for _, r in df[[chrom_col, ctr_col]].dropna().iterrows():
        c = str(r[chrom_col])
        try:
            p = int(float(r[ctr_col]))
        except Exception:
            continue
        sizes[c] = max(sizes[c], p + half)
    # 0 は避ける
    for k, v in list(sizes.items()):
        sizes[k] = max(1, int(v))
    return dict(sizes)

# ───────────────────────────────────────────────────────────────
# 3) 1 行 → PAF レコード化
# ───────────────────────────────────────────────────────────────
def paf_row(
    rid: str,
    c1: str, p1: int, L1: int,
    c2: str, p2: int, L2: int,
    half: int
) -> Optional[str]:
    qs = max(0, p1 - half)
    qe = min(L1, p1 + half)
    ts = max(0, p2 - half)
    te = min(L2, p2 + half)
    if qe <= qs or te <= ts:
        return None
    nmatch = min(qe - qs, te - ts)
    alen   = nmatch
    mapq   = 60
    fields = [
        c1, L1, qs, qe, "+",
        c2, L2, ts, te,
        nmatch, alen, mapq,
        f"id:Z:{rid}"
    ]
    return "\t".join(map(str, fields))

# ───────────────────────────────────────────────────────────────
# 4) families → PAF 変換（1ペア分）
# ───────────────────────────────────────────────────────────────
def families_to_paf(
    df: pd.DataFrame,
    ref_suf: str,
    tgt_suf: str,
    half: int,
    out_path: str
) -> int:
    ref_chr = f"chrom_{ref_suf}"
    ref_ctr = f"center_{ref_suf}"
    tgt_chr = f"chrom_{tgt_suf}"
    tgt_ctr = f"center_{tgt_suf}"
    for need in (ref_chr, ref_ctr, tgt_chr, tgt_ctr):
        if need not in df.columns:
            raise SystemExit(f"[ERR] missing column: {need}")

    ref_sizes = infer_chrom_sizes(df, ref_suf, half)
    tgt_sizes = infer_chrom_sizes(df, tgt_suf, half)
    if not ref_sizes or not tgt_sizes:
        print(f"[WARN] sizes could not be fully inferred for {ref_suf} or {tgt_suf}; "
              f"some rows may be skipped.", file=sys.stderr)

    n = 0
    with open(out_path, "w") as out:
        for _, r in df.iterrows():
            c1, c2 = r.get(ref_chr), r.get(tgt_chr)
            if pd.isna(c1) or pd.isna(c2):
                continue
            c1 = str(c1)
            c2 = str(c2)
            if c1 not in ref_sizes or c2 not in tgt_sizes:
                continue
            try:
                p1 = int(float(r.get(ref_ctr)))
                p2 = int(float(r.get(tgt_ctr)))
            except Exception:
                continue
            rid = str(r.get("ref_anchor_id", "."))
            line = paf_row(rid, c1, p1, ref_sizes[c1], c2, p2, tgt_sizes[c2], half)
            if line:
                out.write(line + "\n")
                n += 1
    return n

# ───────────────────────────────────────────────────────────────
# 5) メイン
# ───────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    df = pd.read_csv(args.families, sep="\t")

    # サンプル候補を列から自動検出
    samples = detect_samples(df.columns)
    if not samples:
        raise SystemExit("[ERR] chrom_* / center_* 形式の列が見つかりませんでした。")

    # 出力ディレクトリ
    outdir = args.outdir or os.path.dirname(os.path.abspath(args.families)) or "."
    os.makedirs(outdir, exist_ok=True)

    # 生成するペア
    pairs: List[Tuple[str, str]] = []
    if args.pairs:
        for token in args.pairs.split(","):
            token = token.strip()
            if not token or ":" not in token: 
                continue
            ref, tgt = token.split(":", 1)
            if ref not in samples or tgt not in samples:
                raise SystemExit(f"[ERR] --pairs の指定 '{ref}:{tgt}' のいずれかが列から検出できません。検出={samples}")
            pairs.append((ref, tgt))
    else:
        # 既定：最左（列順で最初に見つかった）を参照として、それ以外すべてをターゲットに
        ref = samples[0]
        for tgt in samples[1:]:
            pairs.append((ref, tgt))

    # 変換実行
    base = os.path.splitext(os.path.basename(args.families))[0]
    total_links = 0
    outputs = []
    for ref, tgt in pairs:
        safe_ref = re.sub(r"[^A-Za-z0-9._-]+", "_", ref)
        safe_tgt = re.sub(r"[^A-Za-z0-9._-]+", "_", tgt)
        out_paf = os.path.join(outdir, f"{base}__{safe_ref}_vs_{safe_tgt}.paf")
        n = families_to_paf(df, ref, tgt, args.half, out_paf)
        print(f"[OK] wrote {n} links -> {out_paf}")
        total_links += n
        outputs.append(out_paf)

    if total_links == 0:
        print("[WARN] 出力0行でした。列名（chrom_/center_）や座標値を確認してください。", file=sys.stderr)

if __name__ == "__main__":
    main()
