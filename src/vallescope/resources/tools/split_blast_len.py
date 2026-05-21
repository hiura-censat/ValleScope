#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

def split_blast_by_length(infile, ge100_file, lt100_file, threshold=100):
    with open(infile, "r") as fin, \
         open(ge100_file, "w") as fout_ge, \
         open(lt100_file, "w") as fout_lt:

        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                # コメント行は両方に出力
                fout_ge.write(line + "\n")
                fout_lt.write(line + "\n")
                continue

            cols = line.split("\t")
            if len(cols) < 4:
                continue  # フォーマット不正はスキップ

            try:
                length = int(cols[3])
            except ValueError:
                continue

            if length >= threshold:
                fout_ge.write(line + "\n")
            else:
                fout_lt.write(line + "\n")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="BLAST outfmt6 を長さで2分割")
    ap.add_argument("infile", help="入力: BLAST outfmt6 ファイル")
    ap.add_argument("--ge100", default="blast_len_ge100.txt", help="出力: 長さ>=100 のファイル")
    ap.add_argument("--lt100", default="blast_len_lt100.txt", help="出力: 長さ<100 のファイル")
    ap.add_argument("--threshold", type=int, default=100, help="しきい値 (デフォルト100)")
    args = ap.parse_args()

    split_blast_by_length(args.infile, args.ge100, args.lt100, args.threshold)
