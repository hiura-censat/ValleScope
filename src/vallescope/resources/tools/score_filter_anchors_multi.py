#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ValleScope: context-based anchor dictionary + query assignment
MUS/MUM + alpha-separation (token-sequence version)  [THEORY-ORIENTED REVISION]

FULL script (with inversion search in 1st pass; minimal change)

This "FULL" script includes the fixes discussed previously PLUS:
  - 1st pass tries BOTH orientations per (qid, k):
      * forward: token-MUM -> ED
      * inversion: token-MUM_inv -> ED_inv (ED uses L/R swapped)
  - assignments.tsv includes "assign_strand" (the chosen strand for the assignment: '+' or '-')
  - Downstream logic remains as-is; inversion links can contribute to later bundling
    (final bundle extraction already supports inversion bundles by coordinate transform).

NOTE (minimal change request):
  - PAF output is removed (no .bundles.paf is produced).
  - Everything else is left as-is as much as possible (TSV bundles + plots remain).
"""

import subprocess
import tempfile
import shutil

import argparse
import gzip
import hashlib
import os
from collections import Counter, defaultdict
from typing import Dict, IO, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import edlib
from functools import lru_cache


# ----------------- IO helpers -----------------
def _split_csv_arg(x: str) -> List[str]:
    return [s.strip() for s in x.split(",") if s.strip()]


def _derive_sample_name_from_bed(bed_path: str) -> str:
    stem = os.path.basename(bed_path)
    for suf in [".filtered.anchors.bed", ".bed"]:
        if stem.endswith(suf):
            stem = stem[: -len(suf)]
            break
    return stem


def load_bed(path: str) -> pd.DataFrame:
    cols = ["chrom", "start", "end", "anchor_id", "score", "strand"]
    df = pd.read_csv(path, sep="\t", header=None, names=cols)
    df["start"] = df["start"].astype(float)
    df["end"] = df["end"].astype(float)
    df["center"] = (df["start"] + df["end"]) / 2.0
    df["anchor_id"] = df["anchor_id"].astype(str)
    df["chrom"] = df["chrom"].astype(str)
    df["strand"] = df["strand"].astype(str)
    return df


def load_fasta(path: str) -> Dict[str, str]:
    d: Dict[str, str] = {}
    name: Optional[str] = None
    buf: List[str] = []

    fh: IO[str]
    if path.endswith((".gz", ".bgz")):
        fh = gzip.open(path, "rt")
    else:
        fh = open(path, "rt")

    with fh as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    d[name] = "".join(buf)
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())

    if name is not None:
        d[name] = "".join(buf)

    return {str(k): v for k, v in d.items()}


# ----------------- DNA helpers -----------------
def revcomp(s: str) -> str:
    t = str.maketrans("ACGTacgtnN", "TGCAtgcaNn")
    return s.translate(t)[::-1]


# ----------------- MinHash + LSH grouping -----------------
def _hash64_bytes(x: bytes) -> np.uint64:
    h = hashlib.blake2b(x, digest_size=8).digest()
    return np.frombuffer(h, dtype=np.uint64)[0]


def _kmer_hashes_canonical(seq: str, k: int) -> np.ndarray:
    """
    Return unique canonical k-mer hashes (uint64) for seq.
    canonical: min(kmer, revcomp(kmer)) so strand-invariant.
    """
    seq = seq.upper()
    n = len(seq)
    if n < k:
        return np.array([], dtype=np.uint64)

    hs = set()
    for i in range(n - k + 1):
        kmer = seq[i : i + k]
        if "N" in kmer:
            continue
        rc = revcomp(kmer)
        canon = kmer if kmer <= rc else rc
        hs.add(int(_hash64_bytes(canon.encode("ascii"))))
    if not hs:
        return np.array([], dtype=np.uint64)
    return np.array(list(hs), dtype=np.uint64)


def _make_minhash_params(num_perm: int, seed: int = 1):
    rng = np.random.default_rng(seed)
    P = np.uint64(18446744073709551557)  # near 2^64 prime
    a = rng.integers(1, np.iinfo(np.uint64).max, size=num_perm, dtype=np.uint64)
    b = rng.integers(0, np.iinfo(np.uint64).max, size=num_perm, dtype=np.uint64)
    return a, b, P


def minhash_signature_from_hashes(
    hashes: np.ndarray, a: np.ndarray, b: np.ndarray, P: np.uint64
) -> np.ndarray:
    num_perm = a.shape[0]
    if hashes.size == 0:
        return np.full(num_perm, np.uint64(np.iinfo(np.uint64).max), dtype=np.uint64)

    sig = np.full(num_perm, P - np.uint64(1), dtype=np.uint64)
    for x in hashes:
        vals = (a * x + b) % P
        sig = np.minimum(sig, vals)
    return sig


def lsh_band_hashes(sig: np.ndarray, bands: int, rows_per_band: int) -> List[str]:
    if sig.size != bands * rows_per_band:
        raise ValueError(f"signature length ({sig.size}) != bands*rows ({bands*rows_per_band})")
    out = []
    for bi in range(bands):
        chunk = sig[bi * rows_per_band : (bi + 1) * rows_per_band]
        h = hashlib.blake2b(chunk.tobytes(), digest_size=8).hexdigest()
        out.append(h)
    return out


def group_id_from_band_hashes(band_hashes: List[str], group_bands: int) -> str:
    m = int(group_bands)
    if m <= 0:
        raise ValueError("group_bands must be >=1")
    if m > len(band_hashes):
        m = len(band_hashes)
    return "LSH:" + "-".join(band_hashes[:m])


def build_groups_all_samples(
    fastas: List[Dict[str, str]],
    kmer_k: int,
    num_perm: int,
    bands: int,
    rows: int,
    seed: int,
    group_bands: int,
) -> List[Dict[str, str]]:
    if bands * rows != num_perm:
        raise ValueError("bands * rows must equal num_perm")
    if group_bands <= 0:
        raise ValueError("group_bands must be >= 1")
    if group_bands > bands:
        print(f"[WARN] group_bands ({group_bands}) > bands ({bands}); clamping to {bands}")
        group_bands = bands

    a_mh, b_mh, P_mh = _make_minhash_params(num_perm, seed=seed)

    groups: List[Dict[str, str]] = []
    for fa in fastas:
        gmap: Dict[str, str] = {}
        for aid, seq in fa.items():
            hashes = _kmer_hashes_canonical(seq, kmer_k)
            sig = minhash_signature_from_hashes(hashes, a_mh, b_mh, P_mh)
            bh = lsh_band_hashes(sig, bands=bands, rows_per_band=rows)
            gid = group_id_from_band_hashes(bh, group_bands=group_bands)
            gmap[str(aid)] = gid
        groups.append(gmap)
    return groups


def _write_all_anchors_fasta_for_vsearch(
    out_fa: str,
    fastas: List[Dict[str, str]],
    sample_names: List[str],
) -> List[Tuple[int, str, str]]:
    """
    Write a combined FASTA for vsearch clustering.
    Each record name encodes (sample_index, anchor_id):  S{idx}|{anchor_id}

    Returns:
      recs: list of (sample_index, anchor_id, record_name)
    """
    recs: List[Tuple[int, str, str]] = []

    # Deterministic order: sort by (sample_index, anchor_id)
    for si, fa in enumerate(fastas):
        for aid in sorted(fa.keys(), key=lambda x: str(x)):
            rec_name = f"S{si}|{str(aid)}"
            recs.append((si, str(aid), rec_name))

    with open(out_fa, "w") as f:
        for si, aid, rec_name in recs:
            seq = fastas[si].get(aid, "")
            if not seq:
                continue
            f.write(f">{rec_name}\n{seq}\n")

    return recs


def _parse_uc_file(uc_path: str) -> Dict[str, str]:
    """
    Parse vsearch .uc output, return:
      record_name -> cluster_id
    We accept any record type that contains a query label (col 9).
    """
    m: Dict[str, str] = {}
    with open(uc_path, "r") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            cluno = cols[1]
            qlabel = cols[8]  # query label

            # cluster id
            cid = f"VS:{int(cluno):06d}" if cluno.isdigit() else f"VS:{cluno}"

            if qlabel and qlabel != "*":
                m[qlabel] = cid

    return m


def build_groups_all_samples_vsearch(
    fastas: List[Dict[str, str]],
    sample_names: List[str],
    vsearch_path: str = "vsearch",
    identity: float = 0.97,
    threads: int = 1,
    tmp_dir: Optional[str] = None,
    keep_tmp: bool = False,
) -> List[Dict[str, str]]:
    """
    Replace MinHash+LSH grouping with vsearch clustering.

    Returns:
      groups: List[Dict[anchor_id -> group_id]]
    """
    # Locate vsearch
    if os.path.sep in vsearch_path:
        exe = vsearch_path
    else:
        exe = shutil.which(vsearch_path) or vsearch_path

    # quick check
    try:
        subprocess.run(
            [exe, "--version"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )
    except Exception as e:
        raise RuntimeError(f"vsearch not found or not runnable: {vsearch_path} ({e})")

    # Workdir
    if tmp_dir is None:
        td_obj = tempfile.TemporaryDirectory(prefix="vallescope_vsearch_")
        workdir = td_obj.name
    else:
        os.makedirs(tmp_dir, exist_ok=True)
        td_obj = None
        workdir = tmp_dir

    all_fa = os.path.join(workdir, "all_anchors.fa")
    out_uc = os.path.join(workdir, "clusters.uc")
    out_centroids = os.path.join(workdir, "centroids.fa")

    recs = _write_all_anchors_fasta_for_vsearch(all_fa, fastas, sample_names)

    cmd = [
        exe,
        "--cluster_fast",
        all_fa,
        "--id",
        f"{float(identity):.6f}",
        "--strand",
        "both",
        "--uc",
        out_uc,
        "--centroids",
        out_centroids,
        "--threads",
        str(int(threads)),
    ]

    # run
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        raise RuntimeError(
            "vsearch failed.\n"
            f"cmd: {' '.join(cmd)}\n"
            f"stdout:\n{res.stdout}\n"
            f"stderr:\n{res.stderr}\n"
        )

    label_to_cid = _parse_uc_file(out_uc)

    # Build per-sample mapping
    groups: List[Dict[str, str]] = [dict() for _ in range(len(fastas))]
    n_missing = 0
    for si, aid, rec_name in recs:
        cid = label_to_cid.get(rec_name, "")
        if not cid:
            n_missing += 1
            cid = "VS:UNCLUSTERED"
        groups[si][aid] = cid

    if n_missing > 0:
        print(f"[WARN] vsearch: {n_missing} records missing in UC; assigned VS:UNCLUSTERED")

    if (td_obj is not None) and (not keep_tmp):
        td_obj.cleanup()
    else:
        print(f"[info] vsearch tmpdir kept at: {workdir}")
        print(f"[info] vsearch outputs: {out_uc} {out_centroids}")

    return groups


# ----------------- context signatures -----------------
def _bin_gap(x: float, beta: float) -> int:
    return int(np.round(float(x) / float(beta)))


def build_order_by_chrom(bed_df: pd.DataFrame) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    bed_df = bed_df.sort_values(["chrom", "center"])
    for chrom, sub in bed_df.groupby("chrom"):
        out[str(chrom)] = sub["anchor_id"].astype(str).tolist()
    return out


Sig = Tuple[Tuple[Tuple[str, int], ...], Tuple[Tuple[str, int], ...]]  # (Lvec, Rvec)


def build_signatures(
    order_by_chrom: Dict[str, List[str]],
    centers: Dict[str, float],
    groups: Dict[str, str],
    alpha_L: int,
    beta: float,
) -> Dict[str, Dict[int, Sig]]:
    sig_map: Dict[str, Dict[int, Sig]] = {}

    for chrom, order in order_by_chrom.items():
        c = [float(centers[aid]) for aid in order]
        g = [str(groups.get(aid, "")) for aid in order]
        n = len(order)

        for idx, aid in enumerate(order):
            sig_map.setdefault(aid, {})
            for k in range(1, alpha_L + 1):
                if idx - k < 0 or idx + k >= n:
                    continue
                Lvec = []
                for t in range(1, k + 1):
                    gap = c[idx - (t - 1)] - c[idx - t]
                    Lvec.append((g[idx - t], _bin_gap(gap, beta)))
                Rvec = []
                for t in range(1, k + 1):
                    gap = c[idx + t] - c[idx + (t - 1)]
                    Rvec.append((g[idx + t], _bin_gap(gap, beta)))
                sig_map[aid][k] = (tuple(Lvec), tuple(Rvec))
    return sig_map


# ----------------- token sequence helpers -----------------
def _tok_hash(s: str) -> int:
    return int(_hash64_bytes(s.encode("utf-8")))


def build_token_sequence_for_chrom(
    order: List[str],
    centers: Dict[str, float],
    groups: Dict[str, str],
    beta: float,
) -> Tuple[List[int], Dict[str, int], Dict[int, str]]:
    n = len(order)
    T: List[int] = []
    anchor_to_tok: Dict[str, int] = {}
    tok_to_anchor: Dict[int, str] = {}

    for i, aid in enumerate(order):
        gid = str(groups.get(aid, ""))
        T.append(_tok_hash(f"G:{gid}"))
        tok_idx = len(T) - 1
        anchor_to_tok[aid] = tok_idx
        tok_to_anchor[tok_idx] = aid

        if i < n - 1:
            gap = float(centers[order[i + 1]]) - float(centers[aid])
            b = _bin_gap(gap, beta)
            T.append(_tok_hash(f"D:{b}"))

    return T, anchor_to_tok, tok_to_anchor


def build_token_sequences_all_chrom(
    order_by_chrom: Dict[str, List[str]],
    centers: Dict[str, float],
    groups: Dict[str, str],
    beta: float,
):
    T_by: Dict[str, List[int]] = {}
    a2t_by: Dict[str, Dict[str, int]] = {}
    t2a_by: Dict[str, Dict[int, str]] = {}

    for chrom, order in order_by_chrom.items():
        T, a2t, t2a = build_token_sequence_for_chrom(order, centers, groups, beta)
        T_by[chrom] = T
        a2t_by[chrom] = a2t
        t2a_by[chrom] = t2a

    return T_by, a2t_by, t2a_by


def build_inversion_token_sequences_from_raw(
    T_by_raw: Dict[str, List[int]],
    a2t_by_raw: Dict[str, Dict[str, int]],
) -> Tuple[Dict[str, List[int]], Dict[str, Dict[str, int]]]:
    T_inv_by: Dict[str, List[int]] = {}
    a2t_inv_by: Dict[str, Dict[str, int]] = {}

    for chrom, T in T_by_raw.items():
        L = len(T)
        if L == 0:
            T_inv_by[chrom] = []
            a2t_inv_by[chrom] = {}
            continue
        T_inv = list(reversed(T))
        T_inv_by[chrom] = T_inv

        a2t_raw = a2t_by_raw.get(chrom, {})
        m: Dict[str, int] = {}
        for aid, t0 in a2t_raw.items():
            m[str(aid)] = (L - 1) - int(t0)
        a2t_inv_by[chrom] = m

    return T_inv_by, a2t_inv_by


# ----------------- FAST rolling hash over token arrays -----------------
class RollingHash128:
    __slots__ = ("n", "H1", "H2", "P1", "P2")

    def __init__(self, toks: List[int], B1: int, B2: int, C1: int, C2: int):
        n = len(toks)
        self.n = n
        H1 = np.zeros(n + 1, dtype=np.uint64)
        H2 = np.zeros(n + 1, dtype=np.uint64)
        P1 = np.ones(n + 1, dtype=np.uint64)
        P2 = np.ones(n + 1, dtype=np.uint64)

        b1 = np.uint64(B1)
        b2 = np.uint64(B2)
        c1 = np.uint64(C1)
        c2 = np.uint64(C2)

        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            for i, t in enumerate(toks):
                x = np.uint64(int(t) & ((1 << 64) - 1))
                H1[i + 1] = H1[i] * b1 + (x + c1)
                H2[i + 1] = H2[i] * b2 + (x + c2)
                P1[i + 1] = P1[i] * b1
                P2[i + 1] = P2[i] * b2

        self.H1, self.H2, self.P1, self.P2 = H1, H2, P1, P2

    def interval_hash(self, left: int, right: int) -> Tuple[int, int]:
        """Hash of toks[left:right] inclusive (i.e. toks[left:right+1])."""
        if left < 0:
            left = 0
        if right >= self.n:
            right = self.n - 1
        if left > right:
            return (0, 0)
        length = right - left + 1
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            h1 = self.H1[right + 1] - self.H1[left] * self.P1[length]
            h2 = self.H2[right + 1] - self.H2[left] * self.P2[length]
        return (int(h1), int(h2))


def build_rolling_hashers(
    T_by: Dict[str, List[int]],
    seed: int = 12345,
) -> Dict[str, RollingHash128]:
    rng = np.random.default_rng(seed)
    B1 = int(rng.integers(2**32, 2**48) | 1)
    B2 = int(rng.integers(2**32, 2**48) | 1)
    C1 = int(rng.integers(1, 2**32))
    C2 = int(rng.integers(1, 2**32))
    out = {}
    for chrom, T in T_by.items():
        out[chrom] = RollingHash128(T, B1=B1, B2=B2, C1=C1, C2=C2)
    return out


def precompute_anchor_window_hashes(
    order_by_chrom: Dict[str, List[str]],
    a2t_by: Dict[str, Dict[str, int]],
    hasher_by_chrom: Dict[str, RollingHash128],
    alpha_L: int,
) -> Dict[str, Dict[int, Tuple[int, int]]]:
    """
    win_hash[anchor_id][k] = (h1,h2) for the FULL window [t0-2k, t0+2k] inclusive.
    Only stored when full window fits inside chrom token array.
    """
    out: Dict[str, Dict[int, Tuple[int, int]]] = defaultdict(dict)
    for chrom, order in order_by_chrom.items():
        hasher = hasher_by_chrom.get(chrom, None)
        if hasher is None:
            continue
        L = hasher.n
        a2t = a2t_by.get(chrom, {})
        for aid in order:
            t0 = a2t.get(aid, None)
            if t0 is None:
                continue
            for k in range(1, alpha_L + 1):
                left = t0 - 2 * k
                right = t0 + 2 * k
                if left < 0 or right >= L:
                    continue
                out[aid][k] = hasher.interval_hash(left, right)
    return out


# ----------------- MUS helpers (FAST via rolling hashes) -----------------
def compute_unique_substrings_hash_counts_fast(
    hasher: RollingHash128,
    Lmax: int,
) -> Dict[int, Counter]:
    counts_by_L: Dict[int, Counter] = {}
    n = hasher.n
    for L in range(1, Lmax + 1):
        cnt = Counter()
        if n < L:
            counts_by_L[L] = cnt
            continue
        for i in range(0, n - L + 1):
            h = hasher.interval_hash(i, i + L - 1)
            cnt[h] += 1
        counts_by_L[L] = cnt
    return counts_by_L


def compute_mus_intervals_in_window_fast(
    hasher: RollingHash128,
    counts_by_L: Dict[int, Counter],
    win_start: int,
    win_end: int,
    Lmax: int,
) -> List[Tuple[int, int]]:
    mus: List[Tuple[int, int]] = []
    n = hasher.n
    ws = max(0, win_start)
    we = min(n - 1, win_end)

    for s in range(ws, we + 1):
        best = None
        for L in range(1, Lmax + 1):
            e = s + L - 1
            if e > we:
                break
            h = hasher.interval_hash(s, e)
            if counts_by_L[L].get(h, 0) == 1:
                best = (s, e)
                break
        if best is not None:
            mus.append(best)
    return mus


def max_nonoverlap_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: (x[1], x[0]))
    picked: List[Tuple[int, int]] = []
    cur_end = -1
    for s, e in intervals:
        if s > cur_end:
            picked.append((s, e))
            cur_end = e
    return picked


# ----------------- THEORY: compute ref alpha' ONLY until reaching alpha_target -----------------
def compute_ref_alpha_prime_min_k(
    ref_order_by_chrom: Dict[str, List[str]],
    ref_groups: Dict[str, str],
    ref_hasher_by: Dict[str, RollingHash128],
    ref_a2t_by: Dict[str, Dict[str, int]],
    alpha_L: int,
    mus_Lmax_tokens: int,
    alpha_target: int,
) -> Tuple[Dict[str, int], Dict[str, int], Dict[str, str]]:
    k_min_ok: Dict[str, int] = defaultdict(int)
    alpha_prime_at_min: Dict[str, int] = defaultdict(int)
    evidence_at_min: Dict[str, str] = defaultdict(str)

    counts_by_chrom: Dict[str, Dict[int, Counter]] = {}
    for chrom, hasher in ref_hasher_by.items():
        counts_by_chrom[chrom] = compute_unique_substrings_hash_counts_fast(hasher, Lmax=mus_Lmax_tokens)

    for chrom, order in ref_order_by_chrom.items():
        hasher = ref_hasher_by.get(chrom, None)
        if hasher is None or hasher.n == 0:
            continue
        cnt_by_L = counts_by_chrom[chrom]
        a2t = ref_a2t_by.get(chrom, {})

        for rid in order:
            gid = ref_groups.get(rid, "")
            if not gid:
                continue

            t0 = a2t.get(rid, None)
            if t0 is None:
                continue

            found_k = 0
            found_a = 0
            found_ev = ""

            for k in range(1, alpha_L + 1):
                if (t0 - 2 * k) < 0 or (t0 + 2 * k) >= hasher.n:
                    continue

                ws, we = (t0 - 2 * k, t0 + 2 * k)
                mus = compute_mus_intervals_in_window_fast(
                    hasher, cnt_by_L, ws, we, Lmax=mus_Lmax_tokens
                )
                picked = max_nonoverlap_intervals(mus)
                a_prime = len(picked)

                if a_prime >= int(alpha_target):
                    found_k = k
                    found_a = a_prime
                    show = picked[: min(3, len(picked))]
                    found_ev = f"chrom={chrom} win=[{ws},{we}] mus_picked={show} (n={a_prime})"
                    break

            if found_k > 0:
                k_min_ok[rid] = int(found_k)
                alpha_prime_at_min[rid] = int(found_a)
                evidence_at_min[rid] = str(found_ev)
            else:
                k_min_ok[rid] = 0
                alpha_prime_at_min[rid] = 0
                evidence_at_min[rid] = ""

    return k_min_ok, alpha_prime_at_min, evidence_at_min


def assign_unique_ids_reference_mus(
    ref_order_by_chrom: Dict[str, List[str]],
    ref_groups: Dict[str, str],
    ref_k_min_ok: Dict[str, int],
    ref_alpha_prime_at_min: Dict[str, int],
    ref_evidence_at_min: Dict[str, str],
) -> Tuple[Dict[str, str], Dict[str, bool], Dict[str, int], Dict[str, int], Dict[str, str]]:
    uid_map: Dict[str, str] = {}
    amb_map: Dict[str, bool] = {}
    used_k: Dict[str, int] = {}
    alpha_prime_used: Dict[str, int] = {}
    evidence_used: Dict[str, str] = {}

    uid_counter = 1

    for chrom, order in ref_order_by_chrom.items():
        for rid in order:
            gid = ref_groups.get(rid, "")
            if not gid:
                uid_map[rid] = ""
                amb_map[rid] = True
                used_k[rid] = 0
                alpha_prime_used[rid] = 0
                evidence_used[rid] = ""
                continue

            k0 = int(ref_k_min_ok.get(rid, 0))
            if k0 > 0:
                uid = f"anchor{uid_counter:05d}"
                uid_counter += 1
                uid_map[rid] = uid
                amb_map[rid] = False
                used_k[rid] = k0
                alpha_prime_used[rid] = int(ref_alpha_prime_at_min.get(rid, 0))
                evidence_used[rid] = str(ref_evidence_at_min.get(rid, ""))
            else:
                uid_map[rid] = ""
                amb_map[rid] = True
                used_k[rid] = 0
                alpha_prime_used[rid] = 0
                evidence_used[rid] = ""

    return uid_map, amb_map, used_k, alpha_prime_used, evidence_used


# ----------------- bounded edit distance (edlib) -----------------
_EDLIB_CACHE_MAXSIZE = 500_000


@lru_cache(maxsize=_EDLIB_CACHE_MAXSIZE)
def edit_distance_ops_bounded(a, b, bound: int) -> int:
    bound = int(bound)
    INF = bound + 1

    n, m = len(a), len(b)
    if bound < 0:
        return INF
    if abs(n - m) > bound:
        return INF
    if n == 0:
        return 0 if m <= bound else INF
    if m == 0:
        return 0 if n <= bound else INF

    res = edlib.align(a, b, task="distance", k=bound)
    d = int(res["editDistance"])
    return d if d >= 0 else INF


# ----------------- Reference exact index (token-MUM) [UNIQUE KEYS ONLY] -----------------
def build_ref_exact_index_unique(
    ref_order_by_chrom: Dict[str, List[str]],
    ref_groups: Dict[str, str],
    ref_uid: Dict[str, str],
    ref_amb: Dict[str, bool],
    ref_win_hash: Dict[str, Dict[int, Tuple[int, int]]],  # rid->k->hash
    alpha_L: int,
) -> Dict[Tuple[str, int, Tuple[int, int]], str]:
    tmp: Dict[Tuple[str, int, Tuple[int, int]], List[str]] = defaultdict(list)

    for chrom, order in ref_order_by_chrom.items():
        for rid in order:
            if ref_amb.get(rid, True):
                continue
            if not ref_uid.get(rid, ""):
                continue
            gid = ref_groups.get(rid, "")
            if not gid:
                continue
            kh = ref_win_hash.get(rid, {})
            for k in range(1, alpha_L + 1):
                h = kh.get(k, None)
                if h is None:
                    continue
                tmp[(gid, k, h)].append(rid)

    idx: Dict[Tuple[str, int, Tuple[int, int]], str] = {}
    for key, rids in tmp.items():
        if len(rids) == 1:
            idx[key] = rids[0]
    return idx


# ----------------- Conflict-aware query assignment (α' gates candidates) -----------------
def match_query_to_reference_conflict_resolve(
    q_order_by_chrom: Dict[str, List[str]],
    q_groups: Dict[str, str],
    q_sig: Dict[str, Dict[int, Sig]],
    q_win_hash_by_chrom_fwd: Dict[str, Dict[str, Dict[int, Tuple[int, int]]]],  # chrom->qid->k->hash
    q_win_hash_by_chrom_inv: Dict[str, Dict[str, Dict[int, Tuple[int, int]]]],  # chrom->qid->k->hash
    qid_to_chrom: Dict[str, str],
    ref_exact_idx: Dict[Tuple[str, int, Tuple[int, int]], str],
    ref_sig: Dict[str, Dict[int, Sig]],
    ref_groups: Dict[str, str],
    ref_uid: Dict[str, str],
    ref_amb: Dict[str, bool],
    ref_k_min_ok: Dict[str, int],
    alpha_L: int,
    gap_beta: int,
) -> Dict[str, Dict[str, object]]:
    ref_by_group: Dict[str, List[str]] = {}
    for rid, gid in ref_groups.items():
        if ref_amb.get(rid, True):
            continue
        if ref_uid.get(rid, "") == "":
            continue
        ref_by_group.setdefault(gid, []).append(rid)

    result: Dict[str, Dict[str, object]] = {}
    unresolved = set()

    for chrom, order in q_order_by_chrom.items():
        for qid in order:
            unresolved.add(qid)
            result[qid] = {"status": "", "rid": "", "uid": "", "used_k": 0, "n_cand": 0, "method": "", "strand": ""}

    def _alpha_ok(rid: str, k: int) -> bool:
        k0 = int(ref_k_min_ok.get(rid, 0))
        return (k0 > 0) and (int(k) >= k0)

    def _candidates_for_qid(qid: str, k: int, inv: bool) -> Tuple[List[str], str, str]:
        gid = q_groups.get(qid, "")
        qchrom = qid_to_chrom.get(qid, "")
        if not gid or not qchrom:
            return [], "", "+"

        strand = "-" if inv else "+"
        method_suffix = "_inv" if inv else ""

        # token-MUM
        if inv:
            h = q_win_hash_by_chrom_inv.get(qchrom, {}).get(qid, {}).get(k, None)
        else:
            h = q_win_hash_by_chrom_fwd.get(qchrom, {}).get(qid, {}).get(k, None)

        if h is not None:
            rid = ref_exact_idx.get((gid, k, h), "")
            if rid and _alpha_ok(rid, k):
                return [rid], f"token_mum{method_suffix}", strand

        # ED fallback
        if k not in q_sig.get(qid, {}):
            return [], f"ed{method_suffix}", strand

        Lq, Rq = q_sig[qid][k]
        if inv:
            Lq, Rq = Rq, Lq  # swap for inversion

        cand = []
        for rid in ref_by_group.get(gid, []):
            if not _alpha_ok(rid, k):
                continue
            if k not in ref_sig.get(rid, {}):
                continue
            Lr, Rr = ref_sig[rid][k]
            dL = edit_distance_ops_bounded(Lq, Lr, bound=gap_beta)
            if dL > gap_beta:
                continue
            dR = edit_distance_ops_bounded(Rq, Rr, bound=gap_beta - dL)
            if dL + dR <= int(gap_beta):
                cand.append(rid)
        return cand, f"ed{method_suffix}", strand

    def _is_token_mum(method: str) -> bool:
        return method.startswith("token_mum")

    for k in range(1, alpha_L + 1):
        if not unresolved:
            break

        tentative: Dict[str, Dict[str, object]] = {}
        rid_to_qids: Dict[str, List[str]] = defaultdict(list)
        next_unresolved = set()

        for qid in unresolved:
            gid = q_groups.get(qid, "")
            qchrom = qid_to_chrom.get(qid, "")
            if not gid or not qchrom:
                tentative[qid] = {
                    "status": "unmatched",
                    "rid": "",
                    "uid": "",
                    "used_k": k,
                    "n_cand": 0,
                    "method": "",
                    "strand": "",
                }
                next_unresolved.add(qid)
                continue

            c_fwd, m_fwd, s_fwd = _candidates_for_qid(qid, k, inv=False)
            c_inv, m_inv, s_inv = _candidates_for_qid(qid, k, inv=True)

            chosen_cands: List[str] = []
            chosen_method = ""
            chosen_strand = ""

            if len(c_fwd) == 1 and _is_token_mum(m_fwd) and not (
                len(c_inv) == 1 and _is_token_mum(m_inv) and c_inv[0] != c_fwd[0]
            ):
                chosen_cands, chosen_method, chosen_strand = c_fwd, m_fwd, s_fwd
            elif len(c_inv) == 1 and _is_token_mum(m_inv) and not (
                len(c_fwd) == 1 and _is_token_mum(m_fwd) and c_fwd[0] != c_inv[0]
            ):
                chosen_cands, chosen_method, chosen_strand = c_inv, m_inv, s_inv
            else:
                union = list(dict.fromkeys([*c_fwd, *c_inv]))
                chosen_cands = union
                if _is_token_mum(m_fwd) or _is_token_mum(m_inv):
                    chosen_method = "token_mum_mix"
                else:
                    chosen_method = "ed_mix"
                chosen_strand = "+" if (len(c_fwd) >= len(c_inv)) else "-"

            if len(chosen_cands) == 0:
                tentative[qid] = {
                    "status": "unmatched",
                    "rid": "",
                    "uid": "",
                    "used_k": k,
                    "n_cand": 0,
                    "method": chosen_method,
                    "strand": chosen_strand,
                }
                next_unresolved.add(qid)
            elif len(chosen_cands) == 1:
                rid = chosen_cands[0]
                tentative[qid] = {
                    "status": "assign",
                    "rid": rid,
                    "uid": ref_uid.get(rid, ""),
                    "used_k": k,
                    "n_cand": 1,
                    "method": chosen_method,
                    "strand": chosen_strand,
                }
                rid_to_qids[rid].append(qid)
            else:
                tentative[qid] = {
                    "status": "discordant",
                    "rid": "",
                    "uid": "",
                    "used_k": k,
                    "n_cand": len(chosen_cands),
                    "method": chosen_method,
                    "strand": chosen_strand,
                }
                next_unresolved.add(qid)

        conflict_qids = set()
        for rid, qids in rid_to_qids.items():
            if len(qids) > 1:
                conflict_qids.update(qids)

        for qid, info in tentative.items():
            if info["status"] == "assign" and qid in conflict_qids:
                result[qid] = {
                    "status": "discordant",
                    "rid": "",
                    "uid": "",
                    "used_k": k,
                    "n_cand": int(info.get("n_cand", 0)),
                    "method": info.get("method", ""),
                    "strand": info.get("strand", ""),
                }
                next_unresolved.add(qid)
            else:
                result[qid] = info

        unresolved = next_unresolved

    for qid in unresolved:
        if result[qid].get("status", "") in ("discordant",):
            continue
        result[qid] = {"status": "unmatched", "rid": "", "uid": "", "used_k": 0, "n_cand": 0, "method": "", "strand": ""}

    return result


# ----------------- bundle extraction + outputs -----------------
def infer_chrom_lengths_from_bed(bed_df: pd.DataFrame) -> Dict[str, int]:
    d: Dict[str, int] = {}
    for chrom, sub in bed_df.groupby("chrom"):
        mx = float(sub["end"].max())
        d[str(chrom)] = int(np.ceil(mx))
    return d


def extract_parallel_bundles(
    items: List[Tuple[float, float, str, str]],  # (ref_center, qry_center, rid, qid)
    max_gap_ref: float,
    max_gap_qry: float,
    slope_min: float,
    slope_max: float,
    min_links: int,
    min_span: float,
    max_skip: int = 5,
    scale_gap_with_skip: bool = True,
) -> List[Dict[str, object]]:
    if not items:
        return []

    items = sorted(items, key=lambda x: x[0])
    bundles: List[List[Tuple[float, float, str, str]]] = []
    cur: List[Tuple[float, float, str, str]] = [items[0]]
    skip_count = 0

    for (r, q, rid, qid) in items[1:]:
        r0, q0, _, _ = cur[-1]
        dr = float(r - r0)
        dq = float(q - q0)

        gap_mult = (skip_count + 1) if scale_gap_with_skip else 1.0

        ok = True
        if dr <= 0:
            ok = False
        else:
            slope = dq / dr
            if dr > max_gap_ref * gap_mult:
                ok = False
            if abs(dq) > max_gap_qry * gap_mult:
                ok = False
            if not (slope_min <= slope <= slope_max):
                ok = False

        if ok:
            cur.append((r, q, rid, qid))
            skip_count = 0
        else:
            if skip_count < max_skip:
                skip_count += 1
                continue
            bundles.append(cur)
            cur = [(r, q, rid, qid)]
            skip_count = 0

    bundles.append(cur)

    out: List[Dict[str, object]] = []
    for b in bundles:
        if len(b) < min_links:
            continue
        ref_start = min(x[0] for x in b)
        ref_end = max(x[0] for x in b)
        qry_start = min(x[1] for x in b)
        qry_end = max(x[1] for x in b)
        if (ref_end - ref_start) < float(min_span):
            continue

        slopes = []
        for i in range(1, len(b)):
            dr = b[i][0] - b[i - 1][0]
            dq = b[i][1] - b[i - 1][1]
            if dr > 0:
                slopes.append(dq / dr)
        slope_med = float(np.median(slopes)) if slopes else float("nan")

        a = slope_med if np.isfinite(slope_med) else 1.0
        rs = np.array([x[0] for x in b], dtype=float)
        qs = np.array([x[1] for x in b], dtype=float)
        b0 = float(np.median(qs - a * rs))
        resid = qs - (a * rs + b0)
        resid_mad = float(np.median(np.abs(resid - np.median(resid))))

        max_dr = float(np.max(np.diff(rs))) if len(rs) >= 2 else 0.0
        max_dq = float(np.max(np.abs(np.diff(qs)))) if len(qs) >= 2 else 0.0

        out.append(
            {
                "ref_start": float(ref_start),
                "ref_end": float(ref_end),
                "qry_start": float(qry_start),
                "qry_end": float(qry_end),
                "n_links": int(len(b)),
                "slope_med": slope_med,
                "b0_med": b0,
                "resid_mad": resid_mad,
                "max_gap_ref": max_dr,
                "max_gap_qry": max_dq,
                "rids": [x[2] for x in b],
                "qids": [x[3] for x in b],
            }
        )

    return out


# --- inversion-mode conversion helpers ---
def _make_items_inversion_mode(
    items_fwd: List[Tuple[float, float, str, str]],
    qry_chr_len: int,
) -> List[Tuple[float, float, str, str]]:
    qlen = float(qry_chr_len)
    out = []
    for rc, qc, rid, qid in items_fwd:
        out.append((float(rc), float(qlen - float(qc)), str(rid), str(qid)))
    out.sort(key=lambda x: x[0])
    return out


def _convert_bundle_qry_interval_inv_to_raw(
    b: Dict[str, object],
    qry_chr_len: int,
) -> Tuple[float, float]:
    qlen = float(qry_chr_len)
    inv_qs = float(b["qry_start"])
    inv_qe = float(b["qry_end"])
    raw_qs = qlen - inv_qe
    raw_qe = qlen - inv_qs
    if raw_qs > raw_qe:
        raw_qs, raw_qe = raw_qe, raw_qs
    return raw_qs, raw_qe


def _convert_raw_interval_to_inv(raw_qs: float, raw_qe: float, qry_chr_len: int) -> Tuple[float, float]:
    qlen = float(qry_chr_len)
    inv_qs = qlen - float(raw_qe)
    inv_qe = qlen - float(raw_qs)
    if inv_qs > inv_qe:
        inv_qs, inv_qe = inv_qe, inv_qs
    return inv_qs, inv_qe


def write_bundles_tsv(
    out_tsv: str,
    ref_name: str,
    qry_name: str,
    chrom: str,
    bundles: List[Dict[str, object]],
    bundle_id_start: int,
) -> int:
    rows = []
    bid = bundle_id_start
    for b in bundles:
        bundle_id = f"B{bid:06d}"
        bid += 1
        strand = str(b.get("strand", "+"))
        rows.append(
            {
                "bundle_id": bundle_id,
                "ref_name": ref_name,
                "qry_name": qry_name,
                "chrom": chrom,
                "ref_start": int(round(b["ref_start"])),
                "ref_end": int(round(b["ref_end"])),
                "qry_start": int(round(b["qry_start"])),
                "qry_end": int(round(b["qry_end"])),
                "strand": strand,
                "n_links": int(b["n_links"]),
                "ref_span": int(round(b["ref_end"] - b["ref_start"])),
                "qry_span": int(round(b["qry_end"] - b["qry_start"])),
                "slope_med": float(b["slope_med"]),
                "resid_mad": float(b["resid_mad"]),
                "max_gap_ref": float(b["max_gap_ref"]),
                "max_gap_qry": float(b["max_gap_qry"]),
                "ref_anchor_ids": ",".join(map(str, b["rids"])),
                "qry_anchor_ids": ",".join(map(str, b["qids"])),
            }
        )

    if rows:
        df = pd.DataFrame(rows)
        if os.path.exists(out_tsv):
            df.to_csv(out_tsv, sep="\t", index=False, mode="a", header=False)
        else:
            df.to_csv(out_tsv, sep="\t", index=False)
    return bid


def _pair_chrom_key(
    rid: str,
    qid: str,
    ref_chroms: Dict[str, str],
    q_chroms: Dict[str, str],
) -> str:
    rchr = str(ref_chroms.get(rid, "") or "")
    qchr = str(q_chroms.get(qid, "") or "")
    return rchr if rchr else qchr


# ----------------- plotting -----------------
def plot_links_by_chrom(
    ref_name: str,
    query_name: str,
    ref_centers: Dict[str, float],
    ref_chroms: Dict[str, str],
    q_centers: Dict[str, float],
    q_chroms: Dict[str, str],
    pairs: List[Tuple[str, str]],
    out_prefix: str,
    max_lines: int = 20000,
):
    by_chrom: Dict[str, List[Tuple[float, float, str, str]]] = {}
    for rid, qid in pairs:
        chrom = _pair_chrom_key(rid, qid, ref_chroms, q_chroms)
        if not chrom:
            continue
        rc = ref_centers.get(rid, np.nan)
        qc = q_centers.get(qid, np.nan)
        if not np.isfinite(rc) or not np.isfinite(qc):
            continue
        by_chrom.setdefault(chrom, []).append((float(rc), float(qc), rid, qid))

    if not by_chrom:
        print(f"[plot] {query_name}: no assign pairs to plot")
        return

    for chrom, items in by_chrom.items():
        items.sort(key=lambda x: x[0])

        if len(items) > max_lines:
            step = int(np.ceil(len(items) / max_lines))
            items = items[::step]
            print(f"[plot] {query_name} {chrom}: downsampled to {len(items)} lines (step={step})")

        rxs = [x[0] for x in items]
        qxs = [x[1] for x in items]
        xmin = min(min(rxs), min(qxs))
        xmax = max(max(rxs), max(qxs))
        pad = (xmax - xmin) * 0.02 if xmax > xmin else 1.0

        fig = plt.figure(figsize=(14, 4), dpi=200)
        ax = fig.add_subplot(111)

        ax.hlines(0.0, xmin - pad, xmax + pad, linewidth=1.5)
        ax.hlines(1.0, xmin - pad, xmax + pad, linewidth=1.5)

        ax.text(xmin - pad, 0.0, f"{ref_name}", va="bottom", ha="left")
        ax.text(xmin - pad, 1.0, f"{query_name}", va="bottom", ha="left")

        for rc, qc, rid, qid in items:
            ax.plot([rc, qc], [0.0, 1.0], linewidth=0.5, alpha=0.5)

        ax.set_title(f"Anchor links: {ref_name} vs {query_name} ({chrom})  n={len(items)}")
        ax.set_xlim(xmin - pad, xmax + pad)
        ax.set_ylim(-0.2, 1.2)
        ax.set_yticks([0.0, 1.0])
        ax.set_yticklabels(["ref", "query"])
        ax.set_xlabel("center position (bp)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        out_png = f"{out_prefix}.links.png"
        fig.tight_layout()
        fig.savefig(out_png)
        plt.close(fig)
        print(f"[plot] wrote {out_png}")


def plot_bundle_pngs_by_chrom(
    ref_name: str,
    query_name: str,
    chrom: str,
    items: List[Tuple[float, float, str, str]],
    bundles: List[Dict[str, object]],
    out_prefix: str,
    max_lines_bg: int = 20000,
):
    if not bundles:
        return
    if not items:
        print(f"[plot] {query_name} {chrom}: no items to plot")
        return

    bg = list(items)
    if len(bg) > max_lines_bg:
        step = int(np.ceil(len(bg) / max_lines_bg))
        bg = bg[::step]
        print(f"[plot] {query_name} {chrom}: downsampled links to {len(bg)} (step={step})")

    rxs = [x[0] for x in bg]
    qxs = [x[1] for x in bg]
    xmin = min(min(rxs), min(qxs))
    xmax = max(max(rxs), max(qxs))
    pad = (xmax - xmin) * 0.02 if xmax > xmin else 1.0

    fig = plt.figure(figsize=(14, 4), dpi=200)
    ax = fig.add_subplot(111)

    ax.hlines(0.0, xmin - pad, xmax + pad, linewidth=1.5)
    ax.hlines(1.0, xmin - pad, xmax + pad, linewidth=1.5)
    ax.text(xmin - pad, 0.0, f"{ref_name}", va="bottom", ha="left")
    ax.text(xmin - pad, 1.0, f"{query_name}", va="bottom", ha="left")

    band_half_height = 0.07
    legend_added = {"fwd": False, "inv": False}

    for b in bundles:
        rs = float(b["ref_start"])
        re = float(b["ref_end"])
        qs = float(b["qry_start"])
        qe = float(b["qry_end"])
        strand = str(b.get("strand", "+"))

        is_inv = (strand == "-")
        face = "orange" if is_inv else "blue"
        alpha = 0.22 if is_inv else 0.18
        hatch = "///" if is_inv else None
        label = "bundle (+)" if (not is_inv) else "bundle (-, inversion)"

        add_label = False
        if (not is_inv) and (not legend_added["fwd"]):
            add_label = True
            legend_added["fwd"] = True
        if is_inv and (not legend_added["inv"]):
            add_label = True
            legend_added["inv"] = True

        ax.fill_between(
            [rs, re],
            [0.0 - band_half_height, 0.0 - band_half_height],
            [0.0 + band_half_height, 0.0 + band_half_height],
            color=face,
            alpha=alpha,
            linewidth=0.6 if hatch else 0.0,
            edgecolor="black" if hatch else None,
            hatch=hatch,
            label=(label if add_label else None),
        )
        ax.fill_between(
            [qs, qe],
            [1.0 - band_half_height, 1.0 - band_half_height],
            [1.0 + band_half_height, 1.0 + band_half_height],
            color=face,
            alpha=alpha,
            linewidth=0.6 if hatch else 0.0,
            edgecolor="black" if hatch else None,
            hatch=hatch,
        )

        mid_r = (rs + re) / 2.0
        ax.text(mid_r, -0.02, strand, ha="center", va="top", fontsize=7, alpha=0.8)

    for r, q, rid, qid in bg:
        ax.plot([r, q], [0.0, 1.0], linewidth=0.5, alpha=0.5)

    ax.set_title(
        f"Anchor links + bundles: {ref_name} vs {query_name} ({chrom})  links={len(bg)} bundles={len(bundles)}"
    )
    ax.set_xlim(xmin - pad, xmax + pad)
    ax.set_ylim(-0.2, 1.2)
    ax.set_yticks([0.0, 1.0])
    ax.set_yticklabels(["ref", "query"])
    ax.set_xlabel("center position (bp)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if legend_added["fwd"] or legend_added["inv"]:
        ax.legend(loc="upper right", frameon=False, fontsize=8)

    out_png = f"{out_prefix}.bundles.png"
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)
    print(f"[plot] wrote {out_png}")


# ----------------- BETWEEN-BUNDLES refinement -----------------
def _items_from_pairs_same_chrom(
    pairs: List[Tuple[str, str]],
    chrom: str,
    ref_centers: Dict[str, float],
    ref_chroms: Dict[str, str],
    q_centers: Dict[str, float],
    q_chroms: Dict[str, str],
) -> List[Tuple[float, float, str, str]]:
    out: List[Tuple[float, float, str, str]] = []
    for rid, qid in pairs:
        if _pair_chrom_key(rid, qid, ref_chroms, q_chroms) != chrom:
            continue
        rc = ref_centers.get(rid, np.nan)
        qc = q_centers.get(qid, np.nan)
        if not np.isfinite(rc) or not np.isfinite(qc):
            continue
        out.append((float(rc), float(qc), str(rid), str(qid)))
    out.sort(key=lambda x: x[0])
    return out


def refine_between_bundles(
    chrom: str,
    bundles_1st: List[Dict[str, object]],
    qids_candidates: List[str],
    q_groups: Dict[str, str],
    q_sig: Dict[str, Dict[int, Sig]],
    q_win_hash_by_qid: Dict[str, Dict[int, Tuple[int, int]]],  # qid->k->hash
    q_centers: Dict[str, float],
    ref_exact_idx: Dict[Tuple[str, int, Tuple[int, int]], str],
    ref_sig: Dict[str, Dict[int, Sig]],
    ref_groups: Dict[str, str],
    ref_uid: Dict[str, str],
    ref_amb: Dict[str, bool],
    ref_k_min_ok: Dict[str, int],
    ref_centers: Dict[str, float],
    alpha_L: int,
    gap_beta: int,
    used_rids_1st: set,
    inversion_mode: bool = False,
) -> List[Tuple[str, str]]:
    if not bundles_1st:
        return []
    bundles_sorted = sorted(bundles_1st, key=lambda b: float(b["ref_start"]))

    ref_by_group: Dict[str, List[str]] = {}
    for rid, gid in ref_groups.items():
        if ref_amb.get(rid, True):
            continue
        if ref_uid.get(rid, "") == "":
            continue
        if gid:
            ref_by_group.setdefault(gid, []).append(rid)

    def _alpha_ok(rid: str, k: int) -> bool:
        k0 = int(ref_k_min_ok.get(rid, 0))
        return (k0 > 0) and (int(k) >= k0)

    def _rid_in_ref_gap(rid: str, ref_gap: Tuple[float, float]) -> bool:
        c = ref_centers.get(rid, np.nan)
        if not np.isfinite(c):
            return False
        a, b = ref_gap
        return (float(a) < float(c)) and (float(c) < float(b))

    gaps: List[Tuple[Tuple[float, float], Tuple[float, float]]] = []
    for j in range(1, len(bundles_sorted)):
        prev = bundles_sorted[j - 1]
        nxt = bundles_sorted[j]
        ref_gap = (float(prev["ref_end"]), float(nxt["ref_start"]))
        qry_gap = (float(prev["qry_end"]), float(nxt["qry_start"]))
        if ref_gap[1] <= ref_gap[0]:
            continue
        if qry_gap[1] <= qry_gap[0]:
            continue
        gaps.append((ref_gap, qry_gap))

    if not gaps:
        return []

    qids_candidates_set = set(qids_candidates)
    refine_pairs: List[Tuple[str, str]] = []
    used_rids_refine: set = set()

    for (ref_gap, qry_gap) in gaps:
        qids_gap = []
        for qid in qids_candidates_set:
            qc = q_centers.get(qid, np.nan)
            if not np.isfinite(qc):
                continue
            if float(qry_gap[0]) < float(qc) < float(qry_gap[1]):
                qids_gap.append(qid)

        if not qids_gap:
            continue

        unresolved = set(qids_gap)
        gap_result: Dict[str, Dict[str, object]] = {
            qid: {"status": "unmatched", "rid": "", "used_k": 0, "n_cand": 0, "method": ""}
            for qid in qids_gap
        }

        for k in range(1, alpha_L + 1):
            if not unresolved:
                break

            tentative: Dict[str, Dict[str, object]] = {}
            rid_to_qids: Dict[str, List[str]] = defaultdict(list)
            next_unresolved = set()

            for qid in unresolved:
                gid = q_groups.get(qid, "")
                if not gid:
                    tentative[qid] = {"status": "unmatched", "rid": "", "used_k": k, "n_cand": 0, "method": ""}
                    next_unresolved.add(qid)
                    continue

                method = "token_mum_inv" if inversion_mode else "token_mum"
                cand: List[str] = []

                h = q_win_hash_by_qid.get(qid, {}).get(k, None)
                if h is not None:
                    rid = ref_exact_idx.get((gid, k, h), "")
                    if (
                        rid
                        and _alpha_ok(rid, k)
                        and _rid_in_ref_gap(rid, ref_gap)
                        and (rid not in used_rids_1st)
                        and (rid not in used_rids_refine)
                    ):
                        cand = [rid]

                if len(cand) == 0:
                    method = "ed_inv" if inversion_mode else "ed"
                    if k in q_sig.get(qid, {}):
                        Lq, Rq = q_sig[qid][k]
                        if inversion_mode:
                            Lq, Rq = Rq, Lq

                        cc = []
                        for rid in ref_by_group.get(gid, []):
                            if rid in used_rids_1st or rid in used_rids_refine:
                                continue
                            if not _rid_in_ref_gap(rid, ref_gap):
                                continue
                            if not _alpha_ok(rid, k):
                                continue
                            if k not in ref_sig.get(rid, {}):
                                continue
                            Lr, Rr = ref_sig[rid][k]
                            dL = edit_distance_ops_bounded(Lq, Lr, bound=gap_beta)
                            if dL > gap_beta:
                                continue
                            dR = edit_distance_ops_bounded(Rq, Rr, bound=gap_beta - dL)
                            if dL + dR <= int(gap_beta):
                                cc.append(rid)
                        cand = cc

                if len(cand) == 0:
                    tentative[qid] = {"status": "unmatched", "rid": "", "used_k": k, "n_cand": 0, "method": method}
                    next_unresolved.add(qid)
                elif len(cand) == 1:
                    rid = cand[0]
                    tentative[qid] = {"status": "assign", "rid": rid, "used_k": k, "n_cand": 1, "method": method}
                    rid_to_qids[rid].append(qid)
                else:
                    tentative[qid] = {"status": "discordant", "rid": "", "used_k": k, "n_cand": len(cand), "method": method}
                    next_unresolved.add(qid)

            conflict_qids = set()
            for rid, qids in rid_to_qids.items():
                if len(qids) > 1:
                    conflict_qids.update(qids)

            for qid, info in tentative.items():
                if info["status"] == "assign" and qid in conflict_qids:
                    gap_result[qid] = {"status": "discordant", "rid": "", "used_k": info["used_k"], "n_cand": info["n_cand"], "method": info["method"]}
                    next_unresolved.add(qid)
                else:
                    gap_result[qid] = info
                    if info["status"] == "assign":
                        used_rids_refine.add(str(info["rid"]))

            unresolved = next_unresolved

        for qid, info in gap_result.items():
            if info.get("status", "") == "assign" and info.get("rid", ""):
                refine_pairs.append((str(info["rid"]), str(qid)))

    return refine_pairs


# ----------------- 3rd pass: around-bundle refinement + parallel-only keep -----------------
def _in_range(center: float, a: float, b: float) -> bool:
    return (float(a) < float(center)) and (float(center) < float(b))


def _rid_in_ref_range(rid: str, ref_centers: Dict[str, float], ref_range: Tuple[float, float]) -> bool:
    c = ref_centers.get(rid, np.nan)
    if not np.isfinite(c):
        return False
    return _in_range(float(c), float(ref_range[0]), float(ref_range[1]))


def _qid_in_qry_range(qid: str, q_centers_use: Dict[str, float], qry_range: Tuple[float, float]) -> bool:
    c = q_centers_use.get(qid, np.nan)
    if not np.isfinite(c):
        return False
    return _in_range(float(c), float(qry_range[0]), float(qry_range[1]))


def _is_parallel_bundle(
    base_bundle: Dict[str, object],
    cand_bundle: Dict[str, object],
    slope_tol: float = 0.08,
    require_overlap: bool = True,
    overlap_pad_bp: float = 20000.0,
) -> bool:
    s0 = float(base_bundle.get("slope_med", np.nan))
    s1 = float(cand_bundle.get("slope_med", np.nan))
    if not (np.isfinite(s0) and np.isfinite(s1)):
        return False
    if abs(s1 - s0) > float(slope_tol):
        return False

    if require_overlap:
        rs0, re0 = float(base_bundle["ref_start"]), float(base_bundle["ref_end"])
        rs1, re1 = float(cand_bundle["ref_start"]), float(cand_bundle["ref_end"])
        rs0 -= overlap_pad_bp
        re0 += overlap_pad_bp
        if re1 < rs0 or re0 < rs1:
            return False

    return True


def _assign_qids_in_local_ranges(
    chrom: str,
    qids_local: List[str],
    q_groups: Dict[str, str],
    q_sig: Dict[str, Dict[int, Sig]],
    q_win_hash_by_qid_use: Dict[str, Dict[int, Tuple[int, int]]],
    q_centers_use: Dict[str, float],
    ref_exact_idx: Dict[Tuple[str, int, Tuple[int, int]], str],
    ref_sig: Dict[str, Dict[int, Sig]],
    ref_groups: Dict[str, str],
    ref_uid: Dict[str, str],
    ref_amb: Dict[str, bool],
    ref_centers: Dict[str, float],
    ref_k_min_ok: Dict[str, int],
    alpha_L: int,
    gap_beta: int,
    ref_range: Tuple[float, float],
    qry_range: Tuple[float, float],
    used_rids_global: set,
    used_rids_local: set,
    inversion_mode: bool = False,
) -> List[Tuple[str, str]]:
    ref_by_group: Dict[str, List[str]] = defaultdict(list)
    for rid, gid in ref_groups.items():
        if ref_amb.get(rid, True):
            continue
        if ref_uid.get(rid, "") == "":
            continue
        if not gid:
            continue
        if rid in used_rids_global or rid in used_rids_local:
            continue
        if not _rid_in_ref_range(rid, ref_centers, ref_range):
            continue
        ref_by_group[gid].append(rid)

    def _alpha_ok(rid: str, k: int) -> bool:
        k0 = int(ref_k_min_ok.get(rid, 0))
        return (k0 > 0) and (int(k) >= k0)

    qids = [q for q in qids_local if _qid_in_qry_range(q, q_centers_use, qry_range)]
    if not qids:
        return []

    unresolved = set(qids)
    result_pairs: List[Tuple[str, str]] = []

    for k in range(1, alpha_L + 1):
        if not unresolved:
            break

        tentative: Dict[str, Dict[str, object]] = {}
        rid_to_qids: Dict[str, List[str]] = defaultdict(list)
        next_unresolved = set()

        for qid in unresolved:
            gid = q_groups.get(qid, "")
            if not gid:
                tentative[qid] = {"status": "unmatched", "rid": "", "used_k": k, "n_cand": 0, "method": ""}
                next_unresolved.add(qid)
                continue

            cand: List[str] = []
            method = "token_mum"

            h = q_win_hash_by_qid_use.get(qid, {}).get(k, None)
            if h is not None:
                rid = ref_exact_idx.get((gid, k, h), "")
                if rid:
                    if (
                        (rid not in used_rids_global)
                        and (rid not in used_rids_local)
                        and _rid_in_ref_range(rid, ref_centers, ref_range)
                        and _alpha_ok(rid, k)
                    ):
                        cand = [rid]

            if len(cand) == 0:
                method = "ed"
                if k in q_sig.get(qid, {}):
                    Lq, Rq = q_sig[qid][k]
                    if inversion_mode:
                        Lq, Rq = Rq, Lq
                    cc = []
                    for rid in ref_by_group.get(gid, []):
                        if not _alpha_ok(rid, k):
                            continue
                        if k not in ref_sig.get(rid, {}):
                            continue
                        Lr, Rr = ref_sig[rid][k]
                        dL = edit_distance_ops_bounded(Lq, Lr, bound=gap_beta)
                        if dL > gap_beta:
                            continue
                        dR = edit_distance_ops_bounded(Rq, Rr, bound=gap_beta - dL)
                        if dL + dR <= int(gap_beta):
                            cc.append(rid)
                    cand = cc

            if len(cand) == 0:
                tentative[qid] = {"status": "unmatched", "rid": "", "used_k": k, "n_cand": 0, "method": method}
                next_unresolved.add(qid)
            elif len(cand) == 1:
                rid = cand[0]
                tentative[qid] = {"status": "assign", "rid": rid, "used_k": k, "n_cand": 1, "method": method}
                rid_to_qids[rid].append(qid)
            else:
                tentative[qid] = {"status": "discordant", "rid": "", "used_k": k, "n_cand": len(cand), "method": method}
                next_unresolved.add(qid)

        conflict_qids = set()
        for rid, qlist in rid_to_qids.items():
            if len(qlist) > 1:
                conflict_qids.update(qlist)

        for qid, info in tentative.items():
            if info["status"] == "assign" and qid in conflict_qids:
                next_unresolved.add(qid)
                continue
            if info["status"] == "assign" and info["rid"]:
                rid = str(info["rid"])
                used_rids_local.add(rid)
                result_pairs.append((rid, str(qid)))
            else:
                next_unresolved.add(qid)

        unresolved = next_unresolved

    return result_pairs


def third_pass_around_bundles_parallel_only(
    chrom: str,
    base_bundles: List[Dict[str, object]],
    base_items_fwd: List[Tuple[float, float, str, str]],
    base_items_inv: List[Tuple[float, float, str, str]],
    qry_chr_len: int,
    qids_candidates: List[str],
    q_groups: Dict[str, str],
    q_sig: Dict[str, Dict[int, Sig]],
    q_win_hash_fwd_by_qid: Dict[str, Dict[int, Tuple[int, int]]],
    q_win_hash_inv_by_qid: Dict[str, Dict[int, Tuple[int, int]]],
    q_centers: Dict[str, float],
    q_centers_inv: Dict[str, float],
    ref_exact_idx: Dict[Tuple[str, int, Tuple[int, int]], str],
    ref_sig: Dict[str, Dict[int, Sig]],
    ref_groups: Dict[str, str],
    ref_uid: Dict[str, str],
    ref_amb: Dict[str, bool],
    ref_centers: Dict[str, float],
    ref_k_min_ok: Dict[str, int],
    alpha_L: int,
    gap_beta: int,
    window_steps: List[int],
    bundle_max_gap_ref: float,
    bundle_max_gap_qry: float,
    bundle_slope_min: float,
    bundle_slope_max: float,
    bundle_min_links: int,
    bundle_min_span: float,
    used_rids_global: set,
    slope_tol: float = 0.08,
) -> List[Tuple[str, str]]:
    if not base_bundles:
        return []

    def _items_in_ranges(
        items: List[Tuple[float, float, str, str]],
        ref_range: Tuple[float, float],
        qry_range: Tuple[float, float],
    ) -> List[Tuple[float, float, str, str]]:
        out = []
        for rc, qc, rid, qid in items:
            if (float(ref_range[0]) <= float(rc) <= float(ref_range[1])) and (
                float(qry_range[0]) <= float(qc) <= float(qry_range[1])
            ):
                out.append((rc, qc, rid, qid))
        out.sort(key=lambda x: x[0])
        return out

    kept_pairs: List[Tuple[str, str]] = []
    used_rids_local_all: set = set()

    base_bundles = sorted(base_bundles, key=lambda b: (str(b.get("strand", "+")), float(b["ref_start"])))

    for b0 in base_bundles:
        strand0 = str(b0.get("strand", "+"))
        rs0, re0 = float(b0["ref_start"]), float(b0["ref_end"])

        if strand0 == "-":
            raw_qs0, raw_qe0 = float(b0["qry_start"]), float(b0["qry_end"])
            qs0, qe0 = _convert_raw_interval_to_inv(raw_qs0, raw_qe0, qry_chr_len=qry_chr_len)

            base_items = base_items_inv
            q_centers_use = q_centers_inv
            q_win_hash_use = q_win_hash_inv_by_qid
            inversion_mode = True
        else:
            qs0, qe0 = float(b0["qry_start"]), float(b0["qry_end"])
            base_items = base_items_fwd
            q_centers_use = q_centers
            q_win_hash_use = q_win_hash_fwd_by_qid
            inversion_mode = False

        new_pairs_all_windows: List[Tuple[str, str]] = []

        for W in window_steps:
            ref_range = (rs0 - float(W), re0 + float(W))
            qry_range = (qs0 - float(W), qe0 + float(W))

            new_pairs = _assign_qids_in_local_ranges(
                chrom=chrom,
                qids_local=qids_candidates,
                q_groups=q_groups,
                q_sig=q_sig,
                q_win_hash_by_qid_use=q_win_hash_use,
                q_centers_use=q_centers_use,
                ref_exact_idx=ref_exact_idx,
                ref_sig=ref_sig,
                ref_groups=ref_groups,
                ref_uid=ref_uid,
                ref_amb=ref_amb,
                ref_centers=ref_centers,
                ref_k_min_ok=ref_k_min_ok,
                alpha_L=alpha_L,
                gap_beta=gap_beta,
                ref_range=ref_range,
                qry_range=qry_range,
                used_rids_global=used_rids_global,
                used_rids_local=used_rids_local_all,
                inversion_mode=inversion_mode,
            )
            if new_pairs:
                new_pairs_all_windows.extend(new_pairs)

        if not new_pairs_all_windows:
            continue

        Wmax = max(window_steps) if window_steps else 0
        ref_range = (rs0 - float(Wmax), re0 + float(Wmax))
        qry_range = (qs0 - float(Wmax), qe0 + float(Wmax))

        base_local_items = _items_in_ranges(base_items, ref_range, qry_range)

        new_items = []
        for rid, qid in new_pairs_all_windows:
            rc = ref_centers.get(rid, np.nan)
            qc = q_centers_use.get(qid, np.nan)
            if np.isfinite(rc) and np.isfinite(qc):
                new_items.append((float(rc), float(qc), str(rid), str(qid)))

        all_local_items = sorted(base_local_items + new_items, key=lambda x: x[0])
        if not all_local_items:
            continue

        cand_bundles = extract_parallel_bundles(
            items=all_local_items,
            max_gap_ref=bundle_max_gap_ref,
            max_gap_qry=bundle_max_gap_qry,
            slope_min=bundle_slope_min,
            slope_max=bundle_slope_max,
            min_links=bundle_min_links,
            min_span=bundle_min_span,
        )
        if not cand_bundles:
            continue

        parallel_bundles = [bb for bb in cand_bundles if _is_parallel_bundle(b0, bb, slope_tol=slope_tol)]
        if not parallel_bundles:
            continue

        keep_set = set()
        for bb in parallel_bundles:
            for rid, qid in zip(bb["rids"], bb["qids"]):
                keep_set.add((str(rid), str(qid)))

        for rid, qid in new_pairs_all_windows:
            if (str(rid), str(qid)) in keep_set:
                kept_pairs.append((str(rid), str(qid)))

    return list(dict.fromkeys(kept_pairs))


# ----------------- merge parallel bundles (default ON) -----------------
def merge_parallel_bundles(
    bundles: List[Dict[str, object]],
    slope_tol: float = 0.06,
    ref_merge_gap: float = 30000.0,
    require_qry_overlap: bool = True,
    qry_merge_gap: float = 30000.0,
    b0_tol: float = 5000.0,
) -> List[Dict[str, object]]:
    if not bundles:
        return []

    bsorted = sorted(bundles, key=lambda b: float(b["ref_start"]))
    out: List[Dict[str, object]] = []
    cur = dict(bsorted[0])

    seed_slope = float(cur.get("slope_med", np.nan))
    seed_strand = str(cur.get("strand", "+"))

    def _finite_slope(b) -> float:
        return float(b.get("slope_med", np.nan))

    def _qry_close(a, b) -> bool:
        a_qs, a_qe = float(a["qry_start"]), float(a["qry_end"])
        b_qs, b_qe = float(b["qry_start"]), float(b["qry_end"])
        if b_qs > a_qe + float(qry_merge_gap):
            return False
        if a_qs > b_qe + float(qry_merge_gap):
            return False
        return True

    def _can_merge(cur_bundle, next_bundle) -> bool:
        nonlocal seed_slope, seed_strand
        sb = _finite_slope(next_bundle)
        sa = seed_slope

        if str(next_bundle.get("strand", "+")) != str(seed_strand):
            return False

        if not (np.isfinite(sa) and np.isfinite(sb)):
            return False
        if abs(sb - sa) > float(slope_tol):
            return False

        b0a = float(cur_bundle.get("b0_med", np.nan))
        b0b = float(next_bundle.get("b0_med", np.nan))
        if not (np.isfinite(b0a) and np.isfinite(b0b)):
            return False
        if abs(b0b - b0a) > float(b0_tol):
            return False

        if float(next_bundle["ref_start"]) > float(cur_bundle["ref_end"]) + float(ref_merge_gap):
            return False

        if require_qry_overlap:
            if not _qry_close(cur_bundle, next_bundle):
                return False

        return True

    for nb in bsorted[1:]:
        if _can_merge(cur, nb):
            cur["ref_start"] = float(min(float(cur["ref_start"]), float(nb["ref_start"])))
            cur["ref_end"] = float(max(float(cur["ref_end"]), float(nb["ref_end"])))
            cur["qry_start"] = float(min(float(cur["qry_start"]), float(nb["qry_start"])))
            cur["qry_end"] = float(max(float(cur["qry_end"]), float(nb["qry_end"])))
            cur["n_links"] = int(cur.get("n_links", 0)) + int(nb.get("n_links", 0))
            cur["rids"] = list(cur.get("rids", [])) + list(nb.get("rids", []))
            cur["qids"] = list(cur.get("qids", [])) + list(nb.get("qids", []))

            cur_s = float(cur.get("slope_med", np.nan))
            nb_s = float(nb.get("slope_med", np.nan))
            if np.isfinite(cur_s) and np.isfinite(nb_s):
                cur["slope_med"] = float(np.median([cur_s, nb_s]))

            cur["resid_mad"] = float(max(float(cur.get("resid_mad", 0.0)), float(nb.get("resid_mad", 0.0))))
            cur["max_gap_ref"] = float(max(float(cur.get("max_gap_ref", 0.0)), float(nb.get("max_gap_ref", 0.0))))
            cur["max_gap_qry"] = float(max(float(cur.get("max_gap_qry", 0.0)), float(nb.get("max_gap_qry", 0.0))))
            cur["strand"] = seed_strand
        else:
            out.append(cur)
            cur = dict(nb)
            seed_slope = float(cur.get("slope_med", np.nan))
            seed_strand = str(cur.get("strand", "+"))

    out.append(cur)
    return out


# ----------------- main -----------------
def main():
    ap = argparse.ArgumentParser(
        description=(
            "ValleScope (THEORY revision, SPEEDUP): MinHash+LSH grouping + per-k alpha' + "
            "alpha'-gated token-MUM/ED assignment + conflict-aware k-expansion + "
            "between-bundles refinement + 3rd pass around-bundle refinement (parallel-only keep) + "
            "optional merge-parallel bundles + all-links plotting + inversion-mode bundle extraction "
            "+ inversion-aware 3rd pass + (NEW) 1st-pass forward+inversion search."
        )
    )
    ap.add_argument("--beds", required=True, help="comma-separated BED paths (>=2). ref is index 0.")
    ap.add_argument("--fastas", required=True, help="comma-separated FASTA paths (>=2). ref is index 0.")
    ap.add_argument("--names", default="", help="optional comma-separated sample names (same count as beds/fastas)")
    ap.add_argument("--out-prefix", required=True)

    # grouping params
    ap.add_argument("--kmer-k", type=int, default=31, help="k for canonical k-mers")
    ap.add_argument("--num-perm", type=int, default=80, help="MinHash signature length")
    ap.add_argument("--bands", type=int, default=40, help="LSH bands")
    ap.add_argument("--rows", type=int, default=2, help="LSH rows per band")
    ap.add_argument("--seed", type=int, default=42, help="seed for MinHash params")
    ap.add_argument("--group-bands", type=int, default=1, help="Use only first N LSH band hashes to define group_id")

    # context params
    ap.add_argument("--alpha-L", type=int, default=50, help="max context radius k")
    ap.add_argument("--beta", type=float, default=50.0, help="gap bin size in bp")
    ap.add_argument("--gap-beta", type=int, default=20, help="ED threshold for signature fallback")

    # MUS/alpha' params
    ap.add_argument("--alpha-target", type=int, default=1, help="Require alpha' >= this to accept mapping at k")
    ap.add_argument("--mus-Lmax-tokens", type=int, default=50, help="Max token-substring length for MUS search")

    # bundle params
    ap.add_argument("--bundle-min-links", type=int, default=5, help="min links to call a bundle")
    ap.add_argument("--bundle-max-gap-ref", type=float, default=5000.0, help="max gap between adjacent links on ref")
    ap.add_argument("--bundle-max-gap-qry", type=float, default=5000.0, help="max gap between adjacent links on query")
    ap.add_argument("--bundle-slope-min", type=float, default=0.99, help="min allowed local slope dq/dr")
    ap.add_argument("--bundle-slope-max", type=float, default=1.01, help="max allowed local slope dq/dr")
    ap.add_argument("--bundle-min-span", type=float, default=1000.0, help="min span on ref to keep bundle")

    # 3rd pass params
    ap.add_argument("--third-pass-windows", type=str, default="10000,20000,50000,100000", help="comma-separated bp windows")
    ap.add_argument("--third-parallel-slope-tol", type=float, default=0.08, help="slope tolerance for 'parallel-only keep'")

    # merge-parallel bundles
    ap.add_argument("--merge-parallel", dest="merge_parallel", action="store_true", help="merge parallel bundles (default: ON)")
    ap.add_argument("--no-merge-parallel", dest="merge_parallel", action="store_false", help="disable merge-parallel")
    ap.set_defaults(merge_parallel=True)
    ap.add_argument("--merge-slope-tol", type=float, default=0.06, help="merge slope tolerance")
    ap.add_argument("--merge-ref-gap", type=float, default=50000.0, help="merge if ref gap <= this")
    ap.add_argument("--merge-b0-tol", type=float, default=5000.0, help="merge intercept(b0) tolerance in bp")
    ap.add_argument("--merge-require-qry-overlap", dest="merge_require_qry_overlap", action="store_true", help="require qry overlap/adjacency (default: ON)")
    ap.add_argument("--no-merge-require-qry-overlap", dest="merge_require_qry_overlap", action="store_false", help="do not require qry overlap/adjacency")
    ap.set_defaults(merge_require_qry_overlap=True)
    ap.add_argument("--merge-qry-gap", type=float, default=50000.0, help="merge if qry gap <= this when qry-overlap is on")

    # rolling-hash seed
    ap.add_argument("--rh-seed", type=int, default=1337, help="seed for rolling hash bases")

    # ID mismatch check options
    ap.add_argument("--id-mismatch-mode", choices=["error", "warn", "off"], default="error", help="BED/FASTA ID mismatch policy")

    # vsearch grouping params
    ap.add_argument("--grouping", choices=["minhash", "vsearch"], default="minhash", help="grouping method for anchor sequences")
    ap.add_argument("--vsearch-path", default="vsearch", help="path to vsearch executable")
    ap.add_argument("--vsearch-id", type=float, default=0.99, help="vsearch clustering identity threshold")
    ap.add_argument("--vsearch-threads", type=int, default=8, help="vsearch threads")
    ap.add_argument("--vsearch-keep-tmp", action="store_true", help="keep vsearch temp files")
    ap.add_argument("--vsearch-tmpdir", default="", help="optional tmp dir for vsearch outputs")

    args = ap.parse_args()

    pref = args.out_prefix
    outdir = os.path.dirname(pref)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    bed_paths = _split_csv_arg(args.beds)
    fa_paths = _split_csv_arg(args.fastas)
    if len(bed_paths) != len(fa_paths):
        raise ValueError("beds and fastas must have same number of items")
    n = len(bed_paths)
    if n < 2:
        raise ValueError("Need >=2 samples (ref + >=1 query)")

    ref_idx = 0

    # names
    if args.names:
        name_list = _split_csv_arg(args.names)
        if len(name_list) != n:
            raise ValueError("--names must match number of samples")
    else:
        raw = [_derive_sample_name_from_bed(p) for p in bed_paths]
        seen: Dict[str, int] = {}
        name_list = []
        for nm in raw:
            if nm not in seen:
                seen[nm] = 1
                name_list.append(nm)
            else:
                seen[nm] += 1
                name_list.append(f"{nm}_{seen[nm]}")

    ref_name = name_list[ref_idx]

    third_windows = [int(x) for x in _split_csv_arg(args.third_pass_windows)]
    third_windows = [w for w in third_windows if w > 0] or [10000, 20000, 50000, 100000]

    print(f"[info] samples={n} ref={ref_name} names={', '.join(name_list)}")
    print(f"[info] context: alpha_L={args.alpha_L} beta={args.beta} gap_beta={args.gap_beta}")
    print(f"[info] MUS/alpha': alpha_target={args.alpha_target} mus_Lmax_tokens={args.mus_Lmax_tokens}")
    print(f"[info] 3rd-pass windows={third_windows} parallel_slope_tol={args.third_parallel_slope_tol}")
    print(
        f"[info] merge_parallel={args.merge_parallel} merge_slope_tol={args.merge_slope_tol} "
        f"merge_ref_gap={args.merge_ref_gap} merge_require_qry_overlap={args.merge_require_qry_overlap} merge_qry_gap={args.merge_qry_gap}"
    )

    # load
    beds = [load_bed(p) for p in bed_paths]
    fastas = [load_fasta(p) for p in fa_paths]

    # --- ID consistency check ---
    mode = args.id_mismatch_mode
    max_missing_frac = float(0.0)
    max_missing_n = int(0)

    if mode != "off":
        for i in range(n):
            bed_ids = set(beds[i]["anchor_id"].astype(str).tolist())
            fa_ids = set(str(k) for k in fastas[i].keys())

            missing = sorted(bed_ids - fa_ids)
            extra = sorted(fa_ids - bed_ids)

            n_bed = len(bed_ids)
            n_miss = len(missing)
            miss_frac = (n_miss / n_bed) if n_bed else 0.0

            print(f"[idcheck] sample={name_list[i]} BED={n_bed} FASTA={len(fa_ids)} missing={n_miss} extra={len(extra)}")

            ok = (n_miss <= max_missing_n) or (miss_frac <= max_missing_frac)
            if (mode == "error") and (not ok):
                ex = ", ".join(missing[:20])
                raise ValueError(
                    f"[idcheck][FAIL] sample={name_list[i]}\n"
                    f"  BED : {bed_paths[i]}\n"
                    f"  FA  : {fa_paths[i]}\n"
                    f"  Missing in FASTA: {n_miss} ({miss_frac:.4%})\n"
                    f"  Examples: {ex}\n"
                )
            if (mode == "warn") and (n_miss > 0):
                ex = ", ".join(missing[:20])
                print(f"[idcheck][WARN] sample={name_list[i]} missing_in_fasta={n_miss} examples: {ex}")

            if len(extra) > 0:
                ex = ", ".join(extra[:20])
                print(f"[idcheck][WARN] sample={name_list[i]} extra_in_fasta={len(extra)} examples: {ex}")

    # --- end ID consistency check ---
    centers: List[Dict[str, float]] = []
    chroms: List[Dict[str, str]] = []
    strands: List[Dict[str, str]] = []
    for i in range(n):
        di = beds[i].set_index("anchor_id")
        centers.append(di["center"].to_dict())
        chroms.append(di["chrom"].to_dict())
        strands.append(di["strand"].to_dict())

    # 1) grouping for ALL samples
    if args.grouping == "vsearch":
        groups = build_groups_all_samples_vsearch(
            fastas=fastas,
            sample_names=name_list,
            vsearch_path=str(args.vsearch_path),
            identity=float(args.vsearch_id),
            threads=int(args.vsearch_threads),
            tmp_dir=(str(args.vsearch_tmpdir) if str(args.vsearch_tmpdir) else None),
            keep_tmp=bool(args.vsearch_keep_tmp),
        )
    else:
        groups = build_groups_all_samples(
            fastas=fastas,
            kmer_k=int(args.kmer_k),
            num_perm=int(args.num_perm),
            bands=int(args.bands),
            rows=int(args.rows),
            seed=int(args.seed),
            group_bands=int(args.group_bands),
        )

    out_groups_tsv = f"{pref}.anchor_groups.tsv"
    rows_out = []
    for i in range(n):
        nm = name_list[i]
        for aid, gid in groups[i].items():
            rows_out.append(
                {
                    "sample": nm,
                    "anchor_id": aid,
                    "group_id": gid,
                    "chrom": chroms[i].get(aid, ""),
                    "center": centers[i].get(aid, np.nan),
                    "strand": strands[i].get(aid, ""),
                    "len": len(fastas[i].get(aid, "")),
                }
            )
    pd.DataFrame(rows_out).to_csv(out_groups_tsv, sep="\t", index=False)
    print(f"[out] {out_groups_tsv}")

    # 2) reference dictionary
    ref_order_by_chrom = build_order_by_chrom(beds[ref_idx])
    ref_groups = groups[ref_idx]

    ref_T_by, ref_a2t_by, _ = build_token_sequences_all_chrom(
        order_by_chrom=ref_order_by_chrom,
        centers=centers[ref_idx],
        groups=ref_groups,
        beta=float(args.beta),
    )

    ref_hasher_by = build_rolling_hashers(ref_T_by, seed=int(args.rh_seed))
    ref_win_hash = precompute_anchor_window_hashes(
        order_by_chrom=ref_order_by_chrom,
        a2t_by=ref_a2t_by,
        hasher_by_chrom=ref_hasher_by,
        alpha_L=int(args.alpha_L),
    )

    ref_sig = build_signatures(
        order_by_chrom=ref_order_by_chrom,
        centers=centers[ref_idx],
        groups=ref_groups,
        alpha_L=int(args.alpha_L),
        beta=float(args.beta),
    )

    ref_k_min_ok, ref_alpha_prime_at_min, ref_evidence_at_min = compute_ref_alpha_prime_min_k(
        ref_order_by_chrom=ref_order_by_chrom,
        ref_groups=ref_groups,
        ref_hasher_by=ref_hasher_by,
        ref_a2t_by=ref_a2t_by,
        alpha_L=int(args.alpha_L),
        mus_Lmax_tokens=int(args.mus_Lmax_tokens),
        alpha_target=int(args.alpha_target),
    )

    ref_uid, ref_amb, ref_used_k, ref_alpha_prime_used, ref_evidence_used = assign_unique_ids_reference_mus(
        ref_order_by_chrom=ref_order_by_chrom,
        ref_groups=ref_groups,
        ref_k_min_ok=ref_k_min_ok,
        ref_alpha_prime_at_min=ref_alpha_prime_at_min,
        ref_evidence_at_min=ref_evidence_at_min,
    )

    ref_dict_tsv = f"{pref}.ref_dictionary.tsv"
    ref_rows = []
    for chrom, order in ref_order_by_chrom.items():
        for rid in order:
            ref_rows.append(
                {
                    "ref_anchor_id": rid,
                    "uid": ref_uid.get(rid, ""),
                    "ambiguous": bool(ref_amb.get(rid, True)),
                    "group_id": ref_groups.get(rid, ""),
                    "chrom": chroms[ref_idx].get(rid, ""),
                    "center": centers[ref_idx].get(rid, np.nan),
                    "strand": strands[ref_idx].get(rid, ""),
                    "len": len(fastas[ref_idx].get(rid, "")),
                    "used_k": int(ref_used_k.get(rid, 0)),
                    "alpha_prime": int(ref_alpha_prime_used.get(rid, 0)),
                    "evidence": ref_evidence_used.get(rid, ""),
                }
            )
    pd.DataFrame(ref_rows).to_csv(ref_dict_tsv, sep="\t", index=False)
    print(f"[out] {ref_dict_tsv}")

    ref_exact_idx = build_ref_exact_index_unique(
        ref_order_by_chrom=ref_order_by_chrom,
        ref_groups=ref_groups,
        ref_uid=ref_uid,
        ref_amb=ref_amb,
        ref_win_hash=ref_win_hash,
        alpha_L=int(args.alpha_L),
    )
    print(f"[ref] exact_index_keys(unique)={len(ref_exact_idx)}")

    # 3) queries
    assign_rows = []
    for i in range(n):
        if i == ref_idx:
            continue

        nm = name_list[i]
        q_order_by_chrom = build_order_by_chrom(beds[i])
        q_groups = groups[i]

        q_T_by, q_a2t_by, _ = build_token_sequences_all_chrom(
            order_by_chrom=q_order_by_chrom,
            centers=centers[i],
            groups=q_groups,
            beta=float(args.beta),
        )

        q_T_inv_by, q_a2t_inv_by = build_inversion_token_sequences_from_raw(q_T_by, q_a2t_by)

        q_hasher_by = build_rolling_hashers(q_T_by, seed=int(args.rh_seed))
        q_hasher_inv_by = build_rolling_hashers(q_T_inv_by, seed=int(args.rh_seed))

        q_win_hash_fwd = precompute_anchor_window_hashes(
            order_by_chrom=q_order_by_chrom,
            a2t_by=q_a2t_by,
            hasher_by_chrom=q_hasher_by,
            alpha_L=int(args.alpha_L),
        )
        q_win_hash_inv = precompute_anchor_window_hashes(
            order_by_chrom=q_order_by_chrom,
            a2t_by=q_a2t_inv_by,
            hasher_by_chrom=q_hasher_inv_by,
            alpha_L=int(args.alpha_L),
        )

        q_win_hash_by_chrom_fwd: Dict[str, Dict[str, Dict[int, Tuple[int, int]]]] = defaultdict(dict)
        q_win_hash_by_chrom_inv: Dict[str, Dict[str, Dict[int, Tuple[int, int]]]] = defaultdict(dict)
        for chrom, order in q_order_by_chrom.items():
            for qid in order:
                if qid in q_win_hash_fwd:
                    q_win_hash_by_chrom_fwd[chrom][qid] = q_win_hash_fwd[qid]
                if qid in q_win_hash_inv:
                    q_win_hash_by_chrom_inv[chrom][qid] = q_win_hash_inv[qid]

        q_sig = build_signatures(
            order_by_chrom=q_order_by_chrom,
            centers=centers[i],
            groups=q_groups,
            alpha_L=int(args.alpha_L),
            beta=float(args.beta),
        )

        qid_to_chrom: Dict[str, str] = {}
        for ch, order in q_order_by_chrom.items():
            for qid in order:
                qid_to_chrom[qid] = ch

        # ---- 1st pass (NOW: fwd + inv) ----
        q_res_1st = match_query_to_reference_conflict_resolve(
            q_order_by_chrom=q_order_by_chrom,
            q_groups=q_groups,
            q_sig=q_sig,
            q_win_hash_by_chrom_fwd=q_win_hash_by_chrom_fwd,
            q_win_hash_by_chrom_inv=q_win_hash_by_chrom_inv,
            qid_to_chrom=qid_to_chrom,
            ref_exact_idx=ref_exact_idx,
            ref_sig=ref_sig,
            ref_groups=ref_groups,
            ref_uid=ref_uid,
            ref_amb=ref_amb,
            ref_k_min_ok=ref_k_min_ok,
            alpha_L=int(args.alpha_L),
            gap_beta=int(args.gap_beta),
        )

        for qid, r in q_res_1st.items():
            assign_rows.append(
                {
                    "sample": nm,
                    "query_anchor_id": qid,
                    "status": r["status"],
                    "method": r.get("method", ""),
                    "assign_strand": r.get("strand", ""),
                    "assigned_ref_anchor_id": r["rid"],
                    "assigned_uid": r["uid"],
                    "used_k": int(r["used_k"]),
                    "n_cand": int(r["n_cand"]),
                    "query_group_id": q_groups.get(qid, ""),
                    "chrom": chroms[i].get(qid, ""),
                    "center": centers[i].get(qid, np.nan),
                    "strand": strands[i].get(qid, ""),
                    "len": len(fastas[i].get(qid, "")),
                }
            )

        assigned_pairs_1st: List[Tuple[str, str]] = []
        for qid, rr in q_res_1st.items():
            if rr["status"] == "assign" and rr["rid"]:
                assigned_pairs_1st.append((str(rr["rid"]), str(qid)))

        used_rids_1st = set(rid for rid, _ in assigned_pairs_1st)

        plot_links_by_chrom(
            ref_name=ref_name,
            query_name=nm,
            ref_centers=centers[ref_idx],
            ref_chroms=chroms[ref_idx],
            q_centers=centers[i],
            q_chroms=chroms[i],
            pairs=assigned_pairs_1st,
            out_prefix=pref,
            max_lines=20000,
        )

        # ---- 2nd pass ----
        qids_unassigned = [qid for qid, rr in q_res_1st.items() if rr.get("status", "") != "assign"]
        assigned_pairs_refine: List[Tuple[str, str]] = []

        qry_chr_len_map_2pass = infer_chrom_lengths_from_bed(beds[i])
        q_centers_inv_2pass: Dict[str, float] = {}
        for qid, qc in centers[i].items():
            ch = chroms[i].get(qid, "")
            if not ch:
                continue
            qlen_ch = int(qry_chr_len_map_2pass.get(str(ch), 0))
            if qlen_ch <= 0:
                continue
            q_centers_inv_2pass[qid] = float(qlen_ch) - float(qc)

        chrom_keys_2pass = sorted(
            set(_pair_chrom_key(rid, qid, chroms[ref_idx], chroms[i]) for rid, qid in assigned_pairs_1st)
        )
        chrom_keys_2pass = [c for c in chrom_keys_2pass if c]

        for chrom in chrom_keys_2pass:
            items_1st_fwd = _items_from_pairs_same_chrom(
                pairs=assigned_pairs_1st,
                chrom=chrom,
                ref_centers=centers[ref_idx],
                ref_chroms=chroms[ref_idx],
                q_centers=centers[i],
                q_chroms=chroms[i],
            )
            if items_1st_fwd:
                bundles_1st_fwd = extract_parallel_bundles(
                    items=items_1st_fwd,
                    max_gap_ref=float(args.bundle_max_gap_ref),
                    max_gap_qry=float(args.bundle_max_gap_qry),
                    slope_min=float(args.bundle_slope_min),
                    slope_max=float(args.bundle_slope_max),
                    min_links=int(args.bundle_min_links),
                    min_span=float(args.bundle_min_span),
                    max_skip=5,
                    scale_gap_with_skip=True,
                )

                if bundles_1st_fwd:
                    refine_pairs_chr_fwd = refine_between_bundles(
                        chrom=chrom,
                        bundles_1st=bundles_1st_fwd,
                        qids_candidates=qids_unassigned,
                        q_groups=q_groups,
                        q_sig=q_sig,
                        q_win_hash_by_qid=q_win_hash_fwd,
                        q_centers=centers[i],
                        ref_exact_idx=ref_exact_idx,
                        ref_sig=ref_sig,
                        ref_groups=ref_groups,
                        ref_uid=ref_uid,
                        ref_amb=ref_amb,
                        ref_centers=centers[ref_idx],
                        ref_k_min_ok=ref_k_min_ok,
                        alpha_L=int(args.alpha_L),
                        gap_beta=int(args.gap_beta),
                        used_rids_1st=used_rids_1st,
                        inversion_mode=False,
                    )
                    if refine_pairs_chr_fwd:
                        assigned_pairs_refine.extend(refine_pairs_chr_fwd)

            if items_1st_fwd:
                qlen_chr = int(qry_chr_len_map_2pass.get(str(chrom), 0))
                if qlen_chr <= 0:
                    if (beds[i]["chrom"].astype(str) == str(chrom)).any():
                        qlen_chr = int(np.ceil(beds[i].loc[beds[i]["chrom"].astype(str) == str(chrom), "end"].max()))
                if qlen_chr > 0:
                    items_1st_inv = _make_items_inversion_mode(items_fwd=items_1st_fwd, qry_chr_len=qlen_chr)
                    bundles_1st_inv = extract_parallel_bundles(
                        items=items_1st_inv,
                        max_gap_ref=float(args.bundle_max_gap_ref),
                        max_gap_qry=float(args.bundle_max_gap_qry),
                        slope_min=float(args.bundle_slope_min),
                        slope_max=float(args.bundle_slope_max),
                        min_links=int(args.bundle_min_links),
                        min_span=float(args.bundle_min_span),
                        max_skip=5,
                        scale_gap_with_skip=True,
                    )

                    if bundles_1st_inv:
                        refine_pairs_chr_inv = refine_between_bundles(
                            chrom=chrom,
                            bundles_1st=bundles_1st_inv,
                            qids_candidates=qids_unassigned,
                            q_groups=q_groups,
                            q_sig=q_sig,
                            q_win_hash_by_qid=q_win_hash_inv,
                            q_centers=q_centers_inv_2pass,
                            ref_exact_idx=ref_exact_idx,
                            ref_sig=ref_sig,
                            ref_groups=ref_groups,
                            ref_uid=ref_uid,
                            ref_amb=ref_amb,
                            ref_centers=centers[ref_idx],
                            ref_k_min_ok=ref_k_min_ok,
                            alpha_L=int(args.alpha_L),
                            gap_beta=int(args.gap_beta),
                            used_rids_1st=used_rids_1st,
                            inversion_mode=True,
                        )
                        if refine_pairs_chr_inv:
                            assigned_pairs_refine.extend(refine_pairs_chr_inv)

        assigned_pairs_all = list(dict.fromkeys(list(assigned_pairs_1st) + list(assigned_pairs_refine)))
        print(
            f"[refine] sample={nm} pairs_1st={len(assigned_pairs_1st)} "
            f"pairs_refine={len(assigned_pairs_refine)} pairs_all_2pass={len(assigned_pairs_all)}"
        )

        # ---- 3rd pass ----
        used_rids_global = set(rid for rid, _ in assigned_pairs_all)
        assigned_qids_all = set(qid for _, qid in assigned_pairs_all)
        qids_unassigned_after2 = [q for q in qids_unassigned if q not in assigned_qids_all]

        assigned_pairs_3rd: List[Tuple[str, str]] = []

        qry_chr_len_map = infer_chrom_lengths_from_bed(beds[i])

        q_centers_inv: Dict[str, float] = {}
        for qid, qc in centers[i].items():
            ch = chroms[i].get(qid, "")
            if not ch:
                continue
            qlen_ch = int(qry_chr_len_map.get(str(ch), 0))
            if qlen_ch <= 0:
                continue
            q_centers_inv[qid] = float(qlen_ch) - float(qc)

        chrom_keys_3pass = sorted(set(_pair_chrom_key(rid, qid, chroms[ref_idx], chroms[i]) for rid, qid in assigned_pairs_all))
        chrom_keys_3pass = [c for c in chrom_keys_3pass if c]

        for chrom in chrom_keys_3pass:
            base_items_chr_fwd = _items_from_pairs_same_chrom(
                pairs=assigned_pairs_all,
                chrom=chrom,
                ref_centers=centers[ref_idx],
                ref_chroms=chroms[ref_idx],
                q_centers=centers[i],
                q_chroms=chroms[i],
            )
            if not base_items_chr_fwd:
                continue

            qlen_chr = int(qry_chr_len_map.get(chrom, 0))
            if qlen_chr <= 0:
                qlen_chr = int(np.ceil(beds[i].loc[beds[i]["chrom"].astype(str) == chrom, "end"].max())) if (beds[i]["chrom"].astype(str) == chrom).any() else 0
            base_items_chr_inv = _make_items_inversion_mode(base_items_chr_fwd, qry_chr_len=qlen_chr)

            base_bundles_fwd = extract_parallel_bundles(
                items=base_items_chr_fwd,
                max_gap_ref=float(args.bundle_max_gap_ref),
                max_gap_qry=float(args.bundle_max_gap_qry),
                slope_min=float(args.bundle_slope_min),
                slope_max=float(args.bundle_slope_max),
                min_links=int(args.bundle_min_links),
                min_span=float(args.bundle_min_span),
            )
            for b in base_bundles_fwd:
                b["strand"] = "+"

            base_bundles_inv = extract_parallel_bundles(
                items=base_items_chr_inv,
                max_gap_ref=float(args.bundle_max_gap_ref),
                max_gap_qry=float(args.bundle_max_gap_qry),
                slope_min=float(args.bundle_slope_min),
                slope_max=float(args.bundle_slope_max),
                min_links=int(args.bundle_min_links),
                min_span=float(args.bundle_min_span),
            )
            for b in base_bundles_inv:
                raw_qs, raw_qe = _convert_bundle_qry_interval_inv_to_raw(b, qry_chr_len=qlen_chr)
                b["qry_start"] = float(raw_qs)
                b["qry_end"] = float(raw_qe)
                b["strand"] = "-"

            base_bundles_chr = list(base_bundles_fwd) + list(base_bundles_inv)
            if not base_bundles_chr:
                continue

            kept_pairs_chr = third_pass_around_bundles_parallel_only(
                chrom=chrom,
                base_bundles=base_bundles_chr,
                base_items_fwd=base_items_chr_fwd,
                base_items_inv=base_items_chr_inv,
                qry_chr_len=qlen_chr,
                qids_candidates=qids_unassigned_after2,
                q_groups=q_groups,
                q_sig=q_sig,
                q_win_hash_fwd_by_qid=q_win_hash_fwd,
                q_win_hash_inv_by_qid=q_win_hash_inv,
                q_centers=centers[i],
                q_centers_inv=q_centers_inv,
                ref_exact_idx=ref_exact_idx,
                ref_sig=ref_sig,
                ref_groups=ref_groups,
                ref_uid=ref_uid,
                ref_amb=ref_amb,
                ref_centers=centers[ref_idx],
                ref_k_min_ok=ref_k_min_ok,
                alpha_L=int(args.alpha_L),
                gap_beta=int(args.gap_beta),
                window_steps=third_windows,
                bundle_max_gap_ref=float(args.bundle_max_gap_ref),
                bundle_max_gap_qry=float(args.bundle_max_gap_qry),
                bundle_slope_min=float(args.bundle_slope_min),
                bundle_slope_max=float(args.bundle_slope_max),
                bundle_min_links=int(args.bundle_min_links),
                bundle_min_span=float(args.bundle_min_span),
                used_rids_global=used_rids_global,
                slope_tol=float(args.third_parallel_slope_tol),
            )

            if kept_pairs_chr:
                assigned_pairs_3rd.extend(kept_pairs_chr)
                for rid, _ in kept_pairs_chr:
                    used_rids_global.add(str(rid))

        assigned_pairs_3rd = list(dict.fromkeys(assigned_pairs_3rd))
        assigned_pairs_all = list(dict.fromkeys(list(assigned_pairs_all) + list(assigned_pairs_3rd)))

        print(f"[3rd] sample={nm} q_unassigned_after2={len(qids_unassigned_after2)} pairs_3rd_kept={len(assigned_pairs_3rd)} pairs_all_now={len(assigned_pairs_all)}")

        plot_links_by_chrom(
            ref_name=ref_name,
            query_name=f"{nm}.all",
            ref_centers=centers[ref_idx],
            ref_chroms=chroms[ref_idx],
            q_centers=centers[i],
            q_chroms=chroms[i],
            pairs=assigned_pairs_all,
            out_prefix=pref,
            max_lines=20000,
        )

        # ---- bundle outputs using ALL links ----
        # NOTE: PAF output removed. TSV output remains.
        qry_chr_len_map2 = infer_chrom_lengths_from_bed(beds[i])

        out_bundles_tsv = f"{pref}.bundles.tsv"
        if os.path.exists(out_bundles_tsv):
            os.remove(out_bundles_tsv)

        bundle_id_counter = 1

        by_chrom_all: Dict[str, List[Tuple[float, float, str, str]]] = defaultdict(list)
        for rid, qid in assigned_pairs_all:
            chrom_key = _pair_chrom_key(rid, qid, chroms[ref_idx], chroms[i])
            if not chrom_key:
                continue
            rc = centers[ref_idx].get(rid, np.nan)
            qc = centers[i].get(qid, np.nan)
            if not np.isfinite(rc) or not np.isfinite(qc):
                continue
            by_chrom_all[chrom_key].append((float(rc), float(qc), str(rid), str(qid)))

        for chrom, items in by_chrom_all.items():
            chrom = str(chrom)
            items = sorted(items, key=lambda x: x[0])

            bundles_fwd = extract_parallel_bundles(
                items=items,
                max_gap_ref=float(args.bundle_max_gap_ref),
                max_gap_qry=float(args.bundle_max_gap_qry),
                slope_min=float(args.bundle_slope_min),
                slope_max=float(args.bundle_slope_max),
                min_links=int(args.bundle_min_links),
                min_span=float(args.bundle_min_span),
            )
            for b in bundles_fwd:
                b["strand"] = "+"

            qlen_chr = int(qry_chr_len_map2.get(chrom, 0))
            if qlen_chr <= 0:
                qlen_chr = int(np.ceil(max(b["qry_end"] for b in bundles_fwd))) if bundles_fwd else 0
                if qlen_chr <= 0 and (beds[i]["chrom"].astype(str) == chrom).any():
                    qlen_chr = int(np.ceil(beds[i].loc[beds[i]["chrom"].astype(str) == chrom, "end"].max()))

            items_inv = _make_items_inversion_mode(items_fwd=items, qry_chr_len=qlen_chr)
            bundles_inv = extract_parallel_bundles(
                items=items_inv,
                max_gap_ref=float(args.bundle_max_gap_ref),
                max_gap_qry=float(args.bundle_max_gap_qry),
                slope_min=float(args.bundle_slope_min),
                slope_max=float(args.bundle_slope_max),
                min_links=int(args.bundle_min_links),
                min_span=float(args.bundle_min_span),
            )
            for b in bundles_inv:
                raw_qs, raw_qe = _convert_bundle_qry_interval_inv_to_raw(b, qry_chr_len=qlen_chr)
                b["qry_start"] = float(raw_qs)
                b["qry_end"] = float(raw_qe)
                b["strand"] = "-"

            bundles_all = list(bundles_fwd) + list(bundles_inv)
            if not bundles_all:
                continue

            if bool(args.merge_parallel):
                before = len(bundles_all)
                bundles_all = merge_parallel_bundles(
                    bundles=bundles_all,
                    slope_tol=float(args.merge_slope_tol),
                    ref_merge_gap=float(args.merge_ref_gap),
                    require_qry_overlap=bool(args.merge_require_qry_overlap),
                    qry_merge_gap=float(args.merge_qry_gap),
                    b0_tol=float(args.merge_b0_tol),
                )
                after = len(bundles_all)
                if after != before:
                    print(f"[merge] {nm} {chrom}: bundles {before} -> {after}")

            plot_bundle_pngs_by_chrom(
                ref_name=ref_name,
                query_name=nm,
                chrom=str(chrom),
                items=items,
                bundles=bundles_all,
                out_prefix=pref,
                max_lines_bg=20000,
            )

            bundle_id_counter = write_bundles_tsv(
                out_tsv=out_bundles_tsv,
                ref_name=ref_name,
                qry_name=nm,
                chrom=str(chrom),
                bundles=bundles_all,
                bundle_id_start=bundle_id_counter,
            )

        if os.path.exists(out_bundles_tsv):
            print(f"[out] {out_bundles_tsv}")

        st = pd.Series([v["status"] for v in q_res_1st.values()]).value_counts().to_dict()
        st_m = pd.Series([v.get("method", "") for v in q_res_1st.values() if v["status"] == "assign"]).value_counts().to_dict()
        print(f"[assign-1st] sample={nm} {st} assign_method={st_m}")

        rid_counts = Counter([v["rid"] for v in q_res_1st.values() if v["status"] == "assign" and v["rid"]])
        dup = {rid: c for rid, c in rid_counts.items() if c > 1}
        if dup:
            print(f"[WARN] duplicate rid assignments in 1st pass: {len(dup)} rids")
            for rid, c in list(dup.items())[:10]:
                print(f"  rid={rid} count={c}")

    assign_tsv = f"{pref}.assignments.tsv"
    pd.DataFrame(assign_rows).to_csv(assign_tsv, sep="\t", index=False)
    print(f"[out] {assign_tsv}")
    print("[done]")


if __name__ == "__main__":
    main()