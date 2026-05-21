#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# src/vallescope/cli.py

from __future__ import annotations

from typing import List, Optional, Tuple
import os
import gzip
import shutil
import typer
from rich import print
import shlex
import glob
from importlib.resources import files, as_file

from .runner import (
    run_snakemake,
    check_requirements,
)

app = typer.Typer(add_completion=False, help="ValleScope: anchors & plots CLI")

@app.callback()
def main() -> None:
    """
    ValleScope: context-based correspondence assignment for satellite DNA.
    """
    pass

# ---------- helpers ----------
def _stem(p: str) -> str:
    """a.fa.gz -> a / a.fasta -> a"""
    name = os.path.basename(p)
    for suf in (".fa.gz", ".fasta.gz", ".fna.gz", ".fa", ".fasta", ".fna"):
        if name.endswith(suf):
            return name[: -len(suf)]
    return os.path.splitext(name)[0]


def _gzip_if_needed(src: str, dst_gz: str) -> None:
    os.makedirs(os.path.dirname(dst_gz), exist_ok=True)
    if src.endswith(".gz"):
        if not os.path.exists(dst_gz):
            try:
                os.symlink(os.path.abspath(src), dst_gz)
            except Exception:
                shutil.copy2(src, dst_gz)
        return
    with open(src, "rb") as fi, gzip.open(dst_gz, "wb", compresslevel=6) as fo:
        shutil.copyfileobj(fi, fo)


def _symlink_or_copy(src: str, dst: str) -> None:
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if not os.path.exists(dst):
        try:
            os.symlink(os.path.abspath(src), dst)
        except Exception:
            shutil.copy2(src, dst)


def _expand_inputs(values: Optional[List[str]]) -> List[str]:
    """
    Expand CLI option values:
      - repeated flags:   --f a.fa -f b.fa
      - single string:    --f "a.fa b.fa" or --f "a.fa,b.fa"
      - glob patterns:    --f "*.fa"
      - list file:        --f @list.txt   (one path per line)
    """
    if not values:
        return []

    tokens: List[str] = []
    for v in values:
        for t in shlex.split(v):
            if t.startswith("@"):
                path = t[1:]
                with open(path) as fh:
                    for line in fh:
                        s = line.strip()
                        if not s or s.startswith("#"):
                            continue
                        tokens.append(s)
            else:
                tokens.extend([s for s in t.split(",") if s])

    expanded: List[str] = []
    for tok in tokens:
        hits = glob.glob(tok)
        if hits:
            expanded.extend(sorted(hits))
        else:
            expanded.append(tok)
    return expanded


def _real(p: str) -> str:
    return os.path.realpath(os.path.abspath(p))


def _find_index_by_path(paths: List[str], target: str) -> int:
    rt = _real(target)
    for i, p in enumerate(paths):
        if _real(p) == rt:
            return i
    return -1


def _read_fasta_lengths(path: str) -> List[Tuple[str, int]]:
    """
    Return list of (seq_name, length) for a FASTA/FASTA.GZ.
    """
    out: List[Tuple[str, int]] = []
    name: Optional[str] = None
    length = 0

    if path.endswith(".gz"):
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if not line:
                    continue
                if line.startswith(">"):
                    if name is not None:
                        out.append((name, length))
                    name = line[1:].strip().split()[0]
                    length = 0
                else:
                    length += len(line.strip())
    else:
        with open(path, "rt") as fh:
            for line in fh:
                if not line:
                    continue
                if line.startswith(">"):
                    if name is not None:
                        out.append((name, length))
                    name = line[1:].strip().split()[0]
                    length = 0
                else:
                    length += len(line.strip())

    if name is not None:
        out.append((name, length))
    return out


def _write_whole_genome_bed_from_fasta(fasta_path: str, out_bed: str, label: str) -> None:
    """
    Create a BED6 that covers the whole length of each FASTA record:
      chrom  start  end  label  score  strand
    This makes mask_by_bed_labels(--invert) effectively do "no masking".
    """
    os.makedirs(os.path.dirname(out_bed), exist_ok=True)
    recs = _read_fasta_lengths(fasta_path)
    with open(out_bed, "wt") as fo:
        for chrom, seq_len in recs:
            if seq_len <= 0:
                continue
            fo.write(f"{chrom}\t0\t{seq_len}\t{label}\t0\t+\n")


def _load_mask_label(cfg_path: str, default: str = "HSat2") -> str:
    """
    Best-effort YAML parsing to get config.mask_label.label.
    If PyYAML isn't installed or parsing fails, return default.
    """
    try:
        import yaml

        with open(cfg_path, "rt") as fh:
            cfg = yaml.safe_load(fh) or {}
        return str(((cfg.get("mask_label") or {}).get("label")) or default)
    except Exception:
        return default


def _with_pair_prefix(prefix: str, name: str) -> str:
    """
    Pair-mode only: enforce stable ordering in downstream sorted() by adding 00_/01_.
    Keep original name visible.
    """
    if name.startswith("00_") or name.startswith("01_"):
        return name
    return f"{prefix}_{name}"


# ---------- commands ----------
@app.command(help="Run the Snakemake pipeline with given FASTA/BED inputs.")
def run(
    ref: Optional[str] = typer.Option(
        None,
        "--ref",
        help="Pair-mode: reference FASTA (use together with --query). If set, pair-mode is activated and order is forced to ref,query.",
    ),
    query: Optional[str] = typer.Option(
        None,
        "--query",
        help="Pair-mode: query FASTA (use together with --ref).",
    ),
    ref_bed_pair: Optional[str] = typer.Option(
        None,
        "--ref-bed-pair",
        help="Pair-mode: reference BED corresponding to --ref (optional).",
    ),
    query_bed_pair: Optional[str] = typer.Option(
        None,
        "--query-bed-pair",
        help="Pair-mode: query BED corresponding to --query (optional).",
    ),
    fasta: List[str] = typer.Option(
        None,
        "--fasta",
        "-f",
        help="Input FASTA files (same order as --bed when --bed is used). Accepts repeats, globs, comma-separated strings, and @list files.",
        hidden=True,
    ),
    bed: Optional[List[str]] = typer.Option(
        None,
        "--bed",
        "-b",
        help="Matching CENSAT BED files (same order as --fasta). Optional: if omitted, a whole-genome BED will be generated from FASTA (no masking).",
        hidden=True,
    ),
    ref_fasta: Optional[str] = typer.Option(
        None,
        "--ref-fasta",
        help="FASTA path to use as reference (it will be moved to the first position) in multi-mode.",
        hidden=True,
    ),
    ref_bed: Optional[str] = typer.Option(
        None,
        "--ref-bed",
        help="BED path to use as reference (must correspond to --ref-fasta when both are provided) in multi-mode.",
        hidden=True,
    ),
    names: Optional[List[str]] = typer.Option(
        None,
        "--names",
        "-S",
        help="Sample names (default: derived from FASTA basenames).",
        hidden=True,
    ),
    outdir: str = typer.Option(".", "--outdir", help="Work directory."),
    config_file: Optional[str] = typer.Option(
        None,
        "--config",
        "-c",
        help="Path to config.yaml. If omitted, a template will be copied into the work directory.",
    ),
    cores: int = typer.Option(8, "--cores", "-j", help="Number of cores for Snakemake."),
    dryrun: bool = typer.Option(False, "--dry-run", "-n", help="Snakemake dry-run."),
    keep_going: bool = typer.Option(True, "--keep-going/--no-keep-going", help="Keep going after errors."),
    use_conda: bool = typer.Option(False, "--use-conda/--no-conda", help="Use per-rule conda environments."),
    target: str = typer.Option("all", "--target", help="Snakemake target rule (default: all)."),
    alpha: int = typer.Option(
        1,
        "--alpha",
        help="alpha' threshold passed to score_filter_anchors_multi.py as --alpha-target.",
    ),
    beta: int = typer.Option(
        20,
        "--beta",
        help="ED threshold passed to score_filter_anchors_multi.py as --gap-beta.",
    ),
) -> None:
    if not check_requirements(["snakemake", "blastn", "makeblastdb", "samtools", "python"]):
        typer.secho(
            "[ERR] Missing required executables: snakemake, blastn, makeblastdb, samtools, or python.",
            fg="red",
            err=True,
        )
        raise typer.Exit(code=1)

    pair_mode = (ref is not None) or (query is not None) or (ref_bed_pair is not None) or (query_bed_pair is not None)

    if pair_mode:
        if not ref or not query:
            typer.secho("[ERR] Pair-mode requires both --ref and --query.", fg="red", err=True)
            raise typer.Exit(code=2)

        if fasta or bed or ref_fasta or ref_bed:
            typer.secho(
                "[ERR] Do not mix pair-mode (--ref/--query) with multi-mode (--fasta/--bed/--ref-fasta/--ref-bed).",
                fg="red",
                err=True,
            )
            raise typer.Exit(code=2)

        fasta_files = _expand_inputs([ref, query])
        if len(fasta_files) != 2:
            typer.secho(f"[ERR] Pair-mode expects exactly 2 FASTA paths after expansion. Got: {fasta_files}", fg="red", err=True)
            raise typer.Exit(code=2)

        fasta_files = [fasta_files[0], fasta_files[1]]

        bed_files: List[str] = []
        if ref_bed_pair or query_bed_pair:
            if not (ref_bed_pair and query_bed_pair):
                typer.secho("[ERR] Pair-mode: if you provide BEDs, provide both --ref-bed-pair and --query-bed-pair.", fg="red", err=True)
                raise typer.Exit(code=2)
            bed_files = _expand_inputs([ref_bed_pair, query_bed_pair])
            if len(bed_files) != 2:
                typer.secho(f"[ERR] Pair-mode expects exactly 2 BED paths after expansion. Got: {bed_files}", fg="red", err=True)
                raise typer.Exit(code=2)
            bed_files = [bed_files[0], bed_files[1]]

        if names:
            if len(names) != 2:
                typer.secho("[ERR] Pair-mode: --names must have exactly 2 entries (ref,query).", fg="red", err=True)
                raise typer.Exit(code=2)
            raw_names = list(names)
        else:
            raw_names = [_stem(fasta_files[0]), _stem(fasta_files[1])]

        sample_names = [
            _with_pair_prefix("00", raw_names[0]),
            _with_pair_prefix("01", raw_names[1]),
        ]

        print(f"[cyan]Pair-mode[/cyan]: ref={fasta_files[0]} query={fasta_files[1]}")
        print(f"[cyan]Pair-mode names[/cyan]: {sample_names[0]} , {sample_names[1]}")

    else:
        fasta_files = _expand_inputs(fasta)
        bed_files = _expand_inputs(bed) if bed else []

        if not fasta_files:
            typer.secho("[ERR] --fasta must be provided (or use --ref/--query pair-mode).", fg="red", err=True)
            raise typer.Exit(code=2)

        if bed_files and (len(fasta_files) != len(bed_files)):
            typer.secho(
                "[ERR] --fasta and --bed must be provided in equal numbers when --bed is used.",
                fg="red",
                err=True,
            )
            raise typer.Exit(code=2)

        if names:
            if len(names) != len(fasta_files):
                typer.secho("[ERR] The number of --names must match --fasta.", fg="red", err=True)
                raise typer.Exit(code=2)
            sample_names = list(names)
        else:
            sample_names = [_stem(f) for f in fasta_files]

        ref_idx = 0

        if ref_fasta:
            j = _find_index_by_path(fasta_files, ref_fasta)
            if j < 0:
                typer.secho(f"[ERR] --ref-fasta was not found in --fasta inputs: {ref_fasta}", fg="red", err=True)
                raise typer.Exit(code=2)
            ref_idx = j

        if ref_bed and bed_files:
            k = _find_index_by_path(bed_files, ref_bed)
            if k < 0:
                typer.secho(f"[ERR] --ref-bed was not found in --bed inputs: {ref_bed}", fg="red", err=True)
                raise typer.Exit(code=2)
            if ref_fasta and (k != ref_idx):
                typer.secho("[ERR] --ref-fasta and --ref-bed point to different indices in the input lists.", fg="red", err=True)
                raise typer.Exit(code=2)
            ref_idx = k

        if ref_idx != 0:
            fasta_files[0], fasta_files[ref_idx] = fasta_files[ref_idx], fasta_files[0]
            sample_names[0], sample_names[ref_idx] = sample_names[ref_idx], sample_names[0]
            if bed_files:
                bed_files[0], bed_files[ref_idx] = bed_files[ref_idx], bed_files[0]

    work = os.path.abspath(outdir)
    results_root = os.path.join(work, "results")
    data_root = os.path.join(results_root, "assemblies")
    result_data = os.path.join(results_root, "data")

    os.makedirs(results_root, exist_ok=True)
    os.makedirs(data_root, exist_ok=True)
    os.makedirs(result_data, exist_ok=True)

    if config_file:
        cfg_path = os.path.abspath(config_file)
        if not os.path.exists(cfg_path):
            typer.secho(f"[ERR] Config file does not exist: {cfg_path}", fg="red", err=True)
            raise typer.Exit(code=2)
        print(f"[green]Using user-specified config:[/green] {cfg_path}")
    else:
        tpl_res = files("vallescope.resources.templates").joinpath("config.yaml")
        cfg_path = os.path.join(work, "config.yaml")
        with as_file(tpl_res) as tpl_path:
            if not os.path.exists(cfg_path):
                shutil.copy2(tpl_path, cfg_path)
                print(f"[yellow]No config specified → copied template to[/yellow] {cfg_path}")
            else:
                print(f"[cyan]Using existing config:[/cyan] {cfg_path}")

    mask_label = _load_mask_label(cfg_path, default="HSat2")

    staged: List[Tuple[str, str, str, str]] = []
    for idx, f_in in enumerate(fasta_files):
        samp = sample_names[idx]
        base = _stem(f_in)

        dst_fa = os.path.join(result_data, samp, "filfa", f"{base}.t2t.filtered.fasta.gz")
        dst_bed = os.path.join(result_data, samp, "censat", f"{base}.censat.bed")

        _gzip_if_needed(f_in, dst_fa)

        if bed_files:
            b_in = bed_files[idx]
            _symlink_or_copy(b_in, dst_bed)
        else:
            _write_whole_genome_bed_from_fasta(dst_fa, dst_bed, label=mask_label)

        staged.append((samp, base, dst_fa, dst_bed))

    snake_res = files("vallescope.resources.workflows").joinpath("Snakefile")
    with as_file(snake_res) as snake_path:
        ok = run_snakemake(
            snakefile=str(snake_path),
            configfile=cfg_path,
            config={
                "alpha_target": alpha,
                "gap_beta": beta,
            },
            targets=[target],
            cores=cores,
            use_conda=use_conda,
            profile=None,
            keep_going=keep_going,
            dryrun=dryrun,
        )
    raise typer.Exit(code=0 if ok else 1)

if __name__ == "__main__":
    app()