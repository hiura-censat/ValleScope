# src/vallescope/runner.py
from __future__ import annotations
import subprocess
import shutil 
from typing import List, Optional, Dict
from rich import print

def _exists(cmd: str) -> bool:
    return shutil.which(cmd) is not None

def check_requirements(cmds):
    ok = True
    for c in cmds:
        if not _exists(c):
            print(f"[red]missing[/red]: {c}")
            ok = False
    return ok

def run_snakemake(
    snakefile: str,
    configfile: str,
    targets: List[str],
    cores: int = 8,
    use_conda: bool = True,
    profile: Optional[str] = None,
    keep_going: bool = True,
    dryrun: bool = False,
    config: Optional[Dict[str, object]] = None,
) -> bool:
    cmd = ["snakemake", "-s", snakefile, "-c", str(cores), "--printshellcmds"]
    if configfile:
        cmd += ["--configfile", configfile]
    if config:
        cmd += ["--config"]
        for k, v in config.items():
            cmd.append(f"{k}={v}")
    if use_conda:
        cmd += ["--use-conda", "--conda-frontend", "mamba"]
    if keep_going:
        cmd += ["--keep-going"]
    if dryrun:
        cmd += ["-n"]
    if profile:
        cmd += ["--profile", profile]
    if targets:
        cmd += targets
    print(" ".join(cmd))
    return subprocess.call(cmd) == 0

def rule_graph(snakefile: str, configfile: str, out_png: str):
    cmd = ["snakemake", "-s", snakefile, "--rulegraph", "--configfile", configfile]
    dot = subprocess.check_output(cmd, text=True)
    p = subprocess.Popen(["dot", "-Tpng", "-o", out_png], stdin=subprocess.PIPE, text=True)
    p.communicate(input=dot)
    if p.returncode != 0:
        raise RuntimeError("Graph rendering failed")

def unlock_snakemake(workdir: str):
    subprocess.check_call(["snakemake", "--unlock", "--directory", workdir])

