# workflows/3_anchor_and_plot.smk
import os, re, hashlib, shutil, gzip, glob, sys
from snakemake.shell import shell
from pathlib import Path
from importlib.resources import files as _pkg_files
from vallescope.utils import split_name, public_pairs, sha256sum, find_meta

_HERE = Path(__file__).resolve()  # .../resources/workflows/rules/3_anchor_and_plot.smk

try:
    _PKG_TOOLS = Path(_pkg_files("vallescope.resources").joinpath("tools")).resolve()
except Exception:
    _PKG_TOOLS = (_HERE.parents[1]).parent / "tools"

_REPO_TOOLS = Path.cwd() / "workflows" / "tools"

def _tool(name: str) -> str:
    for base in (_PKG_TOOLS, _REPO_TOOLS):
        cand = base / name
        if cand.exists():
            return str(cand)
    return str(_PKG_TOOLS / name)


DATA_ROOT   = config["data_root"]
THREADS     = int(config.get("threads", 8))

MASK_LABEL     = config.get("mask_label", {}).get("label", "HSat2")
MERYL_K        = 31
MAX_KMER_COUNT = 5
MIN_UKF        = 0.90
METRICS_K      = 24

PROJECT_ROOT    = os.path.join(DATA_ROOT, "..")
RESULT_ROOT     = os.path.join(PROJECT_ROOT, "data") 
ANCHOR_ROOT     = os.path.join(PROJECT_ROOT, "anchor") 
SAMPLES, BASES = glob_wildcards(f"{RESULT_ROOT}" + "/{sample}/censat/{base}.censat.bed")
PAIRS = list(zip(SAMPLES, BASES))  

def filfa_path(s, b): return os.path.join(RESULT_ROOT, s, "filfa",  f"{b}.t2t.filtered.fasta.gz")
def censat_path(s,b): return os.path.join(RESULT_ROOT, s, "censat", f"{b}.censat.bed")
def masked_path(s,b): return os.path.join(RESULT_ROOT, s, "masked", f"{b}.masked.fasta")
def local_fa_path(s, b): return os.path.join(ANCHOR_ROOT, s, b, "anchoring_data", "sat_fasta", f"{b}.masked.fasta")
def genid_dir(s, b): return os.path.join(ANCHOR_ROOT, s, b, "anchoring_data", "genmap_index")
def genout_dir(s, b): return os.path.join(ANCHOR_ROOT, s, b, "anchoring_data", "genmap_output")
def genbg_path(s, b): return os.path.join(ANCHOR_ROOT, s, b, "anchoring_data", "genmap_output", "genout.bedgraph")

def blast6_path(s, b): return os.path.join(ANCHOR_ROOT, s, b, "anchoring_data", "blast_output", "blast_self.outfmt6.txt")
def ge100_path(s, b): return os.path.join(ANCHOR_ROOT, s, b, "anchoring_data", "blast_output", f"{b}.blast_out.len_ge100.txt")
def lt100_path(s, b): return os.path.join(ANCHOR_ROOT, s, b, "anchoring_data", "blast_output", f"{b}.blast_out.len_lt100.txt")
def anchor_path(s, b): return os.path.join(ANCHOR_ROOT, s, b, "anchoring_data", "anchors", f"{b}.masked.len100_anchors.bed")
def ancplot_path(s, b): return os.path.join(ANCHOR_ROOT, s, b, "anchoring_data", "anchors", f"{b}.masked.len100.region.png")
def ancprefix_path(s, b): return os.path.join(ANCHOR_ROOT, s, b, "anchoring_data", "anchors", f"{b}.masked.len100")
def bigwig_path(s, b): return os.path.join(ANCHOR_ROOT, s, b, "anchoring_data", "anchors", f"{b}.masked.len100_smooth.bw")
def chrsize_path(s, b): return os.path.join(ANCHOR_ROOT, s, b, "anchoring_data", "anchors", f"{b}.masked.chrom.sizes")
def bedgsort_path(s, b): return os.path.join(ANCHOR_ROOT, s, b, "anchoring_data", "anchors", f"{b}.masked.len100_coverage.sorted.bedGraph")

def uqbigwig_path(s, b): return os.path.join(ANCHOR_ROOT, s, b, "anchor_scoring", f"{b}.masked.anchors.limited.UniqueKmerFrac.bw")
def uqanctsv_path(s, b): return os.path.join(ANCHOR_ROOT, s, b, "anchor_scoring", f"{b}.masked.anchors.limited.UniqueKmerFrac.tsv")
def ancfa_path(s, b): return os.path.join(ANCHOR_ROOT, s, b, "anchor_scoring", f"{b}.masked.anchors.fa")
def filids_path(s, b): return os.path.join(ANCHOR_ROOT, s, b, "anchor_scoring", f"{b}.masked.limited.filtered.ids.txt")
def ancsprefix_path(s, b): return os.path.join(ANCHOR_ROOT, s, b, "anchor_scoring", f"{b}.masked")
def afilprefix_path(s, b): return os.path.join(ANCHOR_ROOT, s, b, "anchor_scoring", f"{b}.filtered.anchors")
def afilfa_path(s, b): return os.path.join(ANCHOR_ROOT, s, b, "anchor_scoring", f"{b}.filtered.anchors.fa")
def afiltsv_path(s, b): return os.path.join(ANCHOR_ROOT, s, b, "anchor_scoring", f"{b}.filtered.anchors.tsv")
def afilbed_path(s, b): return os.path.join(ANCHOR_ROOT, s, b, "anchor_scoring", f"{b}.filtered.anchors.bed")

MSK_TARGETS     = [masked_path(s, b) for s, b in PAIRS]
GBG_TARGETS     = [genbg_path(s, b) for s, b in PAIRS]
ANC_TARGETS     = [ancplot_path(s, b) for s, b in PAIRS]

SUQ_TARGETS     = [uqanctsv_path(s, b) for s, b in PAIRS]
FIDS_TARGETS    = [afilfa_path(s, b) for s, b in PAIRS]

rule all3:
    input: 
        MSK_TARGETS,
        GBG_TARGETS,
        ANC_TARGETS,

        SUQ_TARGETS,
        FIDS_TARGETS

rule stage_masked_ok:
    input:  MSK_TARGETS
    output: ".stages/masked.ok"
    shell:
        r"""
        mkdir -p .stages
        touch {output}
        """

rule stage_mkgen_ok:
    input:  GBG_TARGETS
    output: ".stages/mkgen.ok"
    shell: "touch {output}"

rule stage_anchor_ok:
    input:  ANC_TARGETS
    output: ".stages/anchor.ok"
    shell: "touch {output}"

rule stage_satuniq_ok:
    input:  SUQ_TARGETS
    output: ".stages/satuniq.ok"
    shell: "touch {output}"

rule stage_afilfa_ok:
    input:  FIDS_TARGETS
    output: ".stages/afilfa.ok"
    shell: "touch {output}"

rule mask_by_bed_labels:
    input:
      filfa     = filfa_path("{sample}", "{base}"),
      censat    = censat_path("{sample}", "{base}")
    output:
      masked    = masked_path("{sample}", "{base}") 
    params:
      script    = _tool("mask_by_bed_labels.py"),
      label     = MASK_LABEL
    shell:
        r"""
        set -euo pipefail
        python {params.script} \
          --fasta {input.filfa} \
          --bed {input.censat} \
          --out {output.masked} \
          --labels {params.label} \
          --column 4 \
          --mask-char n \
          --merge-distance 0 \
          --padding 0 \
          --invert
        samtools faidx {output.masked}
        """

# -----------------------------------------------------------------------------
# 1) run genmap
# -----------------------------------------------------------------------------
rule run_genmap:
    input:
        barrier = ".stages/masked.ok",
        masked  = masked_path("{sample}", "{base}")
    output:
        gen_bg   = genbg_path("{sample}", "{base}"),
    params:
        genid_dir    = genid_dir("{sample}", "{base}"),
        genout_dir   = genout_dir("{sample}", "{base}"),
    threads: 1
    shell:
        r"""
        set -euo pipefail
        if [[ -d "{params.genid_dir}" ]]; then
          rm -rf "{params.genid_dir}"
        fi
        mkdir -p "{params.genout_dir}"
        genmap index -F {input.masked} -I {params.genid_dir}
        genmap map -K 6 -E 0 -I {params.genid_dir} -O {params.genout_dir}/genout -t -w -bg -fl
        test -s "{output.gen_bg}"
        """

# -----------------------------------------------------------------------------
# 2) make anchors
# -----------------------------------------------------------------------------
rule make_anchors_from_valleys:
    input:
        barrier     = ".stages/mkgen.ok",
        fa          = masked_path("{sample}", "{base}"),
        gen_bg      = genbg_path("{sample}", "{base}"),
        censat      = censat_path("{sample}", "{base}"),
    output:
        anchors_bed = anchor_path("{sample}", "{base}"),
        plot_png    = ancplot_path("{sample}", "{base}"),
        bw          = bigwig_path("{sample}", "{base}"),
    params:
        script_anchors  = _tool("make_anchors_from_genmap.py"),
        outprefix   = ancprefix_path("{sample}", "{base}"),
        label = MASK_LABEL 
    threads: 2
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {params.outprefix})"
        python3 "{params.script_anchors}" \
          --fasta "{input.fa}" \
          --bedgraph "{input.gen_bg}" \
          --outprefix "{params.outprefix}" \
          --smooth 50 \
          --valley-pct 30 \
          --valley-radius 50 \
          --auto-span \
          --auto-span-factor 1.2 \
          --plot-png "{output.plot_png}" \
          --plot-range auto \
          --skip-bed "{input.censat}" \
          --skip-labels "{params.label}" \
          --focus-labels \
          --skip-col 4 \
          --skip-merge-distance 0 \
          --skip-pad 0 \
          --export-highcov-bed \
          --no-coverage-tsv \
          --write-bigwig \
          --bigwig-kind smooth
        """

rule run_satuniq:
    input:
        barrier     = ".stages/anchor.ok",
        fa          = masked_path("{sample}", "{base}"),
        anchors_bed = anchor_path("{sample}", "{base}"),
    output:
        tsv        = uqanctsv_path("{sample}", "{base}"),
        anchors_fa = ancfa_path("{sample}", "{base}"),
    params:
        script_satuniq = _tool("detect_satuniq.py"),
        k              = MERYL_K,
        maxcnt         = MAX_KMER_COUNT,
        outprefix= ancsprefix_path("{sample}", "{base}"),
    threads: THREADS
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {params.outprefix})"
        python "{params.script_satuniq}" \
          --genome-fa   "{input.fa}" \
          --anchors-bed "{input.anchors_bed}" \
          --k           {params.k} \
          --threads     {threads} \
          --out-prefix  "{params.outprefix}" \
          --max-kmer-count {params.maxcnt} 
        """

rule filter_one_ids:
  input:
    barrier     = ".stages/satuniq.ok",
    bed         = anchor_path("{sample}", "{base}"),
    fasta       = ancfa_path("{sample}", "{base}"),
    tsv         = uqanctsv_path("{sample}", "{base}"),
  output:
    out_fa  = afilfa_path("{sample}", "{base}"),
    out_tsv = afiltsv_path("{sample}", "{base}"),
    out_bed = afilbed_path("{sample}", "{base}")
  params:
    script = _tool("filter_ukf_and_emit_bed_fa.py"),
    prefix = afilprefix_path("{sample}", "{base}"),
  threads: 1
  shell:
    r"""
    set -euo pipefail
    python {params.script} \
      --tsv "{input.tsv}" \
      --bed "{input.bed}" \
      --fa "{input.fasta}" \
      --out-prefix "{params.prefix}"
    """
