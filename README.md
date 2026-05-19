# RecombTracer

**Version:** 0.0.1  
English | [中文](docs/README_zh.md)

RecombTracer is a toolkit for identifying recombination breakpoints and local ancestry in MAGIC (Multi-parent Advanced Generation Inter-Cross) populations. It combines PBWT-style chromosome painting with a Hidden Markov Model (HMM) to infer parental haplotype segments and recombination events.

## Overview

The pipeline consists of two core modules:

- **`recombtracer.core.recombiner`** — PBWT-style haplotype painting that assigns each progeny SNP to the most likely parental origin based on maximal matching segments.
- **`recombtracer.core.hmm`** — HMM refinement that smooths noisy PBWT assignments into contiguous ancestry blocks and calls recombination breakpoints with confidence scores.

A utility module is provided for converting standard VCF files into the NumPy-based inputs required by the analysis pipeline.

---

## Dependencies

- Python ^3.10
- numpy 2.2.6
- pandas 2.3.3
- scipy 1.15.3
- cyvcf2 0.33.0 (for VCF I/O)
- pyyaml >=6.0
- rich 15.0.0
- rich-gradient 0.3.12
- rich-argparse 1.8.0
- loguru 0.7.2

Install via pip (editable install from source):

```bash
pip install -e .
```

Or with Poetry:

```bash
poetry install
```

---

## VCF Conversion (`recombtracer convert-vcf`)

If your data starts from a VCF file (e.g. after variant calling and filtering), use `recombtracer convert-vcf` to extract haplotypes and convert them into the matrix format used by `MagicRecombiner`.

### Input requirements

- The VCF should be indexed (`.tbi` or `.csi`).
- All sites are expected to be homozygous (e.g. `0/0` or `1/1`). Heterozygous sites should be removed beforehand (e.g. `bcftools view -e GT="het"`).
- Missing genotypes (`./.`) are dropped by default to avoid analysis errors.

### Command-line usage

Convert one chromosome at a time to a compressed `.npz` archive:

```bash
recombtracer convert-vcf test/subset_remove_het.vcf.gz \
    --parents 1,2,3,4,5,6,7,8,9,10,11 \
    --progeny 98,99,100,144,145,146 \
    --chrom LG1 \
    --out-dir ./magic_input
```

| Flag | Description |
|------|-------------|
| `--parents` | **Required.** Comma-separated founder/parent sample names. |
| `--progeny` | Comma-separated progeny sample names. If omitted, all non-parent samples in the VCF are used automatically. |
| `--chrom` | **Required.** Chromosome/contig to extract (e.g. `LG1`). |
| `--out-dir` | Output directory for the `.npz` file (default: current directory). |
| `--keep-missing` | Keep SNPs that have missing genotypes in some samples (default: drop them). |

The output file (`{chrom}_magic.npz`) contains:

| Key | Shape / Type | Description |
|-----|--------------|-------------|
| `parent_haps` | `(n_parents, n_snps)` int8 | Parent haplotype matrix. |
| `parent_names` | list of str | Parent sample names. |
| `progeny_{name}` | `(1, n_snps)` int8 | One haplotype per progeny individual (homozygous sites store only one allele). |
| `progeny_names` | list of str | Progeny sample names. |
| `positions` | `(n_snps,)` int32 | 1-based genomic coordinates. |
| `chrom` | str | Chromosome name. |

### Python API

```python
from recombtracer import vcf_to_magic_inputs, load_chromosome_npz

# Convert directly from VCF
data = vcf_to_magic_inputs(
    vcf_path="test/subset_remove_het.vcf.gz",
    parent_samples=["1", "2", "3", "4"],
    progeny_samples=["5", "6", "7", "8"],
    chrom="LG0",
)

# Or load a previously saved .npz file
# data = load_chromosome_npz("magic_input/LG0_magic.npz")

print(data["parent_haps"].shape)      # (n_parents, n_snps)
print(data["progeny_haps"]["5"].shape) # (1, n_snps)
```

---

## Recombination Analysis CLI (`recombtracer run`)

After converting the VCF to `.npz`, run the full PBWT + HMM pipeline from the command line:

```bash
recombtracer run magic_input/LG1_magic.npz \
    --out-dir ./results \
    --smooth-window 5 \
    --save-raw
```

| Flag | Description |
|------|-------------|
| `npz` | **Required.** Input `.npz` file produced by `convert-vcf`. |
| `--out-dir` | Output directory for result CSVs (default: current directory). |
| `--min-match-len` | Minimum length of PBWT match segments in SNPs (default: 2). |
| `--smooth-window` | PBWT median-filter window size (default: 5). |
| `--min-segment-snps` | Minimum SNPs per ancestry segment (default: 5). |
| `--min-segment-bp` | Minimum bp length per ancestry segment (default: 1000). |
| `--min-posterior` | Minimum HMM posterior probability to keep a breakpoint (default: 0.8). |
| `--save-raw` | Also save raw PBWT results alongside HMM-refined results. |
| `--progeny` | Comma-separated list of progeny to analyze (default: all in `.npz`). |
| `--haplotype` | Only analyze a specific haplotype index (default: all). |

#### Why is `--min-match-len` set to 2 by default?

> **Note:** This default is provisional. We will perform extensive benchmarking across diverse MAGIC datasets in upcoming releases to determine the optimal default.

In MAGIC (Multi-parent Advanced Generation Inter-Cross) populations, a default of `2` is intentionally lenient for the following reasons:

1. **Small parental panels** — MAGIC populations typically have only 4–8 founders, so the risk of spurious matches is inherently low.
2. **Short recombination fragments** — Crossover breakpoints can be separated by just a few SNPs; a higher threshold would miss these short segments.
3. **Downstream filtering** — `extract_segments()` enforces `min_segment_snps` and `min_segment_bp`, and the HMM further refines calls, so being permissive at the PBWT stage is safe.

##### Typical values in other PBWT-based tools

| Application | Typical `min_len` | Rationale |
|-------------|-------------------|-----------|
| IBD detection (e.g., P-smoother) | 20–30 sites | Long exact matches are needed to rule out false-positive IBD segments. |
| Genetic genealogy (Syllable-PBWT) | 127–255 sites (~0.4–0.9 Mb) | Searching for distant relatives requires matches spanning large physical distances. |
| Recombination-rate inference | Tens to hundreds of sites | Balances runtime and match quality; overly short seeds explode the number of matches. |
| Phasing / Imputation | >10 sites | Needs dense, informative matches to condition state selection. |

##### Practical guidance

- **High-density WGS data**: Try `--min-match-len 5` or `10` to reduce noise and speed up analysis.
- **Array / capture data (lower density)**: Keep the default (`2` or `3`) to avoid missing short fragments.
- **When in doubt**: Run with the default (`2`), inspect the raw confidence distribution in the output, and tighten the threshold if needed.

### Output files

For each progeny × haplotype you will get:

| File | Description |
|------|-------------|
| `{name}_hap{h}_{chrom}_hmm_viterbi.csv` | Per-SNP Viterbi parent assignment and posterior probabilities. |
| `{name}_hap{h}_{chrom}_hmm_segments.csv` | Smoothed contiguous ancestry segments from HMM. |
| `{name}_hap{h}_{chrom}_hmm_recombinations.csv` | HMM-filtered recombination breakpoints with confidence scores. |
| `{name}_hap{h}_{chrom}_paint.csv` | *(with `--save-raw`)* Raw PBWT per-SNP parent assignments. |
| `{name}_hap{h}_{chrom}_segments_raw.csv` | *(with `--save-raw`)* Raw PBWT ancestry segments. |
| `{name}_hap{h}_{chrom}_recombinations_raw.csv` | *(with `--save-raw`)* Raw PBWT breakpoint calls. |
| `summary_{chrom}.csv` | Summary table: breakpoints and segment counts per individual. |

---

## Full Pipeline CLI (`recombtracer pipeline`)

Run the complete workflow — VCF conversion followed by PBWT + HMM analysis — in a single command:

```bash
recombtracer pipeline test/subset_remove_het.vcf.gz \
    --parents 1,2,3,4,5,6,7,8,9,10,11 \
    --chrom LG1 \
    --out-dir ./results \
    --save-raw
```

This is equivalent to running `convert-vcf` and then `run` with the same arguments. All algorithm and filtering parameters from `run` are also available here.

---

## Recombination Analysis (`recombtracer.core.recombiner`)

Once you have the haplotype matrices, run PBWT painting:

```python
from recombtracer import MagicRecombiner

recombiner = MagicRecombiner(
    parent_haps=data["parent_haps"],
    parent_names=data["parent_names"],
    positions=data["positions"],
    chrom=data["chrom"],
)

# Analyze one progeny
prog_name = data["progeny_names"][0]
prog_hap = data["progeny_haps"][prog_name]  # shape (1, n_snps)

paint_df = recombiner.paint_progeny(
    progeny_haps=prog_hap,
    progeny_name=prog_name,
    smooth_window=5,
)

# Extract contiguous ancestry segments
segments = recombiner.extract_segments(paint_df, min_segment_snps=5)

# Call recombination breakpoints
rec_df = recombiner.call_recombinations(segments)
```

### Output columns

`paint_df` (per-SNP parent assignment):
- `chrom`, `position`, `progeny`, `haplotype`
- `parent` — inferred parental origin
- `confidence` — PBWT match confidence

`segments` (contiguous blocks):
- `chrom`, `start_pos`, `end_pos`, `parent`, `haplotype`
- `score` — mean confidence across the segment

`rec_df` (breakpoint calls):
- `chrom`, `position`, `haplotype`
- `left_parent`, `right_parent`
- `confidence` — average score of the two flanking segments

---

## HMM Refinement (`recombtracer.core.hmm`)

To smooth noisy PBWT calls and obtain posterior probabilities, run the HMM:

```python
from recombtracer import run_hmm_refinement

viterbi_df, segments_df, rec_df = run_hmm_refinement(
    paint_df=paint_df[paint_df["haplotype"] == 0],
    parent_haps=data["parent_haps"],
    parent_names=data["parent_names"],
    progeny_hap=prog_hap[0],  # flatten to (n_snps,)
    progeny_name=prog_name,
    chrom=data["chrom"],
)
```

- `viterbi_df` — per-SNP most likely parent and posterior probabilities.
- `segments_df` — smoothed ancestry segments.
- `rec_df` — HMM-filtered recombination breakpoints.

---

## Complete Workflow Example

### Quick CLI pipeline (one-shot)

```bash
# Full pipeline: VCF → .npz → PBWT + HMM
recombtracer pipeline test/subset_remove_het.vcf.gz \
    --parents 1,2,3,4,5,6,7,8,9,10,11 \
    --chrom LG1 \
    --out-dir ./results \
    --save-raw
```

### Step-by-step CLI

```bash
# 1. Convert VCF to NumPy inputs
recombtracer convert-vcf test/subset_remove_het.vcf.gz \
    --parents 1,2,3,4,5,6,7,8,9,10,11 \
    --chrom LG1 \
    --out-dir ./magic_input

# 2. Run PBWT + HMM analysis
recombtracer run magic_input/LG1_magic.npz \
    --out-dir ./results \
    --save-raw
```

### Python API (for custom analysis)

If you need programmatic control over parameters or downstream analysis:

```python
from recombtracer import (
    load_chromosome_npz,
    MagicRecombiner,
    run_hmm_refinement,
)

data = load_chromosome_npz("magic_input/LG1_magic.npz")

recombiner = MagicRecombiner(
    parent_haps=data["parent_haps"],
    parent_names=data["parent_names"],
    positions=data["positions"],
    chrom=data["chrom"],
)

for prog_name in data["progeny_names"]:
    prog_hap = data["progeny_haps"][prog_name]

    paint_df = recombiner.paint_progeny(prog_hap, progeny_name=prog_name)
    segments = recombiner.extract_segments(paint_df)
    rec_raw = recombiner.call_recombinations(segments)

    viterbi_df, seg_hmm, rec_hmm = run_hmm_refinement(
        paint_df[paint_df["haplotype"] == 0],
        data["parent_haps"],
        data["parent_names"],
        prog_hap[0],
        progeny_name=prog_name,
        chrom=data["chrom"],
    )

    print(f"{prog_name}: {len(rec_raw)} raw breaks, {len(rec_hmm)} HMM breaks")
```

---

## Project Structure

```
.
├── recombtracer/              # Main package
│   ├── __init__.py            # Public API exports
│   ├── cli.py                 # Unified command-line interface
│   ├── core/                  # Core analysis modules
│   │   ├── recombiner.py      # PBWT chromosome painting
│   │   ├── hmm.py             # HMM smoothing & breakpoint calling
│   │   ├── vcf.py             # VCF → NumPy conversion utilities
│   │   ├── convert.py         # CLI handler for convert-vcf
│   │   ├── run.py             # CLI handler for run
│   │   └── pipeline.py        # CLI handler for pipeline
│   ├── config/                # Package configuration
│   │   ├── software.yaml      # Software metadata
│   │   └── default.yaml       # Default analysis parameters
│   └── utils/                 # Utility modules (logging, logo, etc.)
├── test/
│   └── subset_remove_het.vcf.gz   # Example VCF (het sites removed)
├── pyproject.toml             # Poetry-based build configuration
└── README.md
```

---

## License

MIT
