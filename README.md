# RecombTracer

RecombTracer is a toolkit for identifying recombination breakpoints and local ancestry in MAGIC (Multi-parent Advanced Generation Inter-Cross) populations. It combines PBWT-style chromosome painting with a Hidden Markov Model (HMM) to infer parental haplotype segments and recombination events.

## Overview

The pipeline consists of two core modules:

- **`recombtracer.recombiner`** — PBWT-style haplotype painting that assigns each progeny SNP to the most likely parental origin based on maximal matching segments.
- **`recombtracer.hmm`** — HMM refinement that smooths noisy PBWT assignments into contiguous ancestry blocks and calls recombination breakpoints with confidence scores.

A utility module is provided for converting standard VCF files into the NumPy-based inputs required by the analysis pipeline.

---

## Dependencies

- Python >= 3.8
- numpy
- pandas
- scipy
- cyvcf2 (for VCF I/O)

Install via pip:

```bash
pip install numpy pandas scipy cyvcf2
```

---

## VCF Conversion (`recombtracer.vcf`)

If your data starts from a VCF file (e.g. after variant calling and filtering), use `recombtracer.vcf` to extract haplotypes and convert them into the matrix format used by `MagicRecombiner`.

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
from src.vcf_utils import vcf_to_magic_inputs, load_chromosome_npz

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
| `--smooth-window` | PBWT median-filter window size (default: 5). |
| `--min-segment-snps` | Minimum SNPs per ancestry segment (default: 5). |
| `--min-segment-bp` | Minimum bp length per ancestry segment (default: 1000). |
| `--min-posterior` | Minimum HMM posterior probability to keep a breakpoint (default: 0.8). |
| `--save-raw` | Also save raw PBWT results alongside HMM-refined results. |
| `--progeny` | Comma-separated list of progeny to analyze (default: all in `.npz`). |
| `--haplotype` | Only analyze a specific haplotype index (default: all). |

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

## Recombination Analysis (`recombtracer.recombiner`)

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

## HMM Refinement (`recombtracer.hmm`)

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

Run the built-in demo to see the full pipeline on synthetic data:

```bash
recombtracer demo-hmm
recombtracer demo-recombiner
```

---

## Complete Workflow Example

### Quick CLI pipeline

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
├── src/
│   ├── __init__.py           # Public API exports
│   ├── cli.py                # Unified command-line interface
│   ├── recombiner.py         # PBWT chromosome painting
│   ├── hmm.py                # HMM smoothing & breakpoint calling
│   └── vcf.py                # VCF → NumPy conversion utilities
├── test/
│   └── subset_remove_het.vcf.gz   # Example VCF (het sites removed)
├── pbwt/                     # PBWT C implementation
└── README.md
```

---

## License

MIT
