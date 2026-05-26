#!/usr/bin/env python3
"""
Generate synthetic test data for the RecombTracer report module.

Produces:
    test/report_test_data/
        ├── {prog}_hap0_chr1_hmm_recombinations.csv
        ├── {prog}_hap0_chr1_hmm_segments.csv
        ├── {prog}_hap1_chr1_hmm_recombinations.csv
        ├── {prog}_hap1_chr1_hmm_segments.csv
        └── summary_chr1.csv

Usage:
    python test/generate_report_test_data.py
"""

import os
import random
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "report_test_data")
CHROM = "chr1"
CHROM_LENGTH = 50_000_000   # 50 Mb
N_PROGENY = 8
N_RECOMBINATIONS = 5        # expected breakpoints per haplotype
PARENTS = ["ParentA", "ParentB"]
SEED = 42

random.seed(SEED)
np.random.seed(SEED)

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def generate_breakpoints(chrom_length: int, n_expected: int):
    """Return sorted list of breakpoint positions."""
    n = max(0, int(np.random.poisson(n_expected)))
    if n == 0:
        return []
    positions = np.sort(np.random.choice(range(1_000_000, chrom_length - 1_000_000), size=n, replace=False))
    return positions.tolist()


def build_segments(breakpoints, chrom_length, parents):
    """
    Given sorted breakpoints, build contiguous segments alternating parents.
    Returns list of dicts with start, end, parent.
    """
    segments = []
    start = 0
    parent_idx = 0
    for bp in breakpoints:
        segments.append({
            "start": start,
            "end": bp,
            "parent": parents[parent_idx % len(parents)],
        })
        start = bp
        parent_idx += 1
    # Final segment
    segments.append({
        "start": start,
        "end": chrom_length,
        "parent": parents[parent_idx % len(parents)],
    })
    return segments


# ---------------------------------------------------------------------------
# Main generation
# ---------------------------------------------------------------------------
summary_rows = []

for prog_i in range(1, N_PROGENY + 1):
    prog_name = f"Progeny{prog_i:02d}"

    for hap in (0, 1):
        bps = generate_breakpoints(CHROM_LENGTH, N_RECOMBINATIONS)
        segs = build_segments(bps, CHROM_LENGTH, PARENTS)

        # --- recombinations CSV ---
        rec_rows = []
        for pos in bps:
            # Infer left/right parent from segments
            for seg in segs:
                if seg["start"] <= pos < seg["end"]:
                    right_parent = seg["parent"]
                    break
            else:
                right_parent = "Unknown"
            # Left parent is the previous segment's parent
            left_parent = "Unknown"
            for idx, seg in enumerate(segs):
                if seg["start"] <= pos < seg["end"] and idx > 0:
                    left_parent = segs[idx - 1]["parent"]
                    break

            rec_rows.append({
                "chrom": CHROM,
                "position": int(pos),
                "left_parent": left_parent,
                "right_parent": right_parent,
                "confidence": round(random.uniform(0.75, 0.99), 4),
            })

        rec_df = pd.DataFrame(rec_rows)
        rec_path = os.path.join(
            OUTPUT_DIR,
            f"{prog_name}_hap{hap}_{CHROM}_hmm_recombinations.csv",
        )
        rec_df.to_csv(rec_path, index=False)

        # --- segments CSV ---
        seg_rows = []
        for seg in segs:
            n_snps = int((seg["end"] - seg["start"]) / 1000) + random.randint(5, 50)
            seg_rows.append({
                "chrom": CHROM,
                "start": seg["start"],
                "end": seg["end"],
                "parent": seg["parent"],
                "n_snps": n_snps,
                "mean_posterior": round(random.uniform(0.85, 0.999), 4),
            })

        seg_df = pd.DataFrame(seg_rows)
        seg_path = os.path.join(
            OUTPUT_DIR,
            f"{prog_name}_hap{hap}_{CHROM}_hmm_segments.csv",
        )
        seg_df.to_csv(seg_path, index=False)

        # --- summary row ---
        summary_rows.append({
            "progeny": prog_name,
            "chrom": CHROM,
            "haplotype": hap,
            "hmm_breakpoints": len(bps),
        })

# --- summary CSV ---
summary_df = pd.DataFrame(summary_rows)
summary_path = os.path.join(OUTPUT_DIR, f"summary_{CHROM}.csv")
summary_df.to_csv(summary_path, index=False)

print(f"Test data written to: {OUTPUT_DIR}")
print(f"  Progeny      : {N_PROGENY}")
print(f"  Chromosome   : {CHROM} ({CHROM_LENGTH:,} bp)")
print(f"  Files        : {len(list(os.listdir(OUTPUT_DIR)))}")
