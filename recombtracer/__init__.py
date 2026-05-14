"""
RecombTracer — Recombination breakpoint detection in MAGIC populations.

Provides PBWT-style chromosome painting, HMM refinement, and VCF I/O utilities.
"""

__version__ = "0.1.0"

# Core recombination painting
from .recombiner import (
    MagicRecombiner,
    MatchSegment,
    AncestrySegment,
)

# HMM refinement
from .hmm import (
    HMMParams,
    MagicHMM,
    run_hmm_refinement,
)

# VCF I/O utilities
from .vcf import (
    list_vcf_samples,
    extract_chromosome,
    vcf_to_magic_inputs,
    load_chromosome_npz,
    save_chromosome_npz,
)

__all__ = [
    "MagicRecombiner",
    "MatchSegment",
    "AncestrySegment",
    "HMMParams",
    "MagicHMM",
    "run_hmm_refinement",
    "list_vcf_samples",
    "extract_chromosome",
    "vcf_to_magic_inputs",
    "load_chromosome_npz",
    "save_chromosome_npz",
]
