"""
RecombTracer core package.

All core logic (recombination painting, HMM refinement, and VCF processing)
can be imported directly from this package::

    from recombtracer.core import MagicRecombiner, MagicHMM
    from recombtracer.core import vcf_to_magic_inputs, run_hmm_refinement
"""

from .hmm import HMMParams, MagicHMM, run_hmm_refinement
from .recombiner import AncestrySegment, MagicRecombiner, MatchSegment
from .pipeline import handle_convert_vcf, handle_run
from .vcf import (
    extract_homozygous_chromosome,
    list_vcf_samples,
    load_chromosome_npz,
    save_chromosome_npz,
    vcf_to_magic_inputs,
)

__all__ = [
    # High-level pipelines
    "handle_convert_vcf",
    "handle_run",
    # HMM refinement
    "HMMParams",
    "MagicHMM",
    "run_hmm_refinement",
    # Core recombination painting
    "MagicRecombiner",
    "MatchSegment",
    "AncestrySegment",
    # VCF I/O utilities
    "list_vcf_samples",
    "extract_homozygous_chromosome",
    "vcf_to_magic_inputs",
    "load_chromosome_npz",
    "save_chromosome_npz",
]
