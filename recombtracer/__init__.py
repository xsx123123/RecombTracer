"""
RecombTracer — Recombination breakpoint detection in MAGIC populations.

Provides PBWT-style chromosome painting, HMM refinement, and VCF I/O utilities.
"""

import yaml
from importlib import resources

# Load software metadata directly from package-internal YAML.
# This avoids importing utils (which may pull heavy deps like loguru/omegaconf)
# just to read the version number.
_software_conf = yaml.safe_load(
    (resources.files("recombtracer.config") / "software.yaml")
    .read_text(encoding="utf-8")
)
SOFTWARE_INFO = _software_conf.get("software", {})

__version__ = SOFTWARE_INFO.get("version", "unknown")
__app_name__ = SOFTWARE_INFO.get("app_name", "RecombTracer")
__author__ = SOFTWARE_INFO.get("author", "unknown")
__email__ = SOFTWARE_INFO.get("email", "")
__description__ = SOFTWARE_INFO.get("description", "")

# Core recombination painting
from .core.recombiner import (
    MagicRecombiner,
    MatchSegment,
    AncestrySegment,
)

# HMM refinement
from .core.hmm import (
    HMMParams,
    MagicHMM,
    run_hmm_refinement,
)

# VCF I/O utilities
from .core.vcf import (
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
