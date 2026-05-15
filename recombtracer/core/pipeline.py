"""
High-level analysis pipeline that chains ``convert-vcf`` and ``run``.

This module provides the ``pipeline`` CLI workflow: convert a VCF to .npz
and immediately run PBWT + HMM recombination analysis on the result.
"""

import os
import types

from .convert import handle_convert_vcf
from .run import handle_run


def handle_pipeline(args):
    """
    Run the full pipeline: convert VCF → .npz → PBWT + HMM analysis.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments containing both convert-vcf and run options.
    """
    from ..utils import logger_generator

    logger, _ = logger_generator(args.out_dir, log_level=args.log_level)

    # ------------------------------------------------------------------
    # Step 1: Convert VCF to .npz
    # ------------------------------------------------------------------
    logger.info("=" * 50)
    logger.info("STEP 1/2: Converting VCF to .npz")
    logger.info("=" * 50)

    convert_args = types.SimpleNamespace(
        vcf=args.vcf,
        parents=args.parents,
        chrom=args.chrom,
        progeny=args.progeny,
        out_dir=args.out_dir,
        keep_missing=args.keep_missing,
        log_level=args.log_level,
    )
    handle_convert_vcf(convert_args)

    npz_path = os.path.join(args.out_dir, f"{args.chrom}_magic.npz")

    # ------------------------------------------------------------------
    # Step 2: Run PBWT + HMM analysis
    # ------------------------------------------------------------------
    logger.info("=" * 50)
    logger.info("STEP 2/2: Running PBWT + HMM analysis")
    logger.info("=" * 50)

    run_args = types.SimpleNamespace(
        npz=npz_path,
        out_dir=args.out_dir,
        smooth_window=args.smooth_window,
        min_segment_snps=args.min_segment_snps,
        min_segment_bp=args.min_segment_bp,
        min_posterior=args.min_posterior,
        save_raw=args.save_raw,
        haplotype=args.haplotype,
        progeny=args.progeny,
        log_level=args.log_level,
    )
    handle_run(run_args)
