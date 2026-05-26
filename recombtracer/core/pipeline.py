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

    logger.info("=" * 60)
    logger.info(f"RecombTracer Pipeline: {args.chrom}")
    logger.info(f"Input VCF : {args.vcf}")
    logger.info(f"Output Dir: {args.out_dir}")
    logger.info("=" * 60)

    # ------------------------------------------------------------------
    # Step 1: Convert VCF to .npz
    # ------------------------------------------------------------------
    logger.info("[STEP 1/2] Converting VCF to .npz format")
    logger.debug(f"Conversion params: parents={args.parents}, progeny={args.progeny}, keep_missing={args.keep_missing}")

    convert_args = types.SimpleNamespace(
        vcf=args.vcf,
        parents=args.parents,
        chrom=args.chrom,
        progeny=args.progeny,
        out_dir=args.out_dir,
        keep_missing=args.keep_missing,
        log_level=args.log_level,
    )
    handle_convert_vcf(convert_args, logger=logger)

    npz_path = os.path.join(args.out_dir, f"{args.chrom}_magic.npz")
    if not os.path.exists(npz_path):
        logger.error(f"Failed to create .npz file at {npz_path}")
        return

    # ------------------------------------------------------------------
    # Step 2: Run PBWT + HMM analysis
    # ------------------------------------------------------------------
    logger.info("-" * 60)
    logger.info("[STEP 2/2] Running PBWT + HMM recombination analysis")
    logger.debug(f"Analysis params: window={args.smooth_window}, min_snps={args.min_segment_snps}, min_bp={args.min_segment_bp}, posterior={args.min_posterior}")

    run_args = types.SimpleNamespace(
        npz=npz_path,
        out_dir=args.out_dir,
        smooth_window=args.smooth_window,
        min_match_len=args.min_match_len,
        min_segment_snps=args.min_segment_snps,
        min_segment_bp=args.min_segment_bp,
        min_posterior=args.min_posterior,
        save_raw=args.save_raw,
        haplotype=args.haplotype,
        progeny=args.progeny,
        log_level=args.log_level,
    )
    handle_run(run_args, logger=logger)

    logger.info("=" * 60)
    logger.info("Pipeline completed successfully.")
    logger.info("=" * 60)
