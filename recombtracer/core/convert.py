"""
VCF-to-NPZ conversion pipeline for the ``convert-vcf`` CLI command.
"""

import os
import sys

from .vcf import list_vcf_samples, save_chromosome_npz, vcf_to_magic_inputs


def handle_convert_vcf(args, logger=None):
    """
    Pipeline for the ``convert-vcf`` command.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments.
    logger : loguru.Logger, optional
        External logger to use. If None, a new one is created.
    """
    if logger is None:
        from ..utils import logger_generator
        logger, _ = logger_generator(args.out_dir, log_level=args.log_level)

    logger.info(f"VCF Conversion: {args.chrom}")
    logger.debug(f"Input VCF: {args.vcf}")
    logger.debug(f"Output Dir: {args.out_dir}")

    parent_samples = [s.strip() for s in args.parents.split(",")]
    all_samples = list_vcf_samples(args.vcf)

    if args.progeny:
        progeny_samples = [s.strip() for s in args.progeny.split(",")]
    else:
        logger.info("  No progeny specified, using all non-parent samples from VCF")
        progeny_samples = [s for s in all_samples if s not in parent_samples]

    missing_parents = set(parent_samples) - set(all_samples)
    missing_progeny = set(progeny_samples) - set(all_samples)
    if missing_parents:
        logger.error(f"Parent samples not in VCF: {missing_parents}")
        sys.exit(f"Error: parent samples not in VCF: {missing_parents}")
    if missing_progeny:
        logger.error(f"Progeny samples not in VCF: {missing_progeny}")
        sys.exit(f"Error: progeny samples not in VCF: {missing_progeny}")

    logger.info(f"  Chromosome: {args.chrom}")
    logger.info(f"  Parents   ({len(parent_samples)}): {', '.join(parent_samples)}")
    logger.info(f"  Progeny   ({len(progeny_samples)}): {len(progeny_samples)} total")
    logger.debug(f"  Progeny list: {progeny_samples}")

    logger.info("  Extracting genotypes and filtering SNPs...")
    magic_data = vcf_to_magic_inputs(
        args.vcf,
        parent_samples,
        progeny_samples,
        args.chrom,
        skip_missing=not args.keep_missing,
    )

    out_name = f"{args.chrom}_magic.npz"
    out_path = os.path.join(args.out_dir, out_name)
    os.makedirs(args.out_dir, exist_ok=True)

    logger.info(f"  Saving processed data to: {out_path}")
    save_chromosome_npz(
        out_path,
        magic_data["parent_haps"],
        magic_data["parent_names"],
        magic_data["progeny_haps"],
        magic_data["progeny_names"],
        magic_data["positions"],
        magic_data["chrom"],
    )

    logger.info("  Conversion Summary:")
    logger.info(f"    SNPs    : {len(magic_data['positions']):,}")
    logger.info(f"    Parents : {magic_data['parent_haps'].shape[0]} haplotypes")
    logger.info(f"    Progeny : {len(magic_data['progeny_haps'])} individuals")
    logger.debug(f"    Genotype matrix shape: {magic_data['parent_haps'].shape}")
