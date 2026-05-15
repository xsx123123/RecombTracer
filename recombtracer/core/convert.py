"""
VCF-to-NPZ conversion pipeline for the ``convert-vcf`` CLI command.
"""

import os
import sys

from .vcf import list_vcf_samples, save_chromosome_npz, vcf_to_magic_inputs


def handle_convert_vcf(args):
    """
    Pipeline for the ``convert-vcf`` command.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments.
    """
    from ..utils import logger_generator

    logger, _ = logger_generator(args.out_dir, log_level=args.log_level)

    parent_samples = [s.strip() for s in args.parents.split(",")]
    all_samples = list_vcf_samples(args.vcf)

    if args.progeny:
        progeny_samples = [s.strip() for s in args.progeny.split(",")]
    else:
        progeny_samples = [s for s in all_samples if s not in parent_samples]

    missing_parents = set(parent_samples) - set(all_samples)
    missing_progeny = set(progeny_samples) - set(all_samples)
    if missing_parents:
        sys.exit(f"Error: parent samples not in VCF: {missing_parents}")
    if missing_progeny:
        sys.exit(f"Error: progeny samples not in VCF: {missing_progeny}")

    logger.info(f"Chromosome: {args.chrom}")
    logger.info(f"Parents   ({len(parent_samples)}): {', '.join(parent_samples)}")
    logger.info(f"Progeny   ({len(progeny_samples)}): {', '.join(progeny_samples)}")

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

    save_chromosome_npz(
        out_path,
        magic_data["parent_haps"],
        magic_data["parent_names"],
        magic_data["progeny_haps"],
        magic_data["progeny_names"],
        magic_data["positions"],
        magic_data["chrom"],
    )

    logger.info(f"  SNPs    : {len(magic_data['positions']):,}")
    logger.info(f"  Parents : {magic_data['parent_haps'].shape}")
    logger.info(f"  Progeny : {len(magic_data['progeny_haps'])} individuals")
