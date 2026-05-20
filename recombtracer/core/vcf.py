#!/usr/bin/env python3
"""
VCF utilities for MAGIC Recombiner.

Convert VCF files (especially homozygous-only VCFs) into the NumPy arrays
required by MagicRecombiner and MagicHMM.

Dependencies: cyvcf2, numpy
"""

import os
from typing import List, Tuple, Dict, Optional

import numpy as np
from cyvcf2 import VCF

from ..utils.log_utils import logger


def list_vcf_samples(vcf_path: str) -> List[str]:
    """
    Return all sample names in the VCF.
    """
    logger.info(f"Reading sample list from VCF: {vcf_path}")
    vcf = VCF(vcf_path)
    samples = list(vcf.samples)
    logger.info(f"Found {len(samples)} samples in VCF")
    logger.debug(f"Sample names (first 10): {samples[:10]}")
    return samples


def extract_homozygous_chromosome(
    vcf_path: str,
    samples: List[str],
    chrom: str,
    skip_missing: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract haplotypes for a single chromosome from a VCF.

    Since the VCF has had heterozygous sites removed, every genotype is
    either 0/0 or 1/1.  We therefore store only one allele per sample
    (ploidy=1) to save memory and match the haplotype-based API of
    MagicRecombiner.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF (gzip-compressed and indexed are fine).
    samples : list of str
        Samples to extract, in the desired order.
    chrom : str
        Chromosome/contig name (e.g. 'LG1').
    skip_missing : bool
        If True (default), drop any SNP that has a missing genotype
        (./. or -1) in any of the requested samples.

    Returns
    -------
    positions : np.ndarray, shape (n_snps,)
        1-based genomic positions.
    haplotypes : np.ndarray, shape (n_samples, n_snps), dtype=int8
        Allele calls (0 = ref, 1 = alt).
    """
    n_samples = len(samples)
    logger.info(
        f"Extracting chromosome {chrom} from {vcf_path} "
        f"for {n_samples} samples (skip_missing={skip_missing})"
    )
    logger.debug(f"Target samples: {samples}")

    vcf = VCF(vcf_path, samples=samples)

    positions = []
    rows = []
    skipped_missing = 0
    skipped_heterozygous = 0
    n_variants = 0

    for variant in vcf(chrom):
        n_variants += 1
        gts = variant.genotypes
        row = np.empty(n_samples, dtype=np.int8)
        has_missing = False

        for i in range(n_samples):
            a1, a2, _ = gts[i]
            if a1 == -1 or a2 == -1:
                has_missing = True
                if skip_missing:
                    break
                else:
                    row[i] = -1
            elif a1 != a2:
                # Heterozygous site — skip regardless of skip_missing
                skipped_heterozygous += 1
                row = None
                break
            else:
                # Homozygous: store one copy
                row[i] = a1

        if row is None:
            continue
        if has_missing and skip_missing:
            skipped_missing += 1
            continue

        positions.append(variant.POS)
        rows.append(row)

    if not positions:
        logger.error(f"No valid SNPs found on {chrom} for the requested samples.")
        raise ValueError(
            f"No valid SNPs found on {chrom} for the requested samples."
        )

    positions_arr = np.array(positions, dtype=np.int32)
    haps_arr = np.stack(rows, axis=1)

    bp_span = int(positions_arr[-1] - positions_arr[0])
    logger.info(
        f"Extracted {len(positions)} SNPs from {n_variants} variants "
        f"(skipped {skipped_missing} missing, {skipped_heterozygous} heterozygous)"
    )
    logger.info(
        f"Genotype matrix shape: {haps_arr.shape}, "
        f"span: {bp_span:,} bp, density: {len(positions) / max(bp_span, 1) * 1e6:.2f} SNPs/Mb"
    )
    logger.debug(f"Position range: {positions_arr[0]} - {positions_arr[-1]}")
    logger.debug(f"First 5 positions: {positions_arr[:5].tolist()}")
    logger.debug(f"Last 5 positions: {positions_arr[-5:].tolist()}")

    return positions_arr, haps_arr


def vcf_to_magic_inputs(
    vcf_path: str,
    parent_samples: List[str],
    progeny_samples: List[str],
    chrom: str,
    skip_missing: bool = True,
) -> Dict:
    """
    Load one chromosome from a VCF and prepare inputs for MagicRecombiner.
    """
    logger.info(f"Preparing MAGIC inputs for chromosome {chrom}")
    logger.debug(f"  VCF path: {vcf_path}")
    logger.debug(f"  Parents ({len(parent_samples)}): {parent_samples}")
    logger.debug(f"  Progeny ({len(progeny_samples)}): {len(progeny_samples)} total")

    # Merge parent and progeny samples, preserving order and uniqueness
    all_samples = list(dict.fromkeys(parent_samples + progeny_samples))
    logger.debug(f"  Total samples to extract: {len(all_samples)}")

    # Extract chromosome data from VCF
    positions, haps = extract_homozygous_chromosome(
        vcf_path, all_samples, chrom, skip_missing=skip_missing
    )

    # Build index mappings
    parent_idx = [all_samples.index(s) for s in parent_samples]
    progeny_idx = [all_samples.index(s) for s in progeny_samples]

    # Extract parent haplotypes
    parent_haps = haps[parent_idx, :]
    logger.debug(f"  Parent haplotype matrix shape: {parent_haps.shape}")

    # Extract progeny haplotypes
    progeny_haps_dict = {}
    for s, idx in zip(progeny_samples, progeny_idx):
        progeny_haps_dict[s] = haps[idx : idx + 1, :]
    
    logger.info(f"  Successfully prepared data for {len(progeny_samples)} progeny and {len(parent_samples)} parents")
    logger.info(f"  Total SNPs retained: {len(positions):,}")

    return {
        "parent_haps": parent_haps,
        "parent_names": parent_samples,
        "progeny_haps": progeny_haps_dict,
        "progeny_names": progeny_samples,
        "positions": positions,
        "chrom": chrom,
    }


def save_chromosome_npz(
    out_path: str,
    parent_haps: np.ndarray,
    parent_names: List[str],
    progeny_haps: Dict[str, np.ndarray],
    progeny_names: List[str],
    positions: np.ndarray,
    chrom: str,
):
    """
    Save prepared synthetic multi-parental population inputs to an .npz archive.
    """
    logger.info(f"Saving data to .npz: {out_path}")
    logger.debug(f"  Chromosome: {chrom}")
    logger.debug(f"  SNPs: {len(positions):,}")
    logger.debug(f"  Parents: {len(parent_names)}")
    logger.debug(f"  Progeny: {len(progeny_names)}")

    data = {
        "parent_haps": parent_haps,
        "parent_names": np.array(parent_names),
        "progeny_names": np.array(progeny_names),
        "positions": positions,
        "chrom": np.array([chrom]),
    }
    for name, haps in progeny_haps.items():
        safe_name = name.replace("/", "_").replace("\\", "_")
        data[f"progeny_{safe_name}"] = haps

    np.savez_compressed(out_path, **data)
    logger.debug(f"  Compression completed. Total arrays saved: {len(data)}")



def load_chromosome_npz(npz_path: str) -> Dict:
    """
    Load data written by save_chromosome_npz.
    """
    logger.info(f"Loading chromosome data from {npz_path}")

    data = np.load(npz_path, allow_pickle=True)
    logger.debug(f"Available keys in npz: {list(data.keys())}")

    result = {
        "parent_haps": data["parent_haps"],
        "parent_names": data["parent_names"].tolist(),
        "progeny_names": data["progeny_names"].tolist(),
        "positions": data["positions"],
        "chrom": str(data["chrom"][0]),
        "progeny_haps": {},
    }

    for name in result["progeny_names"]:
        safe_name = name.replace("/", "_").replace("\\", "_")
        result["progeny_haps"][name] = data[f"progeny_{safe_name}"]

    n_parents = len(result["parent_names"])
    n_progeny = len(result["progeny_names"])
    n_snps = len(result["positions"])
    logger.info(
        f"Loaded {result['chrom']}: {n_snps} SNPs, "
        f"{n_parents} parents, {n_progeny} progeny"
    )
    logger.debug(f"Parent names: {result['parent_names']}")
    logger.debug(f"Progeny names (first 10): {result['progeny_names'][:10]}")

    return result
