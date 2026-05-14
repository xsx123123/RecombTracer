#!/usr/bin/env python3
"""
VCF utilities for MAGIC Recombiner.

Convert VCF files (especially homozygous-only VCFs) into the NumPy arrays
required by MagicRecombiner and MagicHMM.

Dependencies: cyvcf2, numpy
"""

import numpy as np
from cyvcf2 import VCF
from typing import List, Tuple, Dict, Optional
import os


def list_vcf_samples(vcf_path: str) -> List[str]:
    """
    Return all sample names in the VCF.
    """
    vcf = VCF(vcf_path)
    return list(vcf.samples)


def extract_chromosome(vcf_path: str,samples: List[str],chrom: str,skip_missing: bool = True,) -> Tuple[np.ndarray, np.ndarray]:
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
    # get sample variant data by cyvcf2
    vcf = VCF(vcf_path, samples=samples)

    positions = []
    rows = []
    n_samples = len(samples)

    for variant in vcf(chrom):
        # vcf(chrom) already filters by contig when an index is present

        # get variant genotypes
        gts = variant.genotypes
        # create empty numpy array
        row = np.empty(n_samples, dtype=np.int8)

        for i in range(n_samples):
            a1, a2, _ = gts[i]
            if a1 == -1 or a2 == -1:
                if skip_missing:
                    row = None
                    break
                else:
                    row[i] = -1
            else:
                # Homozygous sites: a1 == a2, store one copy
                row[i] = a1

        if row is None:
            continue

        positions.append(variant.POS)
        rows.append(row)

    if not positions:
        raise ValueError(
            f"No valid SNPs found on {chrom} for the requested samples."
        )

    positions_arr = np.array(positions, dtype=np.int32)
    haps_arr = np.stack(rows, axis=1) 
    return positions_arr, haps_arr


def vcf_to_magic_inputs(vcf_path: str,parent_samples: List[str],progeny_samples: List[str],chrom: str,skip_missing: bool = True) -> Dict:
    """
    Load one chromosome from a VCF and prepare inputs for MagicRecombiner.

    Parameters
    ----------
    vcf_path : str
    parent_samples : list of str
        Founder/parent sample names.
    progeny_samples : list of str
        Progeny sample names.
    chrom : str
        Chromosome to load.
    skip_missing : bool
        Drop SNPs with missing genotypes (recommended).

    Returns
    -------
    dict with keys:
        - 'parent_haps'   : np.ndarray, shape (n_parents, n_snps)
        - 'parent_names'  : list of str
        - 'progeny_haps'  : dict, {name: np.ndarray shape (1, n_snps)}
        - 'progeny_names' : list of str
        - 'positions'     : np.ndarray, shape (n_snps,)
        - 'chrom'         : str
    """
    all_samples = list(dict.fromkeys(parent_samples + progeny_samples))
    positions, haps = extract_chromosome(
        vcf_path, all_samples, chrom, skip_missing=skip_missing
    )

    parent_idx = [all_samples.index(s) for s in parent_samples]
    progeny_idx = [all_samples.index(s) for s in progeny_samples]

    parent_haps = haps[parent_idx, :]
    progeny_haps_dict = {}
    for s, idx in zip(progeny_samples, progeny_idx):
        progeny_haps_dict[s] = haps[idx : idx + 1, :]   # shape (1, n_snps)

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
    Save prepared MAGIC inputs to an .npz archive.

    Because npz stores flat arrays, progeny haplotypes are saved with keys
    like 'progeny_<name>'.
    """
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
    print(f"Saved: {out_path}")


def load_chromosome_npz(npz_path: str) -> Dict:
    """Load data written by save_chromosome_npz."""
    data = np.load(npz_path, allow_pickle=True)
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
    return result


# Command-line entry point moved to recombtracer.cli
