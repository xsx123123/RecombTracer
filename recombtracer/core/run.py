"""
PBWT + HMM analysis pipeline for the ``run`` CLI command.
"""

import os
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

from .hmm import run_hmm_refinement
from .recombiner import MagicRecombiner
from .vcf import load_chromosome_npz


def _analyze_single_haplotype(task):
    """
    Worker function for parallel haplotype analysis.

    Parameters
    ----------
    task : tuple
        (idx, total, prog_name, h, prog_hap, recombiner, parent_haps,
         parent_names, chrom, args, logger)

    Returns
    -------
    dict
        Summary row for this haplotype.
    """
    from ..utils.log_utils import logger

    (
        idx,
        total,
        prog_name,
        h,
        prog_hap,
        recombiner,
        parent_haps,
        parent_names,
        chrom,
        args,
    ) = task

    ploidy = prog_hap.shape[0]
    logger.info(f"[{idx}/{total}] Processing individual: {prog_name} (ploidy={ploidy})")

    if h < 0 or h >= ploidy:
        logger.warning(
            f"  Warning: haplotype {h} out of range for {prog_name} (ploidy={ploidy})"
        )
        return None

    logger.info(f"  Processing haplotype {h} ...")
    logger.debug(f"  Starting PBWT painting for {prog_name}_hap{h}")

    # --- PBWT paint ---
    paint_df = recombiner.paint_progeny(
        prog_hap[h : h + 1],
        progeny_name=prog_name,
        min_match_len=args.min_match_len,
        smooth_window=args.smooth_window,
    )
    paint_hap = paint_df[paint_df["haplotype"] == 0].copy()
    paint_hap["haplotype"] = h

    logger.debug(f"  PBWT paint completed. Extracting raw segments...")
    segments_raw = recombiner.extract_segments(
        paint_hap,
        min_segment_snps=args.min_segment_snps,
        min_segment_bp=args.min_segment_bp,
    )
    rec_raw = recombiner.call_recombinations(segments_raw)
    logger.debug(f"  Raw analysis: {len(segments_raw)} segments, {len(rec_raw)} recombinations")

    # --- HMM refinement ---
    logger.debug(f"  Starting HMM refinement for {prog_name}_hap{h}")
    viterbi_df, seg_hmm, rec_hmm = run_hmm_refinement(
        paint_hap,
        parent_haps,
        parent_names,
        prog_hap[h],
        progeny_name=prog_name,
        chrom=chrom,
    )

    if args.min_posterior > 0 and not rec_hmm.empty:
        count_before = len(rec_hmm)
        rec_hmm = rec_hmm[rec_hmm["confidence"] >= args.min_posterior].reset_index(drop=True)
        logger.debug(f"  HMM confidence filtering (>= {args.min_posterior}): {count_before} -> {len(rec_hmm)}")

    logger.info(f"    Done: raw_rec={len(rec_raw)} HMM_rec={len(rec_hmm)} HMM_segments={len(seg_hmm)}")

    # --- Save results ---
    safe_name = prog_name.replace("/", "_").replace("\\", "_")
    base = f"{safe_name}_hap{h}_{chrom}"

    logger.debug(f"  Saving HMM results to: {args.out_dir}/{base}_hmm_*.csv")
    viterbi_df.to_csv(os.path.join(args.out_dir, f"{base}_hmm_viterbi.csv"), index=False)
    seg_hmm.to_csv(os.path.join(args.out_dir, f"{base}_hmm_segments.csv"), index=False)
    rec_hmm.to_csv(os.path.join(args.out_dir, f"{base}_hmm_recombinations.csv"), index=False)

    if args.save_raw:
        logger.debug(f"  Saving raw results to: {args.out_dir}/{base}_paint.csv etc.")
        paint_hap.to_csv(os.path.join(args.out_dir, f"{base}_paint.csv"), index=False)
        pd.DataFrame(segments_raw).to_csv(
            os.path.join(args.out_dir, f"{base}_segments_raw.csv"), index=False
        )
        rec_raw.to_csv(
            os.path.join(args.out_dir, f"{base}_recombinations_raw.csv"), index=False
        )

    return {
        "progeny": prog_name,
        "haplotype": h,
        "chrom": chrom,
        "raw_breakpoints": len(rec_raw),
        "hmm_breakpoints": len(rec_hmm),
        "hmm_segments": len(seg_hmm),
    }


def handle_run(args, logger=None):
    """
    Pipeline for the ``run`` command.

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

    logger.debug(f"Input arguments: {args}")

    if not os.path.isfile(args.npz):
        logger.error(f"Input .npz file not found: {args.npz}")
        sys.exit(f"Error: file not found: {args.npz}")

    os.makedirs(args.out_dir, exist_ok=True)

    logger.info(f"Loading chromosome data: {args.npz}")
    data = load_chromosome_npz(args.npz)
    chrom = data["chrom"]

    logger.info(f"Analysis parameters: chrom={chrom}, parents={len(data['parent_names'])}, progeny={len(data['progeny_names'])}, SNPs={len(data['positions']):,}")
    logger.debug(f"Parent names: {data['parent_names']}")
    logger.debug(f"Position range: {data['positions'][0]} - {data['positions'][-1]}")

    recombiner = MagicRecombiner(
        parent_haps=data["parent_haps"],
        parent_names=data["parent_names"],
        positions=data["positions"],
        chrom=chrom,
    )

    summary_rows = []
    progeny_names = data["progeny_names"]
    if args.progeny:
        selected = [s.strip() for s in args.progeny.split(",")]
        logger.info(f"Filtering progeny: {len(selected)} requested")
        missing = set(selected) - set(progeny_names)
        if missing:
            logger.error(f"Progeny samples not found in .npz: {missing}")
            sys.exit(f"Error: progeny not in .npz: {missing}")
        progeny_names = selected

    # ------------------------------------------------------------------
    # Build task list
    # ------------------------------------------------------------------
    tasks = []
    total_progeny = len(progeny_names)
    for idx, prog_name in enumerate(progeny_names, 1):
        prog_hap = data["progeny_haps"][prog_name]
        ploidy = prog_hap.shape[0]
        hap_indices = [args.haplotype] if args.haplotype is not None else list(range(ploidy))
        for h in hap_indices:
            tasks.append(
                (
                    idx,
                    total_progeny,
                    prog_name,
                    h,
                    prog_hap,
                    recombiner,
                    data["parent_haps"],
                    data["parent_names"],
                    chrom,
                    args,
                )
            )

    workers = getattr(args, "workers", 1)
    if workers == 0:
        import multiprocessing
        workers = multiprocessing.cpu_count() or 1

    # ------------------------------------------------------------------
    # Execute analysis
    # ------------------------------------------------------------------
    if workers <= 1:
        for task in tasks:
            result = _analyze_single_haplotype(task)
            if result is not None:
                summary_rows.append(result)
    else:
        logger.info(f"Running with {workers} parallel workers ...")
        summary_lock = threading.Lock()
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(_analyze_single_haplotype, t): t for t in tasks}
            for future in as_completed(futures):
                try:
                    result = future.result()
                    if result is not None:
                        with summary_lock:
                            summary_rows.append(result)
                except Exception as exc:
                    task = futures[future]
                    prog_name, h = task[2], task[3]
                    logger.error(f"Exception in {prog_name}_hap{h}: {exc}")

    # ------------------------------------------------------------------
    # Save summary
    # ------------------------------------------------------------------
    summary_df = pd.DataFrame(summary_rows)
    if not summary_df.empty:
        summary_df = summary_df.sort_values(["progeny", "haplotype"]).reset_index(drop=True)
    summary_path = os.path.join(args.out_dir, f"summary_{chrom}.csv")
    summary_df.to_csv(summary_path, index=False)
    logger.info(f"\nFinal Summary for {chrom}:\n{summary_df.to_string(index=False)}")
    logger.info(f"All results saved to directory: {args.out_dir}")
