"""
PBWT + HMM analysis pipeline for the ``run`` CLI command.
"""

import os
import sys

import pandas as pd

from .hmm import run_hmm_refinement
from .recombiner import MagicRecombiner
from .vcf import load_chromosome_npz


def handle_run(args):
    """
    Pipeline for the ``run`` command.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments.
    """
    from ..utils import logger_generator

    logger, _ = logger_generator(args.out_dir, log_level=args.log_level)

    if not os.path.isfile(args.npz):
        sys.exit(f"Error: file not found: {args.npz}")

    os.makedirs(args.out_dir, exist_ok=True)

    logger.info(f"Loading: {args.npz}")
    data = load_chromosome_npz(args.npz)
    chrom = data["chrom"]

    logger.info(f"Chromosome: {chrom}")
    logger.info(f"Parents   : {len(data['parent_names'])}")
    logger.info(f"Progeny   : {len(data['progeny_names'])}")
    logger.info(f"SNPs      : {len(data['positions']):,}")

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
        missing = set(selected) - set(progeny_names)
        if missing:
            sys.exit(f"Error: progeny not in .npz: {missing}")
        progeny_names = selected

    for prog_name in progeny_names:
        prog_hap = data["progeny_haps"][prog_name]
        ploidy = prog_hap.shape[0]

        hap_indices = [args.haplotype] if args.haplotype is not None else list(range(ploidy))
        for h in hap_indices:
            if h < 0 or h >= ploidy:
                logger.warning(
                    f"  Warning: haplotype {h} out of range for {prog_name} (ploidy={ploidy})"
                )
                continue

            logger.info(f"  Processing {prog_name} haplotype {h} ...")

            # --- PBWT paint ---
            paint_df = recombiner.paint_progeny(
                prog_hap,
                progeny_name=prog_name,
                min_match_len=args.min_match_len,
                smooth_window=args.smooth_window,
            )
            paint_hap = paint_df[paint_df["haplotype"] == h]

            segments_raw = recombiner.extract_segments(
                paint_df,
                min_segment_snps=args.min_segment_snps,
                min_segment_bp=args.min_segment_bp,
            )
            rec_raw = recombiner.call_recombinations(segments_raw)

            # --- HMM refinement ---
            viterbi_df, seg_hmm, rec_hmm = run_hmm_refinement(
                paint_hap,
                data["parent_haps"],
                data["parent_names"],
                prog_hap[h],
                progeny_name=prog_name,
                chrom=chrom,
            )

            if args.min_posterior > 0 and not rec_hmm.empty:
                rec_hmm = rec_hmm[rec_hmm["confidence"] >= args.min_posterior].reset_index(drop=True)

            logger.info(f"    Done: raw={len(rec_raw)} HMM={len(rec_hmm)}")

            # --- Save results ---
            safe_name = prog_name.replace("/", "_").replace("\\", "_")
            base = f"{safe_name}_hap{h}_{chrom}"

            viterbi_df.to_csv(os.path.join(args.out_dir, f"{base}_hmm_viterbi.csv"), index=False)
            seg_hmm.to_csv(os.path.join(args.out_dir, f"{base}_hmm_segments.csv"), index=False)
            rec_hmm.to_csv(os.path.join(args.out_dir, f"{base}_hmm_recombinations.csv"), index=False)

            if args.save_raw:
                paint_df.to_csv(os.path.join(args.out_dir, f"{base}_paint.csv"), index=False)
                pd.DataFrame(segments_raw).to_csv(
                    os.path.join(args.out_dir, f"{base}_segments_raw.csv"), index=False
                )
                rec_raw.to_csv(
                    os.path.join(args.out_dir, f"{base}_recombinations_raw.csv"), index=False
                )

            summary_rows.append({
                "progeny": prog_name,
                "haplotype": h,
                "chrom": chrom,
                "raw_breakpoints": len(rec_raw),
                "hmm_breakpoints": len(rec_hmm),
                "hmm_segments": len(seg_hmm),
            })

    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(args.out_dir, f"summary_{chrom}.csv")
    summary_df.to_csv(summary_path, index=False)
    logger.info(f"\n{summary_df.to_string(index=False)}")
    logger.info(f"Results saved to: {args.out_dir}")
