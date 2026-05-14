#!/usr/bin/env python3
"""Command-line interface for RecombTracer."""

import argparse
import sys
import os
from .core.recombiner import demo as demo_recombiner
from .core.hmm import demo_hmm


def main():
    parser = argparse.ArgumentParser(
        prog="recombtracer",
        description="RecombTracer — Recombination detection in MAGIC populations",
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # ------------------------------------------------------------------
    # demo-recombiner
    # ------------------------------------------------------------------
    parser_demo_r = subparsers.add_parser(
        "demo-recombiner",
        help="Run the synthetic MAGIC recombiner demo",
    )

    # ------------------------------------------------------------------
    # demo-hmm
    # ------------------------------------------------------------------
    parser_demo_h = subparsers.add_parser(
        "demo-hmm",
        help="Run the synthetic PBWT + HMM demo",
    )

    # ------------------------------------------------------------------
    # convert-vcf
    # ------------------------------------------------------------------
    parser_vcf = subparsers.add_parser(
        "convert-vcf",
        help="Convert a VCF to MAGIC Recombiner inputs (.npz)",
    )
    parser_vcf.add_argument("vcf", help="Input VCF (.vcf.gz)")
    parser_vcf.add_argument(
        "--parents",
        required=True,
        help="Comma-separated list of parent sample names",
    )
    parser_vcf.add_argument(
        "--progeny",
        default=None,
        help="Comma-separated list of progeny sample names (default: all non-parents)",
    )
    parser_vcf.add_argument(
        "--chrom",
        required=True,
        help="Chromosome to extract (e.g. LG1)",
    )
    parser_vcf.add_argument(
        "--out-dir",
        default=".",
        help="Output directory for .npz files (default: .)",
    )
    parser_vcf.add_argument(
        "--keep-missing",
        action="store_true",
        help="Keep SNPs with missing genotypes (default: drop them)",
    )

    # ------------------------------------------------------------------
    # run
    # ------------------------------------------------------------------
    parser_run = subparsers.add_parser(
        "run",
        help="Run PBWT + HMM recombination analysis from a .npz file",
    )
    parser_run.add_argument("npz", help="Input .npz file (from convert-vcf)")
    parser_run.add_argument(
        "--out-dir",
        default=".",
        help="Output directory for result CSVs (default: .)",
    )
    parser_run.add_argument(
        "--smooth-window",
        type=int,
        default=5,
        help="PBWT median filter window size (default: 5)",
    )
    parser_run.add_argument(
        "--min-segment-snps",
        type=int,
        default=5,
        help="Minimum SNPs per ancestry segment (default: 5)",
    )
    parser_run.add_argument(
        "--min-segment-bp",
        type=int,
        default=1000,
        help="Minimum bp length per ancestry segment (default: 1000)",
    )
    parser_run.add_argument(
        "--min-posterior",
        type=float,
        default=0.8,
        help="Minimum HMM posterior probability at breakpoint (default: 0.8)",
    )
    parser_run.add_argument(
        "--save-raw",
        action="store_true",
        help="Also save raw PBWT results (default: only HMM results)",
    )
    parser_run.add_argument(
        "--haplotype",
        type=int,
        default=None,
        help="Only analyze a specific haplotype index (default: all)",
    )
    parser_run.add_argument(
        "--progeny",
        default=None,
        help="Comma-separated list of progeny to analyze (default: all in .npz)",
    )

    args = parser.parse_args()

    if args.command == "demo-recombiner":
        demo_recombiner()
    elif args.command == "demo-hmm":
        demo_hmm()
    elif args.command == "convert-vcf":
        from .core.vcf import list_vcf_samples, vcf_to_magic_inputs, save_chromosome_npz

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

        print(f"Chromosome: {args.chrom}")
        print(f"Parents   ({len(parent_samples)}): {', '.join(parent_samples)}")
        print(f"Progeny   ({len(progeny_samples)}): {', '.join(progeny_samples)}")

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

        print(f"  SNPs    : {len(magic_data['positions']):,}")
        print(f"  Parents : {magic_data['parent_haps'].shape}")
        print(f"  Progeny : {len(magic_data['progeny_haps'])} individuals")

    elif args.command == "run":
        from .core.vcf import load_chromosome_npz
        from .core.recombiner import MagicRecombiner
        from .core.hmm import run_hmm_refinement
        import pandas as pd

        if not os.path.isfile(args.npz):
            sys.exit(f"Error: file not found: {args.npz}")

        os.makedirs(args.out_dir, exist_ok=True)

        print(f"Loading: {args.npz}")
        data = load_chromosome_npz(args.npz)
        chrom = data["chrom"]

        print(f"Chromosome: {chrom}")
        print(f"Parents   : {len(data['parent_names'])}")
        print(f"Progeny   : {len(data['progeny_names'])}")
        print(f"SNPs      : {len(data['positions']):,}")
        print()

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

            # Determine which haplotypes to process
            hap_indices = [args.haplotype] if args.haplotype is not None else list(range(ploidy))
            for h in hap_indices:
                if h < 0 or h >= ploidy:
                    print(f"  Warning: haplotype {h} out of range for {prog_name} (ploidy={ploidy})")
                    continue

                print(f"  Processing {prog_name} haplotype {h} ...", end=" ", flush=True)

                # --- PBWT paint ---
                paint_df = recombiner.paint_progeny(
                    prog_hap,
                    progeny_name=prog_name,
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

                # Filter HMM recombinations by posterior if requested
                if args.min_posterior > 0 and not rec_hmm.empty:
                    rec_hmm = rec_hmm[rec_hmm["confidence"] >= args.min_posterior].reset_index(drop=True)

                print(f"raw={len(rec_raw)} HMM={len(rec_hmm)}")

                # --- Save per-progeny results ---
                safe_name = prog_name.replace("/", "_").replace("\\", "_")
                base = f"{safe_name}_hap{h}_{chrom}"

                viterbi_df.to_csv(
                    os.path.join(args.out_dir, f"{base}_hmm_viterbi.csv"), index=False
                )
                seg_hmm.to_csv(
                    os.path.join(args.out_dir, f"{base}_hmm_segments.csv"), index=False
                )
                rec_hmm.to_csv(
                    os.path.join(args.out_dir, f"{base}_hmm_recombinations.csv"), index=False
                )

                if args.save_raw:
                    paint_df.to_csv(
                        os.path.join(args.out_dir, f"{base}_paint.csv"), index=False
                    )
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

        # --- Save summary ---
        summary_df = pd.DataFrame(summary_rows)
        summary_path = os.path.join(args.out_dir, f"summary_{chrom}.csv")
        summary_df.to_csv(summary_path, index=False)
        print()
        print(summary_df.to_string(index=False))
        print()
        print(f"Results saved to: {args.out_dir}")
        print(f"Summary saved to: {summary_path}")

    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
