#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Command-line interface for RecombTracer.

This module is intentionally thin: it only defines argument parsers and routes
commands. All business logic lives in ``recombtracer.core.pipeline``.
"""

import argparse
import sys

from rich_argparse import RichHelpFormatter

from . import SOFTWARE_INFO
from .utils import config2logo, show_versions
from .core.convert import handle_convert_vcf
from .core.run import handle_run
from .core.pipeline import handle_pipeline
from .report.generator import ReportGenerator


class _LogoHelpAction(argparse.Action):
    """Custom help action that prints logo before help text."""

    def __call__(self, parser, namespace, values, option_string=None):
        config2logo(SOFTWARE_INFO)
        parser.print_help()
        parser.exit()


# ── Common arguments shared by all subcommands 
_common_parent = argparse.ArgumentParser(add_help=False)
_common_parent.add_argument(
    "-l", "--log-level",
    default="INFO",
    choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    help="Set logging level (default: INFO)",
)


def add_convert_vcf_parser(subparsers):
    """
    Add parser for convert-vcf command.
    """
    parser = subparsers.add_parser(
        "convert-vcf",
        parents=[_common_parent],
        formatter_class=RichHelpFormatter,
        add_help=False,
        help="Convert a VCF to synthetic multi-parental Recombiner inputs (.npz)",
    )
    parser.add_argument(
        "-h", "--help",
        action=_LogoHelpAction,
        nargs=0,
        default=False,
        help="Show this help message and exit",
    )
    parser.add_argument("vcf", help="Input VCF (.vcf.gz)")
    parser_req = parser.add_argument_group("required arguments")
    parser_req.add_argument(
        "--parents",
        required=True,
        help="Comma-separated list of parent sample names",
    )
    parser_req.add_argument(
        "--chrom",
        required=True,
        help="Chromosome to extract (e.g. chr1)",
    )
    parser.add_argument(
        "--progeny",
        default=None,
        help="Comma-separated list of progeny sample names (default: all non-parents)",
    )
    parser.add_argument(
        "--out-dir",
        default=".",
        help="Output directory for .npz files (default: .)",
    )
    parser.add_argument(
        "--keep-missing",
        action="store_true",
        help="Keep SNPs with missing genotypes (default: drop them)",
    )


def add_run_parser(subparsers):
    """
    Add parser for run command.
    """
    parser = subparsers.add_parser(
        "run",
        parents=[_common_parent],
        formatter_class=RichHelpFormatter,
        add_help=False,
        help="Run PBWT + HMM recombination analysis from a .npz file",
    )
    parser.add_argument(
        "-h", "--help",
        action=_LogoHelpAction,
        nargs=0,
        default=False,
        help="Show this help message and exit",
    )
    parser.add_argument("npz", help="Input .npz file (from convert-vcf)")
    parser.add_argument(
        "--out-dir",
        default=".",
        help="Output directory for result CSVs (default: .)",
    )

    group_algo = parser.add_argument_group("algorithm parameters")
    group_algo.add_argument(
        "--min-match-len",
        type=int,
        default=2,
        help="Minimum length of PBWT match segments (default: 2)",
    )
    group_algo.add_argument(
        "--smooth-window",
        type=int,
        default=5,
        help="PBWT median filter window size (default: 5)",
    )
    group_algo.add_argument(
        "--min-segment-snps",
        type=int,
        default=5,
        help="Minimum SNPs per ancestry segment (default: 5)",
    )
    group_algo.add_argument(
        "--min-segment-bp",
        type=int,
        default=1000,
        help="Minimum bp length per ancestry segment (default: 1000)",
    )

    group_algo.add_argument(
        "--min-posterior",
        type=float,
        default=0.8,
        help="Minimum HMM posterior probability at breakpoint (default: 0.8)",
    )

    parser.add_argument(
        "--save-raw",
        action="store_true",
        help="Also save raw PBWT results (default: only HMM results)",
    )

    parser.add_argument(
        "--haplotype",
        type=int,
        default=None,
        help="Only analyze a specific haplotype index (default: all)",
    )

    parser.add_argument(
        "--progeny",
        default=None,
        help="Comma-separated list of progeny to analyze (default: all in .npz)",
    )


def _add_common_run_args(parser):
    """Add algorithm parameters shared by ``run`` and ``pipeline``."""
    group_algo = parser.add_argument_group("algorithm parameters")
    group_algo.add_argument(
        "--min-match-len",
        type=int,
        default=2,
        help="Minimum length of PBWT match segments (default: 2)",
    )
    group_algo.add_argument(
        "--smooth-window",
        type=int,
        default=5,
        help="PBWT median filter window size (default: 5)",
    )
    group_algo.add_argument(
        "--min-segment-snps",
        type=int,
        default=5,
        help="Minimum SNPs per ancestry segment (default: 5)",
    )
    group_algo.add_argument(
        "--min-segment-bp",
        type=int,
        default=1000,
        help="Minimum bp length per ancestry segment (default: 1000)",
    )
    group_algo.add_argument(
        "--min-posterior",
        type=float,
        default=0.8,
        help="Minimum HMM posterior probability at breakpoint (default: 0.8)",
    )
    parser.add_argument(
        "--save-raw",
        action="store_true",
        help="Also save raw PBWT results (default: only HMM results)",
    )
    parser.add_argument(
        "--haplotype",
        type=int,
        default=None,
        help="Only analyze a specific haplotype index (default: all)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=10,
        help="Number of parallel workers for progeny analysis (default: 10, 0=auto-detect CPU count)",
    )


def add_report_parser(subparsers):
    """
    Add parser for report command (generate HTML report from CSV results).
    """
    parser = subparsers.add_parser(
        "report",
        parents=[_common_parent],
        formatter_class=RichHelpFormatter,
        add_help=False,
        help="Generate an interactive HTML report from run results",
    )
    parser.add_argument(
        "-h", "--help",
        action=_LogoHelpAction,
        nargs=0,
        default=False,
        help="Show this help message and exit",
    )
    parser.add_argument("results_dir", help="Directory containing *_hmm_*.csv result files")
    parser.add_argument(
        "--out",
        default="recombtracer_report.html",
        help="Output HTML file path (default: recombtracer_report.html)",
    )
    parser.add_argument(
        "--chrom",
        default=None,
        help="Only report data for this chromosome (default: all found)",
    )
    parser.add_argument(
        "--window-size",
        type=int,
        default=1_000_000,
        help="Window size (bp) for hotspot calculation (default: 1,000,000)",
    )
    parser.add_argument(
        "--step-size",
        type=int,
        default=None,
        help="Step size for sliding windows (default: same as window-size)",
    )
    parser.add_argument(
        "--hotspot-threshold",
        type=float,
        default=2.0,
        help="Standard-deviation multiplier for hotspot threshold (default: 2.0)",
    )


def add_pipeline_parser(subparsers):
    """
    Add parser for pipeline command (convert-vcf + run in one step).
    """
    parser = subparsers.add_parser(
        "pipeline",
        parents=[_common_parent],
        formatter_class=RichHelpFormatter,
        add_help=False,
        help="Run the full pipeline: VCF → .npz → PBWT + HMM analysis",
    )
    parser.add_argument(
        "-h", "--help",
        action=_LogoHelpAction,
        nargs=0,
        default=False,
        help="Show this help message and exit",
    )
    parser.add_argument("vcf", help="Input VCF (.vcf.gz)")

    parser_req = parser.add_argument_group("required arguments")
    parser_req.add_argument(
        "--parents",
        required=True,
        help="Comma-separated list of parent sample names",
    )
    parser_req.add_argument(
        "--chrom",
        required=True,
        help="Chromosome to extract (e.g. chr1)",
    )

    parser.add_argument(
        "--progeny",
        default=None,
        help="Comma-separated list of progeny sample names (default: all non-parents)",
    )
    parser.add_argument(
        "--out-dir",
        default=".",
        help="Output directory for all results (default: .)",
    )
    parser.add_argument(
        "--keep-missing",
        action="store_true",
        help="Keep SNPs with missing genotypes (default: drop them)",
    )

    _add_common_run_args(parser)


def main():
    """
    define main function for recombtracer scripts
    """
    # Create main parser (disable default -h so we can handle it manually)
    parser = argparse.ArgumentParser(
        formatter_class=RichHelpFormatter,
        prog="recombtracer",
        description="Recombination breakpoint detection for multi-parental populations",
        add_help=False,
    )
    # Add version flagß
    parser.add_argument("-v", "--version", action="store_true", help="Show version information")

    # Add help flag manually
    parser.add_argument("-h", "--help", action="store_true", help="Show this help message and exit")

    # Add subparsers
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Add subcommands
    add_convert_vcf_parser(subparsers)
    add_run_parser(subparsers)
    add_pipeline_parser(subparsers)
    add_report_parser(subparsers)

    args = parser.parse_args()

    # Show version and exit
    if args.version:
        show_versions(
            project_name=SOFTWARE_INFO.get("app_name", "RecombTracer"),
            deps=None,
            extras={
                "Version": SOFTWARE_INFO.get("version", "unknown"),
                "Author": SOFTWARE_INFO.get("author", "unknown"),
                "Email": SOFTWARE_INFO.get("email", ""),
                "URL": SOFTWARE_INFO.get("url", ""),
                "Description": SOFTWARE_INFO.get("description", ""),
            },
        )
        return

    # Show logo + help when: no command, -h/--help, or unknown command
    if args.help or not args.command:
        config2logo(SOFTWARE_INFO)
        parser.print_help()
        return

    if args.command in ["run", "convert-vcf", "pipeline"]:
        config2logo(SOFTWARE_INFO)

    # Route commands
    if args.command == "convert-vcf":
        handle_convert_vcf(args)
    elif args.command == "run":
        handle_run(args)
    elif args.command == "pipeline":
        handle_pipeline(args)
    elif args.command == "report":
        gen = ReportGenerator(
            results_dir=args.results_dir,
            out_html=args.out,
            chrom=args.chrom,
            window_size=args.window_size,
            step_size=args.step_size,
            hotspot_threshold_std=args.hotspot_threshold,
        )
        out_path = gen.generate()
        print(f"Report written to: {out_path}")
    else:
        config2logo(SOFTWARE_INFO)
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
