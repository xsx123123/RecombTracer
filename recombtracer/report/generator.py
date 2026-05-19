"""
HTML report generator for RecombTracer.

Loads PBWT + HMM CSV outputs and produces a single interactive HTML file.
"""

import glob
import os
from datetime import datetime
from typing import Dict, List, Optional

import pandas as pd
from jinja2 import Environment, PackageLoader, select_autoescape

from .plots import (
    plot_breakpoint_distribution,
    plot_mosaic,
    plot_parent_contribution_summary,
    plot_recombination_count_per_progeny,
    plot_recombination_landscape,
)
from .stats import compute_recombination_hotspots, summarize_breakpoints


class ReportGenerator:
    """
    Generate an interactive HTML report from RecombTracer CSV results.
    """

    def __init__(
        self,
        results_dir: str,
        out_html: str = "recombtracer_report.html",
        chrom: Optional[str] = None,
        window_size: int = 1_000_000,
        step_size: Optional[int] = None,
        hotspot_threshold_std: float = 2.0,
    ):
        """
        Parameters
        ----------
        results_dir : str
            Directory containing ``*_hmm_recombinations.csv``,
            ``*_hmm_segments.csv``, and ``summary_*.csv`` files.
        out_html : str
            Output HTML file path.
        chrom : str, optional
            If given, only include files matching this chromosome.
        window_size : int
            Window size (bp) for hotspot calculation.
        step_size : int, optional
            Step size for sliding windows. Defaults to non-overlapping (window_size).
        hotspot_threshold_std : float
            Standard-deviation multiplier for hotspot threshold.
        """
        self.results_dir = os.path.abspath(results_dir)
        self.out_html = os.path.abspath(out_html)
        self.chrom = chrom
        self.window_size = window_size
        self.step_size = step_size
        self.hotspot_threshold_std = hotspot_threshold_std

        self._env = Environment(
            loader=PackageLoader("recombtracer.report", "templates"),
            autoescape=select_autoescape(["html", "xml"]),
        )

        # Data containers populated by load_data()
        self.all_recs: Optional[pd.DataFrame] = None
        self.all_segments: Optional[pd.DataFrame] = None
        self.summary_df: Optional[pd.DataFrame] = None
        self.progeny_segments: Dict[str, pd.DataFrame] = {}
        self.chrom_start: int = 0
        self.chrom_end: int = 0
        self.parent_names: List[str] = []
        self.chrom_found: Optional[str] = None

    def load_data(self) -> None:
        """
        Scan ``results_dir`` and load all relevant CSVs.
        """
        if not os.path.isdir(self.results_dir):
            raise FileNotFoundError(f"Results directory not found: {self.results_dir}")

        rec_files = glob.glob(os.path.join(self.results_dir, "*_hmm_recombinations.csv"))
        seg_files = glob.glob(os.path.join(self.results_dir, "*_hmm_segments.csv"))
        summary_files = glob.glob(os.path.join(self.results_dir, "summary_*.csv"))

        if not rec_files and not seg_files:
            raise ValueError(
                f"No *_hmm_recombinations.csv or *_hmm_segments.csv found in {self.results_dir}"
            )

        # --- Load recombination breakpoints ---
        rec_dfs = []
        for f in sorted(rec_files):
            basename = os.path.basename(f)
            # Parse progeny and haplotype from filename: {prog}_hap{h}_{chrom}_hmm_recombinations.csv
            parts = basename.replace("_hmm_recombinations.csv", "").split("_")
            if len(parts) >= 3:
                file_chrom = parts[-1]
                if self.chrom and file_chrom != self.chrom:
                    continue
                self.chrom_found = file_chrom
                hap = parts[-2].replace("hap", "")
                prog = "_".join(parts[:-2])
            else:
                prog = basename
                hap = "0"
                file_chrom = self.chrom or "unknown"
                self.chrom_found = file_chrom

            df = pd.read_csv(f)
            df["progeny"] = prog
            if "haplotype" not in df.columns:
                df["haplotype"] = int(hap) if hap.isdigit() else 0
            rec_dfs.append(df)

        if rec_dfs:
            self.all_recs = pd.concat(rec_dfs, ignore_index=True)
        else:
            self.all_recs = pd.DataFrame()

        # --- Load segments ---
        seg_dfs = []
        progeny_seg_map: Dict[str, List[pd.DataFrame]] = {}
        for f in sorted(seg_files):
            basename = os.path.basename(f)
            parts = basename.replace("_hmm_segments.csv", "").split("_")
            if len(parts) >= 3:
                file_chrom = parts[-1]
                if self.chrom and file_chrom != self.chrom:
                    continue
                hap = parts[-2].replace("hap", "")
                prog = "_".join(parts[:-2])
            else:
                prog = basename
                hap = "0"
                file_chrom = self.chrom or "unknown"

            df = pd.read_csv(f)
            df["progeny"] = prog
            # Ensure haplotype column exists
            if "haplotype" not in df.columns:
                df["haplotype"] = int(hap) if hap.isdigit() else 0
            seg_dfs.append(df)
            progeny_seg_map.setdefault(prog, []).append(df)

        if seg_dfs:
            self.all_segments = pd.concat(seg_dfs, ignore_index=True)
        else:
            self.all_segments = pd.DataFrame()

        # Build per-progeny combined segments dict for mosaic plot
        for prog, dfs in progeny_seg_map.items():
            self.progeny_segments[prog] = pd.concat(dfs, ignore_index=True)

        # --- Load summary ---
        summary_dfs = []
        for f in sorted(summary_files):
            df = pd.read_csv(f)
            if self.chrom and "chrom" in df.columns:
                df = df[df["chrom"] == self.chrom]
            summary_dfs.append(df)
        if summary_dfs:
            self.summary_df = pd.concat(summary_dfs, ignore_index=True)
        else:
            self.summary_df = pd.DataFrame()

        # --- Infer chromosome range and parent names ---
        if self.all_segments is not None and not self.all_segments.empty:
            self.chrom_start = int(self.all_segments["start"].min())
            self.chrom_end = int(self.all_segments["end"].max())
            if "parent" in self.all_segments.columns:
                self.parent_names = sorted(self.all_segments["parent"].unique().tolist())
        elif self.all_recs is not None and not self.all_recs.empty:
            self.chrom_start = int(self.all_recs["position"].min())
            self.chrom_end = int(self.all_recs["position"].max())

        # Fallback: if positions look like 1-based small numbers, pad range
        if self.chrom_end <= self.chrom_start:
            self.chrom_end = self.chrom_start + 1

    def generate(self) -> str:
        """
        Load data, compute statistics, generate plots, render template, and write HTML.

        Returns
        -------
        str
            Path to the generated HTML file.
        """
        self.load_data()

        # Compute hotspot / coldspot stats
        hotspot_df = compute_recombination_hotspots(
            self.all_recs,
            self.chrom_start,
            self.chrom_end,
            window_size=self.window_size,
            step_size=self.step_size,
            hotspot_threshold_std=self.hotspot_threshold_std,
        )

        # Summary stats
        rec_summary = summarize_breakpoints(self.all_recs)
        total_progeny = len(self.progeny_segments)
        total_segments = len(self.all_segments) if self.all_segments is not None else 0
        hotspot_count = int(hotspot_df["is_hotspot"].sum()) if not hotspot_df.empty else 0
        coldspot_count = int(hotspot_df["is_coldspot"].sum()) if not hotspot_df.empty else 0

        # Generate Plotly figures as HTML divs
        fig_mosaic = plot_mosaic(
            self.progeny_segments,
            self.parent_names,
            self.chrom_found or self.chrom or "Chromosome",
            self.chrom_start,
            self.chrom_end,
        )
        fig_landscape = plot_recombination_landscape(
            hotspot_df, self.chrom_found or self.chrom or "Chromosome"
        )
        fig_dist = plot_breakpoint_distribution(
            self.all_recs,
            self.chrom_found or self.chrom or "Chromosome",
            self.chrom_start,
            self.chrom_end,
        )
        fig_contrib = plot_parent_contribution_summary(self.all_segments)
        fig_per_prog = plot_recombination_count_per_progeny(self.summary_df)

        # Convert to HTML div strings (embed plotly.js for offline use)
        div_mosaic = fig_mosaic.to_html(full_html=False, include_plotlyjs=False)
        div_landscape = fig_landscape.to_html(full_html=False, include_plotlyjs=False)
        div_dist = fig_dist.to_html(full_html=False, include_plotlyjs=False)
        div_contrib = fig_contrib.to_html(full_html=False, include_plotlyjs=False)
        div_per_prog = fig_per_prog.to_html(full_html=False, include_plotlyjs=False)

        # Prepare breakpoint table
        all_table_cols = [
            "progeny", "haplotype", "chrom", "position",
            "left_parent", "right_parent", "confidence"
        ]
        if self.all_recs is not None and not self.all_recs.empty:
            table_cols = [c for c in all_table_cols if c in self.all_recs.columns]
            table_data = self.all_recs[table_cols].to_dict(orient="records")
        else:
            table_cols = []
            table_data = []

        table_columns = [{"key": c, "label": c.replace("_", " ").title()} for c in table_cols]

        # Render template
        template = self._env.get_template("report.html.j2")
        html = template.render(
            title=f"RecombTracer Report — {self.chrom_found or self.chrom or 'Genome'}",
            generated_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            chrom=self.chrom_found or self.chrom or "N/A",
            total_progeny=total_progeny,
            total_breakpoints=rec_summary["total_breakpoints"],
            mean_confidence=round(rec_summary["mean_confidence"], 4),
            total_segments=total_segments,
            hotspot_count=hotspot_count,
            coldspot_count=coldspot_count,
            parent_names=self.parent_names,
            div_mosaic=div_mosaic,
            div_landscape=div_landscape,
            div_dist=div_dist,
            div_contrib=div_contrib,
            div_per_prog=div_per_prog,
            table_data=table_data,
            table_columns=table_columns,
        )

        os.makedirs(os.path.dirname(self.out_html) or ".", exist_ok=True)
        with open(self.out_html, "w", encoding="utf-8") as fh:
            fh.write(html)

        return self.out_html
