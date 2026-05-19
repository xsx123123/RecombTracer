"""
RecombTracer HTML report generation module.

Provides interactive HTML reports from PBWT + HMM analysis results,
including mosaic plots, recombination hotspots, and breakpoint summaries.
"""

from .generator import ReportGenerator

__all__ = ["ReportGenerator"]
