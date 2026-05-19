"""
Recombination hotspot and coldspot (conserved region) statistics.
"""

from typing import Optional, Tuple

import numpy as np
import pandas as pd


def compute_recombination_hotspots(
    rec_df: pd.DataFrame,
    chrom_start: int,
    chrom_end: int,
    window_size: int = 1_000_000,
    step_size: Optional[int] = None,
    hotspot_threshold_std: float = 2.0,
    coldspot_max_rate: float = 0.0,
) -> pd.DataFrame:
    """
    Compute recombination density across sliding windows.

    Parameters
    ----------
    rec_df : pd.DataFrame
        Merged recombination breakpoints (HMM). Must contain 'position' column.
    chrom_start : int
        Start coordinate of the chromosome (bp).
    chrom_end : int
        End coordinate of the chromosome (bp).
    window_size : int
        Window size in bp (default 1 Mb).
    step_size : int, optional
        Step size for sliding windows. If None, uses non-overlapping windows
        (step_size = window_size).
    hotspot_threshold_std : float
        Hotspot defined as rate > mean + threshold_std * std.
    coldspot_max_rate : float
        Coldspot (conserved region) defined as rate <= this value.
        Default 0 means absolutely no recombination events in the window.

    Returns
    -------
    pd.DataFrame
        Columns: start, end, center, count, rate_per_mb, is_hotspot, is_coldspot
    """
    if rec_df.empty:
        # Return empty windows with zeros
        step = step_size if step_size is not None else window_size
        n_windows = int(np.ceil((chrom_end - chrom_start) / step))
        rows = []
        for i in range(n_windows):
            s = chrom_start + i * step
            e = min(s + window_size, chrom_end)
            rows.append(
                {
                    "start": s,
                    "end": e,
                    "center": (s + e) // 2,
                    "count": 0,
                    "rate_per_mb": 0.0,
                    "is_hotspot": False,
                    "is_coldspot": True,
                }
            )
        return pd.DataFrame(rows)

    step = step_size if step_size is not None else window_size
    positions = rec_df["position"].values

    rows = []
    n_windows = int(np.ceil((chrom_end - chrom_start) / step))
    for i in range(n_windows):
        s = chrom_start + i * step
        e = min(s + window_size, chrom_end)
        count = int(((positions >= s) & (positions < e)).sum())
        span_mb = max((e - s) / 1_000_000, 1e-6)
        rate = count / span_mb
        rows.append(
            {
                "start": s,
                "end": e,
                "center": (s + e) // 2,
                "count": count,
                "rate_per_mb": rate,
                "is_hotspot": False,
                "is_coldspot": False,
            }
        )

    df = pd.DataFrame(rows)
    mean_rate = df["rate_per_mb"].mean()
    std_rate = df["rate_per_mb"].std()

    if std_rate > 0:
        df["is_hotspot"] = df["rate_per_mb"] > (mean_rate + hotspot_threshold_std * std_rate)
    else:
        df["is_hotspot"] = False

    df["is_coldspot"] = df["rate_per_mb"] <= coldspot_max_rate
    return df


def summarize_breakpoints(rec_df: pd.DataFrame) -> dict:
    """
    Generate text summary statistics from recombination breakpoint DataFrame.
    """
    if rec_df.empty:
        return {
            "total_breakpoints": 0,
            "unique_positions": 0,
            "mean_confidence": 0.0,
            "progeny_count": 0,
        }

    return {
        "total_breakpoints": int(len(rec_df)),
        "unique_positions": int(rec_df["position"].nunique()),
        "mean_confidence": float(rec_df["confidence"].mean()) if "confidence" in rec_df.columns else 0.0,
        "progeny_count": int(rec_df["progeny"].nunique()) if "progeny" in rec_df.columns else 0,
    }
