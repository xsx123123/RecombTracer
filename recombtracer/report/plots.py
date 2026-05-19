"""
Plotly-based interactive visualizations for RecombTracer HTML reports.
"""

import colorsys
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# --- Professional Color Palette ---
PALETTE = {
    "primary": "#0f172a",    # Slate 900
    "secondary": "#64748b",  # Slate 500
    "accent": "#0ea5e9",     # Sky 500
    "hotspot": "#ef4444",    # Red 500
    "conserved": "#3b82f6",  # Blue 500
    "neutral": "#94a3b8",    # Slate 400
    "grid": "#f1f5f9",       # Slate 100
    "background": "#ffffff",
}

DEFAULT_LAYOUT = dict(
    font_family="'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif",
    font_size=12,
    font_color=PALETTE["primary"],
    plot_bgcolor=PALETTE["background"],
    paper_bgcolor=PALETTE["background"],
    hoverlabel=dict(
        bgcolor="white",
        font_size=13,
        font_family="'Inter', sans-serif",
        bordercolor=PALETTE["grid"],
    ),
    margin=dict(l=60, r=40, t=80, b=60),
)


def _apply_style(fig: go.Figure, title: str) -> None:
    """Apply consistent styling to a figure."""
    fig.update_layout(**DEFAULT_LAYOUT)
    fig.update_layout(
        title=dict(
            text=f"<b>{title}</b>",
            x=0.05,
            y=0.95,
            font=dict(size=18, color=PALETTE["primary"]),
        ),
    )
    fig.update_xaxes(
        gridcolor=PALETTE["grid"],
        linecolor=PALETTE["secondary"],
        zeroline=False,
        tickfont=dict(color=PALETTE["secondary"]),
    )
    fig.update_yaxes(
        gridcolor=PALETTE["grid"],
        linecolor=PALETTE["secondary"],
        zeroline=False,
        tickfont=dict(color=PALETTE["secondary"]),
    )


def _generate_parent_colors(parent_names: List[str]) -> Dict[str, str]:
    """
    Assign a distinct, perceptually uniform color to each parent.
    Uses HSL space for good separation.
    """
    n = len(parent_names)
    colors = {}
    for i, name in enumerate(parent_names):
        # Using a more curated set of colors for the first few parents
        standard_colors = [
            "#0ea5e9", "#10b981", "#f59e0b", "#8b5cf6",
            "#f43f5e", "#06b6d4", "#84cc16", "#d946ef"
        ]
        if i < len(standard_colors):
            colors[name] = standard_colors[i]
        else:
            hue = (i * 0.618033988749895) % 1.0
            lightness = 0.55
            saturation = 0.7
            rgb = colorsys.hls_to_rgb(hue, lightness, saturation)
            colors[name] = "#%02x%02x%02x" % tuple(int(c * 255) for c in rgb)
    return colors


def plot_mosaic(
    progeny_segments: Dict[str, pd.DataFrame],
    parent_names: List[str],
    chrom: str,
    chrom_start: int,
    chrom_end: int,
) -> go.Figure:
    """
    Create an interactive mosaic plot showing ancestry segments per progeny / haplotype.
    """
    parent_colors = _generate_parent_colors(parent_names)
    progeny_names = sorted(progeny_segments.keys())

    fig = go.Figure()

    # Build y-axis mapping
    y_ticks = []
    y_ticktext = []
    y_pos = 0
    progeny_y_map = {}

    for prog in progeny_names:
        progeny_y_map[prog] = (y_pos + 0.4, y_pos + 1.4)
        y_ticks.append(y_pos + 0.9)
        y_ticktext.append(f"<b>{prog}</b>")
        y_pos += 2.2

    for prog in progeny_names:
        seg_df = progeny_segments[prog]
        if seg_df.empty:
            continue

        hap0_y, hap1_y = progeny_y_map[prog]

        for _, row in seg_df.iterrows():
            hap = int(row.get("haplotype", 0))
            y_center = hap0_y if hap == 0 else hap1_y
            parent = row.get("parent", "Unknown")
            start = int(row["start"])
            end = int(row["end"])
            width = end - start

            hover_text = (
                f"<b>{prog}</b> (H{hap})<br>"
                f"<span style='color:{parent_colors.get(parent, '#999999')}'>●</span> <b>{parent}</b><br>"
                f"Range: {start:,} - {end:,}<br>"
                f"Size: {width:,} bp"
            )
            if "n_snps" in row:
                hover_text += f"<br>Markers: {int(row['n_snps'])}"
            if "mean_posterior" in row:
                hover_text += f"<br>Confidence: {row['mean_posterior']:.3f}"

            fig.add_trace(
                go.Bar(
                    x=[width],
                    y=[y_center],
                    base=[start],
                    orientation="h",
                    marker=dict(
                        color=parent_colors.get(parent, "#999999"),
                        line=dict(width=0.5, color="white")
                    ),
                    width=0.8,
                    showlegend=False,
                    hovertemplate=hover_text + "<extra></extra>",
                )
            )

    # Add legend
    for parent in parent_names:
        fig.add_trace(
            go.Scatter(
                x=[None], y=[None],
                mode="markers",
                marker=dict(size=10, color=parent_colors[parent], symbol="square"),
                legendgroup=parent,
                showlegend=True,
                name=parent,
            )
        )

    _apply_style(fig, f"Ancestry Mosaic — {chrom}")
    fig.update_layout(
        xaxis=dict(title="Genomic Position (bp)", range=[chrom_start, chrom_end], tickformat=",d"),
        yaxis=dict(
            tickmode="array",
            tickvals=y_ticks,
            ticktext=y_ticktext,
            range=[-0.5, y_pos],
            showgrid=False,
            side="left"
        ),
        height=max(500, len(progeny_names) * 60),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1,
            title=None,
        ),
    )

    return fig


def plot_recombination_landscape(
    hotspot_df: pd.DataFrame,
    chrom: str,
) -> go.Figure:
    """
    Bar plot of recombination rate per window.
    """
    fig = go.Figure()

    types = [
        ("is_hotspot", PALETTE["hotspot"], "Hotspot", "Hotspot"),
        ("is_coldspot", PALETTE["conserved"], "Conserved", "Conserved"),
        (None, PALETTE["neutral"], "Normal", "Normal"),
    ]

    for condition, color, name, label in types:
        if condition:
            mask = hotspot_df[condition]
        else:
            mask = (~hotspot_df["is_hotspot"]) & (~hotspot_df["is_coldspot"])
        
        subset = hotspot_df[mask]
        if subset.empty:
            continue

        fig.add_trace(
            go.Bar(
                x=subset["center"],
                y=subset["rate_per_mb"],
                width=subset["end"] - subset["start"],
                marker_color=color,
                name=name,
                hovertemplate=(
                    f"<b>{label} Region</b><br>"
                    "Range: %{customdata[0]:,} - %{customdata[1]:,}<br>"
                    "Events: %{customdata[2]}<br>"
                    "Rate: %{y:.2f} / Mb<extra></extra>"
                ),
                customdata=np.stack(
                    [subset["start"], subset["end"], subset["count"]], axis=-1
                ),
            )
        )

    _apply_style(fig, f"Recombination Landscape — {chrom}")
    fig.update_layout(
        xaxis=dict(title="Position (bp)", tickformat=",d"),
        yaxis=dict(title="Recombination events / Mb"),
        barmode="overlay",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        height=400,
    )
    return fig


def plot_breakpoint_distribution(
    rec_df: pd.DataFrame,
    chrom: str,
    chrom_start: int,
    chrom_end: int,
    n_bins: int = 50,
) -> go.Figure:
    """
    Histogram of breakpoint positions with KDE overlay.
    """
    fig = make_subplots(specs=[[{"secondary_y": True}]])

    if rec_df.empty:
        _apply_style(fig, f"Breakpoint Distribution — {chrom}")
        return fig

    positions = rec_df["position"].values

    fig.add_trace(
        go.Histogram(
            x=positions,
            nbinsx=n_bins,
            name="Count",
            marker=dict(color=PALETTE["accent"], line=dict(width=0.5, color="white")),
            opacity=0.6,
            hovertemplate="<b>Position:</b> %{x:,}<br><b>Count:</b> %{y}<extra></extra>",
        ),
        secondary_y=False,
    )

    try:
        from scipy.stats import gaussian_kde
        kde = gaussian_kde(positions)
        x_range = np.linspace(chrom_start, chrom_end, 500)
        density = kde(x_range) * len(positions)
        fig.add_trace(
            go.Scatter(
                x=x_range,
                y=density,
                mode="lines",
                name="Density",
                line=dict(color=PALETTE["primary"], width=2.5),
                hovertemplate="<b>Density:</b> %{y:.2f}<extra></extra>",
            ),
            secondary_y=True,
        )
    except Exception:
        pass

    _apply_style(fig, f"Breakpoint Distribution — {chrom}")
    fig.update_layout(
        xaxis=dict(title="Position (bp)", tickformat=",d", range=[chrom_start, chrom_end]),
        yaxis=dict(title="Event Count"),
        yaxis2=dict(title="Density", overlaying="y", side="right", showgrid=False),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        height=400,
    )
    return fig


def plot_parent_contribution_summary(seg_df: pd.DataFrame) -> go.Figure:
    """
    Bar chart of total genomic contribution per parent.
    """
    if seg_df.empty or "parent" not in seg_df.columns:
        fig = go.Figure()
        _apply_style(fig, "Genomic Contribution")
        return fig

    seg_df = seg_df.copy()
    seg_df["span"] = seg_df["end"] - seg_df["start"]
    contrib = seg_df.groupby("parent")["span"].sum().reset_index()
    contrib = contrib.sort_values("span", ascending=True)

    parent_colors = _generate_parent_colors(contrib["parent"].tolist())
    colors = [parent_colors.get(p, PALETTE["neutral"]) for p in contrib["parent"]]

    total_span = contrib["span"].sum()
    contrib["pct"] = contrib["span"] / total_span * 100

    fig = go.Figure(
        go.Bar(
            y=contrib["parent"],
            x=contrib["span"],
            orientation="h",
            marker_color=colors,
            text=[f"  {p:.1f}%" for p in contrib["pct"]],
            textposition="outside",
            hovertemplate="<b>%{y}</b><br>Total span: %{x:,.0f} bp<br>Fraction: %{customdata:.2f}%<extra></extra>",
            customdata=contrib["pct"],
        )
    )

    _apply_style(fig, "Parental Contribution")
    fig.update_layout(
        xaxis=dict(title="Total Span (bp)", tickformat=",d"),
        yaxis=dict(title=None),
        height=max(400, len(contrib) * 40),
    )
    return fig


def plot_recombination_count_per_progeny(summary_df: pd.DataFrame) -> go.Figure:
    """
    Bar chart of HMM breakpoint counts per progeny.
    """
    if summary_df.empty:
        fig = go.Figure()
        _apply_style(fig, "Events per Progeny")
        return fig

    agg = summary_df.groupby("progeny")["hmm_breakpoints"].sum().reset_index()
    agg = agg.sort_values("hmm_breakpoints", ascending=True)

    fig = go.Figure(
        go.Bar(
            y=agg["progeny"],
            x=agg["hmm_breakpoints"],
            orientation="h",
            marker_color=PALETTE["accent"],
            text=agg["hmm_breakpoints"],
            textposition="outside",
            hovertemplate="<b>%{y}</b><br>Breakpoints: %{x}<extra></extra>",
        )
    )

    _apply_style(fig, "Events per Progeny")
    fig.update_layout(
        xaxis=dict(title="Number of Breakpoints"),
        yaxis=dict(title=None),
        height=max(400, len(agg) * 30),
    )
    return fig
