#!/usr/bin/env python3
"""Standalone trajectory plot for email attachment."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


POOL_COLORS = {
    "CEA5": "#0072B2",
    "H460 PCNB": "#E69F00",
    "TIM-3": "#009E73",
    "huAEC": "#CC79A7",
}

ROUND_ORDER = {"R1": 1, "R2": 2, "R3": 3, "R4": 4, "Other": 99}
ROUND_MARKERS = {"R2": "o", "R3": "^", "R4": "s"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--centroids-csv", required=True)
    p.add_argument("--out", required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    cent = pd.read_csv(args.centroids_csv)

    fig, ax = plt.subplots(figsize=(10, 8))

    for pool in sorted(cent["pool"].unique()):
        sub = cent[cent["pool"] == pool].copy()
        sub["r_ord"] = sub["round_label"].map(lambda x: ROUND_ORDER.get(x, 99))
        sub = sub.sort_values("r_ord")
        if len(sub) < 2:
            continue
        color = POOL_COLORS.get(pool, "#6b7280")

        # Trajectory line
        ax.plot(sub["centroid_pc1"], sub["centroid_pc2"],
                color=color, lw=2.5, alpha=0.9)

        # Points with round-specific markers
        for row in sub.itertuples(index=False):
            marker = ROUND_MARKERS.get(row.round_label, "D")
            ax.scatter(row.centroid_pc1, row.centroid_pc2,
                       s=120, color=color, marker=marker,
                       edgecolors="white", linewidths=0.8, zorder=5)
            ax.text(row.centroid_pc1, row.centroid_pc2,
                    f"  {row.round_label}", fontsize=9, ha="left", va="center",
                    fontweight="bold", color="#111827")

    # Target legend
    target_handles = [
        Line2D([0], [0], color=POOL_COLORS[p], lw=2.5, marker="o",
               markersize=8, label=p)
        for p in sorted(cent["pool"].unique())
    ]
    # Round legend
    round_handles = [
        Line2D([0], [0], color="#111827", lw=0, marker=m,
               markersize=8, label=r)
        for r, m in sorted(ROUND_MARKERS.items())
    ]

    leg1 = ax.legend(handles=target_handles, title="Target",
                     loc="upper right", fontsize=9, title_fontsize=9, frameon=True)
    ax.add_artist(leg1)
    ax.legend(handles=round_handles, title="Round",
              loc="lower right", fontsize=9, title_fontsize=9, frameon=True)

    ax.set_title("Target trajectories across panning rounds (weighted centroids)",
                 fontsize=13, fontweight="bold")
    ax.set_xlabel("Centroid PC1", fontsize=10)
    ax.set_ylabel("Centroid PC2", fontsize=10)
    ax.grid(alpha=0.15)
    fig.tight_layout()
    fig.savefig(args.out, dpi=300)
    plt.close(fig)
    print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()
