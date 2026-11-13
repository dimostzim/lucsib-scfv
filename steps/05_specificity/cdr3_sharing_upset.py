#!/usr/bin/env python3
"""UpSet plot of the heavy and light CDR3 pairs shared between barcodes 07-12 (figure 4)."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import ticker
from upsetplot import UpSet, from_indicators

from common import load_samples

BARCODE_COLS = ["barcode07", "barcode08", "barcode09", "barcode10", "barcode11", "barcode12"]
LEFT_GROUP = ["barcode07", "barcode08", "barcode09"]
RIGHT_GROUP = ["barcode10", "barcode11", "barcode12"]

# Intersections with fewer CDR3 pairs are not shown, so the x-axis stays readable.
MIN_SUBSET_SIZE = 10


def _plot_matrix(self, ax):
    """UpSet.plot_matrix rewritten to work with pandas 3 Copy-on-Write."""
    ax = self._reorient(ax)
    data = self.intersections
    n_cats = data.index.nlevels
    inclusion = data.index.to_frame().values

    styles = [
        [
            self.subset_styles[i]
            if inclusion[i, j]
            else {"facecolor": self._other_dots_color, "linewidth": 0}
            for j in range(n_cats)
        ]
        for i in range(len(data))
    ]
    styles = sum(styles, [])
    style_columns = {
        "facecolor": "facecolors",
        "edgecolor": "edgecolors",
        "linewidth": "linewidths",
        "linestyle": "linestyles",
        "hatch": "hatch",
    }
    styles = (
        pd.DataFrame(styles)
        .reindex(columns=style_columns.keys())
        .astype(
            {
                "facecolor": "O",
                "edgecolor": "O",
                "linewidth": float,
                "linestyle": "O",
                "hatch": "O",
            }
        )
    )
    styles["linewidth"] = styles["linewidth"].fillna(1)
    styles["facecolor"] = styles["facecolor"].fillna(self._facecolor)
    styles["edgecolor"] = styles["edgecolor"].fillna(styles["facecolor"])
    styles["linestyle"] = styles["linestyle"].fillna("solid")
    del styles["hatch"]

    x = np.repeat(np.arange(len(data)), n_cats)
    y = np.tile(np.arange(n_cats), len(data))

    dot_size = (self._element_size * 0.35) ** 2 if self._element_size is not None else 200
    ax.scatter(
        *self._swapaxes(x, y),
        s=dot_size,
        zorder=10,
        **styles.rename(columns=style_columns),
    )

    if self._with_lines:
        idx = np.flatnonzero(inclusion)
        line_data = (
            pd.Series(y[idx], index=x[idx]).groupby(level=0).aggregate(["min", "max"])
        )
        colors = pd.Series(
            [
                style.get("edgecolor", style.get("facecolor", self._facecolor))
                for style in self.subset_styles
            ],
            name="color",
        )
        line_data = line_data.join(colors)
        ax.vlines(
            line_data.index.values,
            line_data["min"],
            line_data["max"],
            lw=2,
            colors=line_data["color"],
            zorder=5,
        )

    tick_axis = ax.yaxis
    tick_axis.set_ticks(np.arange(n_cats))
    tick_axis.set_ticklabels(data.index.names, rotation=0 if self._horizontal else -90)
    ax.xaxis.set_visible(False)
    ax.tick_params(axis="both", which="both", length=0)
    if not self._horizontal:
        ax.yaxis.set_ticks_position("top")
    ax.set_frame_on(False)
    ax.set_xlim(-0.5, x[-1] + 0.5, auto=False)
    ax.grid(False)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--cdr3-csv", type=Path, required=True)
    p.add_argument("--samples-csv", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    args = p.parse_args()

    UpSet.plot_matrix = _plot_matrix
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "xtick.labelsize": 9,
            "ytick.labelsize": 10,
            "axes.linewidth": 0.8,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )

    samples = load_samples(args.samples_csv)
    labels = dict(zip(samples["barcode"], samples["sample_label"]))
    presence = pd.read_csv(args.cdr3_csv)[BARCODE_COLS] > 0
    presence = presence[presence.any(axis=1)].rename(columns=labels)

    upset = UpSet(
        from_indicators(presence),
        subset_size="count",
        sort_by="cardinality",
        sort_categories_by="input",
        min_subset_size=MIN_SUBSET_SIZE,
        facecolor="#111827",
        other_dots_color=0.3,
        shading_color=0.08,
        with_lines=True,
        element_size=42,
        intersection_plot_elements=10,
        totals_plot_elements=3,
        show_counts=False,
        show_percentages=False,
    )

    # Highlight any intersection containing at least one barcode from each group.
    first_highlight = True
    for left in LEFT_GROUP:
        for right in RIGHT_GROUP:
            upset.style_subsets(
                present=[labels[left], labels[right]],
                facecolor="#0f766e",
                edgecolor="#0f766e",
                linewidth=0.8,
                label="Cross-group subset (07-09 and 10-12)" if first_highlight else None,
            )
            first_highlight = False

    fig = plt.figure(figsize=(15, 9))
    axes = upset.plot(fig=fig)
    axes["intersections"].set_ylabel("Intersection size (unique CDR3 pairs)")
    axes["intersections"].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:,.0f}"))
    axes["intersections"].spines["top"].set_visible(False)
    axes["intersections"].spines["right"].set_visible(False)
    axes["totals"].set_xlabel("Per-barcode total")
    axes["totals"].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:,.0f}"))
    axes["totals"].spines["top"].set_visible(False)
    axes["totals"].spines["right"].set_visible(False)

    fig.suptitle("Shared scFv CDR3 pairs across barcodes 07-12", fontsize=16, weight="bold", y=0.98)
    fig.text(
        0.5,
        0.945,
        f"UpSet plot with subset size >= {MIN_SUBSET_SIZE} (cross-group subsets highlighted in teal)",
        ha="center",
        va="center",
        fontsize=10,
        color="#374151",
    )
    fig.subplots_adjust(left=0.08, right=0.98, top=0.90, bottom=0.08)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out_dir / "upset_CDR3s_barcodes7_12.png", dpi=600, bbox_inches="tight")


if __name__ == "__main__":
    main()
