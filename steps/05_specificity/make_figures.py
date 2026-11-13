#!/usr/bin/env python3
"""Figures 1 and 2 (PCA, nearest-neighbour similarity, cluster sizes, similarity heatmap); figure 3 (trajectories)."""

from __future__ import annotations

import argparse
from collections.abc import Callable
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import all_rounds_map
import r4_h460pcnb_families
from common import bin_by_cluster, binned_similarity, kmer_matrix

POOL_COLORS = {
    "CEA5": "#0072B2",
    "H460 PCNB": "#E69F00",
    "TIM-3": "#009E73",
    "huAEC": "#CC79A7",
}
ROUND_MARKERS = {"R2": "o", "R3": "^", "R4": "s"}
# Colour range of the heatmaps: these quantiles of the off-diagonal similarities, widened to at least MIN_SPAN.
LOWER_Q = 0.02
UPPER_Q = 0.995
MIN_SPAN = 0.04


def _heatmap_limits(sim: np.ndarray) -> tuple[float, float]:
    vals = np.asarray(sim, dtype=np.float64)[~np.eye(len(sim), dtype=bool)]
    lo = max(0.0, min(float(np.quantile(vals, LOWER_Q)), 1.0))
    hi = max(lo + 1e-6, min(1.0, float(np.quantile(vals, UPPER_Q))))
    if hi - lo < MIN_SPAN:
        center = float(np.median(vals))
        half_span = max(MIN_SPAN / 2.0, (float(vals.max()) - float(vals.min())) / 2.0, 1e-6)
        lo = max(0.0, center - half_span)
        hi = min(1.0, center + half_span)
        if hi - lo < MIN_SPAN:
            if hi >= 1.0:
                lo, hi = max(0.0, 1.0 - MIN_SPAN), 1.0
            else:
                hi = min(1.0, lo + MIN_SPAN)
    return lo, hi


def _cluster_colors(labels: np.ndarray) -> np.ndarray:
    cmap = plt.get_cmap("tab20")
    color_map = {u: cmap(i % 20) for i, u in enumerate(sorted(np.unique(labels).tolist()))}
    return np.array([color_map[x] for x in labels], dtype=object)


def _ecdf(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.sort(np.asarray(values, dtype=np.float64))
    return x, np.arange(1, x.size + 1, dtype=np.float64) / x.size


def _style_pca_panel(ax: plt.Axes, title: str) -> None:
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("PC1", fontsize=8)
    ax.set_ylabel("PC2", fontsize=8)
    ax.grid(alpha=0.15)


def panel_pca_by_cluster(ax: plt.Axes, assign: pd.DataFrame) -> None:
    sizes = 6 + 10 * np.log10(assign[r4_h460pcnb_families.PRESENT].to_numpy(dtype=float) + 1.0)
    ax.scatter(assign["PC1"], assign["PC2"], c=_cluster_colors(assign["cluster_id"].to_numpy()), s=sizes,
               alpha=0.6, edgecolors="none", rasterized=True)
    _style_pca_panel(ax, f"A) PCA colored by cluster (n={len(assign):,})")


def panel_pca_by_target(ax: plt.Axes, assign: pd.DataFrame) -> None:
    """Sequences coloured by the target of their dominant sample, with each target's read-weighted centroid."""
    for pool in sorted(assign["dominant_pool"].unique()):
        sub = assign[assign["dominant_pool"] == pool]
        color = POOL_COLORS[pool]
        ax.scatter(sub["PC1"], sub["PC2"], s=3, alpha=0.25, linewidths=0, color=color,
                   label=f"{pool} (n={len(sub):,})", rasterized=True)
        weights = sub["total_count"].to_numpy(dtype=float)
        cx = float(np.average(sub["PC1"], weights=weights))
        cy = float(np.average(sub["PC2"], weights=weights))
        ax.scatter([cx], [cy], s=100, color=color, edgecolors="white", linewidths=0.8, zorder=6)
        ax.text(cx, cy, f"  {pool}", fontsize=7, ha="left", va="center",
                fontweight="bold", color="#111827", zorder=7)
    _style_pca_panel(ax, f"A) PCA colored by target (n={len(assign):,})")
    ax.legend(loc="best", fontsize=7, frameon=True, markerscale=4)


def panel_ecdf(ax: plt.Axes, nearest: np.ndarray) -> None:
    x, y = _ecdf(nearest)
    ax.plot(x, y, lw=2, color="#be123c")
    med = float(np.median(nearest))
    ax.axvline(med, color="#6b7280", ls="--", lw=0.8, alpha=0.6)
    ax.text(med, 0.05, f" median={med:.3f}", fontsize=7, color="#6b7280")
    ax.set_title("B) Nearest-neighbor similarity ECDF", fontsize=10)
    ax.set_xlabel("Nearest-neighbor cosine similarity", fontsize=8)
    ax.set_ylabel("Cumulative fraction", fontsize=8)
    ax.grid(alpha=0.2)


def panel_rank(ax: plt.Axes, cluster_sizes: pd.DataFrame) -> None:
    vals = cluster_sizes["cluster_size"].sort_values(ascending=False).to_numpy()
    ax.plot(np.arange(1, len(vals) + 1), vals, marker="o", ms=3, lw=1.3, color="#065f46")
    ax.set_yscale("log")
    ax.set_title(f"C) Cluster size rank ({len(cluster_sizes)} clusters)", fontsize=10)
    ax.set_xlabel("Cluster rank", fontsize=8)
    ax.set_ylabel("Cluster size (log scale)", fontsize=8)
    ax.grid(alpha=0.25)


def panel_heatmap(ax: plt.Axes, sim: np.ndarray, bin_clusters: np.ndarray) -> matplotlib.image.AxesImage:
    """Binned similarity heatmap with white lines between clusters and cluster IDs as tick labels."""
    vmin, vmax = _heatmap_limits(sim)
    im = ax.imshow(sim, cmap="inferno", vmin=vmin, vmax=vmax, interpolation="nearest", aspect="auto")
    ax.set_title("D) Similarity heatmap (binned, cluster-labeled)", fontsize=10)
    ax.set_xlabel("Cluster ID", fontsize=8)
    ax.set_ylabel("Cluster ID", fontsize=8)
    changes = np.where(np.diff(bin_clusters) != 0)[0]
    for c in changes:
        ax.axhline(c + 0.5, color="white", lw=0.45, alpha=0.5)
        ax.axvline(c + 0.5, color="white", lw=0.45, alpha=0.5)
    starts = np.concatenate([[0], changes + 1])
    ends = np.concatenate([changes, [bin_clusters.size - 1]])
    tick_pos = list((starts + ends) / 2.0)
    tick_labels = [str(bin_clusters[s]) for s in starts]
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_labels, fontsize=5, rotation=0)
    ax.set_yticks(tick_pos)
    ax.set_yticklabels(tick_labels, fontsize=5)
    return im


def make_cluster_figure(
    panel_a: Callable[[plt.Axes, pd.DataFrame], None],
    assign: pd.DataFrame,
    cluster_sizes: pd.DataFrame,
    sim: np.ndarray,
    bin_clusters: np.ndarray,
    suptitle: str,
    out_path: Path,
) -> None:
    fig, axs = plt.subplots(2, 2, figsize=(15, 12))
    panel_a(axs[0, 0], assign)
    panel_ecdf(axs[0, 1], assign["nearest_similarity"].to_numpy())
    panel_rank(axs[1, 0], cluster_sizes)
    im = panel_heatmap(axs[1, 1], sim, bin_clusters)
    cbar = fig.colorbar(im, ax=axs[1, 1], fraction=0.046, pad=0.04)
    label = f"Cosine similarity (contrast-scaled; off-diagonal p{LOWER_Q * 100:g}-p{UPPER_Q * 100:g})"
    cbar.set_label(label, fontsize=8)
    fig.suptitle(suptitle, y=0.995, fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def make_trajectory_figure(cent: pd.DataFrame, out_path: Path) -> None:
    """Read-weighted centroid of each sample, joined in round order for each target."""
    fig, ax = plt.subplots(figsize=(10, 8))
    pools = sorted(cent["pool"].unique())
    for pool in pools:
        sub = cent[cent["pool"] == pool].sort_values("round_label", key=lambda r: r.str[1:].astype(int))
        if len(sub) < 2:
            continue
        color = POOL_COLORS[pool]
        ax.plot(sub["centroid_pc1"], sub["centroid_pc2"], color=color, lw=2.5, alpha=0.9)
        for row in sub.itertuples(index=False):
            ax.scatter(row.centroid_pc1, row.centroid_pc2, s=120, color=color, marker=ROUND_MARKERS[row.round_label],
                       edgecolors="white", linewidths=0.8, zorder=5)
            ax.text(row.centroid_pc1, row.centroid_pc2, f"  {row.round_label}", fontsize=9, ha="left", va="center",
                    fontweight="bold", color="#111827")

    target_handles = [Line2D([0], [0], color=POOL_COLORS[p], lw=2.5, marker="o", markersize=8, label=p) for p in pools]
    round_handles = [Line2D([0], [0], color="#111827", lw=0, marker=m, markersize=8, label=r)
                     for r, m in sorted(ROUND_MARKERS.items())]
    ax.add_artist(ax.legend(handles=target_handles, title="Target", loc="upper right", fontsize=9, title_fontsize=9,
                            frameon=True))
    ax.legend(handles=round_handles, title="Round", loc="lower right", fontsize=9, title_fontsize=9, frameon=True)

    ax.set_title("Target trajectories across panning rounds (weighted centroids)", fontsize=13, fontweight="bold")
    ax.set_xlabel("Centroid PC1", fontsize=10)
    ax.set_ylabel("Centroid PC2", fontsize=10)
    ax.grid(alpha=0.15)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--r4-dir", type=Path, required=True)
    p.add_argument("--all-rounds-dir", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    args = p.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    r4 = pd.read_csv(args.r4_dir / "allseq_assignments.csv")
    _, r4_bin_clusters = bin_by_cluster(r4["cluster_id"].to_numpy(), r4["PC1"].to_numpy())
    make_cluster_figure(
        panel_pca_by_cluster,
        r4,
        pd.read_csv(args.r4_dir / "allseq_cluster_sizes.csv"),
        pd.read_csv(args.r4_dir / "allseq_similarity_heatmap_binned.csv").to_numpy(dtype=float),
        r4_bin_clusters,
        f"Figure 1: R4 H460 PCNB-specific candidates (n={len(r4):,})",
        args.outdir / "figure1_r4_h460pcnb_specific.png",
    )

    allrounds = pd.read_csv(args.all_rounds_dir / "allrounds_sequence_assignments.csv")
    bins, bin_clusters = bin_by_cluster(allrounds["cluster_id"].to_numpy(), allrounds["PC1"].to_numpy())
    make_cluster_figure(
        panel_pca_by_target,
        allrounds,
        pd.read_csv(args.all_rounds_dir / "allrounds_cluster_sizes.csv"),
        binned_similarity(kmer_matrix(allrounds["nt_scFv"].tolist(), all_rounds_map.K), bins),
        bin_clusters,
        f"Figure 2: All rounds / all targets (n={len(allrounds):,}; R2-R4)",
        args.outdir / "figure2_allrounds_alltargets.png",
    )

    make_trajectory_figure(pd.read_csv(args.all_rounds_dir / "sample_centroids_pca.csv"),
                           args.outdir / "figure3_target_trajectories.png")


if __name__ == "__main__":
    main()
