#!/usr/bin/env python3
"""Generate two parallel 4-panel figures for email reply.

Figure 1: R4 H460 PCNB-specific candidates (8,878 sequences)
Figure 2: All rounds / all targets (182,345 sequences)

Both figures have the same panel layout:
  A) PCA colored by cluster
  B) Nearest-neighbor similarity ECDF
  C) Cluster size rank
  D) Binned similarity heatmap with cluster labels
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from common import kmer_matrix, similarity_heatmap_label, similarity_heatmap_limits


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build two parallel 4-panel email figures.")
    p.add_argument("--allseq-dir", required=True, help="Path to allseq_email run dir")
    p.add_argument("--allrounds-dir", required=True, help="Path to allrounds run dir")
    p.add_argument("--outdir", required=True, help="Output directory for figures")
    return p.parse_args()


POOL_COLORS = {
    "CEA5": "#0072B2",
    "H460 PCNB": "#E69F00",
    "TIM-3": "#009E73",
    "huAEC": "#CC79A7",
}


def ecdf(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.sort(np.asarray(values, dtype=np.float64))
    if x.size == 0:
        return x, x
    y = np.arange(1, x.size + 1, dtype=np.float64) / x.size
    return x, y


def _cluster_colors(labels: np.ndarray) -> np.ndarray:
    uniq = sorted(np.unique(labels).tolist())
    cmap = plt.get_cmap("tab20")
    color_map = {u: cmap(i % 20) for i, u in enumerate(uniq)}
    return np.array([color_map[x] for x in labels], dtype=object)


def _add_cluster_labels(ax: plt.Axes, bin_labels: np.ndarray) -> None:
    if bin_labels.size < 2:
        return
    changes = np.where(np.diff(bin_labels) != 0)[0]
    for c in changes:
        pos = c + 0.5
        ax.axhline(pos, color="white", lw=0.45, alpha=0.5)
        ax.axvline(pos, color="white", lw=0.45, alpha=0.5)
    boundaries = np.concatenate([[-1], changes, [bin_labels.size - 1]])
    tick_pos = []
    tick_labels = []
    for i in range(len(boundaries) - 1):
        start = boundaries[i] + 1
        end = boundaries[i + 1]
        mid = (start + end) / 2.0
        tick_pos.append(mid)
        tick_labels.append(str(bin_labels[int(start)]))
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_labels, fontsize=5, rotation=0)
    ax.set_yticks(tick_pos)
    ax.set_yticklabels(tick_labels, fontsize=5)


def _binned_heatmap_from_assignments(
    assign: pd.DataFrame, n_bins: int, kmer_func, seq_col: str, k: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Recompute binned similarity heatmap from assignment data."""
    seqs = assign[seq_col].astype(str).tolist()
    x = kmer_func(seqs, k=k)
    labels = assign["cluster_id"].to_numpy()
    pc1 = assign["PC1"].to_numpy()

    n = x.shape[0]
    n_bins = max(10, min(n_bins, n))
    order = np.lexsort((pc1, labels))
    chunks = np.array_split(order, n_bins)

    bin_vecs = []
    bin_major_labels = []
    for idx in chunks:
        if idx.size == 0:
            continue
        m = x[idx].mean(axis=0)
        norm = np.linalg.norm(m)
        if norm > 0:
            m = m / norm
        bin_vecs.append(m)
        lab = labels[idx]
        binc = np.bincount(lab.astype(int), minlength=int(labels.max()) + 1)
        bin_major_labels.append(int(np.argmax(binc)))

    b = np.vstack(bin_vecs).astype(np.float32)
    sim = b @ b.T
    sim = np.clip(sim, 0.0, 1.0)
    np.fill_diagonal(sim, 1.0)
    return sim, np.array(bin_major_labels, dtype=np.int32)


def panel_pca(ax, assign, sizes, colors, title, var_pc1=None, var_pc2=None):
    ax.scatter(assign["PC1"], assign["PC2"], c=colors, s=sizes,
               alpha=0.6, edgecolors="none", rasterized=True)
    xlabel = "PC1" + (f" ({var_pc1*100:.1f}% var)" if var_pc1 else "")
    ylabel = "PC2" + (f" ({var_pc2*100:.1f}% var)" if var_pc2 else "")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.grid(alpha=0.15)


def panel_pca_targets(ax, assign, title):
    """PCA panel colored by target with legend and centroid labels."""
    for pool in sorted(assign["dominant_pool"].unique()):
        sub = assign[assign["dominant_pool"] == pool]
        color = POOL_COLORS.get(pool, "#6b7280")
        ax.scatter(
            sub["PC1"], sub["PC2"],
            s=3, alpha=0.25, linewidths=0,
            color=color,
            label=f"{pool} (n={len(sub):,})",
            rasterized=True,
        )
        # Weighted centroid with label
        weights = sub["total_count"].to_numpy(dtype=float)
        cx = float(np.average(sub["PC1"], weights=weights))
        cy = float(np.average(sub["PC2"], weights=weights))
        ax.scatter([cx], [cy], s=100, color=color,
                   edgecolors="white", linewidths=0.8, zorder=6)
        ax.text(cx, cy, f"  {pool}", fontsize=7, ha="left", va="center",
                fontweight="bold", color="#111827", zorder=7)

    ax.set_title(title, fontsize=10)
    ax.set_xlabel("PC1", fontsize=8)
    ax.set_ylabel("PC2", fontsize=8)
    ax.grid(alpha=0.15)
    ax.legend(loc="best", fontsize=7, frameon=True, markerscale=4)


def panel_ecdf(ax, nearest, title):
    finite = nearest[np.isfinite(nearest)]
    x, y = ecdf(finite)
    ax.plot(x, y, lw=2, color="#be123c")
    med = float(np.median(finite))
    ax.axvline(med, color="#6b7280", ls="--", lw=0.8, alpha=0.6)
    ax.text(med, 0.05, f" median={med:.3f}", fontsize=7, color="#6b7280")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Nearest-neighbor cosine similarity", fontsize=8)
    ax.set_ylabel("Cumulative fraction", fontsize=8)
    ax.grid(alpha=0.2)


def panel_rank(ax, cluster_sizes, title):
    vals = cluster_sizes["cluster_size"].sort_values(ascending=False).to_numpy()
    ranks = np.arange(1, len(vals) + 1)
    ax.plot(ranks, vals, marker="o", ms=3, lw=1.3, color="#065f46")
    ax.set_yscale("log")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Cluster rank", fontsize=8)
    ax.set_ylabel("Cluster size (log scale)", fontsize=8)
    ax.grid(alpha=0.25)


def panel_heatmap(ax, sim, bin_labels, title):
    vmin, vmax = similarity_heatmap_limits(sim)
    im = ax.imshow(sim, cmap="inferno", vmin=vmin, vmax=vmax,
                   interpolation="nearest", aspect="auto")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Cluster ID", fontsize=8)
    ax.set_ylabel("Cluster ID", fontsize=8)
    _add_cluster_labels(ax, bin_labels)
    return im


def make_figure(
    assign: pd.DataFrame,
    cluster_sizes: pd.DataFrame,
    sim: np.ndarray,
    bin_labels: np.ndarray,
    sizes: np.ndarray,
    colors: np.ndarray,
    suptitle: str,
    out_path: Path,
    var_pc1: float | None = None,
    var_pc2: float | None = None,
    pca_by_target: bool = False,
) -> None:
    fig, axs = plt.subplots(2, 2, figsize=(15, 12))

    if pca_by_target:
        panel_pca_targets(axs[0, 0], assign,
                          f"A) PCA colored by target (n={len(assign):,})")
    else:
        panel_pca(axs[0, 0], assign, sizes, colors,
                  f"A) PCA colored by cluster (n={len(assign):,})",
                  var_pc1=var_pc1, var_pc2=var_pc2)
    panel_ecdf(axs[0, 1], assign["nearest_similarity"].to_numpy(),
               "B) Nearest-neighbor similarity ECDF")
    panel_rank(axs[1, 0], cluster_sizes,
               f"C) Cluster size rank ({len(cluster_sizes)} clusters)")
    im = panel_heatmap(axs[1, 1], sim, bin_labels,
                       "D) Similarity heatmap (binned, cluster-labeled)")
    cbar = fig.colorbar(im, ax=axs[1, 1], fraction=0.046, pad=0.04)
    cbar.set_label(similarity_heatmap_label(), fontsize=8)

    fig.suptitle(suptitle, y=0.995, fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Wrote: {out_path}")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    allseq_dir = Path(args.allseq_dir)
    allrounds_dir = Path(args.allrounds_dir)

    # --- Figure 1: R4 H460 PCNB-specific ---
    allseq = pd.read_csv(allseq_dir / "allseq_assignments.csv")
    allseq_sizes = pd.read_csv(allseq_dir / "allseq_cluster_sizes.csv")
    allseq_sim = pd.read_csv(allseq_dir / "allseq_similarity_heatmap_binned.csv")
    sim1 = allseq_sim.to_numpy(dtype=float)

    # Reconstruct bin_labels
    labels1 = allseq["cluster_id"].to_numpy()
    pc1_1 = allseq["PC1"].to_numpy()
    order1 = np.lexsort((pc1_1, labels1))
    chunks1 = np.array_split(order1, sim1.shape[0])
    bl1 = []
    for idx in chunks1:
        if idx.size == 0:
            continue
        lab = labels1[idx]
        binc = np.bincount(lab.astype(int), minlength=int(labels1.max()) + 1)
        bl1.append(int(np.argmax(binc)))
    bin_labels1 = np.array(bl1, dtype=np.int32)

    count_col = "barcode09" if "barcode09" in allseq.columns else "total_count"
    counts1 = allseq[count_col].to_numpy(dtype=float)
    sizes1 = 6 + 10 * np.log10(counts1 + 1.0)
    colors1 = _cluster_colors(labels1)

    make_figure(
        assign=allseq,
        cluster_sizes=allseq_sizes,
        sim=sim1,
        bin_labels=bin_labels1,
        sizes=sizes1,
        colors=colors1,
        suptitle="Figure 1: R4 H460 PCNB-specific candidates (n=8,878)",
        out_path=outdir / "figure1_r4_h460pcnb_specific.png",
    )

    # --- Figure 2: All rounds / all targets ---
    allrounds = pd.read_csv(allrounds_dir / "allrounds_sequence_assignments.csv")
    allrounds_sizes = pd.read_csv(allrounds_dir / "allrounds_cluster_sizes.csv")

    # Color by target instead of cluster
    pools = allrounds["dominant_pool"].to_numpy()
    colors2 = np.array([POOL_COLORS.get(p, "#6b7280") for p in pools], dtype=object)
    counts2 = allrounds["total_count"].to_numpy(dtype=float)
    sizes2 = 3 + 6 * np.log10(counts2 + 1.0)

    # Recompute binned heatmap for allrounds
    sim2, bin_labels2 = _binned_heatmap_from_assignments(
        assign=allrounds,
        n_bins=140,
        kmer_func=lambda seqs, k: kmer_matrix(seqs, k=k, l2_normalize=True),
        seq_col="nt_scFv",
        k=4,
    )

    make_figure(
        assign=allrounds,
        cluster_sizes=allrounds_sizes,
        sim=sim2,
        bin_labels=bin_labels2,
        sizes=sizes2,
        colors=colors2,
        suptitle="Figure 2: All rounds / all targets (n=182,345; R2-R4)",
        out_path=outdir / "figure2_allrounds_alltargets.png",
        pca_by_target=True,
    )


if __name__ == "__main__":
    main()
