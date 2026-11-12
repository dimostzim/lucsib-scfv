#!/usr/bin/env python3
"""All-round/all-target sequence analysis and email-ready plots."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from sklearn.cluster import MiniBatchKMeans

from common import (
    AA_ALPHABET,
    average_linkage,
    ensure_dir,
    kmer_matrix,
    kmer_matrix_custom,
    pca_from_features,
    similarity_heatmap_label,
    similarity_heatmap_limits,
    topk_neighbors_cosine,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="All-sequence analysis across all mapped rounds/targets.")
    p.add_argument("--input-csv", default="r_output/tables/filt_df.csv", help="Input filt_df CSV")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--sequence-col", default="nt_scFv", help="Sequence column")
    p.add_argument(
        "--barcode-map",
        default="/homes/dtzim01/lucsib/barcode_map.csv",
        help="Barcode map CSV with columns Sample,Pool,Round,Barcode",
    )
    p.add_argument("--kmer", type=int, default=4, help="k-mer size for scalable full-set embedding")
    p.add_argument(
        "--alphabet-type",
        choices=["nt", "aa"],
        default="nt",
        help="Feature alphabet: nt for nucleotide, aa for amino-acid.",
    )
    p.add_argument("--n-clusters", type=int, default=24, help="MiniBatchKMeans clusters")
    p.add_argument("--seed", type=int, default=42, help="Random seed")
    return p.parse_args()


def _barcode_name_from_value(v: object) -> str:
    s = str(v).strip().lower()
    digits = "".join(ch for ch in s if ch.isdigit())
    if not digits:
        return s
    return f"barcode{int(digits):02d}"


def load_barcode_map(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"Sample", "Pool", "Round", "Barcode"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"barcode_map missing columns: {sorted(missing)}")

    out = df.copy()
    out["barcode"] = out["Barcode"].map(_barcode_name_from_value)
    out["sample_id"] = out["Sample"].map(lambda x: f"{int(x):02d}")
    out["round_label"] = out["Round"].map(lambda x: f"R{str(x).strip()}")
    out["sample_label"] = (
        out["sample_id"] + " " + out["Pool"].astype(str).str.strip() + " " + out["round_label"]
    )
    out = out.sort_values(["Pool", "Round", "Sample"]).reset_index(drop=True)
    return out[["Sample", "Pool", "Round", "barcode", "sample_id", "round_label", "sample_label"]]


def load_all_sequences(input_csv: Path, sequence_col: str, barcode_cols: list[str]) -> pd.DataFrame:
    usecols = [sequence_col] + barcode_cols
    df = pd.read_csv(input_csv, usecols=usecols)
    df = df[df[sequence_col].notna() & (df[sequence_col] != "")].copy()

    for c in barcode_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(np.int64)

    agg = (
        df.groupby(sequence_col, as_index=False)[barcode_cols]
        .sum()
        .sort_values(barcode_cols, ascending=False)
        .reset_index(drop=True)
    )
    agg["total_count"] = agg[barcode_cols].sum(axis=1)
    agg = agg[agg["total_count"] > 0].copy().reset_index(drop=True)
    agg["seq_id"] = [f"seq_{i + 1:06d}" for i in range(len(agg))]
    agg["seq_len"] = agg[sequence_col].str.len()
    return agg


def annotate_dominant_sample(agg: pd.DataFrame, meta: pd.DataFrame, barcode_cols: list[str]) -> pd.DataFrame:
    vals = agg[barcode_cols].to_numpy(dtype=np.int64)
    dom_idx = vals.argmax(axis=1)
    dom_count = vals[np.arange(vals.shape[0]), dom_idx]
    dom_barcode = np.array(barcode_cols, dtype=object)[dom_idx]

    map_pool = dict(zip(meta["barcode"], meta["Pool"]))
    map_round = dict(zip(meta["barcode"], meta["round_label"]))
    map_sample = dict(zip(meta["barcode"], meta["sample_label"]))

    agg = agg.copy()
    agg["dominant_barcode"] = dom_barcode
    agg["dominant_count"] = dom_count
    agg["dominant_pool"] = [map_pool.get(x, "Unknown") for x in dom_barcode]
    agg["dominant_round"] = [map_round.get(x, "Other") for x in dom_barcode]
    agg["dominant_sample"] = [map_sample.get(x, x) for x in dom_barcode]
    agg["n_positive_barcodes"] = (vals > 0).sum(axis=1)
    return agg


def _pool_colors(pools: list[str]) -> dict[str, tuple[float, float, float, float]]:
    fixed = {
        "CEA5": "#0072B2",
        "H460 PCNB": "#E69F00",
        "TIM-3": "#009E73",
        "huAEC": "#CC79A7",
    }
    palette = plt.get_cmap("tab10")
    out = {}
    extra_i = 0
    for pool in sorted(set(pools)):
        if pool in fixed:
            out[pool] = fixed[pool]
        else:
            out[pool] = palette(extra_i % 10)
            extra_i += 1
    return out


def _round_markers(rounds: list[str]) -> dict[str, str]:
    default = {"R1": "o", "R2": "o", "R3": "^", "R4": "s", "Other": "X"}
    out = {}
    for r in sorted(set(rounds)):
        out[r] = default.get(r, "D")
    return out


def _weighted_centroids_pca(
    pca: np.ndarray,
    agg: pd.DataFrame,
    meta: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for r in meta.itertuples(index=False):
        col = r.barcode
        w = agg[col].to_numpy(dtype=float)
        total = float(w.sum())
        if total <= 0:
            continue
        x = float(np.average(pca[:, 0], weights=w))
        y = float(np.average(pca[:, 1], weights=w))
        rows.append(
            {
                "barcode": r.barcode,
                "sample_label": r.sample_label,
                "pool": r.Pool,
                "round_label": r.round_label,
                "total_reads": int(total),
                "centroid_pc1": x,
                "centroid_pc2": y,
            }
        )
    return pd.DataFrame(rows)


def _weighted_centroids_feature(
    x: np.ndarray,
    agg: pd.DataFrame,
    meta: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for r in meta.itertuples(index=False):
        col = r.barcode
        w = agg[col].to_numpy(dtype=float)
        total = float(w.sum())
        if total <= 0:
            continue
        c = np.average(x, axis=0, weights=w)
        n = float(np.linalg.norm(c))
        if n > 0:
            c = c / n
        rows.append(
            {
                "barcode": r.barcode,
                "sample_label": r.sample_label,
                "pool": r.Pool,
                "round_label": r.round_label,
                "total_reads": int(total),
                "centroid_vec": c,
            }
        )
    out = pd.DataFrame(rows)
    if len(out) == 0:
        return out
    mat = np.vstack(out["centroid_vec"].to_numpy())
    out = out.drop(columns=["centroid_vec"]).reset_index(drop=True)
    out.attrs["centroid_matrix"] = mat
    return out


def _plot_pca_all(
    ax: plt.Axes,
    df: pd.DataFrame,
    var: np.ndarray,
    color_map: dict[str, tuple[float, float, float, float]],
    marker_map: dict[str, str],
) -> None:
    for pool in sorted(df["dominant_pool"].unique()):
        for rnd in sorted(df["dominant_round"].unique()):
            sub = df[(df["dominant_pool"] == pool) & (df["dominant_round"] == rnd)]
            if len(sub) == 0:
                continue
            ax.scatter(
                sub["PC1"],
                sub["PC2"],
                s=4,
                alpha=0.26,
                linewidths=0,
                color=color_map[pool],
                marker=marker_map[rnd],
                rasterized=True,
            )

    ax.set_title("A) All sequences PCA (color=target, marker=round)")
    ax.set_xlabel(f"PC1 ({var[0] * 100:.2f}% var)")
    ax.set_ylabel(f"PC2 ({var[1] * 100:.2f}% var)")
    ax.grid(alpha=0.15)

    pool_handles = [
        Line2D([0], [0], marker="o", linestyle="None", color=color_map[p], label=p, markersize=6)
        for p in sorted(color_map)
    ]
    round_handles = [
        Line2D([0], [0], marker=marker_map[r], linestyle="None", color="#111827", label=r, markersize=6)
        for r in sorted(marker_map)
    ]

    leg1 = ax.legend(handles=pool_handles, title="Target", loc="upper right", fontsize=8, title_fontsize=8, frameon=True)
    ax.add_artist(leg1)
    ax.legend(handles=round_handles, title="Round", loc="lower right", fontsize=8, title_fontsize=8, frameon=True)


def _plot_pca_r4(
    ax: plt.Axes,
    df: pd.DataFrame,
    var: np.ndarray,
    color_map: dict[str, tuple[float, float, float, float]],
) -> None:
    r4 = df[df["dominant_round"] == "R4"].copy()
    for pool in sorted(r4["dominant_pool"].unique()):
        sub = r4[r4["dominant_pool"] == pool]
        if len(sub) == 0:
            continue
        ax.scatter(
            sub["PC1"],
            sub["PC2"],
            s=6,
            alpha=0.35,
            linewidths=0,
            color=color_map[pool],
            marker="o",
            rasterized=True,
            label=pool,
        )
    ax.set_title("B) R4-only PCA view by target")
    ax.set_xlabel(f"PC1 ({var[0] * 100:.2f}% var)")
    ax.set_ylabel(f"PC2 ({var[1] * 100:.2f}% var)")
    ax.grid(alpha=0.15)
    ax.legend(loc="best", fontsize=8, frameon=True)


def _plot_trajectories(
    ax: plt.Axes,
    centroids: pd.DataFrame,
    color_map: dict[str, tuple[float, float, float, float]],
    marker_map: dict[str, str],
) -> None:
    round_order = {"R1": 1, "R2": 2, "R3": 3, "R4": 4, "Other": 99}

    for pool in sorted(centroids["pool"].unique()):
        sub = centroids[centroids["pool"] == pool].copy()
        sub["round_ord"] = sub["round_label"].map(lambda x: round_order.get(x, 99))
        sub = sub.sort_values("round_ord")
        if len(sub) == 0:
            continue

        # trajectory line
        ax.plot(
            sub["centroid_pc1"],
            sub["centroid_pc2"],
            color=color_map[pool],
            lw=2.0,
            alpha=0.85,
        )
        # points
        for row in sub.itertuples(index=False):
            ax.scatter(
                row.centroid_pc1,
                row.centroid_pc2,
                color=color_map[pool],
                marker=marker_map.get(row.round_label, "D"),
                s=80,
                edgecolor="white",
                linewidth=0.6,
                zorder=5,
            )
            ax.text(
                row.centroid_pc1,
                row.centroid_pc2,
                f" {row.round_label}",
                fontsize=8,
                va="center",
                ha="left",
                color="#111827",
            )

    ax.set_title("C) Target trajectories across rounds (weighted centroids)")
    ax.set_xlabel("Centroid PC1")
    ax.set_ylabel("Centroid PC2")
    ax.grid(alpha=0.15)


def _plot_sample_heatmap(
    ax: plt.Axes,
    sim: np.ndarray,
    labels: list[str],
    pools: list[str],
) -> None:
    dist = 1.0 - sim
    np.fill_diagonal(dist, 0.0)
    link = average_linkage(dist)
    order = link.leaf_order

    sim_ord = sim[np.ix_(order, order)]
    labels_ord = [labels[i] for i in order]
    pools_ord = [pools[i] for i in order]

    vmin, vmax = similarity_heatmap_limits(sim_ord)
    im = ax.imshow(sim_ord, cmap="inferno", vmin=vmin, vmax=vmax, interpolation="nearest", aspect="auto")
    ax.set_title("D) Hierarchical sample similarity heatmap")
    ax.set_xticks(np.arange(len(labels_ord)))
    ax.set_yticks(np.arange(len(labels_ord)))
    ax.set_xticklabels(labels_ord, rotation=90, fontsize=7)
    ax.set_yticklabels(labels_ord, fontsize=7)

    changes = np.where(np.array(pools_ord[:-1]) != np.array(pools_ord[1:]))[0]
    for c in changes:
        pos = c + 0.5
        ax.axhline(pos, color="white", lw=0.6, alpha=0.65)
        ax.axvline(pos, color="white", lw=0.6, alpha=0.65)

    return im, link, order


def main() -> None:
    args = parse_args()
    outdir = ensure_dir(Path(args.outdir))
    plot_dir = ensure_dir(outdir / "plots")

    map_path = Path(args.barcode_map)
    meta = load_barcode_map(map_path)
    barcode_cols = meta["barcode"].tolist()

    agg = load_all_sequences(
        input_csv=Path(args.input_csv),
        sequence_col=args.sequence_col,
        barcode_cols=barcode_cols,
    )
    agg = annotate_dominant_sample(agg=agg, meta=meta, barcode_cols=barcode_cols)

    seqs = agg[args.sequence_col].astype(str).tolist()
    if args.alphabet_type == "nt":
        x = kmer_matrix(seqs, k=args.kmer, l2_normalize=True)
    else:
        x = kmer_matrix_custom(seqs, k=args.kmer, alphabet=AA_ALPHABET, l2_normalize=True)
    pcs, var = pca_from_features(x, n_components=2)

    km = MiniBatchKMeans(
        n_clusters=args.n_clusters,
        random_state=args.seed,
        n_init=10,
        batch_size=4096,
        max_iter=300,
    )
    cluster_labels = km.fit_predict(x).astype(np.int32)

    out = agg[["seq_id", args.sequence_col, "seq_len", "total_count"] + barcode_cols].copy()
    out["dominant_barcode"] = agg["dominant_barcode"].to_numpy()
    out["dominant_pool"] = agg["dominant_pool"].to_numpy()
    out["dominant_round"] = agg["dominant_round"].to_numpy()
    out["dominant_sample"] = agg["dominant_sample"].to_numpy()
    out["n_positive_barcodes"] = agg["n_positive_barcodes"].to_numpy()
    out["cluster_id"] = cluster_labels
    out["PC1"] = pcs[:, 0]
    out["PC2"] = pcs[:, 1]

    print("Computing nearest-neighbor similarities...")
    _, nearest = topk_neighbors_cosine(x, top_k=1, chunk_size=512)
    out["nearest_similarity"] = nearest

    out.to_csv(outdir / "allrounds_sequence_assignments.csv", index=False)

    cluster_sizes = (
        out.groupby("cluster_id", as_index=False)
        .size()
        .rename(columns={"size": "cluster_size"})
        .sort_values("cluster_size", ascending=False)
        .reset_index(drop=True)
    )
    cluster_sizes.to_csv(outdir / "allrounds_cluster_sizes.csv", index=False)

    cluster_comp = (
        out.groupby(["cluster_id", "dominant_pool", "dominant_round"], as_index=False)
        .size()
        .rename(columns={"size": "n_sequences"})
        .sort_values(["cluster_id", "n_sequences"], ascending=[True, False])
        .reset_index(drop=True)
    )
    cluster_comp.to_csv(outdir / "allrounds_cluster_composition.csv", index=False)

    cent_pca = _weighted_centroids_pca(pca=pcs, agg=agg, meta=meta)
    cent_pca.to_csv(outdir / "sample_centroids_pca.csv", index=False)

    cent_feat = _weighted_centroids_feature(x=x, agg=agg, meta=meta)
    cent_feat.to_csv(outdir / "sample_centroids_feature_meta.csv", index=False)
    cent_mat = cent_feat.attrs["centroid_matrix"]
    sim = cent_mat @ cent_mat.T
    sim = np.clip(sim, 0.0, 1.0).astype(np.float32)
    pd.DataFrame(sim, index=cent_feat["sample_label"], columns=cent_feat["sample_label"]).to_csv(
        outdir / "sample_similarity_matrix.csv"
    )

    # Pairwise Jaccard overlap of unique sequences between samples.
    jacc_labels = cent_feat["sample_label"].tolist()
    jacc_barcodes = cent_feat["barcode"].tolist()
    n_samples = len(jacc_barcodes)
    jacc_mat = np.zeros((n_samples, n_samples), dtype=np.float64)
    presence = {}
    for i, bc in enumerate(jacc_barcodes):
        presence[i] = set(np.where(agg[bc].to_numpy() > 0)[0])
    for i in range(n_samples):
        jacc_mat[i, i] = 1.0
        for j in range(i + 1, n_samples):
            inter = len(presence[i] & presence[j])
            union = len(presence[i] | presence[j])
            jval = inter / union if union > 0 else 0.0
            jacc_mat[i, j] = jval
            jacc_mat[j, i] = jval
    pd.DataFrame(jacc_mat, index=jacc_labels, columns=jacc_labels).to_csv(
        outdir / "sample_jaccard_matrix.csv"
    )

    # Plot panel functions on individual files + 4-panel combined.
    color_map = _pool_colors(out["dominant_pool"].tolist())
    marker_map = _round_markers(out["dominant_round"].tolist())

    fig, ax = plt.subplots(figsize=(10, 8))
    _plot_pca_all(ax=ax, df=out, var=var, color_map=color_map, marker_map=marker_map)
    fig.tight_layout()
    fig.savefig(plot_dir / "plot_01_pca_all_target_round.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 8))
    _plot_pca_r4(ax=ax, df=out, var=var, color_map=color_map)
    fig.tight_layout()
    fig.savefig(plot_dir / "plot_02_pca_r4_targets.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 8))
    _plot_trajectories(ax=ax, centroids=cent_pca, color_map=color_map, marker_map=marker_map)
    fig.tight_layout()
    fig.savefig(plot_dir / "plot_03_target_round_trajectories.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 9))
    im, link, order = _plot_sample_heatmap(
        ax=ax,
        sim=sim,
        labels=cent_feat["sample_label"].tolist(),
        pools=cent_feat["pool"].tolist(),
    )
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(similarity_heatmap_label("Cosine similarity of weighted sample centroids"))
    fig.tight_layout()
    fig.savefig(plot_dir / "plot_04_sample_similarity_heatmap.png", dpi=300)
    plt.close(fig)

    fig, axs = plt.subplots(2, 2, figsize=(15, 12))
    _plot_pca_all(ax=axs[0, 0], df=out, var=var, color_map=color_map, marker_map=marker_map)
    _plot_pca_r4(ax=axs[0, 1], df=out, var=var, color_map=color_map)
    _plot_trajectories(ax=axs[1, 0], centroids=cent_pca, color_map=color_map, marker_map=marker_map)
    im2, _, _ = _plot_sample_heatmap(
        ax=axs[1, 1],
        sim=sim,
        labels=cent_feat["sample_label"].tolist(),
        pools=cent_feat["pool"].tolist(),
    )
    cbar2 = fig.colorbar(im2, ax=axs[1, 1], fraction=0.046, pad=0.04)
    cbar2.set_label(similarity_heatmap_label())
    fig.suptitle("All rounds/targets sequence analysis", y=0.995, fontsize=14)
    fig.tight_layout()
    fig.savefig(plot_dir / "plot_00_allrounds_summary_4panel.png", dpi=300)
    plt.close(fig)

    # Save linkage/order for reproducibility.
    pd.DataFrame(link.linkage, columns=["left", "right", "distance", "cluster_size"]).to_csv(
        outdir / "sample_similarity_linkage.csv", index=False
    )
    pd.DataFrame({"leaf_order_index": order, "sample_label": [cent_feat["sample_label"].tolist()[i] for i in order]}).to_csv(
        outdir / "sample_similarity_leaf_order.csv", index=False
    )

    summary = [
        "pipeline=allrounds_targets",
        f"input_csv={Path(args.input_csv).resolve()}",
        f"barcode_map={map_path.resolve()}",
        f"barcode_cols={','.join(barcode_cols)}",
        f"n_sequences={len(out)}",
        f"total_reads={int(out['total_count'].sum())}",
        f"kmer={args.kmer}",
        f"alphabet_type={args.alphabet_type}",
        f"n_clusters={args.n_clusters}",
        f"pca_var_pc1={float(var[0]):.6f}",
        f"pca_var_pc2={float(var[1]):.6f}",
        f"dominant_round_r2={int((out['dominant_round'] == 'R2').sum())}",
        f"dominant_round_r3={int((out['dominant_round'] == 'R3').sum())}",
        f"dominant_round_r4={int((out['dominant_round'] == 'R4').sum())}",
        f"outputs_dir={outdir.resolve()}",
    ]
    (outdir / "summary_allrounds.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")

    params = vars(args)
    params["effective_n_sequences"] = int(len(out))
    params["effective_barcodes"] = barcode_cols
    (outdir / "params_allrounds.json").write_text(json.dumps(params, indent=2) + "\n", encoding="utf-8")

    print(f"All-round/target analysis complete: {outdir}")


if __name__ == "__main__":
    main()
