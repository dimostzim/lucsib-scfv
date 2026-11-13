#!/usr/bin/env python3
"""Sequence clusters over all rounds and targets, read-weighted sample centroids on the PCA map and sample overlap."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import count_cluster_sizes, find_clusters, kmer_matrix, load_samples, nearest_similarity, pca_from_features

K = 4
N_CLUSTERS = 24
BATCH_SIZE = 4096


def load_all_sequences(input_csv: Path, barcodes: list[str]) -> pd.DataFrame:
    """Reads per unique sequence and barcode, for sequences with reads in at least one of the barcodes."""
    df = pd.read_csv(input_csv, usecols=["nt_scFv"] + barcodes)
    df = df[df["nt_scFv"].notna() & (df["nt_scFv"] != "")]
    agg = (
        df.groupby("nt_scFv", as_index=False)[barcodes]
        .sum()
        .sort_values(barcodes, ascending=False)
        .reset_index(drop=True)
    )
    agg["total_count"] = agg[barcodes].sum(axis=1)
    agg = agg[agg["total_count"] > 0].reset_index(drop=True)
    agg["seq_id"] = [f"seq_{i + 1:06d}" for i in range(len(agg))]
    agg["seq_len"] = agg["nt_scFv"].str.len()
    return agg


def weighted_centroids(pcs: np.ndarray, seqs: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    """Read-weighted mean PC1 and PC2 of each sample with reads."""
    rows = []
    for s in samples.itertuples(index=False):
        w = seqs[s.barcode].to_numpy(dtype=float)
        total = float(w.sum())
        if total <= 0:
            continue
        rows.append({
            "barcode": s.barcode,
            "sample_label": s.sample_label,
            "pool": s.Pool,
            "round_label": s.round_label,
            "total_reads": int(total),
            "centroid_pc1": float(np.average(pcs[:, 0], weights=w)),
            "centroid_pc2": float(np.average(pcs[:, 1], weights=w)),
        })
    return pd.DataFrame(rows)


def jaccard_matrix(seqs: pd.DataFrame, barcodes: list[str]) -> np.ndarray:
    """Shared sequences divided by the sequences in either sample, for each pair of samples."""
    present = (seqs[barcodes] > 0).to_numpy(dtype=np.int64)
    shared = present.T @ present
    n = np.diag(shared)
    return shared / (n[:, None] + n[None, :] - shared)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input-csv", type=Path, required=True)
    p.add_argument("--samples-csv", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    args = p.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    samples = load_samples(args.samples_csv)
    barcodes = samples["barcode"].tolist()
    seqs = load_all_sequences(args.input_csv, barcodes)

    x = kmer_matrix(seqs["nt_scFv"].tolist(), K)
    pcs, var = pca_from_features(x)

    counts = seqs[barcodes].to_numpy(dtype=np.int64)
    dominant = samples.iloc[counts.argmax(axis=1)]
    out = seqs[["seq_id", "nt_scFv", "seq_len", "total_count"] + barcodes].copy()
    out["dominant_barcode"] = dominant["barcode"].to_numpy()
    out["dominant_pool"] = dominant["Pool"].to_numpy()
    out["dominant_round"] = dominant["round_label"].to_numpy()
    out["dominant_sample"] = dominant["sample_label"].to_numpy()
    out["n_positive_barcodes"] = (counts > 0).sum(axis=1)
    out["cluster_id"] = find_clusters(x, N_CLUSTERS, BATCH_SIZE)
    out["PC1"] = pcs[:, 0]
    out["PC2"] = pcs[:, 1]
    out["nearest_similarity"] = nearest_similarity(x)
    out.to_csv(args.outdir / "allrounds_sequence_assignments.csv", index=False)

    count_cluster_sizes(out["cluster_id"].to_numpy()).to_csv(args.outdir / "allrounds_cluster_sizes.csv", index=False)
    (
        out.groupby(["cluster_id", "dominant_pool", "dominant_round"], as_index=False)
        .size()
        .rename(columns={"size": "n_sequences"})
        .sort_values(["cluster_id", "n_sequences"], ascending=[True, False])
        .reset_index(drop=True)
        .to_csv(args.outdir / "allrounds_cluster_composition.csv", index=False)
    )

    centroids = weighted_centroids(pcs, seqs, samples)
    centroids.to_csv(args.outdir / "sample_centroids_pca.csv", index=False)
    labels = centroids["sample_label"].tolist()
    pd.DataFrame(jaccard_matrix(seqs, centroids["barcode"].tolist()), index=labels, columns=labels).to_csv(
        args.outdir / "sample_jaccard_matrix.csv"
    )

    summary = [
        "pipeline=allrounds_targets",
        f"barcode_cols={','.join(barcodes)}",
        f"n_sequences={len(out)}",
        f"total_reads={int(out['total_count'].sum())}",
        f"kmer={K}",
        "alphabet_type=nt",
        f"n_clusters={N_CLUSTERS}",
        f"pca_var_pc1={float(var[0]):.6f}",
        f"pca_var_pc2={float(var[1]):.6f}",
    ]
    summary += [f"dominant_round_{r.lower()}={int((out['dominant_round'] == r).sum())}" for r in ("R2", "R3", "R4")]
    (args.outdir / "summary_allrounds.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
