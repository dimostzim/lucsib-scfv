#!/usr/bin/env python3
"""Clusters of the scFv DNA sequences with reads in H460 PCNB R4 (barcode09) and none in huAEC R4 (barcode12)."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import (
    bin_by_cluster,
    binned_similarity,
    count_cluster_sizes,
    find_clusters,
    kmer_matrix,
    nearest_similarity,
    pca_from_features,
)

PRESENT = "barcode09"
ABSENT = "barcode12"
K = 5
N_CLUSTERS = 20
BATCH_SIZE = 2048


def load_candidates(input_csv: Path) -> pd.DataFrame:
    """Unique sequences with reads in PRESENT and none in ABSENT, most reads first."""
    df = pd.read_csv(input_csv, usecols=["nt_scFv", PRESENT, ABSENT])
    df = df[df["nt_scFv"].notna() & (df["nt_scFv"] != "")]
    agg = (
        df.groupby("nt_scFv", as_index=False)[[PRESENT, ABSENT]]
        .sum()
        .sort_values([PRESENT, "nt_scFv"], ascending=[False, True])
        .reset_index(drop=True)
    )
    cand = agg[(agg[PRESENT] > 0) & (agg[ABSENT] == 0)].reset_index(drop=True)
    cand["seq_id"] = [f"seq_{i + 1:05d}" for i in range(len(cand))]
    cand["seq_len"] = cand["nt_scFv"].str.len()
    return cand


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input-csv", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    args = p.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    cand = load_candidates(args.input_csv)
    counts = cand[PRESENT].to_numpy()
    lengths = cand["seq_len"].to_numpy()
    candidates_summary = [
        f"filter={PRESENT}>0 and {ABSENT}==0",
        f"n_candidates={len(cand)}",
        f"total_present_reads={int(counts.sum())}",
        f"median_present_reads={float(np.median(counts)):.3f}",
        f"max_present_reads={int(counts.max())}",
        f"len_min={int(lengths.min())}",
        f"len_median={float(np.median(lengths)):.3f}",
        f"len_max={int(lengths.max())}",
    ]
    (args.outdir / "candidates_summary.txt").write_text("\n".join(candidates_summary) + "\n", encoding="utf-8")

    x = kmer_matrix(cand["nt_scFv"].tolist(), K)
    pcs, var = pca_from_features(x)
    cluster_ids = find_clusters(x, N_CLUSTERS, BATCH_SIZE)
    nearest = nearest_similarity(x)
    bins, _ = bin_by_cluster(cluster_ids, pcs[:, 0])
    sim = binned_similarity(x, bins)

    out = cand[["seq_id", "nt_scFv", "seq_len", PRESENT]].copy()
    out["cluster_id"] = cluster_ids
    out["PC1"] = pcs[:, 0]
    out["PC2"] = pcs[:, 1]
    out["nearest_similarity"] = nearest
    out.to_csv(args.outdir / "allseq_assignments.csv", index=False)
    count_cluster_sizes(cluster_ids).to_csv(args.outdir / "allseq_cluster_sizes.csv", index=False)
    pd.DataFrame(sim).to_csv(args.outdir / "allseq_similarity_heatmap_binned.csv", index=False)

    summary = [
        "pipeline=r4_h460pcnb_families",
        f"n_candidates={len(cand)}",
        f"kmer={K}",
        f"n_clusters={N_CLUSTERS}",
        f"heatmap_bins={sim.shape[0]}",
        f"pca_var_pc1={float(var[0]):.6f}",
        f"pca_var_pc2={float(var[1]):.6f}",
        f"nearest_similarity_median={float(np.median(nearest)):.6f}",
        f"nearest_similarity_q90={float(np.quantile(nearest, 0.9)):.6f}",
    ]
    (args.outdir / "summary_allseq.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
