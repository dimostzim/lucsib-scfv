#!/usr/bin/env python3
"""Functions shared by the specificity scripts: sample sheet, k-mer profiles, PCA, clusters and binned similarity."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import MiniBatchKMeans

DNA_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}
HEATMAP_BINS = 140


def load_samples(path: Path) -> pd.DataFrame:
    """Sample sheet sorted by target, round and sample, with count-table column and label (e.g. '09 H460 PCNB R4')."""
    samples = pd.read_csv(path)
    samples["barcode"] = samples["Barcode"].str.extract(r"(\d+)$", expand=False).astype(int).map("barcode{:02d}".format)
    samples["round_label"] = "R" + samples["Round"].astype(str)
    sample = samples["Sample"].map("{:02d}".format)
    samples["sample_label"] = sample + " " + samples["Pool"] + " " + samples["round_label"]
    return samples.sort_values(["Pool", "Round", "Sample"]).reset_index(drop=True)


def _kmer_index(kmer: str) -> int | None:
    idx = 0
    for c in kmer:
        v = DNA_INDEX.get(c)
        if v is None:
            return None
        idx = (idx << 2) | v
    return idx


def kmer_matrix(seqs: list[str], k: int) -> np.ndarray:
    """L2-normalised k-mer counts, one row per sequence; k-mers with letters other than ACGT are not counted."""
    x = np.zeros((len(seqs), 4**k), dtype=np.float32)
    for i, seq in enumerate(seqs):
        s = seq.upper()
        vec = x[i]
        for j in range(len(s) - k + 1):
            idx = _kmer_index(s[j : j + k])
            if idx is not None:
                vec[idx] += 1.0
        norm = float(np.linalg.norm(vec))
        if norm > 0:
            vec /= norm
    return x


def pca_from_features(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Scores on the first two principal components and the fraction of variance each explains."""
    xc = x - x.mean(axis=0, keepdims=True)
    cov = (xc.T @ xc) / (len(x) - 1)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    pcs = xc @ eigvecs[:, :2]
    var = eigvals[:2] / float(eigvals.sum())
    return pcs.astype(np.float32), var.astype(np.float32)


def find_clusters(x: np.ndarray, n_clusters: int, batch_size: int) -> np.ndarray:
    km = MiniBatchKMeans(n_clusters=n_clusters, random_state=42, n_init=10, batch_size=batch_size, max_iter=300)
    return km.fit_predict(x).astype(np.int32)


def count_cluster_sizes(cluster_ids: np.ndarray) -> pd.DataFrame:
    return (
        pd.DataFrame({"cluster_id": cluster_ids})
        .groupby("cluster_id", as_index=False)
        .size()
        .rename(columns={"size": "cluster_size"})
        .sort_values("cluster_size", ascending=False)
        .reset_index(drop=True)
    )


def nearest_similarity(x: np.ndarray, chunk_size: int = 512) -> np.ndarray:
    """Cosine similarity of each row of x (L2-normalised) to its most similar other row."""
    nearest = np.empty(len(x), dtype=np.float32)
    for start in range(0, len(x), chunk_size):
        sims = x[start : start + chunk_size] @ x.T
        rows = np.arange(len(sims))
        sims[rows, start + rows] = -np.inf
        nearest[start : start + len(sims)] = sims.max(axis=1)
    return nearest


def bin_by_cluster(cluster_ids: np.ndarray, pc1: np.ndarray) -> tuple[list[np.ndarray], np.ndarray]:
    """Sequences ordered by cluster, then PC1, split into equal bins; returns the bins and each bin's main cluster."""
    bins = np.array_split(np.lexsort((pc1, cluster_ids)), min(HEATMAP_BINS, len(cluster_ids)))
    majority = np.array([np.bincount(cluster_ids[b]).argmax() for b in bins], dtype=np.int32)
    return bins, majority


def binned_similarity(x: np.ndarray, bins: list[np.ndarray]) -> np.ndarray:
    """Cosine similarity between the mean k-mer profiles of the bins."""
    means = []
    for b in bins:
        m = x[b].mean(axis=0)
        norm = np.linalg.norm(m)
        means.append(m / norm if norm > 0 else m)
    v = np.vstack(means).astype(np.float32)
    sim = np.clip(v @ v.T, 0.0, 1.0)
    np.fill_diagonal(sim, 1.0)
    return sim
