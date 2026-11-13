#!/usr/bin/env python3
"""Copy the final tables and figures from the work folder into one results folder.

Large tables are gzip-compressed. Figures are numbered as in docs/results.md, which becomes results/README.md.
Fails if any expected output is missing or empty.
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
from pathlib import Path

SAMPLES = ["02", "03", "05", "06", "07", "08", "09", "10", "11", "12"]

# (source relative to the work folder, destination relative to the results folder)
FILES = [
    ("01_qc/pipeline_stats.tsv", "01_qc/reads_per_qc_stage.tsv"),
    ("01_qc/3_multiqc/multiqc_report.html", "01_qc/multiqc_report.html"),
    ("03_annotate/8.Attrition/attrition_table.tsv", "03_annotate/attrition.tsv"),
    ("04_count_tables/r_output/tables/attrition_merged.csv", "04_count_tables/attrition_with_r_filters.csv"),
    ("04_count_tables/r_output/tables/VH_VL.csv", "04_count_tables/counts_heavy_light_pairs.csv.gz"),
    ("04_count_tables/r_output/tables/CDR3s.csv", "04_count_tables/counts_cdr3_pairs.csv.gz"),
    ("04_count_tables/r_output/tables/HCDR3.csv", "04_count_tables/counts_heavy_cdr3.csv.gz"),
    ("04_count_tables/count_tables.html", "04_count_tables/count_tables_notebook.html"),
    ("05_specificity/r4_h460pcnb_families/candidates_summary.txt", "05_specificity/r4_h460pcnb_families/candidates_summary.txt"),
    ("05_specificity/r4_h460pcnb_families/allseq_assignments.csv", "05_specificity/r4_h460pcnb_families/assignments.csv.gz"),
    ("05_specificity/r4_h460pcnb_families/allseq_cluster_sizes.csv", "05_specificity/r4_h460pcnb_families/cluster_sizes.csv"),
    ("05_specificity/r4_h460pcnb_families/allseq_similarity_heatmap_binned.csv", "05_specificity/r4_h460pcnb_families/similarity_heatmap_binned.csv"),
    ("05_specificity/r4_h460pcnb_families/summary_allseq.txt", "05_specificity/r4_h460pcnb_families/summary.txt"),
    ("05_specificity/all_rounds/allrounds_sequence_assignments.csv", "05_specificity/all_rounds/assignments.csv.gz"),
    ("05_specificity/all_rounds/allrounds_cluster_sizes.csv", "05_specificity/all_rounds/cluster_sizes.csv"),
    ("05_specificity/all_rounds/allrounds_cluster_composition.csv", "05_specificity/all_rounds/cluster_composition.csv"),
    ("05_specificity/all_rounds/sample_jaccard_matrix.csv", "05_specificity/all_rounds/sample_overlap_jaccard.csv"),
    ("05_specificity/all_rounds/sample_centroids_pca.csv", "05_specificity/all_rounds/sample_centroids_pca.csv"),
    ("05_specificity/all_rounds/summary_allrounds.txt", "05_specificity/all_rounds/summary.txt"),
    ("06_tim3_binders/matching/tim3_sanger_vs_sanger.tsv", "06_tim3_binders/labels/sanger_vs_sanger.tsv"),
    ("06_tim3_binders/matching/tim3_cutoff_sensitivity.tsv", "06_tim3_binders/labels/cutoff_sensitivity.tsv"),
    ("06_tim3_binders/matching/tim3_ambiguous_full_context_minreads10.tsv", "06_tim3_binders/labels/ambiguous_full_context_minreads10.tsv"),
    ("06_tim3_binders/matching/tim3_labels_with_full_context_minreads10.tsv", "06_tim3_binders/labels/labels_with_full_context.tsv.gz"),
    ("06_tim3_binders/report/TIM3_LABEL_REPORT.md", "06_tim3_binders/label_report.md"),
    ("06_tim3_binders/report/label_read_summary.tsv", "06_tim3_binders/labels/label_read_summary.tsv"),
    ("06_tim3_binders/report/top_clones_by_label.tsv", "06_tim3_binders/labels/top_clones_by_label.tsv"),
    ("06_tim3_binders/report/sanger_conflict_audit.tsv", "06_tim3_binders/labels/sanger_conflict_audit.tsv"),
    ("06_tim3_binders/report/label_matrix_for_plotting.tsv", "06_tim3_binders/labels/label_matrix_for_plotting.tsv"),
    ("05_specificity/figures/figure1_r4_h460pcnb_specific.png", "figures/01_H460PCNB_R4_sequence_families.png"),
    ("05_specificity/figures/figure2_allrounds_alltargets.png", "figures/02_all_rounds_all_targets.png"),
    ("05_specificity/figures/figure3_target_trajectories.png", "figures/03_target_trajectories.png"),
    ("05_specificity/figures/upset_CDR3s_barcodes7_12.png", "figures/04_CDR3_sharing_barcodes07-12.png"),
    ("04_count_tables/r_output/barplots/scfv_counts.png", "figures/05_reads_and_distinct_scFvs_per_barcode.png"),
    ("04_count_tables/r_output/barplots/binned_counts_histogram.png", "figures/06_scFv_read_count_distribution.png"),
    ("04_count_tables/r_output/diversity/Shannon_diversity_per_sample.png", "figures/07_diversity_per_barcode.png"),
    ("04_count_tables/r_output/spaces/scFv_sequence_space_UMAP-KDE.png", "figures/08_sequence_map_all_barcodes.png"),
    ("04_count_tables/r_output/spaces/UMAP_raster_faceted.png", "figures/09_sequence_map_per_barcode.png"),
    ("04_count_tables/r_output/cooccurence/Cooccurence_faceted_numbers.png", "figures/10_sharing_between_barcodes.png"),
    ("04_count_tables/r_output/cooccurence/Cooccurence_heatmap_clustered_VH_VL.pdf", "figures/13_barcodes_clustered_by_sharing.pdf"),
]
FILES += [
    (f"04_count_tables/r_output/chords/chord_barcode{s}.pdf", f"figures/11_gene_pairing/gene_pairing_barcode{s}.pdf")
    for s in SAMPLES
]
FILES += [
    (f"04_count_tables/r_output/alluvial/alluvial_top10_barcode{s}.svg", f"figures/12_top10_scFvs/top10_scFvs_barcode{s}.svg")
    for s in SAMPLES
]
# whole folders copied as they are
FOLDERS = [
    ("06_tim3_binders/representative_sets", "06_tim3_binders/representative_sets"),
    ("06_tim3_binders/report/fasta", "06_tim3_binders/fasta"),
]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--work-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--repo-dir", type=Path, required=True, help="pipeline repository")
    args = p.parse_args()

    missing = [src for src, _ in FILES if not (args.work_dir / src).is_file() or (args.work_dir / src).stat().st_size == 0]
    missing += [src for src, _ in FOLDERS if not (args.work_dir / src).is_dir()]
    if missing:
        raise SystemExit("Missing outputs:\n" + "\n".join(missing))

    for src, dst in FILES:
        target = args.output_dir / dst
        target.parent.mkdir(parents=True, exist_ok=True)
        if dst.endswith(".gz") and not src.endswith(".gz"):
            with (args.work_dir / src).open("rb") as fin, gzip.open(target, "wb", compresslevel=6) as fout:
                shutil.copyfileobj(fin, fout, length=1 << 24)
        else:
            shutil.copyfile(args.work_dir / src, target)
    for src, dst in FOLDERS:
        shutil.copytree(args.work_dir / src, args.output_dir / dst, dirs_exist_ok=True)
    shutil.copyfile(args.repo_dir / "docs" / "results.md", args.output_dir / "README.md")

    # provenance: pipeline commit, settings, input checksums and software versions
    prov = args.output_dir / "provenance"
    prov.mkdir(parents=True, exist_ok=True)
    commit = subprocess.run(["git", "-C", str(args.repo_dir), "rev-parse", "HEAD"], check=True,
                            capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(args.repo_dir), "status", "--porcelain"], check=True,
                           capture_output=True, text=True).stdout.strip()
    (prov / "pipeline_commit.txt").write_text(f"{commit}{' (with uncommitted changes)' if dirty else ''}\n")
    for name in ("settings.env", "inputs_sha256.tsv"):
        shutil.copyfile(args.repo_dir / "config" / name, prov / name)
    for listing in sorted((args.work_dir / "tools").glob("*/*.packages.txt")):
        shutil.copyfile(listing, prov / listing.name)
    print(f"Collected {len(FILES)} files and {len(FOLDERS)} folders into {args.output_dir}")


if __name__ == "__main__":
    main()
