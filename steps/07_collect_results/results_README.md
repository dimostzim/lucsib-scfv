# LUCSIB scFv sequencing results

These are the results of the lucsib-scfv pipeline on the Nanopore reads of the phage-display panning rounds.
Samples are labelled by barcode: barcodes 2 and 3 are TIM-3 rounds 3 and 4, barcodes 5 and 6 are CEA5 rounds 3 and 4,
barcodes 7 to 9 are H460 PCNB rounds 2 to 4, and barcodes 10 to 12 are huAEC rounds 2 to 4.

A "distinct scFv" is a unique pair of heavy and light chain protein sequences. Large tables are gzip-compressed (`.csv.gz`, `.tsv.gz`).

## Figures (`figures/`)

1. `01_H460PCNB_R4_sequence_families.png`: scFv DNA sequences found in H460 PCNB round 4 and not in huAEC round 4, grouped into 20 families by sequence similarity.
2. `02_all_rounds_all_targets.png`: the same analysis on all scFv DNA sequences from rounds 2 to 4 of the four targets.
3. `03_target_trajectories.png`: the average position of each target at each round on the map from figure 2.
4. `04_CDR3_sharing_barcodes07-12.png`: how many heavy and light chain CDR3 pairs are found in only one H460 PCNB or huAEC round and how many are shared between rounds; teal bars are pairs found in both targets.
5. `05_reads_and_distinct_scFvs_per_barcode.png`: reads and distinct scFvs per barcode.
6. `06_scFv_read_count_distribution.png`: how many distinct scFvs were read once, 2 to 10 times, and so on, across all barcodes.
7. `07_diversity_per_barcode.png`: effective number of scFvs per barcode (how many equally common scFvs would give the same Shannon diversity).
8. `08_sequence_map_all_barcodes.png`: each dot is one scFv, placed so that similar protein sequences sit close together, coloured by barcode.
9. `09_sequence_map_per_barcode.png`: the same map split by barcode.
10. `10_sharing_between_barcodes.png`: number of scFvs found in both barcodes of each pair, counted as identical heavy and light chains, identical heavy and light CDR3s, or identical heavy CDR3 alone.
11. `11_gene_pairing/`: one file per barcode; which heavy chain V genes pair with which light chain V genes, weighted by reads (left) and counting each gene pair once (right).
12. `12_top10_scFvs/`: one file per barcode; the 10 most abundant scFvs of that barcode followed across all barcodes, as a percentage of each barcode's reads.
13. `13_barcodes_clustered_by_sharing.pdf`: barcodes clustered by how many scFvs with identical heavy and light chains they share (log scale).

## Tables

- `01_qc/reads_per_qc_stage.tsv`: reads per barcode before QC, after cleaning, adapter trimming and length/quality filtering. `multiqc_report.html` has the read length and quality plots.
- `03_annotate/attrition.tsv`: sequences kept at each step of the scFv annotation.
- `04_count_tables/`: read counts per barcode for each distinct scFv (`counts_heavy_light_pairs`), each heavy and light CDR3 pair (`counts_cdr3_pairs`) and each heavy CDR3 (`counts_heavy_cdr3`); `attrition_with_r_filters.csv` continues the attrition table through the filters of the count-table notebook; `count_tables_notebook.html` is the rendered notebook.
- `05_specificity/r4_h460pcnb_families/`: the H460 PCNB round 4 candidates (present in H460 PCNB round 4, absent from huAEC round 4) with their family, map position and nearest-neighbour similarity; family sizes; family similarity matrix.
- `05_specificity/all_rounds/`: the same for all sequences of rounds 2 to 4, plus each sample's average map position and the overlap (Jaccard index) between samples.
- `06_tim3_binders/`: TIM-3 scFvs labelled as binder or non-binder by matching to the 30 Sanger-sequenced clones. `label_report.md` explains the labels; `labels/` holds the label tables; `fasta/` holds the labelled sequences; `representative_sets/` holds non-redundant strict-binder sets at 95, 96 and 97% protein identity and 99% DNA identity, each as a table of all candidates with their cluster and a FASTA of the representatives.
