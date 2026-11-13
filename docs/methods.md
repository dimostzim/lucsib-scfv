# Methods

This document describes what each step of the pipeline does and with which settings. Versions are pinned in
`envs/` (conda packages, base images by digest) and `config/settings.env` (seq2scfv and dog reference commits).
Some dependencies still float: the Miniconda installer, apt packages and CRAN packages inside the seq2scfv image
(the CRAN packages are not used by the pipeline), and the dependencies of the pinned conda packages. Each run
records what it used: step 07 copies the package lists of both conda environments and both Docker images into
`results/provenance/`.

## Samples and inputs

Ten barcoded libraries from one Nanopore run (flow cell FBD53635, run `20251008_Salibas_run1`, basecalling model
`dna_r10.4.1_e8.2_400bps_sup`) were sequenced. `samples.csv` maps barcodes to target and panning round:

| Barcodes | Target | Rounds |
| --- | --- | --- |
| 2, 3 | TIM-3 | 3, 4 |
| 5, 6 | CEA5 | 3, 4 |
| 7, 8, 9 | H460 PCNB | 2, 3, 4 |
| 10, 11, 12 | huAEC | 2, 3, 4 |

`TIM-3_sanger_binding.tsv` holds 30 Sanger-sequenced TIM-3 clones, each labelled as binding (`yes`) or not
binding (`no`). All inputs are in a restricted Zenodo record (https://doi.org/10.5281/zenodo.24169294); step 00
downloads them, unless they are already in the work folder, and checks their SHA-256 checksums against
`config/inputs_sha256.tsv`.

## Step 01: read QC

Each barcode is processed in the QC Docker image (seqkit, Porechop, Filtlong, NanoPlot, MultiQC), several barcodes
in parallel. `seqkit sana` removes malformed records (for example a sequence and quality string of different
length). Porechop removes adapters and, with `--discard_middle`, discards reads with an adapter in the middle,
which are chimeric. Filtlong keeps reads of at least 1,000 bases (`--min_length 1000`) and then the best 90% of
bases by quality (`--keep_percent 90`). NanoPlot summarises the reads after each of the three stages, and MultiQC
combines the 30 summaries into one report.

## Step 02: preprocessing

Reads are converted from FASTQ to FASTA. To orient them, IgBLAST 1.20.0 is run on every read against the dog
V, D and J germline databases (`-organism dog -ig_seqtype Ig -domain_system imgt`, one alignment per gene,
AIRR output, extended alignments at both ends, `dog_gl.aux`). Reads that IgBLAST reports in reverse orientation are
reverse-complemented; all other reads, including reads without a V gene call, are kept as they are.

## Step 03: scFv annotation

seq2scfv (Mikolaj Kocikowski's fork, pinned commit) runs in its Docker image:

1. Identical reads are merged across barcodes using a sequence checksum (SEGUID), keeping the read-to-barcode
   correspondence.
2. Each unique sequence is annotated with IgBLAST against the dog germlines, iteratively to find both domains.
3. Each sequence is split into heavy chain, linker and light chain; in-frame sequences with one heavy and one light
   domain are kept.
4. The linker is inferred: 10% of the inferred linkers are sampled (`seqkit sample -p 0.1`), aligned with Clustal
   Omega, and their consensus (`EGKSSGASGESKVDDA`) becomes the reference; every linker is scored against it.
5. Sequences are flagged (linker identity below 90%, more than 2 mismatches).
6. Each sequence is counted per barcode.

The attrition table lists, in order: reads over all barcodes, unique read sequences, IgBLAST domain hits (one row
per domain, so this can exceed the number of sequences), sequences without one paired heavy and light domain,
out-of-frame scFvs, scFvs with non-standard amino acids, and validated scFvs.

## Step 04: count tables and overview figures

`count_tables.Rmd` is adapted from the R notebook in Mikolaj Kocikowski's `covid-scfv-publication`. It keeps
scFvs with the heavy chain first, at least 75% linker identity, at most 2 linker mismatches, at most 10 extra linker
amino acids, and non-empty heavy and light CDR3s. Sequences are then collapsed three ways and their read counts
summed per barcode: identical heavy and light chain protein sequences (a "distinct scFv"), identical heavy and light
CDR3 pairs, and identical heavy CDR3.

Figures:

- reads and distinct scFvs per barcode; distribution of read counts per distinct scFv;
- effective Shannon diversity per barcode, the exponential of the Shannon entropy of distinct scFv frequencies;
- sequence map: each distinct scFv (heavy chain, consensus linker, light chain protein) is described by its
  3-letter amino-acid motif frequencies and placed in two dimensions with UMAP (cosine distance, 15 neighbours,
  random seed 42), shown for all barcodes and per barcode;
- number of scFvs shared between each pair of barcodes, for each of the three collapses;
- heavy and light chain V gene pairing per barcode, weighted by reads and counting each pair once;
- the 10 most abundant scFvs of each barcode followed across all barcodes, as a percentage of each barcode's reads;
- barcodes clustered by the log number of shared scFvs (Euclidean distance, Ward linkage).

## Step 05: specificity analyses

All analyses use the per-read filtered table from step 04 and unique scFv DNA sequences.

- **H460 PCNB round 4 candidates:** sequences with reads in H460 PCNB round 4 (barcode 9) and none in huAEC round 4
  (barcode 12). Each is described by its 5-letter DNA motif frequencies (L2-normalised), projected on two principal
  components, and grouped into 20 families with mini-batch k-means (random seed 42). For each sequence the cosine
  similarity to its nearest neighbour is reported, and family similarity is shown as a binned heatmap (140 bins).
- **All rounds and targets:** the same with 4-letter motifs and 24 families on all sequences from rounds 2 to 4.
  Each sample's read-weighted average position (centroid) is computed, and the overlap of sequence sets between
  samples is measured by the Jaccard index. Joining each target's centroids by round gives the trajectories.
- **CDR3 sharing:** an UpSet plot of heavy and light CDR3 pairs over barcodes 7 to 12, showing intersections with
  at least 10 pairs.

## Step 06: TIM-3 binder labels

1. TIM-3 reads (barcodes 2 and 3) are collapsed by exact scFv DNA sequence.
2. Each sequence is compared with the 30 Sanger clones: the 5 `yes` and the 5 `no` clones sharing the most 15-letter
   k-mers with it (at least 3) are aligned with edlib (at most 100 edits), and the best `yes` and the best `no`
   match are kept. When two clones align equally well, the one sharing more k-mers, then the lower clone number, is
   reported.
3. First-pass label: a `yes` or `no` match counts if it has at least 99% identity and covers at least 90% of the
   scFv and 95% of the shorter of the two sequences. With no counting match the label is `unknown`; with one, it
   is `bind` or `no_bind`; with both, it is `ambiguous` (two counting matches can differ by at most 1 percentage
   point of identity, so the code's "higher identity wins" branch does not occur).
4. For ambiguous sequences with at least 10 reads, the competing `yes` and `no` Sanger clones are aligned to the full
   Nanopore reads. A sequence is relabelled `bind_context` when most of its reads match the `yes` clone better by
   more than 0.5 percentage points of identity and the `yes` clone aligns to the full reads with at least 91%
   identity; all other labels stay as they are.
5. Representative sets: strict binders (`bind` only) are ranked by TIM-3 reads and selected greedily if their global
   identity to every higher-ranked selected sequence is below a threshold (95, 96 and 97% for protein sequences,
   99% for DNA), up to 100 representatives; every candidate is then assigned to its nearest representative.

**What `bind_context` measures.** In this dataset all ambiguous sequences tested in step 4 pair Sanger clone 30
(`yes`) with clone 17 or 9 (`no`). These three Sanger clones have identical scFv sequences; clones 30 and 17 differ
only downstream of the scFv (from clone-30 position 849 on, in the tag and gene III region, including single-base
gaps and `N` calls in the Sanger reads). The `bind_context` call therefore rests on sequence outside the scFv and is
not evidence that the scFv binds TIM-3. `expanded_bind_positive` (`bind` plus `bind_context`) inherits this.
The representative sets use `bind` only.

## Changes from the February to June 2026 analysis

The pipeline reproduces the earlier analysis, with these differences:

- barcode 5 now uses the complete read file; the earlier analysis used a copy cut off after 363,179 of 383,375
  reads;
- step 01 runs barcodes in parallel; the MultiQC report now contains all three QC stages (sample names collided
  before, so it showed one); MultiQC is version 1.25.1 (the earlier report used 1.33, which cannot be installed
  with the other pinned QC tools); MultiQC builds only the report;
- step 03 no longer annotates the unoriented temporary FASTA copies left by preprocessing; these doubled the early
  attrition counts (but not the per-barcode counts); the first attrition row now counts all reads (it was one too
  low and labelled "Merged libraries");
- step 04 no longer draws figures that depended on COVID-project sample names; the top-10 plots use each barcode's
  own read total; the shared-scFv numbers are no longer shown plus one; the clustered heatmap is written correctly
  and clusters by Euclidean distance; the UMAP has a fixed random seed;
- step 06 includes the two steps (final labels with full-read context, and the labelled table with sequences)
  whose code was not kept in the earlier analysis. They were rebuilt from the saved outputs. The 91% condition on the
  `yes` clone is a rule that reproduces the saved labels exactly (422 `bind_context` clones, 52,406 reads); the
  original threshold was not recorded, and other thresholds between 90.75% and 91.02% give the same result on the
  earlier data;
- step 06 resolves ties between equally good Sanger matches by shared k-mers and clone number; earlier the choice
  depended on Python's random set ordering and could change between runs. On the earlier data the representative
  sets and the label counts are reproduced exactly; the `no` clone reported for the context-tested sequences is now
  clone 9 instead of the identical clone 17, and the cutoff-sensitivity table differs in 15 of 225 rows (all at a
  margin of 0) by one clone;
- the Python environment pins scikit-learn 1.8.0, which reproduces the earlier sequence families exactly
  (version 1.9 groups about 10% of sequences differently).
