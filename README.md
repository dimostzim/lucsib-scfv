# lucsib-scfv

Analysis of Nanopore sequencing of phage-display scFv libraries from the LUCSIB project, from raw reads to final
tables and figures. Ten barcoded samples cover four selection targets (TIM-3, CEA5, H460 PCNB, huAEC) over
panning rounds 2 to 4. scFvs are annotated against dog germline genes.

This repository holds code only. The inputs (reads, sample sheet, TIM-3 Sanger table) are in a restricted Zenodo
record and are downloaded by the pipeline; the results are added to the same record.

The annotation uses Mikolaj Kocikowski's updated fork of seq2scfv and his dog IgBLAST references, and step 04 is an
adapted version of his R notebook from `covid-scfv-publication`. Both are downloaded at pinned commits
(`config/settings.env`).

## Requirements

- Linux with Docker (the user must be able to run `docker`), conda, git and Python 3.
- About 250 GB of disk for the work folder and 64 CPU threads (fewer works, more slowly).
- Internet access to Zenodo, GitHub, conda-forge, bioconda and Docker Hub. Proxy variables in the environment are
  passed to the Docker builds.
- A Zenodo access token that can read the input record.

## Running

```bash
git clone <this repository> lucsib-scfv
cd lucsib-scfv
cp .env.example .env && chmod 600 .env    # then put the Zenodo token in .env
# edit WORK_DIR and CPUS in config/settings.env
bash run_all.sh                           # all steps; or e.g. --from 04 --to 07
```

`run_all.sh` prints `PROGRESS` lines and stops at the first error. Each step writes to its own folder in the work
folder, and step 07 collects the final tables and figures into `<WORK_DIR>/results` (see
`steps/07_collect_results/results_README.md` for what each file is).

## Steps

| Step | Folder | What it does | Main outputs |
| --- | --- | --- | --- |
| 00 | `steps/00_download` | downloads the inputs from Zenodo and checks their checksums (`config/inputs_sha256.tsv`); downloads seq2scfv and the dog references at pinned commits; builds the QC and seq2scfv Docker images; creates the R and Python conda environments (`envs/`) | `inputs/`, `tools/` |
| 01 | `steps/01_qc` | removes malformed reads, trims adapters and drops chimeric reads (Porechop), keeps reads of at least 1,000 bases and the best 90% by quality (Filtlong), NanoPlot and MultiQC reports | `01_qc/2_filtered/reads`, read counts per stage |
| 02 | `steps/02_preprocess` | converts reads to FASTA and orients them using IgBLAST | `02_preprocess/barcodeXX.fasta` |
| 03 | `steps/03_annotate` | seq2scfv: merges identical reads, IgBLAST V(D)J calls, splits each read into heavy chain, linker and light chain, scores the linker, flags problems, counts each sequence per barcode | `03_annotate/7.Counts/…_counts.tsv`, attrition |
| 04 | `steps/04_count_tables` | R notebook: filters the scFvs, collapses identical heavy and light chains, CDR3 pairs and heavy CDR3s, and draws the overview figures | count tables, per-read filtered table, figures 5 to 13 |
| 05 | `steps/05_specificity` | sequence families of H460 PCNB round 4 candidates and of all sequences, sample trajectories, CDR3 sharing | family tables, figures 1 to 4 |
| 06 | `steps/06_tim3_binders` | matches TIM-3 scFvs to 30 Sanger-sequenced bind/no-bind clones, labels them, and picks non-redundant representative sets | labels, report, representative sets |
| 07 | `steps/07_collect_results` | copies the final tables and figures into one folder | `results/` |

`docs/methods.md` describes the methods and settings of each step.

## Repository layout

```text
run_all.sh                 runs every step
config/settings.env        work folder, CPUs, Zenodo record, pinned tool versions, image names
config/inputs_sha256.tsv   checksums of the Zenodo input files
envs/                      Dockerfiles (QC, seq2scfv) and conda environments (R, Python)
steps/                     one folder per step
tools/zenodo_upload.py     maintainer tool to create the Zenodo draft and upload files
docs/methods.md            methods
```
