# lucsib-scfv

This pipeline analyses the Nanopore sequencing of the LUCSIB phage-display scFv libraries, from the raw reads to the
final tables and figures.

An scFv (single-chain variable fragment) is the heavy chain and light chain variable domains of an antibody joined
by a short linker. The scFv libraries were panned against four targets: TIM-3, CEA5, H460 PCNB and huAEC. In each
panning round, the phage that bind the target are kept and amplified, so binders should become more common in later
rounds. The scFv plasmids of rounds 2 to 4 were sequenced in one Nanopore run, with one barcode for each target and
round.

The pipeline finds the scFv in each read, splits it into heavy chain, linker and light chain, annotates it against
dog germline genes, counts how often each scFv was read in each sample, and compares the samples. For TIM-3 it also
uses 30 Sanger-sequenced clones with known binding results to label the scFvs found by Nanopore as likely binders or
non-binders.

This repository holds only code. The input files and the results are in a restricted Zenodo record,
https://doi.org/10.5281/zenodo.24169294. The scFv annotation uses Mikolaj Kocikowski's updated fork of seq2scfv and
his dog germline references, and step 04 is adapted from his R notebook in `covid-scfv-publication`; both are
downloaded at fixed versions.

## Input files

The pipeline needs 12 files, all in the Zenodo record:

| File | What it is |
| --- | --- |
| `barcode02.fastq.gz` to `barcode12.fastq.gz` (10 files) | The Nanopore reads of each sample: the DNA sequence and per-base quality of every read, gzip-compressed. Each sample has 282,314 to 568,247 reads, and each file is 0.7 to 1.4 GB (9.9 GB in total). |
| `samples.csv` | The target and panning round of each barcode (table below). |
| `TIM-3_sanger_binding.tsv` | The Sanger sequences of 30 TIM-3 clones, each with its binding result. Of the 30 clones, 19 bind TIM-3 and 11 do not. |

| Barcodes | Target | Panning rounds |
| --- | --- | --- |
| 2, 3 | TIM-3 | 3, 4 |
| 5, 6 | CEA5 | 3, 4 |
| 7, 8, 9 | H460 PCNB | 2, 3, 4 |
| 10, 11, 12 | huAEC | 2, 3, 4 |

`config/inputs_sha256.tsv` lists the 12 files with their checksums, which step 00 checks.

## Results

Step 07 copies the final tables and figures into `<WORK_DIR>/results/`; the same files are in the results archive on
Zenodo. `results/README.md` describes every file. The main results are:

- 13 numbered figures, for example the sequence families of H460 PCNB round 4, how often each scFv was read, the
  diversity of each sample, and how many scFvs the samples share.
- Count tables with the number of reads of each distinct scFv in each sample. A distinct scFv is a unique pair of
  heavy and light chain protein sequences. Each chain has three antigen-binding loops, the complementarity-determining
  regions (CDRs), and the CDR3 is the most variable of them. There are also tables for pairs of heavy and light chain
  CDR3s and for heavy chain CDR3s alone.
- Read counts after each QC stage and through the scFv annotation.
- The TIM-3 scFvs labelled by matching them to the Sanger clones, and non-redundant sets of likely binders.
- The pipeline version, settings and software versions of the run (`provenance/`).

## Analysis steps

`run_all.sh` runs the steps in order. Each step writes to its own folder in the work folder.

- **00 Inputs and tools.** Gets the input files, downloads seq2scfv and the dog germline references at fixed
  versions, and installs the software (two Docker images and two conda environments).
- **01 Read QC.** Removes malformed reads, trims the sequencing adapters, discards chimeric reads (reads with an
  adapter in the middle, which come from two molecules joined together), and keeps reads of at least 1,000 bases and
  the best 90% by quality (seqkit, Porechop and Filtlong, with NanoPlot and MultiQC reports).
- **02 Read orientation.** Converts the reads to FASTA and turns each read so that the scFv runs in the same
  direction, using IgBLAST against the dog germline genes.
- **03 scFv annotation.** seq2scfv finds the heavy chain, linker and light chain in each read, assigns germline genes
  and CDRs, keeps reads with one complete in-frame scFv, and counts each scFv in each sample.
- **04 Count tables.** An R notebook filters the scFvs, counts each distinct scFv, CDR3 pair and heavy chain CDR3 in
  each sample, and draws figures 5 to 13.
- **05 Sequence families.** Groups scFv DNA sequences into families by sequence similarity, for H460 PCNB round 4
  and for all rounds and targets, follows each target across rounds, and compares CDR3 sharing between the H460 PCNB
  and huAEC rounds (figures 1 to 4).
- **06 TIM-3 binders.** Matches each TIM-3 scFv to the 30 Sanger clones, labels it from the clones it matches best,
  and picks non-redundant sets of likely binders.
- **07 Results.** Copies the final tables and figures into `results/`.

`docs/methods.md` describes the methods and settings of each step, and how this pipeline differs from the earlier
analysis of February to June 2026.

## Running the pipeline

### Requirements

- Linux with Docker (your user must be able to run `docker`), conda, git and Python 3.8 or newer.
- About 250 GB of disk, 64 CPU threads (fewer also works, more slowly) and at least 64 GB of memory.
- Internet access to Zenodo, GitHub, conda-forge, bioconda and Docker Hub. Proxy variables in the environment are
  passed to the Docker builds.

A full run took 6 h 09 min with 64 threads. About 3 h 40 min of it is the scFv annotation (step 03).

### Getting the input files

Step 00 gets the input files in one of two ways:

- **Zenodo access.** The files in the record are restricted. Ask the record owner (Dimosthenis Tzimotoudis) to share
  the record with your Zenodo account. Then create a personal access token in your Zenodo account settings
  (Applications, Personal access tokens) and put it in `.env`. Step 00 downloads the files.
- **Your own copies.** Put the 12 files, with the same names, in `<WORK_DIR>/inputs/downloads/`. Step 00 checks their
  checksums and needs no token. If a file is missing or its checksum differs, step 00 downloads it, which needs the
  token, and stops with a message if there is none.

### Commands

```bash
git clone git@github.com:dimostzim/lucsib-scfv.git
cd lucsib-scfv
cp .env.example .env && chmod 600 .env    # the Zenodo token goes here (not needed with your own copies)
# set WORK_DIR and CPUS in config/settings.env
bash run_all.sh                           # all steps
bash run_all.sh --from 04 --to 07         # only steps 04 to 07
```

`run_all.sh` prints `PROGRESS` lines and stops at the first error. Rerunning a step first deletes that step's
folder, so its outputs are never mixed with older or partial files. Step 00 keeps input files that are already
downloaded and checked, and reuses conda environments whose definition in `envs/` has not changed.

## Repository layout

```text
run_all.sh                 runs every step
config/settings.env        work folder, CPUs, Zenodo record, pinned tool versions, image names
config/inputs_sha256.tsv   names and checksums of the 12 input files
envs/                      Dockerfiles (QC, seq2scfv) and conda environments (R, Python)
steps/                     one folder per step
tools/zenodo_upload.py     maintainer tool to create the Zenodo draft and upload files
docs/methods.md            methods
```
