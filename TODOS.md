# README improvements — pre/post publication

## Blocking

- [x] **Glimmer training file — add a download source.**
  File is bundled in `resources/training-file_refseq.icm`. `setup_databases.sh` copies it
  to `DB_ROOT` automatically (step 5/5). README manual-setup table and the stale "future
  release" note updated accordingly.

## Should fix

- [x] **Inconsistency between copy instruction and run command.**
  Dropped the copy suggestion — `config.yml` is the live config, users edit it directly.
  Run commands already referenced `config.yml` and are unchanged.

- [x] **Add disk space estimate to the Database setup section.**
  PHROGs ~155 MB, Pfam-A ~2.1 GB, ECOD ~4 GB → ~6.3 GB total. Added as a callout at
  the top of the Database setup section.

- [x] **Pin HH-suite version in `environment.yml` and the Software versions table.**
  Pinned to `3.3.0` in `environment.yml`, `Dockerfile`, `phage_annotation.def`, and the
  README Software versions table.

## Nice to have

- [x] **Add a Quick Start section** near the top — three commands (clone, setup DBs, run)
  for users who already know Nextflow and just want the gist. Added after the title blurb,
  before Overview.

- [x] **Add typical runtime / memory guidance.** Added "Resource requirements" section with
  per-instance RAM for HHsuite steps, total CPU-hour scaling table (10/100/1k genomes),
  peak RAM table for BUILD_ANNOTATION, and suggested job allocations.

- [] **Clean everything up and create V2.** Make sure everything is publication-ready, finish all TODO's, start versioning again and replace the `main` branch with contents of `kc`