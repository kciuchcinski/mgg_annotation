# Equivalence with the original MGG annotation pipeline

This folder holds everything needed to show, stage by stage, that this Nextflow pipeline
reproduces the original Snakemake pipeline (`mgg_annotation`, branch `prodigal-gv`, commit
`66ba474`), and to explain every remaining difference.

## What is expected to differ, and why

These changes are **intentional** and are not counted as failures:

| Change | Effect on outputs |
|---|---|
| ALANDB search and the MGG PHROGs table removed | No `ALAN*` annotations (the MGG table never affected the original's output) |
| N best hits per database (`N_FUNCTIONS_PER_DB`) | `PFAM` → `PFAM1, PFAM2`, `ECOD` → `ECOD1, ECOD2` |
| E-value and bit score columns | `*_bits`, `*_evalue` replace the original `*_params` strings |
| Hit filter: prob ≥ 0.7, qcov ≥ 0.1, tcov ≥ 0.1 | Weak PFAM/ECOD hits (and PHROGs hits below the thresholds) are no longer reported |
| Real A3M MSAs (`mmseqs result2msa --msa-format-mode 5`) | Slightly different HH-suite scores; borderline hits can flip |
| `report.tsv` holds the filtered hits | The original `report.tsv` (best hits per cluster) is not written |
| GenBank: `/seq` holds the ORF nucleotides, header `Prophage <id>` | Qualifier values differ; the original also wrote the last CDS of each contig twice |
| One genome per input file; length cutoff `<` instead of `≤` | Multi-record FASTA files are skipped (listed in `dropped_inputs.tsv`) instead of being split |
| MMseqs2 13.45111 → 16.747c6 | Clusters can differ slightly |

## Runs

| Run | Pipeline | Where | Config | Purpose |
|---|---|---|---|---|
| 1 | original | HPC | `hpc/original_config.yml` | Reference |
| 2(i) | new, hit filter off | HPC (Apptainer) | `hpc/new_filters_off.yml` | Isolates the refactor from the hit filter |
| 2(ii) | new, defaults | HPC (Apptainer) | `hpc/new_defaults.yml` | Shows the effect of the intentional changes |
| 3 | new, defaults | MacBook (Docker, arm64) | same as 2(ii), DB paths changed | Platform check: Mac must equal HPC |

All runs use the **same input genomes and the same database files** (the original's PHROGs,
Pfam 32 and ECOD_F70_20200207 copies), and `SEARCH_TOOL: hhblits` for the original so that
the number of search iterations matches (PHROGs 1, Pfam 2, ECOD 2).

Inputs: `test/data` (3 genomes pass the filter; `phage4` is 420 bp and is dropped) and a real set
of 10–50 phages.

### On the HPC

```bash
# 1. Edit the /EDIT/ paths in benchmark/hpc/*.yml (same DB paths in all three files)

# 2. Reference: the original pipeline (checks out 66ba474)
bash benchmark/hpc/run_reference.sh original /path/to/mgg_annotation 16

# 3. The new pipeline, both configs (needs phage_annotation.sif, see README)
bash benchmark/hpc/run_reference.sh new /path/to/mgg_annotation_kc 16
```

Each output directory gets a `PROVENANCE.txt` (commit, config, database checksums, tool versions).
Check that the database checksums are identical across the three runs before comparing.

### Compare

Copy the three output directories back, then (from the repository root):

```bash
run() { docker run --rm -u $(id -u):$(id -g) -v "$PWD":"$PWD" -w "$PWD" phage_annotation:latest "$@"; }

run python benchmark/compare_runs.py --a out/original --a-layout orig --a-label original \
    --b out/new_filters_off --b-layout new --b-label "new (filter off)" --out benchmark/results/orig_vs_filters_off.md
run python benchmark/compare_runs.py --a out/original --a-layout orig --a-label original \
    --b out/new_defaults --b-layout new --b-label "new (defaults)" --out benchmark/results/orig_vs_defaults.md
run python benchmark/compare_runs.py --a out/new_defaults --a-layout new --a-label HPC \
    --b out/mac_defaults --b-layout new --b-label Mac --out benchmark/results/hpc_vs_mac.md
```

## Acceptance criteria (proposal — to be agreed with the PI)

| Stage | Original vs new (filter off) | HPC vs Mac |
|---|---|---|
| 1. Input & ORFs | identical | identical |
| 2. Clustering | ≥ 95% of proteins in identical clusters | identical |
| 3. Hits | ≥ 95% same top hit per database | identical top hits; scores equal up to rounding |
| 4. Annotation | ≥ 95% agreement per column; every mismatch has a cause | identical |
| 5. GenBank | identical CDS (after removing the original's duplicated last CDS) | identical |

The thresholds are command-line options of `compare_runs.py` (`--min-clustering`, `--min-hits`,
`--min-annotation`). If clustering agreement is low, the next step is to feed the original's
`raw_PCs.tsv` into the new pipeline so the search and annotation stages can be compared on
identical clusters.

## Results

### Done on the MacBook (2026-10-17, Apple Silicon, 11 cores, Docker Desktop 8 GB)

- **Glimmer built from source (arm64) vs bioconda's linux-64 build:** byte-identical `.predict`
  output on `test/data` (phage1–3).
- **PI setup path** (fresh copy of the repository, `setup_databases.sh` with stock macOS tools,
  `docker build`, `nextflow run main.nf -params-file test/config.yml`, no `-profile`):
  84/84 tasks succeeded in 9 min 27 s; largest task peak memory 1.9 GB (`ENRICH_MSA_PHROGS`);
  `COLLECT_HITS` / `BUILD_ANNOTATION` ~0.35 GB.
- **`compare_runs.py`** checked on a run against itself (all PASS), against a copy with planted
  ORF and annotation changes (FAIL / REVIEW as expected), and against a copy rearranged into the
  original's layout and string formats (all PASS).
- **ECOD known issue confirmed:** with ECOD v294 F40 from `setup_databases.sh`, `ECOD1_function`
  holds raw IDs (e.g. `4676609|P44225_nD2|4038.1.1.10`) and `ECOD2_function` is always empty.

### Pending: HPC runs 1, 2(i), 2(ii) and the comparison reports
