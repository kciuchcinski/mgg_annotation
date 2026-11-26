#!/usr/bin/env bash
# ============================================================================
# Equivalence runs on the HPC: ORIGINAL (Snakemake) vs NEW (Nextflow).
#
# Usage (from any directory, after editing the three *.yml files):
#   bash run_reference.sh original  /path/to/mgg_annotation      [CORES]
#   bash run_reference.sh new       /path/to/mgg_annotation_kc   [CORES]
#
#   original : runs mgg_annotation @ 66ba474 (branch prodigal-gv) with original_config.yml
#   new      : runs the Nextflow pipeline twice (new_filters_off.yml, new_defaults.yml)
#              with the apptainer profile (build phage_annotation.sif first, see README)
#
# Each run writes PROVENANCE.txt into its OUTPUT_DIR: commit, tool versions and
# checksums of the database index files, so the comparison is reproducible.
# Copy the OUTPUT_DIRs back and run benchmark/compare_runs.py (see EQUIVALENCE.md).
# ============================================================================
set -euo pipefail

MODE="${1:?usage: run_reference.sh original|new <repo_dir> [cores]}"
REPO="$(cd "${2:?repo dir required}" && pwd)"
CORES="${3:-16}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

yaml_get() {  # yaml_get <file> <key> : top-level or one-level-nested "key: value"
    grep -E "^\s*$2:" "$1" | head -1 | sed -E "s/^\s*$2:\s*//; s/\s*#.*$//; s/^['\"]//; s/['\"]$//"
}

checksum_db() {  # checksum_db <prefix> : md5 of ffindex files + sizes of ffdata files
    local prefix="$1"
    for f in "${prefix}"_*.ffindex; do [[ -f "$f" ]] && md5sum "$f"; done
    for f in "${prefix}"_*.ffdata;  do [[ -f "$f" ]] && echo "size $(stat -c %s "$f")  $f"; done
}

provenance() {  # provenance <config> <outdir> <repo>
    local cfg="$1" out="$2" repo="$3"
    mkdir -p "$out"
    {
        echo "date: $(date -Iseconds)"
        echo "host: $(hostname)"
        echo "repo: $repo"
        echo "commit: $(git -C "$repo" rev-parse HEAD) ($(git -C "$repo" rev-parse --abbrev-ref HEAD))"
        echo "dirty: $(git -C "$repo" status --porcelain | wc -l) changed files"
        echo "config: $cfg"
        echo "----- config -----"
        cat "$cfg"
        echo "----- databases -----"
        for key in PHROGS PFAM ECOD ALANDB; do
            local p; p="$(yaml_get "$cfg" "$key" || true)"
            [[ -n "$p" ]] && { echo "[$key] $p"; checksum_db "$p"; }
        done
        local icm; icm="$(yaml_get "$cfg" GLIMMER_TRAIN || true)"
        [[ -n "$icm" && -f "$icm" ]] && md5sum "$icm"
        local tbl; tbl="$(yaml_get "$cfg" PHROGS_TABLE || true)"
        [[ -n "$tbl" ]] && ( cd "$repo" && md5sum "$tbl" )
    } > "$out/PROVENANCE.txt"
}

case "$MODE" in
original)
    cfg="$HERE/original_config.yml"
    out="$(yaml_get "$cfg" OUTPUT_DIR)"
    git -C "$REPO" checkout --quiet 66ba474
    provenance "$cfg" "$out" "$REPO"
    cd "$REPO"
    snakemake --snakefile ANNOTATION --configfile "$cfg" --use-conda --cores "$CORES" \
        2>&1 | tee "$out/snakemake.log"
    # record the exact tool versions conda resolved for each rule environment
    for env in .snakemake/conda/*/; do
        [[ -d "$env" ]] || continue
        echo "===== $env" >> "$out/PROVENANCE.txt"
        conda list -p "$env" >> "$out/PROVENANCE.txt" 2>/dev/null || true
    done
    ;;
new)
    for name in new_filters_off new_defaults; do
        cfg="$HERE/$name.yml"
        out="$(yaml_get "$cfg" OUTPUT_DIR)"
        provenance "$cfg" "$out" "$REPO"
        ( cd "$REPO" && nextflow run main.nf -profile apptainer -params-file "$cfg" \
              -work-dir "$out/work" -with-trace "$out/trace.tsv" -ansi-log false ) \
            2>&1 | tee "$out/nextflow.log"
        echo "----- container -----" >> "$out/PROVENANCE.txt"
        apptainer exec "$REPO/phage_annotation.sif" bash -c \
            'mmseqs version; hhblits -h | head -1; prodigal-gv -v 2>&1 | head -1; micromamba list' \
            >> "$out/PROVENANCE.txt" 2>&1 || true
    done
    ;;
*)
    echo "unknown mode: $MODE (use original|new)"; exit 1 ;;
esac
