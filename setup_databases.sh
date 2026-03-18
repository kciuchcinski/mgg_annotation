#!/usr/bin/env bash
# ============================================================================
# setup_databases.sh — Download and prepare databases for the MGG Annotation
#                       Pipeline.
#
# Usage:
#   bash setup_databases.sh [DB_ROOT]
#
#   DB_ROOT  Directory where databases will be stored (default: ./databases)
#
# The script downloads each database, extracts it, and renames the directory
# to the version-agnostic path expected by the pipeline defaults. If a
# database directory already exists it is skipped (delete it to re-download).
# ============================================================================
set -euo pipefail

DB_ROOT="${1:-./databases}"
mkdir -p "$DB_ROOT"
DB_ROOT="$(cd "$DB_ROOT" && pwd)"   # resolve to absolute path

echo "==> Database root: $DB_ROOT"
echo ""

# ---------------------------------------------------------------------------
# Helper: download, extract, and rename an archive into a target directory.
#   download_db <url> <target_dir> <glob_prefix>
#
#   url          — download URL
#   target_dir   — final directory name inside DB_ROOT (e.g. "phrogs")
#   glob_prefix  — prefix of the ffdata/ffindex files to locate after
#                  extraction, used to find the versioned directory
# ---------------------------------------------------------------------------
download_db() {
    local url="$1"
    local target_dir="$2"
    local glob_prefix="$3"
    local archive="${url##*/}"

    if [[ -d "$DB_ROOT/$target_dir" ]]; then
        echo "  [skip] $DB_ROOT/$target_dir already exists"
        return 0
    fi

    echo "  Downloading $archive ..."
    wget -q --show-progress -P "$DB_ROOT" "$url"

    echo "  Extracting ..."
    tar xzf "$DB_ROOT/$archive" -C "$DB_ROOT"
    rm -f "$DB_ROOT/$archive"

    # Find the extracted directory by looking for the ffdata files.
    # The archive may unpack into a versioned directory (e.g. pfamA_35.0/).
    local ffdata_file
    ffdata_file="$(find "$DB_ROOT" -maxdepth 2 -name "${glob_prefix}*_a3m.ffdata" -print -quit)"

    if [[ -z "$ffdata_file" ]]; then
        echo "  [error] Could not find ${glob_prefix}*_a3m.ffdata after extraction."
        echo "          Please check the archive contents and set up this database manually."
        return 1
    fi

    local extracted
    extracted="$(dirname "$ffdata_file")"

    if [[ "$extracted" == "$DB_ROOT" ]]; then
        # Files were extracted directly into DB_ROOT (no subdirectory).
        # Move them into the target directory.
        echo "  Moving files into $target_dir/"
        mkdir -p "$DB_ROOT/$target_dir"
        mv "$DB_ROOT"/${glob_prefix}* "$DB_ROOT/$target_dir/"
    elif [[ "$extracted" != "$DB_ROOT/$target_dir" ]]; then
        echo "  Renaming $(basename "$extracted") -> $target_dir"
        mv "$extracted" "$DB_ROOT/$target_dir"
    fi

    # Rename HHsuite file prefixes to match the target directory name if needed.
    # e.g. pfamA_35.0_a3m.ffdata -> pfam_a3m.ffdata
    local actual_prefix
    actual_prefix="$(basename "${ffdata_file%%_a3m.ffdata}")"
    if [[ "$actual_prefix" != "$target_dir" ]]; then
        echo "  Renaming file prefix: $actual_prefix -> $target_dir"
        for f in "$DB_ROOT/$target_dir/${actual_prefix}"*; do
            local suffix="${f#*"$actual_prefix"}"
            mv "$f" "$DB_ROOT/$target_dir/${target_dir}${suffix}"
        done
    fi

    echo "  Done."
}

# ============================= PHROGs v4 ====================================
echo "[1/4] PHROGs v4 (HHsuite database)"
download_db \
    "https://phrogs.lmge.uca.fr/downloads_from_website/phrogs_hhsuite_db.tar.gz" \
    "phrogs" \
    "phrogs"
echo ""

# ============================= PHROGs metadata ==============================
echo "[2/4] PHROGs annotation table"
TABLES_DIR="$DB_ROOT/tables"
mkdir -p "$TABLES_DIR"
if [[ -f "$TABLES_DIR/phrog_annot_v4.tsv" ]]; then
    echo "  [skip] $TABLES_DIR/phrog_annot_v4.tsv already exists"
else
    echo "  Downloading ..."
    wget -q --show-progress -O "$TABLES_DIR/phrog_annot_v4.tsv" \
        "https://phrogs.lmge.uca.fr/phrog_table/phrogs_table_almostfinal_plusGO_wNA_utf8.tsv"
    echo "  Done."
fi
echo ""

# ============================= Pfam-A =======================================
echo "[3/4] Pfam-A (HHsuite database)"
download_db \
    "https://wwwuser.gwdguser.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_35.0.tar.gz" \
    "pfam" \
    "pfam"
echo ""

# ============================= ECOD =========================================
echo "[4/4] ECOD (HHsuite database)"
download_db \
    "http://prodata.swmed.edu/ecod/distributions/ecod.v294.F40.hhm_db.tar.gz" \
    "ecod" \
    "ecod"
echo ""

# ============================= Glimmer ICM ==================================
# TODO: Add download URL for the Glimmer training file when available.
#       For now, ensure training-file_refseq.icm is placed manually in DB_ROOT.
if [[ -f "$DB_ROOT/training-file_refseq.icm" ]]; then
    echo "[info] Glimmer training file found: $DB_ROOT/training-file_refseq.icm"
else
    echo "[warn] Glimmer training file not found at $DB_ROOT/training-file_refseq.icm"
    echo "       Please place training-file_refseq.icm in $DB_ROOT manually."
    # Placeholder for future automated download:
    # wget -q --show-progress -O "$DB_ROOT/training-file_refseq.icm" \
    #     "https://PLACEHOLDER_URL/training-file_refseq.icm"
fi
echo ""

# ============================= Summary ======================================
echo "============================================"
echo " Database setup complete."
echo " DB_ROOT: $DB_ROOT"
echo ""
echo " Expected layout:"
echo "   $DB_ROOT/"
echo "   ├── training-file_refseq.icm"
echo "   ├── phrogs/"
echo "   │   └── phrogs{_a3m.ffdata, _a3m.ffindex, ...}"
echo "   ├── pfam/"
echo "   │   └── pfam{_a3m.ffdata, _a3m.ffindex, ...}"
echo "   ├── ecod/"
echo "   │   └── ecod{_a3m.ffdata, _a3m.ffindex, ...}"
echo "   └── tables/"
echo "       └── phrog_annot_v4.tsv"
echo ""
echo " Set DB_ROOT in config.yml to: $DB_ROOT"
echo "============================================"