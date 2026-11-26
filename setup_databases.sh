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
#
# Works with the stock tools on macOS and Linux (curl, or wget as fallback).
# Interrupted downloads resume when the script is re-run.
# Needs ~30 GB free: ~6.5 GB of downloads, ~25 GB once unpacked.
# Set SKIP_SPACE_CHECK=1 to skip the free-space check.
# ============================================================================
set -euo pipefail

DB_ROOT="${1:-./databases}"
mkdir -p "$DB_ROOT"
DB_ROOT="$(cd "$DB_ROOT" && pwd)"   # resolve to absolute path

DOWNLOAD_DIR="$DB_ROOT/.downloads"
mkdir -p "$DOWNLOAD_DIR"

echo "==> Database root: $DB_ROOT"
echo ""

# ---------------------------------------------------------------------------
# Free-space check (skipped when all databases are already present)
# ---------------------------------------------------------------------------
REQUIRED_GB=30
if [[ "${SKIP_SPACE_CHECK:-0}" != "1" ]] && \
   ! { [[ -d "$DB_ROOT/phrogs" ]] && [[ -d "$DB_ROOT/pfam" ]] && [[ -d "$DB_ROOT/ecod" ]]; }; then
    avail_kb="$(df -Pk "$DB_ROOT" | awk 'NR==2 {print $4}')"
    avail_gb=$(( avail_kb / 1024 / 1024 ))
    if (( avail_gb < REQUIRED_GB )); then
        echo "[error] Only ${avail_gb} GB free in $DB_ROOT; about ${REQUIRED_GB} GB is needed."
        echo "        Free up space, choose another DB_ROOT, or set SKIP_SPACE_CHECK=1."
        exit 1
    fi
    echo "==> Free space: ${avail_gb} GB (need ~${REQUIRED_GB} GB)"
    echo ""
fi

# ---------------------------------------------------------------------------
# Helper: resumable download with curl (default on macOS) or wget.
#   fetch <url> <dest> [insecure]
#
#   insecure — pass "insecure" to skip TLS certificate checks (needed for the
#              PHROGs server, whose certificate has expired)
# ---------------------------------------------------------------------------
fetch() {
    local url="$1"
    local dest="$2"
    local insecure="${3:-}"

    if [[ -f "$dest" ]]; then
        echo "  [skip] already downloaded: $(basename "$dest")"
        return 0
    fi

    if command -v curl >/dev/null 2>&1; then
        local opts=(-fL --retry 3 --retry-delay 5 -C - -o "$dest.part")
        [[ "$insecure" == "insecure" ]] && opts+=(-k)
        curl "${opts[@]}" "$url"
    elif command -v wget >/dev/null 2>&1; then
        local opts=(-c --progress=bar:force:noscroll -O "$dest.part")
        [[ "$insecure" == "insecure" ]] && opts+=(--no-check-certificate)
        wget "${opts[@]}" "$url"
    else
        echo "  [error] Neither curl nor wget is installed."
        return 1
    fi
    mv "$dest.part" "$dest"
}

# ---------------------------------------------------------------------------
# Helper: download, extract, and rename an archive into a target directory.
#   download_db <url> <target_dir> <glob_prefix> [insecure]
#
#   url          — download URL
#   target_dir   — final directory name inside DB_ROOT (e.g. "phrogs")
#   glob_prefix  — prefix of the ffdata/ffindex files to locate after
#                  extraction, used to find the versioned directory
#   insecure     — passed to fetch (skip TLS certificate checks)
#
# The archive is extracted into a staging directory and only moved to
# DB_ROOT/target_dir once complete, so an interrupted run never leaves a
# half-extracted database behind.
# ---------------------------------------------------------------------------
download_db() {
    local url="$1"
    local target_dir="$2"
    local glob_prefix="$3"
    local insecure="${4:-}"
    local archive="${url##*/}"
    local staging="$DB_ROOT/.staging_$target_dir"

    if [[ -d "$DB_ROOT/$target_dir" ]]; then
        echo "  [skip] $DB_ROOT/$target_dir already exists"
        return 0
    fi

    echo "  Downloading $archive ..."
    fetch "$url" "$DOWNLOAD_DIR/$archive" "$insecure"

    echo "  Extracting ..."
    rm -rf "$staging"
    mkdir -p "$staging"
    tar xzf "$DOWNLOAD_DIR/$archive" -C "$staging"

    # Find the extracted directory by looking for the ffdata files.
    # The archive may unpack into a versioned directory (e.g. pfamA_35.0/).
    local ffdata_file
    ffdata_file="$(find "$staging" -maxdepth 3 -name "${glob_prefix}*_a3m.ffdata" -print -quit)"

    if [[ -z "$ffdata_file" ]]; then
        echo "  [error] Could not find ${glob_prefix}*_a3m.ffdata after extraction."
        echo "          Please check the archive contents and set up this database manually."
        return 1
    fi

    local extracted
    extracted="$(dirname "$ffdata_file")"

    echo "  Moving $(basename "$extracted") -> $target_dir"
    mv "$extracted" "$DB_ROOT/$target_dir"
    rm -rf "$staging"
    rm -f "$DOWNLOAD_DIR/$archive"

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
echo "[1/5] PHROGs v4 (HHsuite database)"
download_db \
    "https://phrogs.lmge.uca.fr/downloads_from_website/phrogs_hhsuite_db.tar.gz" \
    "phrogs" \
    "phrogs" \
    insecure
echo ""

# ============================= PHROGs metadata ==============================
echo "[2/5] PHROGs annotation table"
TABLES_DIR="$DB_ROOT/tables"
mkdir -p "$TABLES_DIR"
if [[ -f "$TABLES_DIR/phrog_annot_v4.tsv" ]]; then
    echo "  [skip] $TABLES_DIR/phrog_annot_v4.tsv already exists"
else
    echo "  Downloading ..."
    fetch "https://phrogs.lmge.uca.fr/downloads_from_website/phrog_annot_v4.tsv" \
        "$TABLES_DIR/phrog_annot_v4.tsv" insecure
    echo "  Done."
fi
echo ""

# ============================= Pfam-A =======================================
echo "[3/5] Pfam-A (HHsuite database)"
download_db \
    "https://wwwuser.gwdguser.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_35.0.tar.gz" \
    "pfam" \
    "pfam"
echo ""

# ============================= ECOD =========================================
echo "[4/5] ECOD (HHsuite database)"
download_db \
    "http://prodata.swmed.edu/ecod/distributions/ecod.v294.F40.hhm_db.tar.gz" \
    "ecod" \
    "ecod"
echo ""

# ============================= Glimmer ICM ==================================
echo "[5/5] Glimmer training file"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ICM_SRC="$SCRIPT_DIR/resources/training-file_refseq.icm"
ICM_DST="$DB_ROOT/training-file_refseq.icm"

if [[ -f "$ICM_DST" ]]; then
    echo "  [skip] $ICM_DST already exists"
elif [[ -f "$ICM_SRC" ]]; then
    echo "  Copying from repository resources/ ..."
    cp "$ICM_SRC" "$ICM_DST"
    echo "  Done."
else
    echo "  [error] resources/training-file_refseq.icm not found in the repository."
    echo "          Please re-clone the repository or place the file manually at:"
    echo "          $ICM_DST"
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
rmdir "$DOWNLOAD_DIR" 2>/dev/null || true
echo " Set DB_ROOT in config.yml to: $DB_ROOT"
echo "============================================"