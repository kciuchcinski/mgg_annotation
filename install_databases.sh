#!/usr/bin/env bash
set -euo pipefail

if [ $# -ne 1 ]; then
  echo "Usage: $0 /absolute/path/to/db_root"
  exit 1
fi

DB_ROOT=$(readlink -f "$1")

echo "Installing databases into: $DB_ROOT"

# Example structure: adapt URLs/paths to your real sources
mkdir -p \
  "$DB_ROOT/hhsuite/phrogs" \
  "$DB_ROOT/hhsuite/pfam" \
  "$DB_ROOT/hhsuite/ecod" \
  "$DB_ROOT/hhsuite/alandb" \
  "$DB_ROOT/tables"

# Ideally, the whole bundle would be downloaded from a single source like zenodo, then unpacked.

# GLIMMER model
# cp or wget/scp the ICM file
# cp /some/source/training-file_refseq.icm "$DB_ROOT/glimmer/training-file_refseq.icm"

# PHROGS
# wget ... -O "$DB_ROOT/hhsuite/phrogs/phrogs_hhsuite.tar.gz"
# tar -C "$DB_ROOT/hhsuite/phrogs" -xzf "$DB_ROOT/hhsuite/phrogs/phrogs_hhsuite.tar.gz"

# PFAM, ECOD, ALANDB similarly ...

# Metadata tables
# cp dependencies/tables/phrog_annot_v4.tsv "$DB_ROOT/metadata/"
# cp dependencies/tables/v3_phrogs-table-rafal-3_12.csv "$DB_ROOT/metadata/"
# cp dependencies/tables/alan_annot.tsv "$DB_ROOT/metadata/"

echo "Done. Use:  --DB_ROOT $DB_ROOT"
