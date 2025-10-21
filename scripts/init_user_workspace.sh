#!/usr/bin/env bash
set -euo pipefail
mkdir -p data/user_fastqs ref/user_ref
touch data/user_fastqs/.gitkeep ref/user_ref/.gitkeep
echo "Created: data/user_fastqs/ and ref/user_ref/"
echo "Drop your FASTQs in data/user_fastqs/"
echo "Put your reference FASTA in ref/user_ref/ (indexes optional; pipeline auto-indexes if missing)"
