#!/usr/bin/env bash
set -euo pipefail

# Initialise scaffold directories for user data/references in a fresh clone.
# Creates folders and .gitkeep placeholders so Git tracks the structure.

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"

mkdir -p "$repo_root/data/user_fastqs" "$repo_root/ref/user_ref"
touch "$repo_root/data/user_fastqs/.gitkeep" "$repo_root/ref/user_ref/.gitkeep"

printf 'Created: %s and %s\n' "data/user_fastqs/" "ref/user_ref/"
printf 'Drop your FASTQs in data/user_fastqs/\n'
printf 'Put your reference FASTA in ref/user_ref/ (indexes optional; auto-indexed if missing)\n'
