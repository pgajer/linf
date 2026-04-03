#!/bin/zsh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PAPER_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
MANUSCRIPT_DIR="$PAPER_DIR/manuscript"
PANDOC="${PANDOC:-/Users/pgajer/bin/pandoc}"
INPUT_MD="gut_application_paper.md"
OUTPUT_TEX="gut_application_paper.tex"

cd "$MANUSCRIPT_DIR"

"$PANDOC" "$INPUT_MD" \
  --standalone \
  --from markdown+yaml_metadata_block+implicit_figures \
  --to latex \
  --citeproc \
  -V colorlinks=true \
  -V linkcolor=blue \
  -V urlcolor=blue \
  -o "$OUTPUT_TEX"

print "$MANUSCRIPT_DIR/$OUTPUT_TEX"
