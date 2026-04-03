#!/bin/zsh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PAPER_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
MANUSCRIPT_DIR="$PAPER_DIR/manuscript"
BUILD_DIR="$PAPER_DIR/build"
PANDOC="${PANDOC:-/Users/pgajer/bin/pandoc}"
PDF_ENGINE="${PDF_ENGINE:-/Library/TeX/texbin/xelatex}"
INPUT_MD="gut_application_paper.md"
OUTPUT_PDF="$BUILD_DIR/gut_application_paper.pdf"

mkdir -p "$BUILD_DIR"
cd "$MANUSCRIPT_DIR"

"$PANDOC" "$INPUT_MD" \
  --standalone \
  --from markdown+yaml_metadata_block+implicit_figures \
  --citeproc \
  --pdf-engine="$PDF_ENGINE" \
  -V colorlinks=true \
  -V linkcolor=blue \
  -V urlcolor=blue \
  -o "$OUTPUT_PDF"

print "$OUTPUT_PDF"
