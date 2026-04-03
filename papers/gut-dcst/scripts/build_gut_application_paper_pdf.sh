#!/bin/zsh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PAPER_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
MANUSCRIPT_DIR="$PAPER_DIR/manuscript"
BUILD_DIR="$PAPER_DIR/build"
LATEXMK="${LATEXMK:-/Library/TeX/texbin/latexmk}"
INPUT_TEX="gut_application_paper.tex"
OUTPUT_PDF="$BUILD_DIR/gut_application_paper.pdf"

mkdir -p "$BUILD_DIR"
cd "$MANUSCRIPT_DIR"

"$LATEXMK" \
  -xelatex \
  -interaction=nonstopmode \
  -halt-on-error \
  -file-line-error \
  -outdir="$BUILD_DIR" \
  "$INPUT_TEX"

print "$OUTPUT_PDF"
