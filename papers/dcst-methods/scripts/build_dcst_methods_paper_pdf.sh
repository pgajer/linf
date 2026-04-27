#!/bin/zsh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PAPER_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
MANUSCRIPT_DIR="$PAPER_DIR/manuscript"
BUILD_DIR="$PAPER_DIR/build"
LATEXMK="${LATEXMK:-/Library/TeX/texbin/latexmk}"
INPUT_TEX="dcst_methods_paper.tex"
OUTPUT_PDF="$BUILD_DIR/dcst_methods_paper.pdf"

"$SCRIPT_DIR/write_manuscript_build_info.sh" >/dev/null

cd "$MANUSCRIPT_DIR"

"$LATEXMK" \
  -pdf \
  -interaction=nonstopmode \
  -halt-on-error \
  -file-line-error \
  -outdir="$BUILD_DIR" \
  "$INPUT_TEX"

printf '%s\n' "$OUTPUT_PDF"
