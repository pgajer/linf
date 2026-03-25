#!/bin/zsh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PANDOC="${PANDOC:-/Users/pgajer/bin/pandoc}"
PDF_ENGINE="${PDF_ENGINE:-/Library/TeX/texbin/xelatex}"
INPUT_MD="$SCRIPT_DIR/gut_application_paper_draft.md"
OUTPUT_PDF="$SCRIPT_DIR/gut_application_paper_draft.pdf"

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
