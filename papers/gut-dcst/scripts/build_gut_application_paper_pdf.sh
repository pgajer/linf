#!/bin/zsh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
. "$SCRIPT_DIR/setup_tool_paths.sh"
PAPER_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_DIR="$(cd "$PAPER_DIR/../.." && pwd)"
MANUSCRIPT_DIR="$PAPER_DIR/manuscript"
BUILD_DIR="$PAPER_DIR/build"
LATEXMK="${LATEXMK:-/Library/TeX/texbin/latexmk}"
INPUT_TEX="gut_application_paper.tex"
OUTPUT_PDF="$BUILD_DIR/gut_application_paper.pdf"
BUILD_INFO_TEX="$BUILD_DIR/manuscript_build_info.tex"
INCLUDE_BUILD_STAMP="${PAPER_INCLUDE_BUILD_STAMP:-1}"

escape_for_tex() {
  local s="$1"
  s=${s//_/\\_}
  printf '%s\n' "$s"
}

mkdir -p "$BUILD_DIR"

python3 "$SCRIPT_DIR/build_validation_assets.py"
python3 "$SCRIPT_DIR/build_ibd_results_overview_figure.py"

GIT_VERSION="$(git -C "$REPO_DIR" describe --tags --always --dirty 2>/dev/null || printf '%s\n' 'unversioned')"
GIT_BUILD_NUMBER="$(git -C "$REPO_DIR" rev-list --count HEAD 2>/dev/null || printf '%s\n' 'NA')"
GIT_COMMIT="$(git -C "$REPO_DIR" rev-parse --short HEAD 2>/dev/null || printf '%s\n' 'unknown')"
BUILD_DATETIME="$(date '+%Y-%m-%d %H:%M:%S %Z')"

cat > "$BUILD_INFO_TEX" <<EOF
\renewcommand{\manuscriptversion}{$(escape_for_tex "$GIT_VERSION")}
\renewcommand{\manuscriptbuildnumber}{$(escape_for_tex "$GIT_BUILD_NUMBER")}
\renewcommand{\manuscriptcommit}{$(escape_for_tex "$GIT_COMMIT")}
\renewcommand{\manuscriptbuilddatetime}{$(escape_for_tex "$BUILD_DATETIME")}
EOF

if [[ "$INCLUDE_BUILD_STAMP" == "1" ]]; then
  printf '%s\n' '\showbuildstamptrue' >> "$BUILD_INFO_TEX"
fi

cd "$MANUSCRIPT_DIR"

"$LATEXMK" \
  -xelatex \
  -interaction=nonstopmode \
  -halt-on-error \
  -file-line-error \
  -outdir="$BUILD_DIR" \
  "$INPUT_TEX"

printf '%s\n' "$OUTPUT_PDF"
