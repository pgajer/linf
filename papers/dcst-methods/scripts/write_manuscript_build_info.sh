#!/bin/zsh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PAPER_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_DIR="$(cd "$PAPER_DIR/../.." && pwd)"
BUILD_DIR="$PAPER_DIR/build"
BUILD_INFO_TEX="$BUILD_DIR/manuscript_build_info.tex"
INCLUDE_BUILD_STAMP="${PAPER_INCLUDE_BUILD_STAMP:-1}"

escape_for_tex() {
  local s="$1"
  s=${s//_/\\_}
  printf '%s\n' "$s"
}

mkdir -p "$BUILD_DIR"

GIT_VERSION="$(git -C "$REPO_DIR" describe --tags --always 2>/dev/null || printf '%s\n' 'unversioned')"
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

printf '%s\n' "$BUILD_INFO_TEX"
