#!/bin/zsh
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ./archive_gut_application_paper_snapshot.sh <milestone-label>

Examples:
  ./archive_gut_application_paper_snapshot.sh draft-1
  ./archive_gut_application_paper_snapshot.sh submission-1
  ./archive_gut_application_paper_snapshot.sh resubmission-1
EOF
}

if [[ "${1:-}" == "" || "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  usage
  exit 0
fi

LABEL_RAW="$1"
LABEL="$(print -r -- "$LABEL_RAW" | tr '[:upper:]' '[:lower:]' | tr ' ' '-' | tr -cs 'a-z0-9._-' '-')"
LABEL="${LABEL#-}"
LABEL="${LABEL%-}"

if [[ -z "$LABEL" ]]; then
  print -u2 "Milestone label resolved to empty after sanitization."
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PAPER_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
MANUSCRIPT_DIR="$PAPER_DIR/manuscript"
BUILD_DIR="$PAPER_DIR/build"
ARCHIVE_DIR="$PAPER_DIR/archive"
DATE_STAMP="$(date +%F)"
SNAPSHOT_DIR="$ARCHIVE_DIR/${DATE_STAMP}-${LABEL}"
REPO_ROOT="$(git -C "$PAPER_DIR" rev-parse --show-toplevel)"
CURRENT_COMMIT="$(git -C "$PAPER_DIR" rev-parse --short HEAD)"
CURRENT_BRANCH="$(git -C "$PAPER_DIR" rev-parse --abbrev-ref HEAD)"

if [[ -e "$SNAPSHOT_DIR" ]]; then
  print -u2 "Snapshot directory already exists: $SNAPSHOT_DIR"
  exit 1
fi

mkdir -p "$SNAPSHOT_DIR"

cp "$MANUSCRIPT_DIR/gut_application_paper.tex" \
  "$SNAPSHOT_DIR/gut_application_paper.tex"
cp "$MANUSCRIPT_DIR/references.bib" \
  "$SNAPSHOT_DIR/references.bib"

if [[ -f "$BUILD_DIR/gut_application_paper.pdf" ]]; then
  cp "$BUILD_DIR/gut_application_paper.pdf" \
    "$SNAPSHOT_DIR/gut_application_paper.pdf"
  PDF_STATUS="Included current built PDF."
else
  PDF_STATUS="No built PDF was present at snapshot time."
fi

cat > "$SNAPSHOT_DIR/README.md" <<EOF
# Gut DCST Milestone Snapshot

- Milestone: \`$LABEL\`
- Date: \`$DATE_STAMP\`
- Git branch: \`$CURRENT_BRANCH\`
- Git commit: \`$CURRENT_COMMIT\`
- Repo root: \`$REPO_ROOT\`

$PDF_STATUS

Files in this snapshot are frozen milestone copies. Continue normal manuscript
editing in:
\`$MANUSCRIPT_DIR/gut_application_paper.tex\`
EOF

print "$SNAPSHOT_DIR"
