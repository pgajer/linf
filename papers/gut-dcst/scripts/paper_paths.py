from __future__ import annotations

from pathlib import Path


SCRIPTS_DIR = Path(__file__).resolve().parent
PAPER_ROOT = SCRIPTS_DIR.parent
MANUSCRIPT_DIR = PAPER_ROOT / "manuscript"
ASSETS_DIR = PAPER_ROOT / "assets"
FIGURES_DIR = ASSETS_DIR / "figures"
TABLES_DIR = ASSETS_DIR / "tables"
SUPPLEMENT_DIR = ASSETS_DIR / "supplement"
BUILD_DIR = PAPER_ROOT / "build"
NOTES_DIR = PAPER_ROOT / "notes"
ARCHIVE_DIR = PAPER_ROOT / "archive"
PHASE1_ARCHIVE_DIR = ARCHIVE_DIR / "2026-03-24-phase1"

CURRENT_PROJECTS_ROOT = PAPER_ROOT.parents[2]
LINF_ROOT = CURRENT_PROJECTS_ROOT / "linf"
GUT_MICROBIOME_ROOT = CURRENT_PROJECTS_ROOT / "gut_microbiome"

DCST_ANALYSIS_DIR = GUT_MICROBIOME_ROOT / "outputs" / "dcst_analysis"
DCST_VALIDATION_DIR = GUT_MICROBIOME_ROOT / "outputs" / "dcst_validation"
