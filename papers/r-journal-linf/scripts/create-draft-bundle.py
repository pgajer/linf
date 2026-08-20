#!/usr/bin/env python3
"""Create a minimal, reproducible R Journal review bundle."""

from __future__ import annotations

import subprocess
import tempfile
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "output" / "linf-r-journal-draft.zip"

required = [
    "linf.Rmd",
    "linf.pdf",
    "linf.html",
    "linf.R",
    "linf.tex",
    "RJournal.sty",
    "RJreferences.bib",
    "_Rpackages.txt",
    "README.md",
    "motivation-letter/motivation-letter.md",
    "motivation-letter/motivation-letter.pdf",
    "scripts/render-paper.R",
]
missing = [name for name in required if not (ROOT / name).is_file()]
if missing:
    raise SystemExit("Missing bundle input(s): " + ", ".join(missing))

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
with zipfile.ZipFile(OUTPUT, "w", compression=zipfile.ZIP_DEFLATED) as archive:
    for name in required:
        archive.write(ROOT / name, arcname=name)

with tempfile.TemporaryDirectory(prefix="linf-rjournal-bundle-") as tmp:
    extracted = Path(tmp)
    with zipfile.ZipFile(OUTPUT) as archive:
        archive.extractall(extracted)
    result = subprocess.run(
        ["Rscript", "-e", "rmarkdown::render('linf.Rmd', output_format='rjtools::rjournal_pdf_article', clean=TRUE)"],
        cwd=extracted,
        text=True,
        capture_output=True,
    )
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise SystemExit("The extracted review bundle did not compile")
    if not (extracted / "linf.pdf").is_file():
        raise SystemExit("The extracted review bundle did not produce linf.pdf")

print(f"Created and clean-built {OUTPUT} with {len(required)} files.")
