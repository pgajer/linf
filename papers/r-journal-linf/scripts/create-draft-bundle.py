#!/usr/bin/env python3
"""Create and independently clean-build the R Journal review bundle."""

from __future__ import annotations

import hashlib
import os
import subprocess
import tempfile
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = ROOT.parents[1]
OUTPUT = ROOT / "output" / "linf-r-journal-draft.zip"

required = [
    "Makefile",
    "linf.Rmd",
    "linf.pdf",
    "linf.html",
    "linf.R",
    "linf.tex",
    "RJournal.sty",
    "RJreferences.bib",
    "_Rpackages.txt",
    "README.md",
    "ACTION_PLAN.md",
    "READINESS.md",
    "citation_verification.html",
    "motivation-letter/motivation-letter.md",
    "motivation-letter/motivation-letter.pdf",
    "scripts/render-paper.R",
    "scripts/check-citation-verification.py",
    "scripts/run-rjtools-checks.R",
    "scripts/readiness-scan.py",
]
missing = [name for name in required if not (ROOT / name).is_file()]
if missing:
    raise SystemExit("Missing bundle input(s): " + ", ".join(missing))

published_pdf = ROOT / "output" / "pdf" / "linf-r-journal.pdf"
if not published_pdf.is_file():
    raise SystemExit("Missing current rendered PDF: output/pdf/linf-r-journal.pdf")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


if sha256(ROOT / "linf.pdf") != sha256(published_pdf):
    raise SystemExit(
        "linf.pdf is not the current output/pdf/linf-r-journal.pdf; rerun make audit"
    )

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
with tempfile.TemporaryDirectory(prefix="linf-rjournal-package-") as package_tmp:
    build = subprocess.run(
        ["R", "CMD", "build", str(PACKAGE_ROOT)],
        cwd=package_tmp,
        text=True,
        capture_output=True,
    )
    if build.returncode != 0:
        print(build.stdout)
        print(build.stderr)
        raise SystemExit("Could not build the exact linf source tarball")
    tarballs = sorted(Path(package_tmp).glob("linf_*.tar.gz"))
    if len(tarballs) != 1:
        raise SystemExit("Expected exactly one linf source tarball")
    package_tarball = tarballs[0]

    with zipfile.ZipFile(OUTPUT, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name in required:
            archive.write(ROOT / name, arcname=name)
        archive.write(package_tarball, arcname=f"package/{package_tarball.name}")

with tempfile.TemporaryDirectory(prefix="linf-rjournal-bundle-") as tmp:
    extracted = Path(tmp)
    with zipfile.ZipFile(OUTPUT) as archive:
        archive.extractall(extracted)

    if sha256(extracted / "linf.pdf") != sha256(published_pdf):
        raise SystemExit("The archived PDF differs from the current rendered PDF")

    clean_env = os.environ.copy()
    for name in ("R_LIBS", "R_LIBS_USER", "R_PROFILE", "R_PROFILE_USER"):
        clean_env.pop(name, None)
    clean_env["R_ENVIRON_USER"] = "/dev/null"
    clean_env["R_PROFILE_USER"] = "/dev/null"
    clean_env["R_MAKEVARS_USER"] = "/dev/null"
    review_library = extracted / "review-library"
    result = subprocess.run(
        ["make", "audit", f"LIBRARY={review_library}"],
        cwd=extracted,
        text=True,
        capture_output=True,
        env=clean_env,
    )
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise SystemExit("The extracted review bundle did not pass its clean audit")
    if not (extracted / "output" / "pdf" / "linf-r-journal.pdf").is_file():
        raise SystemExit("The extracted review bundle did not reproduce the PDF")

print(
    f"Created and independently clean-built {OUTPUT} "
    f"with {len(required) + 1} files."
)
