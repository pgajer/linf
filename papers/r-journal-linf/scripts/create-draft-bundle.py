#!/usr/bin/env python3
"""Create and independently rebuild the checksum-pinned R Journal bundle."""

import hashlib
import json
import os
import subprocess
import tempfile
import time
import zipfile
from pathlib import Path

from article_assets import latex_figure_files

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "output" / "linf-r-journal-draft.zip"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


subprocess.run(["python3", "scripts/fetch-package.py"], cwd=ROOT, check=True)
subprocess.run(
    ["python3", "scripts/artifact-manifest.py", "verify"], cwd=ROOT, check=True
)
pin = json.loads((ROOT / "package-source.json").read_text())
manifest = json.loads((ROOT / "build" / "artifact-manifest.json").read_text())
required = [
    "Makefile", "linf.Rmd", "linf.pdf", "linf.html", "linf.R", "linf.tex", "RJwrapper.tex",
    "RJournal.sty", "RJreferences.bib", "_Rpackages.txt", "README.md",
    "READINESS.md", "citation_verification.html",
    "package-source.json", f"package/{pin['package']}_{pin['version']}.tar.gz",
    "motivation-letter/motivation-letter.md",
    "motivation-letter/motivation-letter.pdf",
    "build/artifact-manifest.json", "build/render-info.txt",
    "build/benchmark-results.csv", "build/session-info.txt",
    "build/benchmark-environment.json",
]
required += sorted(
    str(path.relative_to(ROOT))
    for path in (ROOT / "scripts").glob("*")
    if path.suffix in (".py", ".R")
)
required += latex_figure_files(ROOT)
missing = [name for name in required if not (ROOT / name).is_file()]
if missing:
    raise SystemExit("Missing bundle input(s): " + ", ".join(missing))
checksums = {name: sha256(ROOT / name) for name in required}
checksum_text = "".join(f"{digest}  {name}\n" for name, digest in checksums.items())

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
with tempfile.TemporaryDirectory(prefix="linf-rjournal-bundle-") as tmp:
    staging = Path(tmp)
    candidate = staging / OUTPUT.name
    with zipfile.ZipFile(candidate, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name in required:
            archive.write(ROOT / name, arcname=name)
        archive.writestr("SHA256SUMS", checksum_text)

    # Include spaces in the extraction path to exercise the documented commands.
    extracted = staging / "review copy"
    with zipfile.ZipFile(candidate) as archive:
        archive.extractall(extracted)
    for name, expected in checksums.items():
        if sha256(extracted / name) != expected:
            raise SystemExit(f"Archived file differs from its source: {name}")
    if latex_figure_files(extracted) != latex_figure_files(ROOT):
        raise SystemExit("The archive does not contain all referenced TeX figures")

    # Force regeneration: an old PDF, HTML file, or plot cannot satisfy this test.
    for name in manifest["outputs"]:
        path = extracted / name
        if not path.resolve().is_relative_to(extracted.resolve()):
            raise SystemExit(f"Output path escapes the extracted archive: {name}")
        if path.is_file():
            path.unlink()

    clean_env = os.environ.copy()
    for name in (
        "R_LIBS", "R_LIBS_USER", "R_LIBS_SITE", "R_PROFILE", "R_PROFILE_USER",
        "R_ENVIRON", "R_ENVIRON_USER",
    ):
        clean_env.pop(name, None)
    clean_env["R_ENVIRON_USER"] = "/dev/null"
    clean_env["R_PROFILE_USER"] = "/dev/null"
    clean_env["R_MAKEVARS_USER"] = "/dev/null"
    # Dependency installation is a prerequisite, not part of analysis timing.
    prerequisites = subprocess.run(
        ["make", "dependencies", f"LIBRARY={extracted / 'review-library'}"],
        cwd=extracted, text=True, capture_output=True, env=clean_env,
    )
    if prerequisites.returncode != 0:
        print(prerequisites.stdout)
        print(prerequisites.stderr)
        raise SystemExit("Could not prepare the clean review library dependencies")
    started = time.monotonic()
    result = subprocess.run(
        ["make", "audit", f"LIBRARY={extracted / 'review-library'}"],
        cwd=extracted, text=True, capture_output=True, env=clean_env,
    )
    elapsed = time.monotonic() - started
    (ROOT / "build" / "bundle-check.log").write_text(result.stdout + result.stderr)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise SystemExit("The extracted bundle failed its clean audit")
    for name in manifest["outputs"]:
        if not (extracted / name).is_file():
            raise SystemExit(f"The extracted bundle did not reproduce {name}")
    if elapsed >= 600:
        raise SystemExit("The extracted bundle took at least ten minutes to reproduce")
    if candidate.stat().st_size > 10_000_000:
        raise SystemExit("The bundle exceeds the journal's 10 MB guidance")
    OUTPUT.write_bytes(candidate.read_bytes())

print(
    f"Created and independently rebuilt {OUTPUT} with {len(required) + 1} files "
    f"in {elapsed:.1f} seconds; SHA-256 {sha256(OUTPUT)}."
)
