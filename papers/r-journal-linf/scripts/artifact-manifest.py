#!/usr/bin/env python3
"""Record or verify the correspondence between article sources and outputs."""

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

from article_assets import latex_figure_files

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "build" / "artifact-manifest.json"


def hashes(names: list[str]) -> dict[str, str]:
    return {
        name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        for name in names
    }


mode = sys.argv[1] if len(sys.argv) == 2 else ""
if mode == "record":
    pin = json.loads((ROOT / "package-source.json").read_text())
    inputs = [
        "Makefile", "linf.Rmd", "RJreferences.bib", "RJournal.sty", "package-source.json",
        "citation_verification.html",
        "_Rpackages.txt", "motivation-letter/motivation-letter.md",
        f"package/{pin['package']}_{pin['version']}.tar.gz",
    ]
    inputs += sorted(
        str(path.relative_to(ROOT))
        for path in (ROOT / "scripts").glob("*")
        if path.suffix in (".py", ".R")
    )
    outputs = [
        "linf.pdf", "linf.html", "linf.tex", "linf.R", "RJwrapper.tex",
        "output/pdf/linf-r-journal.pdf", "output/html/linf-r-journal.html",
        "motivation-letter/motivation-letter.pdf",
        "build/benchmark-results.csv", "build/session-info.txt",
        "build/benchmark-environment.json",
    ]
    outputs += latex_figure_files(ROOT)
    MANIFEST.write_text(json.dumps({
        "recorded_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": hashes(inputs),
        "outputs": hashes(outputs),
    }, indent=2) + "\n")
elif mode == "verify":
    manifest = json.loads(MANIFEST.read_text())
    for name in latex_figure_files(ROOT):
        if name not in manifest["outputs"]:
            raise SystemExit(f"Figure omitted from artifact manifest: {name}")
    for group in ("inputs", "outputs"):
        for name, expected in manifest[group].items():
            path = ROOT / name
            if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                raise SystemExit(f"Source/output drift: {name}; rerun make audit")
else:
    raise SystemExit("Usage: artifact-manifest.py record|verify")

print(f"Artifact manifest {mode} passed.")
