#!/usr/bin/env python3
"""Scan the R Journal source and rendered artifacts for release defects."""

from __future__ import annotations

import re
import subprocess
from html.parser import HTMLParser
from pathlib import Path

from article_assets import latex_figure_files


ROOT = Path(__file__).resolve().parents[1]
TEXT_FILES = [
    ROOT / "linf.Rmd",
    ROOT / "RJreferences.bib",
    ROOT / "README.md",
    ROOT / "READINESS.md",
    ROOT / "motivation-letter" / "motivation-letter.md",
    ROOT / "_Rpackages.txt",
]
TEXT_FILES.extend(sorted((ROOT / "scripts").glob("*.R")))
TEXT_FILES.extend(
    path
    for path in sorted((ROOT / "scripts").glob("*.py"))
    if path.name != Path(__file__).name
)

PATTERNS = {
    "private absolute path": re.compile(r"/Users/|/home/|[A-Za-z]:\\\\"),
    "prompt trace": re.compile(r"(?i)ignore previous|system prompt|chatgpt said|as an ai"),
    "unresolved marker": re.compile(r"(?i)\b(?:TODO|FIXME|TBD|XXX)\b|\[\s*insert\b"),
    "template residue": re.compile(r"(?i)quokka|bilby|penguins|ToOoOlTiPs"),
}

errors: list[str] = []
for path in TEXT_FILES:
    if not path.exists():
        errors.append(f"missing required source: {path.relative_to(ROOT)}")
        continue
    text = path.read_text(encoding="utf-8", errors="replace")
    for label, pattern in PATTERNS.items():
        for match in pattern.finditer(text):
            line = text.count("\n", 0, match.start()) + 1
            errors.append(
                f"{path.relative_to(ROOT)}:{line}: {label}: {match.group(0)!r}"
            )

pdf = ROOT / "output" / "pdf" / "linf-r-journal.pdf"
html = ROOT / "output" / "html" / "linf-r-journal.html"
for path in (pdf, html):
    if not path.exists() or path.stat().st_size == 0:
        errors.append(f"missing or empty rendered artifact: {path.relative_to(ROOT)}")

if pdf.exists():
    info = subprocess.run(["pdfinfo", str(pdf)], capture_output=True, text=True)
    if info.returncode != 0:
        errors.append("pdfinfo could not read the rendered PDF")
    else:
        match = re.search(r"^Pages:\s+(\d+)", info.stdout, re.MULTILINE)
        if not match:
            errors.append("could not determine PDF page count")
        elif int(match.group(1)) > 20:
            errors.append(f"rendered PDF exceeds 20 pages: {match.group(1)}")

    fonts = subprocess.run(["pdffonts", str(pdf)], capture_output=True, text=True)
    if fonts.returncode != 0:
        errors.append("pdffonts could not read the rendered PDF")
    else:
        for line in fonts.stdout.splitlines()[2:]:
            fields = line.split()
            if len(fields) >= 8 and fields[-5] == "no":
                errors.append("unembedded PDF font: " + line.strip())

    text_result = subprocess.run(["pdftotext", str(pdf), "-"], capture_output=True, text=True)
    if text_result.returncode != 0:
        errors.append("pdftotext could not read the rendered PDF")
    else:
        for pattern, label in [
            (r"\?\?", "unresolved cross-reference marker"),
            (r"\balt=", "figure-alt placeholder instead of a rendered image"),
            (r"ToOoOlTiPs|Quietest Quokka|Bounciest Bilby", "template residue"),
        ]:
            if re.search(pattern, text_result.stdout):
                errors.append(f"rendered PDF contains {label}")

tex = ROOT / "linf.tex"
if not tex.exists():
    errors.append("missing generated LaTeX source: linf.tex")
else:
    tex_text = tex.read_text(encoding="utf-8", errors="replace")
    if tex_text.count(r"\includegraphics") < 2:
        errors.append("generated LaTeX contains fewer than two plot inclusions")
    try:
        latex_figure_files(ROOT)
    except ValueError as error:
        errors.append(str(error))
    for command in ("begin", "end", "caption"):
        if rf"\textbackslash {command}" in tex_text:
            errors.append(f"generated LaTeX contains escaped {command} markup")
    if tex_text.count(r"\begin{table}") != 2 or tex_text.count(r"\end{table}") != 2:
        errors.append("generated LaTeX must contain two complete table environments")
    for label in ("eq:linf-normalization", "tab:interface", "tab:benchmark-table"):
        if rf"\label{{{label}}}" not in tex_text:
            errors.append(f"missing generated equation/table label: {label}")

for log in sorted(ROOT.glob("*.log")) + sorted((ROOT / "build").glob("*.log")):
    text = log.read_text(encoding="utf-8", errors="replace")
    for pattern, label in [
        (r"undefined references", "undefined references"),
        (r"Citation [`'][^`']+['`].*undefined", "undefined citation"),
        (r"Overfull \\hbox", "overfull box"),
    ]:
        if re.search(pattern, text, re.IGNORECASE):
            errors.append(f"{log.relative_to(ROOT)}: {label}")

if html.exists():
    html_text = html.read_text(encoding="utf-8", errors="replace")
    if re.search(r"citeproc-not-found|data-cites=\"\"", html_text):
        errors.append("rendered HTML contains unresolved citations")

    class AccessibilityParser(HTMLParser):
        def __init__(self):
            super().__init__()
            self.tables = []
            self.alt_images = 0

        def handle_starttag(self, tag, attrs):
            values = dict(attrs)
            if tag == "table":
                self.tables.append(values.get("aria-label", ""))
            if tag == "img" and values.get("alt", "").strip():
                self.alt_images += 1

    accessibility = AccessibilityParser()
    accessibility.feed(html_text)
    if len(accessibility.tables) != 2 or not all(accessibility.tables):
        errors.append("the two HTML tables must have nonempty accessibility labels")
    if accessibility.alt_images < 2:
        errors.append("the HTML figures must have nonempty alternative text")

if errors:
    print("Readiness scan failed:")
    for error in errors:
        print("-", error)
    raise SystemExit(1)

print("R Journal readiness scan passed.")
