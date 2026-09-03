"""Resolve and validate figure files referenced by the generated article TeX."""

import re
from pathlib import Path


def latex_figure_files(root: Path) -> list[str]:
    root = root.resolve()
    tex = (root / "linf.tex").read_text(encoding="utf-8")
    references = re.findall(
        r"\\includegraphics\*?(?:\[[^\]]*\])?\s*\{([^}]+)\}", tex
    )
    if not references:
        raise ValueError("The article TeX has no figure references")
    files = set()
    for reference in references:
        path = (root / reference).resolve()
        if not path.is_relative_to(root):
            raise ValueError(f"Figure reference escapes the article directory: {reference}")
        candidates = [path] if path.suffix else [
            path.with_suffix(suffix) for suffix in (".pdf", ".png", ".jpg", ".jpeg")
        ]
        found = next((p for p in candidates if p.is_file() and p.stat().st_size), None)
        if found is None:
            raise ValueError(f"Missing TeX figure asset: {reference}")
        files.add(str(found.relative_to(root)))
    return sorted(files)
