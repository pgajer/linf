#!/usr/bin/env python3
"""Check Pandoc citation keys, BibTeX entries, and HTML evidence rows."""

from __future__ import annotations

import argparse
import re
from collections import Counter, defaultdict
from html.parser import HTMLParser
from pathlib import Path


ALLOWED = {
    "verified",
    "metadata-only",
    "unsupported",
    "wrong-source",
    "needs-human-review",
}


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def citation_keys(text: str) -> set[str]:
    # The negative lookbehind excludes email addresses and R Markdown's
    # escaped cross-reference form, \@ref(...). The project keys do not use
    # periods, so sentence-final punctuation cannot become part of a key.
    return set(re.findall(r"(?<![A-Za-z0-9_.\\])@([A-Za-z][A-Za-z0-9:_-]*)", text))


def bib_keys(text: str) -> set[str]:
    return set(re.findall(r"@\s*[A-Za-z]+\s*[({]\s*([^,\s]+)\s*,", text))


class EvidenceParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.entries: list[dict[str, object]] = []
        self.stack: list[dict[str, object]] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        values = {key: value or "" for key, value in attrs}
        if "data-citation-key" in values:
            entry: dict[str, object] = {
                "key": values["data-citation-key"].strip(),
                "status": values.get("data-status", "").strip(),
                "links": [],
                "tag": tag,
            }
            self.entries.append(entry)
            self.stack.append(entry)
        if self.stack and tag == "a" and "data-source-link" in values:
            self.stack[-1]["links"].append(values.get("href", "").strip())

    def handle_endtag(self, tag: str) -> None:
        if self.stack and self.stack[-1]["tag"] == tag:
            self.stack.pop()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--bib", required=True, type=Path)
    parser.add_argument("--html", required=True, type=Path)
    args = parser.parse_args()

    missing = [path for path in (args.source, args.bib, args.html) if not path.exists()]
    if missing:
        raise SystemExit("Missing citation input: " + ", ".join(map(str, missing)))

    cited = citation_keys(read(args.source))
    bibliography = bib_keys(read(args.bib))
    evidence_parser = EvidenceParser()
    evidence_parser.feed(read(args.html))
    entries = evidence_parser.entries

    errors: list[str] = []
    missing_bib = sorted(cited - bibliography)
    if missing_bib:
        errors.append("missing BibTeX keys: " + ", ".join(missing_bib))

    keys = [str(entry["key"]) for entry in entries]
    counts = Counter(keys)
    duplicates = sorted(key for key, count in counts.items() if count != 1)
    if duplicates:
        errors.append("verification keys not unique: " + ", ".join(duplicates))
    missing_evidence = sorted(cited - set(keys))
    if missing_evidence:
        errors.append("missing verification rows: " + ", ".join(missing_evidence))
    extra_evidence = sorted(set(keys) - cited)
    if extra_evidence:
        errors.append("uncited verification rows: " + ", ".join(extra_evidence))

    statuses: dict[str, list[str]] = defaultdict(list)
    for entry in entries:
        key = str(entry["key"])
        status = str(entry["status"])
        if status not in ALLOWED:
            errors.append(f"invalid status for {key}: {status or '<blank>'}")
        elif status != "verified":
            statuses[status].append(key)
        if key in cited and not any(str(link).strip() for link in entry["links"]):
            errors.append(f"missing source link for {key}")
    for status, status_keys in sorted(statuses.items()):
        errors.append(f"non-passing status {status}: " + ", ".join(sorted(status_keys)))

    if not cited:
        errors.append("no manuscript citations found")

    if errors:
        print("Citation verification failed:")
        for error in errors:
            print("-", error)
        return 1

    print(f"Citation verification passed for {len(cited)} cited sources.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
