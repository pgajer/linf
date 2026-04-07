from __future__ import annotations

import re


ITALIC_TAXA = (
    "Faecalibacterium prausnitzii",
    "Bacteroides vulgatus",
    "Escherichia-Shigella",
    "Bacteroides_H_857956",
    "Phocaeicola_A",
    "Bacteroides_H",
    "Pseudomonas_E",
    "Prevotella_9",
    "Prevotella_7",
    "Prevotella9",
    "Prevotella7",
    "F. prausnitzii",
    "B. vulgatus",
    "Acinetobacter",
    "Agathobacter",
    "Akkermansia",
    "Alistipes",
    "Bacteroides",
    "Bifidobacterium",
    "Blautia",
    "Corynebacterium",
    "Dialister",
    "Escherichia",
    "Faecalibacterium",
    "Lactobacillus",
    "Morganella",
    "Neisseria",
    "Parabacteroides",
    "Paraprevotella",
    "Phocaeicola",
    "Prevotella",
    "Proteus",
    "Pseudomonas",
    "Rothia",
    "Ruminococcus",
    "Segatella",
    "Shigella",
    "Staphylococcus",
    "Streptococcus",
    "Subdoligranulum",
    "Sutterella",
)

_TAXON_PATTERN = re.compile(
    r"(?<![A-Za-z])(?:"
    + "|".join(re.escape(label) for label in sorted(ITALIC_TAXA, key=len, reverse=True))
    + r")(?![A-Za-z])"
)


def _escape_mathtext(text: str) -> str:
    return text.replace("\\", r"\\").replace("_", r"\_").replace(" ", r"\ ")


def _to_mathtext(label: str) -> str:
    if label == "Escherichia-Shigella":
        return r"$\it{Escherichia}$-$\it{Shigella}$"
    return fr"$\it{{{_escape_mathtext(label)}}}$"


def italicize_taxa_mpl(text: str) -> str:
    if not text:
        return text
    return _TAXON_PATTERN.sub(lambda match: _to_mathtext(match.group(0)), text)
