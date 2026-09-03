#!/usr/bin/env python3
"""Obtain the exact CRAN source used by the article, with checksum validation."""

import hashlib
import io
import json
import re
import tarfile
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
pin = json.loads((ROOT / "package-source.json").read_text())
name, version = pin["package"], pin["version"]
target = ROOT / "package" / f"{name}_{version}.tar.gz"


def validate(data: bytes) -> None:
    if hashlib.sha256(data).hexdigest() != pin["sha256"]:
        raise SystemExit("Source checksum mismatch; refusing to use an unverified package")
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as archive:
        description = archive.extractfile(f"{name}/DESCRIPTION").read().decode()
    for key, expected in (("Package", name), ("Version", version)):
        match = re.search(rf"^{key}:\s*(.+)$", description, re.MULTILINE)
        if not match or match.group(1).strip() != expected:
            raise SystemExit(f"Bundled DESCRIPTION does not match {key}: {expected}")


if target.exists():
    validate(target.read_bytes())
else:
    failures = []
    for url in pin["urls"]:
        try:
            with urlopen(url, timeout=60) as response:
                content = response.read()
        except (URLError, TimeoutError) as error:
            failures.append(f"{url}: {error}")
            continue
        validate(content)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(content)
        break
    else:
        raise SystemExit("Could not obtain the pinned CRAN source:\n" + "\n".join(failures))

print(f"Verified {target.relative_to(ROOT)} (SHA-256 {pin['sha256']})")
