#!/usr/bin/env python3
"""Download the pinned public source files and verify their SHA-256 hashes."""
import argparse, csv, hashlib, urllib.request
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument('directory',type=Path); a=p.parse_args()
a.directory.mkdir(parents=True,exist_ok=True)
manifest=Path(__file__).resolve().parents[1]/'inst/DATA_MANIFEST.csv'
for row in csv.DictReader(manifest.open()):
 if row['kind']!='upstream': continue
 target=a.directory/row['file']
 if target.exists(): data=target.read_bytes()
 else:
  with urllib.request.urlopen(row['url']) as response: data=response.read()
 if hashlib.sha256(data).hexdigest()!=row['sha256']:
  raise SystemExit(f'Checksum mismatch: {target}; no existing file was replaced')
 if not target.exists(): target.write_bytes(data)
 print(f'Verified {target.name}')
