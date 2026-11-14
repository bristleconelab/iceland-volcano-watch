"""Download the Zenodo "Hengill seismic catalog" and the COSEISMIQ catalogues it was derived from.

  Zenodo 10.5281/zenodo.16311162 (University of Campania, 2025, CC BY 4.0): hen_conv.dat, 5 whitespace columns,
      described as "time (sec), latitude, longitude, event magnitude, depth".
  COSEISMIQ (Grigoli et al. 2022, Scientific Data 9:220; figshare collection 10.6084/m9.figshare.c.5725592, CC0):
      catalogue_HQDD.csv (high-quality, double-difference relocated), catalogue_HQ.csv, and the HQDD legend.

Files are stored under data/raw/hengill/ with a manifest (URL, bytes, sha256). They are only parsed as text.
"""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "raw" / "hengill"
UA = {"User-Agent": "iceland-volcano-watch/0.1 (research script)"}
FILES = {
    "hen_conv.dat": "https://zenodo.org/api/records/16311162/files/hen_conv.dat/content",
    "catalogue_HQDD.csv": "https://ndownloader.figshare.com/files/34138983",
    "catalogue_HQ.csv": "https://ndownloader.figshare.com/files/34138974",
    "catalogue_HQDD_legend.txt": "https://ndownloader.figshare.com/files/34138977",
}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    man = []
    for name, url in FILES.items():
        p = OUT / name
        if not p.exists():
            r = requests.get(url, headers=UA, timeout=300)
            r.raise_for_status()
            p.write_bytes(r.content)
        b = p.read_bytes()
        man.append(dict(file=name, url=url, bytes=len(b), sha256=hashlib.sha256(b).hexdigest(),
                        checked_utc=datetime.now(timezone.utc).isoformat(timespec="seconds")))
        print(name, len(b), man[-1]["sha256"][:16])
    (OUT / "manifest.json").write_text(json.dumps(man, indent=1))


if __name__ == "__main__":
    main()
