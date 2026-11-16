"""Fetch LiCSAR short (6-12 day) unwrapped interferograms over Svartsengi inside the inflation intervals.

Pairs are taken from the LiCSAR index of each frame when both dates lie inside one window of
config/svartsengi.json -> insar_short_pairs and the temporal baseline is <= max_days. Window reads over HTTP use
licsar_annual.fetch (nothing downloaded is executed; GeoTIFFs are read as data). Unwrapped phase is converted to LOS
mm with the LiCSBAS sign convention (positive = towards the satellite).

Output: data/raw/licsar_short/<frame>/{los_mm,cc}/<pair>.tif and manifest.json with the HTTP outcome of every file.
    python -I scripts/licsar_short_pairs.py
"""
from __future__ import annotations

import datetime as dt
import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from licsar_annual import INDEX, PUBLIC, fetch, pdate  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT / "config" / "svartsengi.json").read_text())["insar_short_pairs"]
OUT = ROOT / "data" / "raw" / "licsar_short"


def select_pairs(pairs, windows, max_days):
    """Pairs whose two dates are inside one window and at most max_days apart; returns {pair: interval}."""
    out = {}
    for p in pairs:
        a, b = map(pdate, p.split("_"))
        if not 0 < (b - a).days <= max_days:
            continue
        for w in windows:
            if dt.date.fromisoformat(w["start"]) <= a and b <= dt.date.fromisoformat(w["end"]):
                out[p] = w["interval"]
                break
    return out


def main():
    bbox = tuple(CFG["bbox"])
    manifest = {}
    for frame in CFG["frames"]:
        t = int(frame[:3])
        html = requests.get(f"{INDEX}/{t}/{frame}/interferograms/", timeout=180).text
        pairs = sorted(set(re.findall(r'href="(\d{8}_\d{8})/?"', html)))
        sel = select_pairs(pairs, CFG["windows"], CFG["max_days"])
        print(f"{frame}: {len(pairs)} pairs listed, {len(sel)} short pairs in the windows", flush=True)
        out = OUT / frame
        base = f"{PUBLIC}/{t}/{frame}"
        jobs = [("meta", f"{base}/metadata/{frame}.geo.{k}.tif", out / "meta" / f"{k}.tif") for k in ("U", "E", "N")]
        for p in sel:
            jobs.append(("unw", f"{base}/interferograms/{p}/{p}.geo.unw.tif", out / "los_mm" / f"{p}.tif"))
            jobs.append(("cc", f"{base}/interferograms/{p}/{p}.geo.cc.tif", out / "cc" / f"{p}.tif"))
        res = {}
        with ThreadPoolExecutor(6) as ex:
            futs = {ex.submit(fetch, k, u, o, bbox): f"{o.parent.name}/{o.name}" for k, u, o in jobs}
            for f in as_completed(futs):
                res[futs[f]] = f.result()
        bad = {k: v for k, v in res.items() if v.startswith("error")}
        print(f"{frame}: {len(res) - len(bad)} ok/cached, {len(bad)} errors", flush=True)
        manifest[frame] = dict(listed=len(pairs), selected=sel, errors=bad)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "manifest.json").write_text(json.dumps(dict(bbox=bbox, windows=CFG["windows"], frames=manifest), indent=1))


if __name__ == "__main__":
    main()
