"""Download LiCSAR summer-to-summer (~1 year) unwrapped interferograms over the Reykjanes Peninsula and Hengill.

Why annual summer pairs: Iceland has snow from October to May, so winter pairs decorrelate. Pairs whose two dates are
both in June-September and ~1 year apart keep coherence on lava and moss and give one ~annual displacement each.

The window read approach (HTTP range reads of strip GeoTIFFs, nothing executed) follows the sibling insar-watch
`licsar_fetch.py`. Unwrapped phase is converted to LOS displacement with the LiCSBAS convention
    d_LOS [mm] = -unw * lambda / (4 pi) * 1000,   lambda = 55.465763 mm (Sentinel-1 C band)
so that positive values mean motion towards the satellite.

Example:
    python -I scripts/licsar_annual.py --frame 016A_02504_162118 --bbox -22.75 63.80 -21.10 64.16 \
        --years 2015-2025 --per-year 3 --out data/raw/licsar/016A_02504_162118
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import rasterio
import requests
from rasterio.windows import Window, from_bounds

INDEX = "https://gws-access.jasmin.ac.uk/public/nceo_geohazards/LiCSAR_products"
PUBLIC = "https://gws-access.jasmin.ac.uk/public/nceo_geohazards/LiCSAR_products.public"
GDAL_ENV = dict(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
                GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES", GDAL_HTTP_MULTIRANGE="YES", GDAL_HTTP_MAX_RETRY="4",
                GDAL_HTTP_RETRY_DELAY="3", GDAL_HTTP_TIMEOUT="120", VSI_CACHE="TRUE")


def pdate(s: str) -> dt.date:
    return dt.datetime.strptime(s, "%Y%m%d").date()


def summer_pairs(pairs, years, per_year, months=(6, 9), bt=(330, 400)):
    by_year = {}
    for p in pairs:
        a, b = map(pdate, p.split("_"))
        if not (months[0] <= a.month <= months[1] and months[0] <= b.month <= months[1]):
            continue
        d = (b - a).days
        if not (bt[0] <= d <= bt[1]) or a.year not in years:
            continue
        by_year.setdefault(a.year, []).append((abs(d - 365), abs(a.timetuple().tm_yday - 213), p))
    sel = []
    for y, lst in sorted(by_year.items()):
        lst.sort()  # closest to one year, then closest to 1 August
        # keep distinct primaries when possible
        seen, k = set(), []
        for _, _, p in lst:
            if p[:8] in seen:
                continue
            seen.add(p[:8])
            k.append(p)
            if len(k) == per_year:
                break
        sel += k
    return sel


def read_window(url, bbox):
    with rasterio.Env(**GDAL_ENV):
        with rasterio.open("/vsicurl/" + url) as ds:
            w = from_bounds(*bbox, transform=ds.transform).round_offsets().round_lengths()
            w = w.intersection(Window(0, 0, ds.width, ds.height))
            return ds.read(1, window=w), ds.window_transform(w), ds.crs


def write_tif(path: Path, arr, tr, crs):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".part.tif")
    with rasterio.open(tmp, "w", driver="GTiff", height=arr.shape[0], width=arr.shape[1], count=1, dtype=arr.dtype,
                       crs=crs, transform=tr, compress="deflate", tiled=True, blockxsize=256, blockysize=256) as dst:
        dst.write(arr, 1)
    tmp.replace(path)


def fetch(kind, url, out, bbox):
    if out.exists():
        return "cached"
    err = ""
    for k in range(3):
        try:
            arr, tr, crs = read_window(url, bbox)
            break
        except Exception as e:  # network hiccup
            err = str(e)
            time.sleep(5 * (k + 1))
    else:
        return "error: " + err[:150]
    arr = arr.astype(np.float32)
    if kind == "unw":
        arr[arr == 0] = np.nan
        arr = -arr * 55.465763 / (4 * np.pi)  # mm, positive towards satellite
    elif kind == "cc":
        if np.nanmax(arr) > 1.5:
            arr = arr / 255.0
        arr[arr == 0] = np.nan
    write_tif(out, arr, tr, crs)
    return "ok"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frame", required=True)
    ap.add_argument("--bbox", nargs=4, type=float, required=True)
    ap.add_argument("--years", default="2015-2025")
    ap.add_argument("--per-year", type=int, default=3)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    y0, y1 = map(int, a.years.split("-"))
    t, frame, out = int(a.frame[:3]), a.frame, Path(a.out)
    html = requests.get(f"{INDEX}/{t}/{frame}/interferograms/", timeout=180).text
    pairs = sorted(set(re.findall(r'href="(\d{8}_\d{8})/?"', html)))
    sel = summer_pairs(pairs, range(y0, y1 + 1), a.per_year)
    print(f"{frame}: {len(pairs)} pairs listed, {len(sel)} summer annual pairs selected", flush=True)
    base = f"{PUBLIC}/{t}/{frame}"
    jobs = [("meta", f"{base}/metadata/{frame}.geo.{k}.tif", out / "meta" / f"{k}.tif") for k in ("U", "E", "N", "hgt")]
    for p in sel:
        jobs.append(("unw", f"{base}/interferograms/{p}/{p}.geo.unw.tif", out / "los_mm" / f"{p}.tif"))
        jobs.append(("cc", f"{base}/interferograms/{p}/{p}.geo.cc.tif", out / "cc" / f"{p}.tif"))
    res = {}
    with ThreadPoolExecutor(6) as ex:
        futs = {ex.submit(fetch, k, u, o, tuple(a.bbox)): o.name + ":" + k for k, u, o in jobs}
        for f in as_completed(futs):
            res[futs[f]] = f.result()
    bad = {k: v for k, v in res.items() if v.startswith("error")}
    print(f"done: {len(res) - len(bad)} ok/cached, {len(bad)} errors {list(bad.items())[:5]}")
    (out / "manifest.json").write_text(json.dumps(dict(frame=frame, bbox=a.bbox, pairs=sel, errors=bad), indent=1))


if __name__ == "__main__":
    main()
