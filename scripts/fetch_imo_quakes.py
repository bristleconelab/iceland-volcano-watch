"""Download the Icelandic Met Office (IMO) earthquake catalogue for a polygon through the public Quakes API.

API: https://api.vedur.is/quakes/ (OpenAPI spec at /quakes/openapi.json, licence CC BY 4.0, spec version 2026-08-06).
The /events endpoint returns CSV with columns event_id,time,latitude,longitude,depth,magnitude,evaluation_mode,type.
Two catalogue systems are served: 'sil' (legacy SIL system, reviewed "manual" events back to 1991) and 'seiscomp'
(primary system since Feb 2026). Requests are split into windows so that each response stays small; a window is split
again when /events/count says it holds more than --max-per-request events.

The responses are parsed as CSV text only; nothing from the server is executed.

Example:
    python -I scripts/fetch_imo_quakes.py --system sil --start 2020-01-01 --end 2026-11-12 \
        --bbox -22.9 63.75 -21.0 64.15 --out data/raw/imo/sil_reykjanes.csv
"""
from __future__ import annotations

import argparse
import datetime as dt
import io
import sys
import time
from pathlib import Path

import pandas as pd
import requests

API = "https://api.vedur.is/quakes"
UA = {"User-Agent": "iceland-volcano-watch/0.1 (research script)"}


def polygon_wkt(w: float, s: float, e: float, n: float) -> str:
    return f"POLYGON(({w} {s}, {w} {n}, {e} {n}, {e} {s}, {w} {s}))"


def get(path: str, params: dict, tries: int = 4) -> requests.Response:
    for k in range(tries):
        try:
            r = requests.get(f"{API}{path}", params=params, headers=UA, timeout=300)
            if r.status_code == 200:
                return r
            err = f"HTTP {r.status_code}: {r.text[:200]}"
        except requests.RequestException as ex:
            err = str(ex)
        time.sleep(5 * (k + 1))
    raise RuntimeError(f"{path} {params.get('start_time')}..{params.get('end_time')}: {err}")


def fetch_window(base: dict, t0: dt.datetime, t1: dt.datetime, max_n: int, out: list, log) -> None:
    p = dict(base, start_time=t0.isoformat(), end_time=t1.isoformat())
    n = int(get("/events/count", p).text.strip())
    if n == 0:
        return
    if n > max_n and (t1 - t0) > dt.timedelta(hours=6):
        mid = t0 + (t1 - t0) / 2
        fetch_window(base, t0, mid, max_n, out, log)
        fetch_window(base, mid, t1, max_n, out, log)
        return
    r = get("/events", dict(p, format="csv"))
    df = pd.read_csv(io.StringIO(r.text))
    if len(df) != n:
        log(f"  warning: count {n} but {len(df)} rows for {t0}..{t1}")
    out.append(df)
    log(f"  {t0:%Y-%m-%d %H:%M} .. {t1:%Y-%m-%d %H:%M}: {len(df)} events")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--system", choices=["sil", "seiscomp"], default="sil")
    ap.add_argument("--start", required=True)
    ap.add_argument("--end", required=True)
    ap.add_argument("--bbox", nargs=4, type=float, required=True, metavar=("W", "S", "E", "N"))
    ap.add_argument("--mode", choices=["manual", "automatic", "all"], default="manual")
    ap.add_argument("--window-days", type=int, default=31)
    ap.add_argument("--max-per-request", type=int, default=20000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    base = dict(polygon=polygon_wkt(*a.bbox), system=a.system, depth_min=0, depth_max=50, size_min=-3, size_max=10)
    if a.mode != "all":
        base["evaluation_mode"] = a.mode
    t, end = dt.datetime.fromisoformat(a.start), dt.datetime.fromisoformat(a.end)
    frames: list[pd.DataFrame] = []
    log = lambda s: print(s, flush=True)  # noqa: E731
    while t < end:
        t1 = min(t + dt.timedelta(days=a.window_days), end)
        fetch_window(base, t, t1, a.max_per_request, frames, log)
        t = t1
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if len(df):
        df = df.drop_duplicates("event_id").sort_values("time").reset_index(drop=True)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"wrote {len(df)} events to {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
