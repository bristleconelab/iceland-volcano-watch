"""Why do the two tracks disagree? Decompose ascending + descending into east and up, and fit each component alone.

For the windows with both tracks over (almost) the same dates (I3 Apr-Jul 2025, I4 Aug-Nov 2025), the chained LOS of
both tracks is block-medianed on the same 0.01 deg grid; at blocks with both tracks the 2 x 2 system
[E_a U_a; E_d U_d] [e; u] = [los_a - N_a n; los_d - N_d n] is solved with n = 0 (north motion neglected; its LOS
weight is |N| ~ 0.12). A source is then fitted to the up field alone and to the east field alone (offset + plane each,
fixed position and depth from the joint fit, grid over depth for each component). If the source shape were right,
dV(up) and dV(east) would agree. Their ratio measures the model's horizontal-to-vertical misfit.
    python -I scripts/track_hv_check.py -> docs/data/track_hv_check.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import svartsengi_gnss as G  # noqa: E402
import svartsengi_joint as J  # noqa: E402
from ivwlib import mogi, to_local_km  # noqa: E402
from okada import rect_source  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SRC = dict(lon=-22.456, lat=63.872)


def cell_index(df):
    """0.01-deg cell index; the two frames' pixel grids are offset by 1/3 pixel, so exact coordinates never match."""
    out = df.copy()
    out["ix"] = np.floor((out.lon + 22.75) / 0.01 + 1e-6).astype(int)
    out["iy"] = np.floor((out.lat - 63.80) / 0.01 + 1e-6).astype(int)
    return out


def decompose(wa, wd):
    m = cell_index(wa).merge(cell_index(wd).drop(columns=["lon", "lat"]), on=["ix", "iy"], suffixes=("_a", "_d"))
    m = m.rename(columns={"lon": "lon", "lat": "lat"})
    det = m.ue_a * m.uu_d - m.uu_a * m.ue_d
    e = (m.los_a * m.uu_d - m.uu_a * m.los_d) / det
    u = (m.ue_a * m.los_d - m.los_a * m.ue_d) / det
    return pd.DataFrame(dict(lon=m.lon, lat=m.lat, e=e, u=u)).dropna()


def fit_component(df, comp, kind, shape=None):
    x, y = to_local_km(df.lon.values, df.lat.values, SRC["lon"], SRC["lat"])
    best = None
    depths = np.arange(2.5, 7.01, 0.25)
    for z in depths:
        if kind == "mogi":
            ge, _, gz = mogi(x, y, z, 1e6)
            g = (ge if comp == "e" else gz) * 1e3
        else:
            L, W, st = shape
            ge, _, gz = rect_source(x, y, 0.0, 0.0, z, L, W, st, 0.0, opening=1e3 / (L * W))
            g = ge if comp == "e" else gz
        A = np.c_[g, np.ones_like(x), x, y]
        m, *_ = np.linalg.lstsq(A, df[comp].values, rcond=None)
        r = df[comp].values - A @ m
        c = float(r @ r)
        if best is None or c < best[0]:
            best = (c, z, m[0])
    return best[1], best[2]


def main():
    rows = []
    for w in J.CFG["windows"]:
        names = {f: J.chain(f, w["start"], w["end"]) for f in J.CFG["frames"]}
        if any(len(n) < 2 for n in names.values()):
            continue
        blk = {}
        for f, n in names.items():
            los, tr = J.chained_los(f, n)
            blk[f] = J.blocks(los, tr, f)
        fa, fd = J.CFG["frames"]
        df = decompose(blk[fa], blk[fd])
        for kind, shape in (("mogi", None), ("sill", (2.0, 2.0, 0.0)), ("sill", (3.0, 3.0, 0.0)), ("sill", (4.0, 4.0, 0.0)),
                            ("sill", (6.0, 2.0, 45.0))):
            zu, dvu = fit_component(df, "u", kind, shape)
            ze, dve = fit_component(df, "e", kind, shape)
            rows.append(dict(window=f"{w['interval']} {w['start']}", n_blocks=len(df),
                             asc=f"{names[fa][0][:8]}-{names[fa][-1][9:]}", desc=f"{names[fd][0][:8]}-{names[fd][-1][9:]}",
                             model=kind if shape is None else f"sill {shape[0]:g}x{shape[1]:g} N{shape[2]:g}E", depth_up_km=zu, dv_up=round(dvu, 2), depth_east_km=ze, dv_east=round(dve, 2),
                             ratio_east_over_up=round(dve / dvu, 3)))
            print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(ROOT / "docs" / "data" / "track_hv_check.csv", index=False)


if __name__ == "__main__":
    main()
