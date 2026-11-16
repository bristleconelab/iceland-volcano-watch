"""East/up shape test done consistently: the model goes through the same decomposition as the data.

track_hv_check.py decomposes ascending + descending LOS into east and up assuming zero north motion, and then fits a
source's pure east and up components with the depth re-fitted for each. An independent review (2026-11-16) noted two
biases: north motion is neglected in the data but not in the model, and refitting depth per component absorbs part of
the very mismatch being tested. Here, for each window with both tracks:
1. At the blocks seen by both tracks, the model's (E, N, U) displacement per Mm3 is projected onto each track's LOS
   with that track's unit vectors and decomposed with the same 2 x 2 system (north set to zero) as the data.
2. With the geometry fixed at the joint GNSS + InSAR best fit of each source, the volume is fitted to the decomposed
   east field alone and to the up field alone (each with offset + plane). A correct shape gives east/up = 1.
Sources: point (joint fit), line of point sources L = 4 and 6 km (elongated_source.csv), 6 x 2 km sill
(joint_sill.csv), pressurised prolate CDM (cdm_source.csv). Bootstrap over blocks (200 resamples) gives a 1-sigma
spread of the ratio (blocks are spatially correlated, so this is a lower bound).

    python -I scripts/east_up_consistent.py -> docs/data/east_up_consistent.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cdm_source as C  # noqa: E402
import elongated_source as E  # noqa: E402
import svartsengi_gnss as G  # noqa: E402
import svartsengi_joint as J  # noqa: E402
import track_hv_check as H  # noqa: E402
from ivwlib import to_local_km  # noqa: E402
from okada import rect_source  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs" / "data"
RNG = np.random.default_rng(7)


def merged(names):
    blk = {}
    for f, nm in names.items():
        los, tr = J.chained_los(f, nm)
        blk[f] = H.cell_index(J.blocks(los, tr, f))
    fa, fd = J.CFG["frames"]
    m = blk[fa].merge(blk[fd].drop(columns=["lon", "lat"]), on=["ix", "iy"], suffixes=("_a", "_d"))
    return m


def decompose(m, la, ld):
    det = m.ue_a * m.uu_d - m.uu_a * m.ue_d
    e = (la * m.uu_d - m.uu_a * ld) / det
    u = (m.ue_a * ld - la * m.ue_d) / det
    return np.asarray(e), np.asarray(u)


def model_enu(kind, p, x, y):
    if kind in ("point", "line"):
        return E.line_greens(x, y, p["depth_km"], p.get("L", 0.0), p.get("strike", 0))
    if kind == "sill":
        L, W, st = p["shape"]
        return rect_source(x, y, 0.0, 0.0, p["depth_km"], L, W, st, 0.0, opening=1e3 / (L * W))
    return C.cdm_greens(x, y, p)


def ratio(m, de, du, ge, gu, x, y, idx=None):
    if idx is None:
        idx = np.arange(len(x))
    out = []
    for g, dobs in ((ge, de), (gu, du)):
        A = np.c_[g[idx], np.ones(idx.size), x[idx], y[idx]]
        c, *_ = np.linalg.lstsq(A, dobs[idx], rcond=None)
        out.append(c[0])
    return out[0] / out[1], out[0], out[1]


def sources(window):
    es = pd.read_csv(DATA / "elongated_source.csv").query("window == @window")
    sl = pd.read_csv(DATA / "joint_sill.csv").set_index("window").loc[window]
    cd = pd.read_csv(DATA / "cdm_source.csv").set_index("window").loc[window]
    out = []
    for L, kind in ((0.0, "point"), (4.0, "line"), (6.0, "line")):
        r = es[es.L_km == L].iloc[0]
        out.append((f"{kind} L={L:g}", kind, dict(lon=r.lon, lat=r.lat, depth_km=r.depth_km, L=L, strike=r.strike)))
    out.append(("sill 6x2", "sill", dict(lon=sl.sill_lon, lat=sl.sill_lat, depth_km=sl.sill_depth_km,
                                         shape=(sl.sill_L_km, sl.sill_W_km, sl.sill_strike))))
    out.append((f"prolate CDM A={cd.cavity_aspect:g}", "cdm", dict(lon=cd.lon, lat=cd.lat, depth_km=cd.depth_km,
                                                                    ax=cd.ax_km, ay=cd.ay_az_km, strike=cd.strike)))
    return out


def main():
    rows = []
    for w in J.CFG["windows"]:
        names = {f: J.chain(f, w["start"], w["end"]) for f in J.CFG["frames"]}
        if any(len(n) < 2 for n in names.values()):
            continue
        key = f"{w['interval']} {w['start']}"
        m = merged(names)
        de, du = decompose(m, m.los_a.values, m.los_d.values)
        ok = np.isfinite(de) & np.isfinite(du)
        m, de, du = m[ok].reset_index(drop=True), de[ok], du[ok]
        for label, kind, p in sources(key):
            x, y = to_local_km(m.lon.values, m.lat.values, p["lon"], p["lat"])
            x, y = np.asarray(x), np.asarray(y)
            ge3, gn3, gu3 = model_enu(kind, p, x, y)
            la = ge3 * m.ue_a + gn3 * m.un_a + gu3 * m.uu_a
            ld = ge3 * m.ue_d + gn3 * m.un_d + gu3 * m.uu_d
            ge, gu = decompose(m, la.values, ld.values)
            r0, ve, vu = ratio(m, de, du, ge, gu, x, y)
            rb = [ratio(m, de, du, ge, gu, x, y, RNG.integers(0, len(x), len(x)))[0] for _ in range(200)]
            r_naive = ratio(m, de, du, ge3, gu3, x, y)[0]     # pure model E/U, as in track_hv_check (fixed depth)
            rows.append(dict(window=key, source=label, n_blocks=len(x), dv_east=round(ve, 2), dv_up=round(vu, 2),
                             east_over_up=round(r0, 3), boot_sigma=round(float(np.std(rb)), 3),
                             east_over_up_naive=round(r_naive, 3)))
            print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(DATA / "east_up_consistent.csv", index=False)


if __name__ == "__main__":
    main()
