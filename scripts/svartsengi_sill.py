"""Horizontal-sill alternative to the Mogi source in the joint GPS + InSAR fit (outline item 3 and the track mismatch).

Same data and linear nuisance terms as svartsengi_joint.py (GPS translation; per-frame InSAR offset and ramp). The
source is a horizontal rectangular tensile crack (Okada 1985, scripts/okada.py): centroid (lon, lat, depth), length L
along strike, width W, uniform opening. Its volume change is dV = L * W * opening. Grid: depth 2.5-6 km, L x W in
{2, 3, 4, 5, 6, 8} km squares, 2:1 rectangles up to 10 x 5 km with strikes 0, 45, 90, 135 deg and 3:1
rectangles (6 x 2, 9 x 3 km) with strikes 30, 45, 60 deg; the centroid is searched on a
3 x 3 grid around the best Mogi position (speed; the Mogi position search is the full grid).

For each window the best Mogi (from svartsengi_joint) and best sill are compared by chi2, and the sill is refitted with
each track alone to see whether it removes the ascending/descending dV disagreement found for the Mogi source.
    python -I scripts/svartsengi_sill.py -> docs/data/joint_sill.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import svartsengi_gnss as G  # noqa: E402
import svartsengi_joint as J  # noqa: E402
from ivwlib import load_events, linear_velocity  # noqa: E402
from okada import rect_source  # noqa: E402

DEPTHS = np.arange(2.5, 6.01, 0.25)
SHAPES = ([(a, a, 0.0) for a in (2, 3, 4, 5, 6, 8)] + [(2 * a, a, s) for a in (1.5, 2, 3, 4, 5) for s in (0, 45, 90, 135)]
          + [(3 * a, a, s) for a in (2, 3) for s in (30, 45, 60)])


def design_sill(geo, obs, ins, lon, lat, depth, L, W, strike):
    n = len(obs)
    sx, sy = geo.src(lon, lat)
    op = 1e3 / (L * W)  # mm of opening that gives 1 Mm3 (L, W in km)
    ge, gn_, gz = rect_source(geo.xg, geo.yg, sx, sy, depth, L, W, strike, 0.0, opening=op)
    frames = list(ins)
    A = [np.zeros((3 * n, 3 + 3 * len(frames)))]
    A[0][:, 0] = np.r_[ge, gn_, gz]
    A[0][:n, 1] = 1.0
    A[0][n:2 * n, 2] = 1.0
    for k, f in enumerate(frames):
        d = ins[f]
        x, y = geo.xi[f]
        ue, un, uz = rect_source(x, y, sx, sy, depth, L, W, strike, 0.0, opening=op)
        B = np.zeros((len(d), A[0].shape[1]))
        B[:, 0] = ue * d.ue.values + un * d.un.values + uz * d.uu.values
        B[:, 3 + 3 * k], B[:, 4 + 3 * k], B[:, 5 + 3 * k] = 1.0, x, y
        A.append(B)
    return np.vstack(A)


def solve_sill(geo, obs, ins, sig, *p):
    dg = G.data_vector(obs)
    d = np.r_[dg, *[ins[f].los.values for f in ins]]
    s = np.r_[G.sigmas(dg, len(obs)), *[np.full(len(ins[f]), sig[f]) for f in ins]]
    A = design_sill(geo, obs, ins, *p)
    m, *_ = np.linalg.lstsq(A / s[:, None], d / s, rcond=None)
    r = (A @ m - d) / s
    return m, float(np.sum(r ** 2))


def best_sill(obs, ins, sig, near):
    """Centroid searched on a 3 x 3 grid (0.008 deg lon, 0.004 deg lat steps) around the best Mogi position `near`."""
    geo = J.Geometry(obs, ins, G.PARKS["lon"], G.PARKS["lat"])
    best = (np.inf, None)
    lons = near["lon"] + np.array([-0.008, 0.0, 0.008])
    lats = near["lat"] + np.array([-0.004, 0.0, 0.004])
    for z in DEPTHS:
        for L, W, st in SHAPES:
            for lo in lons:
                for la in lats:
                    m, c2 = solve_sill(geo, obs, ins, sig, lo, la, z, L, W, st)
                    if m[0] > 0 and c2 < best[0]:
                        best = (c2, dict(lon=lo, lat=la, depth_km=z, L=L, W=W, strike=st, dv_mm3=m[0]))
    return best


def best_mogi(obs, ins, sig):
    geo = J.Geometry(obs, ins, G.PARKS["lon"], G.PARKS["lat"])
    best = (np.inf, None)
    for z in G.DEPTHS:
        for lo in G.LONS[::2]:
            for la in G.LATS[::2]:
                m, c2, *_ = J.solve(geo, obs, ins, sig, lo, la, z)
                if m[0] > 0 and c2 < best[0]:
                    best = (c2, dict(lon=lo, lat=la, depth_km=z, dv_mm3=m[0]))
    return best


def main():
    gn, _ = G.load_all(load_events())
    v, _, _ = linear_velocity(gn["SENG"], *G.PRE_UNREST)
    bg = {"SENG": float(v.u)}
    rows = []
    for w in J.CFG["windows"]:
        chains, ins = {}, {}
        for f in J.CFG["frames"]:
            names = J.chain(f, w["start"], w["end"])
            if len(names) >= 2:
                los, tr = J.chained_los(f, names)
                df = J.blocks(los, tr, f)
                if len(df) >= 50:
                    chains[f], ins[f] = names, df
        if not ins:
            continue
        t0 = min(J.pdate(c[0].split("_")[0]) for c in chains.values()).strftime("%Y-%m-%d")
        t1 = max(J.pdate(c[-1].split("_")[1]) for c in chains.values()).strftime("%Y-%m-%d")
        obs = G.interval_obs(gn, t0, t1, bg)
        sig = {f: J.insar_sigma(d, G.PARKS)[0] for f, d in ins.items()}
        ndata = 3 * len(obs) + sum(len(d) for d in ins.values())
        c_m, pm = best_mogi(obs, ins, sig)
        c_s, ps = best_sill(obs, ins, sig, pm)
        row = dict(window=f"{w['interval']} {w['start']}", n_data=ndata,
                   mogi_depth_km=pm["depth_km"], mogi_dv_mm3=round(pm["dv_mm3"], 2), mogi_chi2=round(c_m, 1),
                   mogi_chi2_dof=round(c_m / (ndata - 6 - 3 * len(ins)), 2),
                   sill_depth_km=ps["depth_km"], sill_L_km=ps["L"], sill_W_km=ps["W"], sill_strike=ps["strike"],
                   sill_lon=round(ps["lon"], 3), sill_lat=round(ps["lat"], 3), sill_dv_mm3=round(ps["dv_mm3"], 2),
                   sill_chi2=round(c_s, 1), sill_chi2_dof=round(c_s / (ndata - 9 - 3 * len(ins)), 2))
        if len(ins) == 2:
            for f in ins:
                _, mt = best_mogi(obs, {f: ins[f]}, {f: sig[f]})
                _, pt = best_sill(obs, {f: ins[f]}, {f: sig[f]}, mt)
                row[f"sill_dv_{f[:4]}"] = round(pt["dv_mm3"], 2)
                row[f"mogi_dv_{f[:4]}"] = round(mt["dv_mm3"], 2)
        rows.append(row)
        print(row, flush=True)
    pd.DataFrame(rows).to_csv(Path(__file__).resolve().parents[1] / "docs" / "data" / "joint_sill.csv", index=False)


if __name__ == "__main__":
    main()
