"""Is the west-south-west residual uplift a second source? (open item after the second review)

Every single-source model leaves an uplift lobe of 34-46 mm west-south-west of the source in the up field (peak
~215 mm) and a first ring harmonic towards WNW. Here a second point source is added to the 6 km line of point sources
(geometry fixed at each window's joint best fit, elongated_source.csv). The second source is searched on a grid
22.60-22.46 W, 63.835-63.885 N (0.01 x 0.005 deg) and 1.5-5 km depth; both volumes and the nuisance terms are linear
(GNSS + InSAR, full InSAR covariance, as elongated_source.py). Reported: delta-chi2 against the line alone, the second
source, its share of the total volume, and station volumes on 21 Jul 2026 with the window's volume share held fixed.
    python -I scripts/two_source.py -> docs/data/two_source.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import elongated_source as E  # noqa: E402
import joint_covariance as JC  # noqa: E402
import network_volume as N  # noqa: E402
import svartsengi_gnss as G  # noqa: E402
import svartsengi_joint as J  # noqa: E402
from ivwlib import linear_velocity, load_events, mogi, to_local_km  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs" / "data"


def second_columns(geo, obs, ins, lon, lat, z):
    n = len(obs)
    sx, sy = geo.src(lon, lat)
    gx, gy, gz = mogi(geo.xg - sx, geo.yg - sy, z, 1.0)
    cols = [np.r_[gx, gy, gz] * 1e9]
    for f, d in ins.items():
        x, y = geo.xi[f]
        ux, uy, uz = mogi(x - sx, y - sy, z, 1.0)
        cols.append((ux * d.ue.values + uy * d.un.values + uz * d.uu.values) * 1e9)
    return np.concatenate(cols)


def main():
    gn, _ = G.load_all(load_events())
    v, _, _ = linear_velocity(gn["SENG"], *G.PRE_UNREST)
    bg = {"SENG": float(v.u)}
    obs_imo = N.displacements(N.load_clean(), N.background(), "2026-07-21")
    es = pd.read_csv(DATA / "elongated_source.csv")
    rows = []
    for w in J.CFG["windows"]:
        wd = E.window_data(w, gn, bg)
        if wd is None:
            continue
        key = f"{w['interval']} {w['start']}"
        obs, ins, geo, Ls = wd["obs"], wd["ins"], wd["geo"], wd["Ls"]
        d = JC.data_vec(obs, ins)
        r6 = es[(es.window == key) & (es.L_km == 6.0)].iloc[0]
        A1 = E.design_line(geo, obs, ins, r6.lon, r6.lat, r6.depth_km, 6.0, r6.strike)
        _, c2_one, _ = JC.whitened_solve(A1, d, obs, ins, Ls)
        best = (np.inf, None)
        for lo in np.arange(-22.60, -22.4599, 0.01):
            for la in np.arange(63.835, 63.8851, 0.005):
                for z in np.arange(1.5, 5.01, 0.5):
                    A = np.c_[A1[:, :1], second_columns(geo, obs, ins, lo, la, z), A1[:, 1:]]
                    m, c2, _ = JC.whitened_solve(A, d, obs, ins, Ls)
                    if m[0] > 0 and m[1] > 0 and c2 < best[0]:
                        best = (c2, dict(lon=float(lo), lat=float(la), depth_km=float(z), dv_main=float(m[0]),
                                         dv_second=float(m[1])))
        c2, p = best
        share = p["dv_second"] / (p["dv_main"] + p["dv_second"])
        r = dict(window=key, chi2_line=round(c2_one, 1), chi2_two=round(c2, 1), dchi2=round(c2_one - c2, 1),
                 second_lon=round(p["lon"], 3), second_lat=round(p["lat"], 3), second_depth_km=p["depth_km"],
                 dv_main=round(p["dv_main"], 2), dv_second=round(p["dv_second"], 2), second_share=round(share, 3))
        main_p = dict(lon=r6.lon, lat=r6.lat, depth_km=r6.depth_km, L=6.0, strike=r6.strike)
        for s in ("SENG", "AUSV", "GEVK"):
            o = obs_imo[obs_imo.site == s].iloc[0]
            g_main = E.gz_at((o.lon, o.lat), main_p)
            x2, y2 = to_local_km(o.lon, o.lat, p["lon"], p["lat"])
            g_sec = float(mogi(np.atleast_1d(x2), np.atleast_1d(y2), p["depth_km"], 1e6)[2][0] * 1e3)
            r[f"{s.lower()}_alone_0721"] = round(float(o.du / ((1 - share) * g_main + share * g_sec)), 2)
        r["seng_over_ausv"] = round(r["seng_alone_0721"] / r["ausv_alone_0721"], 3)
        rows.append(r)
        print(r, flush=True)
    pd.DataFrame(rows).to_csv(DATA / "two_source.csv", index=False)


if __name__ == "__main__":
    main()
