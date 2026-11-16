"""Robustness of the best line source (elongated_source.py) to a wider search.

elongated_source.py searches the line centre on a 3 x 3 grid (+-0.008 deg lon, +-0.004 deg lat) around the best point
source and strikes N0E-N165E. Here, for the three windows with both tracks, the centre is searched on a 7 x 7 grid of
+-0.03 deg lon and +-0.012 deg lat (about +-1.4 and +-1.3 km), strikes N30E-N105E (15 deg), depth 3-5 km, L in
{0, 3, 4, 6} km. If the best centre, strike and volumes do not move, the narrow search did not bias the result.
    python -I scripts/elongated_centre_check.py -> docs/data/elongated_centre_check.csv
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
from ivwlib import linear_velocity, load_events  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def main():
    gn, _ = G.load_all(load_events())
    v, _, _ = linear_velocity(gn["SENG"], *G.PRE_UNREST)
    bg = {"SENG": float(v.u)}
    obs_imo = N.displacements(N.load_clean(), N.background(), "2026-07-21")
    rows = []
    for w in J.CFG["windows"]:
        if w["interval"] not in ("I2", "I3") and w["start"] != "2025-08-12":
            continue
        wd = E.window_data(w, gn, bg)
        obs, ins, geo, Ls, pm = wd["obs"], wd["ins"], wd["geo"], wd["Ls"], wd["pm"]
        d = JC.data_vec(obs, ins)
        best = {}
        for L in (0.0, 3.0, 4.0, 6.0):
            for st in ((0,) if L == 0 else range(30, 106, 15)):
                for z in np.arange(3.0, 5.01, 0.25):
                    for lo in pm["lon"] + np.linspace(-0.03, 0.03, 7):
                        for la in pm["lat"] + np.linspace(-0.012, 0.012, 7):
                            m, c2, _ = JC.whitened_solve(E.design_line(geo, obs, ins, lo, la, z, L, st), d, obs, ins, Ls)
                            if m[0] > 0 and c2 < best.get(L, (np.inf,))[0]:
                                best[L] = (c2, dict(L=L, strike=st, depth_km=float(z), lon=float(lo), lat=float(la),
                                                    dv=float(m[0])))
        npar0 = 6 + 3 * len(ins)
        for L, (c2, p) in sorted(best.items()):
            r = dict(window=f"{w['interval']} {w['start']}", L_km=L, strike=p["strike"], depth_km=p["depth_km"],
                     lon=round(p["lon"], 4), lat=round(p["lat"], 4), dv_mm3=round(p["dv"], 2), chi2=round(c2, 1),
                     chi2_dof=round(c2 / (len(d) - npar0 - (2 if L else 0)), 3))
            for s in ("SENG", "AUSV", "GEVK"):
                o = obs_imo[obs_imo.site == s].iloc[0]
                r[f"{s.lower()}_alone_0721"] = round(o.du / E.gz_at((o.lon, o.lat), p), 2)
            r["all5_3d_0721"], r["all5_3d_chi2_dof"] = [round(x, 2) for x in E.fit_network(obs_imo, p, "3d")]
            rows.append(r)
            print(r, flush=True)
    pd.DataFrame(rows).to_csv(ROOT / "docs" / "data" / "elongated_centre_check.csv", index=False)


if __name__ == "__main__":
    main()
