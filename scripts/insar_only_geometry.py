"""Station test with a geometry chosen by the radar alone (second-round review, item 7b: circularity).

The joint fits that pick the source geometry include the GNSS stations whose volumes are then compared. Here the
geometry is fitted to the InSAR alone (same chains, blocks, covariance and per-track nuisance terms as
elongated_source.py, no GNSS rows) for a point source and for lines of point sources (L = 2-8 km, strike N0E-N165E),
on the same depth and centre grids. The single-station and five-station volumes from 22 Jul 2025 to 21 Jul 2026
(vertical, as network_volume.py) are then computed with each InSAR-only geometry.
    python -I scripts/insar_only_geometry.py -> docs/data/insar_only_geometry.csv
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
        wd = E.window_data(w, gn, bg)
        if wd is None:
            continue
        ins, Ls, pm = wd["ins"], wd["Ls"], wd["pm"]
        empty = wd["obs"].iloc[0:0]
        geo = J.Geometry(empty, ins, G.PARKS["lon"], G.PARKS["lat"])
        d = JC.data_vec(empty, ins)
        lons = pm["lon"] + np.linspace(-0.016, 0.016, 5)
        lats = pm["lat"] + np.linspace(-0.008, 0.008, 5)
        best = {}
        for L in (0.0, 2.0, 4.0, 6.0, 8.0):
            for st in ((0,) if L == 0 else E.STRIKES):
                for z in E.DEPTHS:
                    for lo in lons:
                        for la in lats:
                            A = np.delete(E.design_line(geo, empty, ins, lo, la, z, L, st), [1, 2], axis=1)
                            m, c2, _ = JC.whitened_solve(A, d, empty, ins, Ls)
                            if m[0] > 0 and c2 < best.get(L, (np.inf,))[0]:
                                best[L] = (c2, dict(L=L, strike=int(st), depth_km=float(z), lon=float(lo),
                                                    lat=float(la), dv=float(m[0])))
        for L, (c2, p) in sorted(best.items()):
            r = dict(window=f"{w['interval']} {w['start']}", L_km=L, strike=p["strike"], depth_km=p["depth_km"],
                     lon=round(p["lon"], 3), lat=round(p["lat"], 3), dv_mm3=round(p["dv"], 2), chi2_insar=round(c2, 1))
            for s in ("SENG", "AUSV", "GEVK"):
                o = obs_imo[obs_imo.site == s].iloc[0]
                r[f"{s.lower()}_alone_0721"] = round(float(o.du / E.gz_at((o.lon, o.lat), p)), 2)
            r["seng_over_ausv"] = round(r["seng_alone_0721"] / r["ausv_alone_0721"], 3)
            r["all5_vertical_0721"], r["all5_vertical_chi2_dof"] = [round(x, 2) for x in E.fit_network(obs_imo, p,
                                                                                                        "vertical")]
            rows.append(r)
            print(r, flush=True)
    pd.DataFrame(rows).to_csv(ROOT / "docs" / "data" / "insar_only_geometry.csv", index=False)


if __name__ == "__main__":
    main()
