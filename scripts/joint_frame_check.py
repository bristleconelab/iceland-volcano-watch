"""Robustness of the joint fit: refit every window with each LiCSAR frame alone (plus GNSS), and with GNSS weights x10.

The spread of depth and dV across these variants is a more honest uncertainty than the formal interval, which treats
the spatially correlated InSAR blocks as independent.  python -I scripts/joint_frame_check.py -> docs/data/joint_frame_check.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import svartsengi_gnss as G  # noqa: E402
import svartsengi_joint as J  # noqa: E402
from ivwlib import load_events, linear_velocity  # noqa: E402


def best_fit(obs, ins, sig):
    geo = J.Geometry(obs, ins, G.PARKS["lon"], G.PARKS["lat"])
    best = (np.inf, None)
    for z in G.DEPTHS:
        for lo in G.LONS[::2]:
            for la in G.LATS[::2]:
                m, c2, *_ = J.solve(geo, obs, ins, sig, lo, la, z)
                if m[0] > 0 and c2 < best[0]:
                    best = (c2, (lo, la, z, m[0]))
    return best[1]


def main():
    gn, _ = G.load_all(load_events())
    v, _, _ = linear_velocity(gn["SENG"], *G.PRE_UNREST)
    bg = {"SENG": float(v.u)}
    rows = []
    for w in J.CFG["windows"]:
        chains, ins = {}, {}
        for f in J.CFG["frames"]:
            names = J.chain(f, w["start"], w["end"])
            if len(names) < 2:
                continue
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
        variants = {"both frames": ins}
        if len(ins) == 2:
            variants.update({f"{f[:4]} only": {f: ins[f]} for f in ins})
        for name, sub in variants.items():
            for gw in (1.0, 0.1):  # 0.1 = GNSS sigmas divided by 10 (GNSS weight x100 in chi2)
                old = (G.SIG_H, G.SIG_U)
                G.SIG_H, G.SIG_U = old[0] * gw, old[1] * gw
                try:
                    lo, la, z, dv = best_fit(obs, sub, {f: sig[f] for f in sub})
                finally:
                    G.SIG_H, G.SIG_U = old
                rows.append(dict(window=f"{w['interval']} {w['start']}", variant=name,
                                 gnss_weight="x1" if gw == 1.0 else "x100", lon=round(lo, 3), lat=round(la, 3),
                                 depth_km=z, dv_mm3=round(dv, 2)))
                print(rows[-1])
    pd.DataFrame(rows).to_csv(Path(__file__).resolve().parents[1] / "docs" / "data" / "joint_frame_check.csv",
                              index=False)


if __name__ == "__main__":
    main()
