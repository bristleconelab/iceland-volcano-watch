"""Which part of the source model controls the single-station volumes: position, depth or shape?

Second-round review (2026-11-16, M1): the station spread shrinks with the pressurised prolate cavity, but a point
source at the same position and depth may do the same. Here the volume from 22 Jul 2025 to 21 Jul 2026 (vertical
motion, as network_volume.py) is computed for SENG, AUSV and GEVK alone and for all five stations with:
IMO's point source (4 km); a point source at the joint-fit position at 3.25, 3.5, 4.0 and 4.25 km; a compact CDM
(0.75 x 0.75 km) and the elongated pressurised CDM of cdm_source.py (I3 geometry: 2.62 x 0.75 km, N60E) at 3.25 and
3.5 km. If the ratio SENG/AUSV follows depth and position and hardly changes with shape, the shape is not the lever.
    python -I scripts/station_volume_geometry.py        -> docs/data/station_volume_geometry.csv
    python -I scripts/station_volume_geometry.py radar  -> docs/data/station_volume_geometry_radar_chi2.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cdm_source as C  # noqa: E402
import network_volume as N  # noqa: E402
import svartsengi_gnss as G  # noqa: E402
from ivwlib import to_local_km  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
JOINT = dict(lon=-22.456, lat=63.872)


def volumes(obs, gz_fn, p):
    out = {}
    for s in ("SENG", "AUSV", "GEVK"):
        o = obs[obs.site == s].iloc[0]
        out[f"{s.lower()}_alone"] = round(float(o.du / gz_fn((o.lon, o.lat), p)), 2)
    x, y = to_local_km(obs.lon.values, obs.lat.values, p["lon"], p["lat"])
    if "ax" in p:
        gz = C.cdm_greens(np.asarray(x), np.asarray(y), p)[2]
    else:
        gz = N.mogi(np.asarray(x), np.asarray(y), p["depth_km"], 1e6)[2] * 1e3
    d = obs.du.to_numpy()
    sig = np.sqrt(G.SIG_U ** 2 + (G.SIG_REL * d) ** 2)
    m = float(np.sum(gz * d / sig ** 2) / np.sum(gz ** 2 / sig ** 2))
    out["all5_vertical"] = round(m, 2)
    out["all5_chi2_dof"] = round(float(np.sum(((gz * m - d) / sig) ** 2) / (len(d) - 1)), 2)
    out["seng_over_ausv"] = round(out["seng_alone"] / out["ausv_alone"], 3)
    return out


def point_gz(lonlat, p):
    x, y = to_local_km(*lonlat, p["lon"], p["lat"])
    return float(N.mogi(np.atleast_1d(x), np.atleast_1d(y), p["depth_km"], 1e6)[2][0] * 1e3)


def main():
    obs = N.displacements(N.load_clean(), N.background(), "2026-07-21")
    rows = [dict(model="IMO point source", lon=G.PARKS["lon"], lat=G.PARKS["lat"], depth_km=4.0,
                 **volumes(obs, point_gz, dict(G.PARKS)))]
    for z in (4.25, 4.0, 3.5, 3.25):
        p = dict(JOINT, depth_km=z)
        rows.append(dict(model="point, joint position", lon=p["lon"], lat=p["lat"], depth_km=z,
                         **volumes(obs, point_gz, p)))
    for label, ax_ in (("compact CDM 0.75 x 0.75 km", 0.75), ("pressurised prolate CDM 2.62 x 0.75 km N60E", 2.62)):
        for z in (3.5, 3.25):
            p = dict(JOINT, depth_km=z, ax=ax_, ay=0.75, strike=60.0)
            rows.append(dict(model=label, lon=p["lon"], lat=p["lat"], depth_km=z, **volumes(obs, C.gz_at, p)))
    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "docs" / "data" / "station_volume_geometry.csv", index=False)
    print(df.to_string(index=False))


if __name__ == "__main__" and len(sys.argv) == 1:
    main()


def radar_chi2():
    """Joint GNSS + InSAR chi2 (full InSAR covariance) at the joint position for a point source at several depths and
    for the pressurised prolate: does the radar accept the shallow depth only for an elongated body?"""
    import elongated_source as E
    import joint_covariance as JC
    import svartsengi_joint as J
    from ivwlib import linear_velocity, load_events
    gn, _ = G.load_all(load_events())
    v, _, _ = linear_velocity(gn["SENG"], *G.PRE_UNREST)
    bg = {"SENG": float(v.u)}
    rows = []
    for w in J.CFG["windows"]:
        if w["interval"] not in ("I2", "I3") and w["start"] != "2025-08-12":
            continue
        wd = E.window_data(w, gn, bg)
        obs, ins, geo, Ls = wd["obs"], wd["ins"], wd["geo"], wd["Ls"]
        d = JC.data_vec(obs, ins)
        best_pt = None
        for z in (3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5):
            m, c2, _ = JC.whitened_solve(J.design(geo, obs, ins, JOINT["lon"], JOINT["lat"], z), d, obs, ins, Ls)
            rows.append(dict(window=f"{w['interval']} {w['start']}", model="point", depth_km=z, chi2=round(c2, 1),
                             dv_mm3=round(float(m[0]), 2)))
        for z in (3.25, 3.5):
            p = dict(depth_km=z, ax=2.62, ay=0.75, strike=60.0)
            m, c2, _ = JC.whitened_solve(C.design(geo, obs, ins, JOINT["lon"], JOINT["lat"], p), d, obs, ins, Ls)
            rows.append(dict(window=f"{w['interval']} {w['start']}", model="prolate CDM 2.62 x 0.75 km N60E",
                             depth_km=z, chi2=round(c2, 1), dv_mm3=round(float(m[0]), 2)))
        print(pd.DataFrame(rows).tail(9).to_string(index=False), flush=True)
    pd.DataFrame(rows).to_csv(ROOT / "docs" / "data" / "station_volume_geometry_radar_chi2.csv", index=False)


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "radar":
    radar_chi2()
