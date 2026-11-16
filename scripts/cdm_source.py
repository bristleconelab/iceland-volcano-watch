"""A pressurised prolate source (compound dislocation model, Nikkhoo et al. 2017) for the Svartsengi inflation.

The line of point sources in elongated_source.py is a stand-in. Here the source is a horizontal CDM (scripts/cdm.py)
with semi-axes ax (along strike) >= ay = az (prolate), which represents a uniformly pressurised horizontal prolate
cavity. The linear parameter is reported as the cavity volume change dV, the quantity IMO and Parks et al. (2025)
report for their point source: dV = DV * (dV/DV of the pressurised prolate cavity whose potency ratios match the CDM),
with dV/DV = 0.556 for a sphere (point-source convention) and 0.56-0.61 for prolate cavities (cdm.pressurised_ellipsoid).

Per radar window (same data, covariance and nuisance terms as joint_covariance.py / elongated_source.py):
1. Grid: depth 2.5-5.5 km, strike N30E-N90E (15 deg), centre on a 3 x 3 grid around the best point source.
   Mode "pressurised" (default, the result): the CDM's potency ratios are those of a uniformly pressurised prolate
   cavity of aspect A in {1, 1.5, 2, 3, 4, 6, 10, 50}, i.e. ax = k(A) ay with k from cdm.pressurised_ellipsoid
   (k <= 3.5), ay = az in {0.25, 0.5, 0.75, 1.0, 1.25} km. Mode "free" (kinematic comparison only): ax 0.5-3 km and
   ay = az 0.25-1 km independently; shapes with ax/ay > 3.5 are not pressurised cavities.
2. East/up test (track_hv_check.py) for the best CDM in the windows with both tracks.
3. Volume from 22 Jul 2025 to 21 Jul 2026 from SENG, AUSV, GEVK alone and from all five stations with each window's
   best geometry, against IMO's 28.3 Mm3.

    python -I scripts/cdm_source.py [pressurised|free] -> docs/data/cdm_source.csv | cdm_source_free.csv
Retrospective only; no forecasts.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cdm  # noqa: E402
import elongated_source as E  # noqa: E402
import joint_covariance as JC  # noqa: E402
import network_volume as N  # noqa: E402
import svartsengi_gnss as G  # noqa: E402
import svartsengi_joint as J  # noqa: E402
import track_hv_check as H  # noqa: E402
from ivwlib import linear_velocity, load_events, to_local_km  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs" / "data"
NU = 0.25
DEPTHS = np.arange(2.5, 5.51, 0.25)
STRIKES = np.arange(30, 91, 15)
AX = (0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
AY = (0.25, 0.5, 0.75, 1.0)
ASPECTS = (1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0, 50.0)
AY_P = (0.25, 0.5, 0.75, 1.0, 1.25)


def shapes(mode):
    """(ax, ay, cavity aspect or nan) triples searched in each mode."""
    if mode == "free":
        return [(ax, ay, np.nan) for ax in AX for ay in AY if ay <= ax]
    out = []
    for A in ASPECTS:
        _, _, pot = cdm.pressurised_ellipsoid(A, 1.0, 1.0, NU)
        k = (pot.sum() / pot[0] - 1) / 2
        out += [(round(k * ay, 4), ay, A) for ay in AY_P]
    return out
_RATIO_CACHE = {}


def cavity_ratio(ax, ay):
    """dV/DV of the pressurised prolate cavity whose potency fraction along its long axis equals the CDM's."""
    key = (ax, ay)
    if key in _RATIO_CACHE:
        return _RATIO_CACHE[key]
    f_cdm = ay * ay / (ay * ay + 2 * ax * ay)           # CDM potency fraction along x (uniform opening)
    best = None
    for asp in np.r_[1.0, np.linspace(1.01, 30, 600)]:
        dV, DV, pot = cdm.pressurised_ellipsoid(asp, 1.0, 1.0, NU)
        d = abs(pot[0] / pot.sum() - f_cdm)
        if best is None or d < best[0]:
            best = (d, dV / DV, asp)
    _RATIO_CACHE[key] = (best[1], best[2])
    return _RATIO_CACHE[key]


def cdm_greens(x, y, p):
    """Displacement (mm) per Mm3 of cavity volume change; x, y in km relative to the CDM centre."""
    ue, un, uv, DV = cdm.cdm(x, y, 0.0, 0.0, p["depth_km"], 0.0, 0.0, p["strike"] - 90.0, p["ax"], p["ay"], p["ay"],
                             1.0, NU)
    r, _ = cavity_ratio(p["ax"], p["ay"])
    scale = 1e3 / (DV * r)          # opening 1 km -> potency DV km3; 1 Mm3 = 1e-3 km3; output mm
    return ue * scale, un * scale, uv * scale


def design(geo, obs, ins, lon, lat, p):
    n = len(obs)
    sx, sy = geo.src(lon, lat)
    ge, gn, gz = cdm_greens(geo.xg - sx, geo.yg - sy, p)
    frames = list(ins)
    A = [np.zeros((3 * n, 3 + 3 * len(frames)))]
    A[0][:, 0] = np.r_[ge, gn, gz]
    A[0][:n, 1] = 1.0
    A[0][n:2 * n, 2] = 1.0
    for k, f in enumerate(frames):
        d = ins[f]
        x, y = geo.xi[f]
        ue, un, uz = cdm_greens(x - sx, y - sy, p)
        B = np.zeros((len(d), A[0].shape[1]))
        B[:, 0] = ue * d.ue.values + un * d.un.values + uz * d.uu.values
        B[:, 3 + 3 * k], B[:, 4 + 3 * k], B[:, 5 + 3 * k] = 1.0, x, y
        A.append(B)
    return np.vstack(A)


def gz_at(lonlat, p):
    x, y = to_local_km(*lonlat, p["lon"], p["lat"])
    return float(cdm_greens(np.atleast_1d(x), np.atleast_1d(y), p)[2][0])


def fit_network(obs, p, variant):
    x, y = to_local_km(obs.lon.values, obs.lat.values, p["lon"], p["lat"])
    gx, gy, gz = cdm_greens(np.asarray(x), np.asarray(y), p)
    n = len(obs)
    if variant == "vertical":
        d = obs.du.to_numpy()
        sig = np.sqrt(G.SIG_U ** 2 + (G.SIG_REL * d) ** 2)
        A = gz[:, None]
    else:
        d = np.r_[obs.de, obs.dn, obs.du]
        sig = G.sigmas(d, n)
        A = np.c_[np.r_[gx, gy, gz], np.r_[np.ones(n), np.zeros(2 * n)], np.r_[np.zeros(n), np.ones(n), np.zeros(n)]]
    m, *_ = np.linalg.lstsq(A / sig[:, None], d / sig, rcond=None)
    r = (A @ m - d) / sig
    return float(m[0]), float(r @ r / max(len(d) - A.shape[1], 1))


def hv_ratio(names, p):
    blk = {}
    for f, nm in names.items():
        los, tr = J.chained_los(f, nm)
        blk[f] = J.blocks(los, tr, f)
    fa, fd = J.CFG["frames"]
    df = H.decompose(blk[fa], blk[fd])
    x, y = to_local_km(df.lon.values, df.lat.values, p["lon"], p["lat"])
    out = {}
    for comp in ("e", "u"):
        best = None
        for z in np.arange(2.5, 7.01, 0.25):
            q = dict(p, depth_km=float(z))
            ge, _, gz = cdm_greens(np.asarray(x), np.asarray(y), q)
            A = np.c_[ge if comp == "e" else gz, np.ones_like(x), x, y]
            m, *_ = np.linalg.lstsq(A, df[comp].values, rcond=None)
            c = float(np.sum((df[comp].values - A @ m) ** 2))
            if best is None or c < best[0]:
                best = (c, z, m[0])
        out[comp] = best
    return out["e"][2] / out["u"][2]


def main(mode="pressurised"):
    gn, _ = G.load_all(load_events())
    v, _, _ = linear_velocity(gn["SENG"], *G.PRE_UNREST)
    bg = {"SENG": float(v.u)}
    obs_imo = N.displacements(N.load_clean(), N.background(), "2026-07-21")
    rows = []
    for w in J.CFG["windows"]:
        wd = E.window_data(w, gn, bg)
        if wd is None:
            continue
        obs, ins, geo, Ls, pm = wd["obs"], wd["ins"], wd["geo"], wd["Ls"], wd["pm"]
        d = JC.data_vec(obs, ins)
        lons = pm["lon"] + np.array([-0.008, 0.0, 0.008])
        lats = pm["lat"] + np.array([-0.004, 0.0, 0.004])
        best = (np.inf, None)
        for ax, ay, asp in shapes(mode):
            if True:
                for st in STRIKES:
                    for z in DEPTHS:
                        if z - ay <= 0.2:
                            continue
                        p = dict(ax=ax, ay=ay, strike=float(st), depth_km=float(z), aspect=asp)
                        for lo in lons:
                            for la in lats:
                                A = design(geo, obs, ins, lo, la, p)
                                m, c2, sdv = JC.whitened_solve(A, d, obs, ins, Ls)
                                if m[0] > 0 and c2 < best[0]:
                                    best = (c2, dict(p, lon=float(lo), lat=float(la), dv=float(m[0])))
        c2, p = best
        npar = 6 + 3 * len(ins) + 3          # + strike and two shape parameters
        A0 = J.design(geo, obs, ins, pm["lon"], pm["lat"], pm["depth_km"])
        _, c2_point, _ = JC.whitened_solve(A0, d, obs, ins, Ls)
        r = dict(window=f"{w['interval']} {w['start']}", mode=mode, cavity_aspect=p["aspect"], ax_km=p["ax"],
                 ay_az_km=p["ay"], strike=p["strike"],
                 depth_km=p["depth_km"], lon=round(p["lon"], 3), lat=round(p["lat"], 3), dv_mm3=round(p["dv"], 2),
                 cavity_ratio=round(cavity_ratio(p["ax"], p["ay"])[0], 4),
                 equiv_cavity_aspect=round(cavity_ratio(p["ax"], p["ay"])[1], 2),
                 chi2=round(c2, 1), chi2_dof=round(c2 / (len(d) - npar), 3), chi2_point=round(c2_point, 1),
                 dchi2_vs_point=round(c2_point - c2, 1), point_depth_km=pm["depth_km"])
        if len(ins) == 2:
            r["east_over_up"] = round(hv_ratio(wd["names"], p), 3)
        for s in ("SENG", "AUSV", "GEVK"):
            o = obs_imo[obs_imo.site == s].iloc[0]
            r[f"{s.lower()}_alone_0721"] = round(o.du / gz_at((o.lon, o.lat), p), 2)
        for var in ("3d", "vertical"):
            dv, cc = fit_network(obs_imo, p, var)
            r[f"all5_{var}_0721"], r[f"all5_{var}_chi2_dof"] = round(dv, 2), round(cc, 2)
        rows.append(r)
        print(r, flush=True)
    pd.DataFrame(rows).to_csv(DATA / ("cdm_source.csv" if mode == "pressurised" else "cdm_source_free.csv"),
                              index=False)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "pressurised")
