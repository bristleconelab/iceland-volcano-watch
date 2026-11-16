"""Realistic uncertainties for the joint GPS + InSAR source: spatially correlated InSAR noise (outline item 2/3).

The joint fit (svartsengi_joint.py) treats the ~1,000-1,900 InSAR blocks as independent, so its formal depth interval
collapses to one grid step. Here, for every window and frame:
1. Residuals of the best Mogi fit at blocks more than 12 km from the source (deformation-free) give an empirical
   covariance C(r) in 1-km distance bins; C(r) = s2 exp(-r / lam) is fitted for r > 0, and the nugget
   n2 = max(var - s2, 0.05 var) is added on the diagonal.
2. The fits are redone with the full covariance (Cholesky whitening), GPS errors unchanged.
3. Depth profile (lon, lat re-optimised on a 3 x 3 grid around the best position at each depth): the 1-sigma interval is
   chi2 <= min + max(chi2/dof, 1). dV is reported at the best depth and as the range over that interval.
4. The same for the best sill shape from svartsengi_sill (shape and strike fixed, depth profiled).
5. The SENG volume scale relative to the GPS-only source (1.29 km from SENG, 3.875 km) is converted into the
   21 Jul 2026 volume (29.1 Mm3 with the GPS-only source) for the ends of each interval.

    python -I scripts/joint_covariance.py -> docs/data/joint_covariance.csv, docs/data/insar_covariance.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parent))
import svartsengi_gnss as G  # noqa: E402
import svartsengi_joint as J  # noqa: E402
import svartsengi_sill as S  # noqa: E402
from ivwlib import load_events, linear_velocity, to_local_km  # noqa: E402
from okada import rect_source  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs" / "data"
SENG = (-22.43374, 63.88224)
V_GPS_ONLY_0721 = 29.1  # Mm3, docs/data/imo_volume_comparison.csv (this_work, 2026-07-21)


def exp_cov(r, s2, lam):
    return s2 * np.exp(-r / lam)


def empirical_covariance(x, y, v, rmax=20.0, dr=1.0):
    """Binned covariance of zero-mean v as a function of separation; returns bin centres, values, counts, variance."""
    v = v - v.mean()
    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    r = np.hypot(dx, dy)
    prod = v[:, None] * v[None, :]
    iu = np.triu_indices(len(v), k=1)
    r, prod = r[iu], prod[iu]
    edges = np.arange(0.0, rmax + dr, dr)
    k = np.digitize(r, edges) - 1
    ok = (k >= 0) & (k < len(edges) - 1)
    cnt = np.bincount(k[ok], minlength=len(edges) - 1)
    sm = np.bincount(k[ok], weights=prod[ok], minlength=len(edges) - 1)
    c = np.where(cnt > 0, sm / np.maximum(cnt, 1), np.nan)
    return (edges[:-1] + dr / 2), c, cnt, float(np.var(v))


def fit_covariance(x, y, v):
    rc, c, cnt, var = empirical_covariance(x, y, v)
    ok = np.isfinite(c) & (cnt >= 30)
    (s2, lam), _ = curve_fit(exp_cov, rc[ok], c[ok], p0=[max(c[ok][0], 1e-3), 3.0],
                             bounds=([0.0, 0.1], [10 * var, 100.0]))
    n2 = max(var - s2, 0.05 * var)
    return dict(s2=float(s2), lam_km=float(lam), nugget=float(n2), var=var, bins_km=rc[ok].round(1).tolist(),
                cov=np.round(c[ok], 3).tolist())


def whitener(ins, geo, covs):
    """Lower Cholesky factors of each frame's covariance matrix."""
    out = {}
    for f, d in ins.items():
        x, y = geo.xi[f]
        r = np.hypot(x[:, None] - x[None, :], y[:, None] - y[None, :])
        C = exp_cov(r, covs[f]["s2"], covs[f]["lam_km"]) + covs[f]["nugget"] * np.eye(len(d))
        out[f] = np.linalg.cholesky(C)
    return out


def whitened_solve(A, d, obs, ins, Ls):
    ng = 3 * len(obs)
    sg = G.sigmas(d[:ng], len(obs))
    Aw, dw = [A[:ng] / sg[:, None]], [d[:ng] / sg]
    i0 = ng
    for f, dd in ins.items():
        n = len(dd)
        Aw.append(solve_triangular(Ls[f], A[i0:i0 + n], lower=True))
        dw.append(solve_triangular(Ls[f], d[i0:i0 + n], lower=True))
        i0 += n
    Aw, dw = np.vstack(Aw), np.concatenate(dw)
    m, *_ = np.linalg.lstsq(Aw, dw, rcond=None)
    r = Aw @ m - dw
    cov_m = np.linalg.inv(Aw.T @ Aw)
    return m, float(r @ r), float(np.sqrt(cov_m[0, 0]))


def data_vec(obs, ins):
    return np.r_[G.data_vector(obs), *[ins[f].los.values for f in ins]]


def profile(kind, obs, ins, geo, Ls, best, shape=None):
    d = data_vec(obs, ins)
    lons = best["lon"] + np.array([-0.008, 0.0, 0.008])
    lats = best["lat"] + np.array([-0.004, 0.0, 0.004])
    depths = G.DEPTHS if kind == "mogi" else S.DEPTHS
    rows = []
    for z in depths:
        bz = (np.inf, None)
        for lo in lons:
            for la in lats:
                A = (J.design(geo, obs, ins, lo, la, z) if kind == "mogi"
                     else S.design_sill(geo, obs, ins, lo, la, z, *shape))
                m, c2, sdv = whitened_solve(A, d, obs, ins, Ls)
                if m[0] > 0 and c2 < bz[0]:
                    bz = (c2, dict(depth_km=float(z), lon=lo, lat=la, dv=float(m[0]), sdv=sdv, chi2=c2))
        if bz[1]:
            rows.append(bz[1])
    return pd.DataFrame(rows)


def seng_gz_shape(kind, p, shape=None):
    if kind == "mogi":
        return G.seng_gz(SENG, dict(lon=p["lon"], lat=p["lat"], depth_km=p["depth_km"]))
    x, y = to_local_km(*SENG, G.PARKS["lon"], G.PARKS["lat"])
    sx, sy = to_local_km(p["lon"], p["lat"], G.PARKS["lon"], G.PARKS["lat"])
    L, W, st = shape
    return float(rect_source(x, y, sx, sy, p["depth_km"], L, W, st, 0.0, opening=1e3 / (L * W))[2])


def main():
    gn, _ = G.load_all(load_events())
    v, _, _ = linear_velocity(gn["SENG"], *G.PRE_UNREST)
    bg = {"SENG": float(v.u)}
    g0 = G.seng_gz(SENG, dict(lon=-22.446, lat=63.872, depth_km=3.875))
    sill_csv = pd.read_csv(DATA / "joint_sill.csv").set_index("window")
    rows, covjson = [], {}
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
        key = f"{w['interval']} {w['start']}"
        t0 = min(J.pdate(c[0].split("_")[0]) for c in chains.values()).strftime("%Y-%m-%d")
        t1 = max(J.pdate(c[-1].split("_")[1]) for c in chains.values()).strftime("%Y-%m-%d")
        obs = G.interval_obs(gn, t0, t1, bg)
        sig = {f: J.insar_sigma(d, G.PARKS)[0] for f, d in ins.items()}
        geo = J.Geometry(obs, ins, G.PARKS["lon"], G.PARKS["lat"])
        _, pm = S.best_mogi(obs, ins, sig)
        # residual covariance from the diagonal-weight Mogi fit, far field only
        _, _, _, _, pred = J.solve(geo, obs, ins, sig, pm["lon"], pm["lat"], pm["depth_km"])
        covs, i0 = {}, 3 * len(obs)
        for f, d in ins.items():
            res = d.los.values - pred[i0:i0 + len(d)]
            i0 += len(d)
            x, y = geo.xi[f]
            sx, sy = geo.src(pm["lon"], pm["lat"])
            far = np.hypot(x - sx, y - sy) > J.FAR_KM
            covs[f] = fit_covariance(x[far], y[far], res[far])
        covjson[key] = covs
        Ls = whitener(ins, geo, covs)
        ndof = 3 * len(obs) + sum(len(d) for d in ins.values())
        srow = sill_csv.loc[key]
        shape = (float(srow.sill_L_km), float(srow.sill_W_km), float(srow.sill_strike))
        for kind, best, shp, npar in (("mogi", pm, None, 6 + 3 * len(ins)),
                                      ("sill", dict(lon=float(srow.sill_lon), lat=float(srow.sill_lat)), shape,
                                       6 + 3 * len(ins))):
            prof = profile(kind, obs, ins, geo, Ls, best, shp)
            b = prof.loc[prof.chi2.idxmin()]
            scale = max(b.chi2 / (ndof - npar), 1.0)
            inside = prof[prof.chi2 <= b.chi2 + scale]
            vols = [V_GPS_ONLY_0721 * g0 / seng_gz_shape(kind, r, shp) for _, r in inside.iterrows()]
            rows.append(dict(window=key, model=kind, shape="" if shp is None else f"{shp[0]:g}x{shp[1]:g} km N{shp[2]:g}E",
                             depth_km=b.depth_km, depth_lo_km=inside.depth_km.min(), depth_hi_km=inside.depth_km.max(),
                             dv_mm3=round(b.dv, 2), dv_formal_1sigma=round(b.sdv * np.sqrt(scale), 2),
                             dv_lo_over_depth=round(inside.dv.min(), 2), dv_hi_over_depth=round(inside.dv.max(), 2),
                             chi2_dof_fullcov=round(b.chi2 / (ndof - npar), 2),
                             seng_vol_0721_lo=round(min(vols), 1), seng_vol_0721_hi=round(max(vols), 1)))
            print(rows[-1], flush=True)
        print(key, {f: (round(c["lam_km"], 2), round(np.sqrt(c["s2"]), 2), round(np.sqrt(c["nugget"]), 2))
                    for f, c in covs.items()}, flush=True)
    pd.DataFrame(rows).to_csv(DATA / "joint_covariance.csv", index=False)
    (DATA / "insar_covariance.json").write_text(json.dumps(covjson, indent=1))


if __name__ == "__main__":
    main()
