"""Does a horizontally elongated source remove the station dependence of the Svartsengi volume?

A point source fits neither the five GNSS stations (network_volume.py, chi2/dof 14-22) nor makes SENG agree with
the rest: SENG, on the north-east axis of the uplift, gives 1.23 times IMO's volume and the multi-station fits
0.69-0.99 times. Parks et al. (2025) note the inflation is elongated NE-SW. Here the source is a horizontal line of
N = 9 equal centres of dilatation (Mogi sources) of total volume change dV, length L and strike s, at depth z: a
simple stand-in for a horizontally elongated pressurised body (L = 0 is the point source). It is not the exact
Yang et al. (1988) spheroid; it keeps each element compact, so it does not inherit a flat crack's low horizontal
motion.

Per radar window (same chains, GNSS dates, correlated InSAR covariance and nuisance terms as joint_covariance.py):
1. Grid: depth 2-6 km (0.25), L in {0, 2, 3, 4, 6, 8} km, strike N0E-N165E (15 deg), centre on a 5 x 5 grid
   (+-0.016 deg lon, +-0.008 deg lat, about +-0.8 / 0.9 km) around the best point source; dV and nuisances linear. chi2 with the full InSAR covariance.
2. Station consistency: with each window's best geometry, the volume implied by SENG's vertical alone over the window
   against the joint dV; and on 21 Jul 2026 (from the 22 Jul 2025 post-eruption level, as network_volume.py) SENG
   alone, AUSV alone, GEVK alone and all five stations, against IMO's 28.3 Mm3.
3. East/up test (as track_hv_check.py) for the best line in the windows with both tracks.

    python -I scripts/elongated_source.py -> docs/data/elongated_source.csv, docs/data/elongated_source_imo.csv
Retrospective only; no forecasts.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import joint_covariance as JC  # noqa: E402
import network_volume as N  # noqa: E402
import svartsengi_gnss as G  # noqa: E402
import svartsengi_joint as J  # noqa: E402
import svartsengi_sill as S  # noqa: E402
import track_hv_check as H  # noqa: E402
from ivwlib import linear_velocity, load_events, mogi, to_local_km  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs" / "data"
NPTS = 9
DEPTHS = np.arange(2.0, 6.01, 0.25)
LENGTHS = (0.0, 2.0, 3.0, 4.0, 6.0, 8.0)
STRIKES = np.arange(0, 180, 15)
SENG = JC.SENG


def line_greens(x, y, depth, L, strike):
    """Displacement (mm) per Mm3 of a horizontal line of NPTS Mogi sources centred at (0, 0); x, y in km."""
    if L == 0:
        ge, gn, gz = mogi(x, y, depth, 1.0)
        return ge * 1e9, gn * 1e9, gz * 1e9
    t = np.linspace(-L / 2, L / 2, NPTS)
    a = np.radians(strike)
    ge = gn = gz = 0.0
    for ti in t:
        e, n, z = mogi(x - ti * np.sin(a), y - ti * np.cos(a), depth, 1.0 / NPTS)
        ge, gn, gz = ge + e, gn + n, gz + z
    return ge * 1e9, gn * 1e9, gz * 1e9


def design_line(geo, obs, ins, lon, lat, depth, L, strike):
    n = len(obs)
    sx, sy = geo.src(lon, lat)
    ge, gn, gz = line_greens(geo.xg - sx, geo.yg - sy, depth, L, strike)
    frames = list(ins)
    A = [np.zeros((3 * n, 3 + 3 * len(frames)))]
    A[0][:, 0] = np.r_[ge, gn, gz]
    A[0][:n, 1] = 1.0
    A[0][n:2 * n, 2] = 1.0
    for k, f in enumerate(frames):
        d = ins[f]
        x, y = geo.xi[f]
        ue, un, uz = line_greens(x - sx, y - sy, depth, L, strike)
        B = np.zeros((len(d), A[0].shape[1]))
        B[:, 0] = ue * d.ue.values + un * d.un.values + uz * d.uu.values
        B[:, 3 + 3 * k], B[:, 4 + 3 * k], B[:, 5 + 3 * k] = 1.0, x, y
        A.append(B)
    return np.vstack(A)


def window_data(w, gn, bg):
    chains, ins = {}, {}
    for f in J.CFG["frames"]:
        names = J.chain(f, w["start"], w["end"])
        if len(names) >= 2:
            los, tr = J.chained_los(f, names)
            df = J.blocks(los, tr, f)
            if len(df) >= 50:
                chains[f], ins[f] = names, df
    if not ins:
        return None
    t0 = min(J.pdate(c[0].split("_")[0]) for c in chains.values()).strftime("%Y-%m-%d")
    t1 = max(J.pdate(c[-1].split("_")[1]) for c in chains.values()).strftime("%Y-%m-%d")
    obs = G.interval_obs(gn, t0, t1, bg)
    sig = {f: J.insar_sigma(d, G.PARKS)[0] for f, d in ins.items()}
    geo = J.Geometry(obs, ins, G.PARKS["lon"], G.PARKS["lat"])
    _, pm = S.best_mogi(obs, ins, sig)
    _, _, _, _, pred = J.solve(geo, obs, ins, sig, pm["lon"], pm["lat"], pm["depth_km"])
    covs, i0 = {}, 3 * len(obs)
    for f, d in ins.items():
        res = d.los.values - pred[i0:i0 + len(d)]
        i0 += len(d)
        x, y = geo.xi[f]
        sx, sy = geo.src(pm["lon"], pm["lat"])
        far = np.hypot(x - sx, y - sy) > J.FAR_KM
        covs[f] = JC.fit_covariance(x[far], y[far], res[far])
    return dict(names=chains, ins=ins, obs=obs, geo=geo, pm=pm, Ls=JC.whitener(ins, geo, covs), t0=t0, t1=t1)


def grid(wd):
    obs, ins, geo, Ls, pm = wd["obs"], wd["ins"], wd["geo"], wd["Ls"], wd["pm"]
    d = JC.data_vec(obs, ins)
    lons = pm["lon"] + np.linspace(-0.016, 0.016, 5)
    lats = pm["lat"] + np.linspace(-0.008, 0.008, 5)
    best = {}
    for L in LENGTHS:
        for st in (STRIKES if L > 0 else [0]):
            for z in DEPTHS:
                for lo in lons:
                    for la in lats:
                        A = design_line(geo, obs, ins, lo, la, z, L, st)
                        m, c2, sdv = JC.whitened_solve(A, d, obs, ins, Ls)
                        if m[0] > 0 and c2 < best.get(L, (np.inf,))[0]:
                            best[L] = (c2, dict(L=L, strike=int(st), depth_km=float(z), lon=float(lo), lat=float(la),
                                                dv=float(m[0]), sdv=sdv))
    return best


def gz_at(lonlat, p):
    x, y = to_local_km(*lonlat, p["lon"], p["lat"])
    return float(line_greens(np.atleast_1d(x), np.atleast_1d(y), p["depth_km"], p["L"], p["strike"])[2][0])


def fit_network(obs, p, variant):
    x, y = to_local_km(obs.lon.values, obs.lat.values, p["lon"], p["lat"])
    gx, gy, gz = line_greens(x, y, p["depth_km"], p["L"], p["strike"])
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
            ge, _, gz = line_greens(x, y, z, p["L"], p["strike"])
            A = np.c_[ge if comp == "e" else gz, np.ones_like(x), x, y]
            m, *_ = np.linalg.lstsq(A, df[comp].values, rcond=None)
            c = float(np.sum((df[comp].values - A @ m) ** 2))
            if best is None or c < best[0]:
                best = (c, z, m[0])
        out[comp] = best
    return out["e"][2] / out["u"][2], len(df)


def main():
    gn, _ = G.load_all(load_events())
    v, _, _ = linear_velocity(gn["SENG"], *G.PRE_UNREST)
    bg = {"SENG": float(v.u)}
    gclean, nbg = N.load_clean(), N.background()
    imo = pd.read_csv(DATA / "imo_volume_comparison.csv").iloc[0]
    obs_imo = N.displacements(gclean, nbg, imo.date)
    rows, imo_rows = [], []
    for w in J.CFG["windows"]:
        wd = window_data(w, gn, bg)
        if wd is None:
            continue
        key = f"{w['interval']} {w['start']}"
        ndof = 3 * len(wd["obs"]) + sum(len(d) for d in wd["ins"].values())
        best = grid(wd)
        seng_row = wd["obs"][wd["obs"].site == "SENG"].iloc[0]
        for L, (c2, p) in sorted(best.items()):
            npar = 6 + 3 * len(wd["ins"]) + (2 if L > 0 else 0)
            r = dict(window=key, gnss=f"{wd['t0']}..{wd['t1']}", L_km=L, strike=p["strike"], depth_km=p["depth_km"],
                     lon=round(p["lon"], 3), lat=round(p["lat"], 3), dv_mm3=round(p["dv"], 2), chi2=round(c2, 1),
                     chi2_dof=round(c2 / (ndof - npar), 3),
                     seng_implied_mm3=round(seng_row.du / gz_at(SENG, p), 2))
            r["seng_over_joint"] = round(r["seng_implied_mm3"] / p["dv"], 3)
            if len(wd["ins"]) == 2:
                r["east_over_up"], r["hv_blocks"] = hv_ratio(wd["names"], p)
                r["east_over_up"] = round(r["east_over_up"], 3)
            rows.append(r)
            ir = dict(window=key, L_km=L, strike=p["strike"], depth_km=p["depth_km"], imo_mm3=imo.imo_mm3)
            for s in ("SENG", "AUSV", "GEVK"):
                o = obs_imo[obs_imo.site == s].iloc[0]
                ir[f"{s.lower()}_alone"] = round(o.du / gz_at((o.lon, o.lat), p), 2)
            for var in ("3d", "vertical"):
                dv, cc = fit_network(obs_imo, p, var)
                ir[f"all5_{var}"], ir[f"all5_{var}_chi2_dof"] = round(dv, 2), round(cc, 2)
            imo_rows.append(ir)
            print(r, "|", {k: ir[k] for k in ("seng_alone", "ausv_alone", "gevk_alone", "all5_3d", "all5_3d_chi2_dof")},
                  flush=True)
    pd.DataFrame(rows).to_csv(DATA / "elongated_source.csv", index=False)
    pd.DataFrame(imo_rows).to_csv(DATA / "elongated_source_imo.csv", index=False)


if __name__ == "__main__":
    main()
