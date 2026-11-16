"""Joint GNSS + InSAR fit of the Svartsengi inflation source (outline item 2).

For every inflation window in config/svartsengi.json -> insar_short_pairs and every LiCSAR frame:
1. Chain consecutive 6/12-day unwrapped interferograms from the earliest to the latest reachable date (fewest pairs for
   the longest span). Summing a chain cancels the atmosphere of the intermediate epochs; only the two end epochs remain.
2. Block-median the chained LOS displacement on 10 x 10 pixels (0.01 deg, about 0.49 km E x 1.1 km N); keep blocks
   with >= 50 % valid pixels.
3. Take the GNSS displacements over exactly the chain dates (same cleaning and SENG background as svartsengi_gnss.py).
4. Grid-search a Mogi point source (lon, lat, depth); for each position solve linearly for dV, a common GNSS
   translation and, per frame, an InSAR offset and plane (x, y ramp). InSAR noise sigma per frame: robust std (1.4826
   MAD) of the plane-removed data more than 12 km from Parks' source. GNSS errors as in svartsengi_gnss.py.

Reported per window: source, chi2/dof for GNSS alone and for InSAR alone at the joint optimum, the GNSS-only optimum on
the same dates for comparison, and the number of InSAR blocks. LOS convention (LiCSBAS): positive = towards the
satellite; E, N, U files are the ground-to-satellite unit vector.

    python -I scripts/svartsengi_joint.py
Outputs: docs/data/joint_fits.csv, docs/data/joint_chains.json, docs/img/joint_fit_<window>.png
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

sys.path.insert(0, str(Path(__file__).resolve().parent))
import svartsengi_gnss as G  # noqa: E402
from ivwlib import load_events, linear_velocity, mogi, to_local_km  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT / "config" / "svartsengi.json").read_text())["insar_short_pairs"]
RAW = ROOT / "data" / "raw" / "licsar_short"
DATA, IMG = ROOT / "docs" / "data", ROOT / "docs" / "img"
BLOCK = 10
FAR_KM = 12.0


def pdate(s):
    return pd.Timestamp(f"{s[:4]}-{s[4:6]}-{s[6:8]}")


def chain(frame, start, end):
    """Longest-span path of consecutive pairs inside [start, end]; ties -> fewer pairs. Returns list of pair names."""
    pairs = []
    for f in sorted((RAW / frame / "los_mm").glob("*.tif")):
        a, b = f.stem.split("_")
        ta, tb = pdate(a), pdate(b)
        if pd.Timestamp(start) <= ta and tb <= pd.Timestamp(end):
            pairs.append((ta, tb, f.stem))
    if not pairs:
        return []
    dates = sorted({p[0] for p in pairs} | {p[1] for p in pairs})
    out = {d: [] for d in dates}
    for ta, tb, name in pairs:
        out[ta].append((tb, name))
    best = {}  # date -> (span_end, -npairs, path) for the best path starting at date
    for d in reversed(dates):
        cand = [(d, 0, [])]
        for tb, name in out[d]:
            e, negn, path = best[tb]
            cand.append((e, negn - 1, [name] + path))
        best[d] = max(cand, key=lambda c: ((c[0] - d).days, c[1]))
    start_d = max(dates, key=lambda d: ((best[d][0] - d).days, best[d][1]))
    return best[start_d][2]


def read(path):
    with rasterio.open(path) as r:
        return r.read(1).astype(float), r.transform


def chained_los(frame, names):
    tot = None
    for n in names:
        a, tr = read(RAW / frame / "los_mm" / f"{n}.tif")
        tot = a if tot is None else tot + a
    return tot, tr


def blocks(arr, tr, frame):
    h, w = arr.shape
    hb, wb = h // BLOCK, w // BLOCK
    a = arr[:hb * BLOCK, :wb * BLOCK].reshape(hb, BLOCK, wb, BLOCK)
    valid = np.isfinite(a).mean(axis=(1, 3))
    med = np.nanmedian(np.where(np.isfinite(a), a, np.nan).transpose(0, 2, 1, 3).reshape(hb, wb, -1), axis=2)
    rows, cols = np.mgrid[0:hb, 0:wb]
    xs, ys = rasterio.transform.xy(tr, rows.ravel() * BLOCK + BLOCK / 2 - 0.5, cols.ravel() * BLOCK + BLOCK / 2 - 0.5)
    lon, lat = np.asarray(xs), np.asarray(ys)
    unit = {}
    for k in "ENU":
        u, _ = read(RAW / frame / "meta" / f"{k}.tif")
        unit[k] = np.nanmedian(u[:hb * BLOCK, :wb * BLOCK].reshape(hb, BLOCK, wb, BLOCK).transpose(0, 2, 1, 3)
                               .reshape(hb, wb, -1), axis=2).ravel()
    keep = (valid.ravel() >= 0.5) & np.isfinite(med.ravel()) & np.isfinite(unit["U"])
    return pd.DataFrame(dict(lon=lon[keep], lat=lat[keep], los=med.ravel()[keep],
                             ue=unit["E"][keep], un=unit["N"][keep], uu=unit["U"][keep]))


def insar_sigma(df, src):
    x, y = to_local_km(df.lon.values, df.lat.values, src["lon"], src["lat"])
    far = np.hypot(x, y) > FAR_KM
    A = np.c_[np.ones(far.sum()), x[far], y[far]]
    m, *_ = np.linalg.lstsq(A, df.los.values[far], rcond=None)
    r = df.los.values[far] - A @ m
    return 1.4826 * float(np.median(np.abs(r - np.median(r)))), int(far.sum())


class Geometry:
    """Local km coordinates of all GNSS sites and InSAR blocks about one fixed origin, computed once (pyproj geodesics).
    A trial source at (lon, lat) is placed at its own local coordinates; within ~30 km the shift error of this
    approximation is < 0.1 % of the distances (checked in tests)."""

    def __init__(self, obs, ins, lon0, lat0):
        self.lon0, self.lat0 = lon0, lat0
        self.xg, self.yg = to_local_km(obs.lon.values, obs.lat.values, lon0, lat0)
        self.xi = {f: to_local_km(d.lon.values, d.lat.values, lon0, lat0) for f, d in ins.items()}

    def src(self, lon, lat):
        return to_local_km(lon, lat, self.lon0, self.lat0)


def design(geo, obs, ins, lon, lat, depth):
    """Columns: dV (Mm3) | GNSS te, tn | per frame: offset, ramp x, ramp y. Rows: GNSS e,n,u then InSAR blocks."""
    n = len(obs)
    sx, sy = geo.src(lon, lat)
    gx, gy, gz = mogi(geo.xg - sx, geo.yg - sy, depth, 1.0)
    frames = list(ins)
    ncol = 3 + 3 * len(frames)
    rows = [np.zeros((3 * n, ncol))]
    rows[0][:, 0] = np.r_[gx, gy, gz] * 1e9          # mm per Mm3
    rows[0][:n, 1] = 1.0
    rows[0][n:2 * n, 2] = 1.0
    for k, f in enumerate(frames):
        d = ins[f]
        x, y = geo.xi[f]
        ux, uy, uz = mogi(x - sx, y - sy, depth, 1.0)
        B = np.zeros((len(d), ncol))
        B[:, 0] = (ux * d.ue.values + uy * d.un.values + uz * d.uu.values) * 1e9
        B[:, 3 + 3 * k] = 1.0
        B[:, 4 + 3 * k] = x
        B[:, 5 + 3 * k] = y
        rows.append(B)
    return np.vstack(rows)


def solve(geo, obs, ins, sig_ins, lon, lat, depth):
    dg = G.data_vector(obs)
    sg = G.sigmas(dg, len(obs))
    d = np.r_[dg, *[ins[f].los.values for f in ins]]
    s = np.r_[sg, *[np.full(len(ins[f]), sig_ins[f]) for f in ins]]
    A = design(geo, obs, ins, lon, lat, depth)
    m, *_ = np.linalg.lstsq(A / s[:, None], d / s, rcond=None)
    r = (A @ m - d) / s
    ng = len(dg)
    return m, float(np.sum(r ** 2)), float(np.sum(r[:ng] ** 2)), float(np.sum(r[ng:] ** 2)), A @ m


def depth_interval(profile, scale):
    """Depths whose profile chi2 lies within `scale` (= max(chi2/dof, 1)) of the minimum, i.e. a 1-sigma interval
    with errors rescaled by sqrt(chi2/dof). Formal only: InSAR blocks are spatially correlated, which this ignores."""
    ok = profile <= np.nanmin(profile) + scale
    return float(G.DEPTHS[ok].min()), float(G.DEPTHS[ok].max())


def fit_window(win, gn, bg):
    chains, ins, info = {}, {}, {}
    for f in CFG["frames"]:
        names = chain(f, win["start"], win["end"])
        if len(names) < 2:
            continue
        los, tr = chained_los(f, names)
        df = blocks(los, tr, f)
        if len(df) < 50:
            continue
        chains[f], ins[f] = names, df
    if not ins:
        return None, None
    t0 = min(pdate(c[0].split("_")[0]) for c in chains.values())
    t1 = max(pdate(c[-1].split("_")[1]) for c in chains.values())
    obs = G.interval_obs(gn, t0.strftime("%Y-%m-%d"), t1.strftime("%Y-%m-%d"), bg)
    parks = {"lon": CFG.get("ref_lon", G.PARKS["lon"]), "lat": CFG.get("ref_lat", G.PARKS["lat"])}
    sig = {}
    for f, df in ins.items():
        sig[f], nfar = insar_sigma(df, parks)
        info[f] = dict(pairs=len(chains[f]), first=chains[f][0], last=chains[f][-1], blocks=len(df), far_blocks=nfar,
                       sigma_mm=round(sig[f], 2))
    geo = Geometry(obs, ins, G.PARKS["lon"], G.PARKS["lat"])
    best, prof = (np.inf, None), np.full(len(G.DEPTHS), np.inf)
    for k, z in enumerate(G.DEPTHS):
        for lo in G.LONS[::2]:
            for la in G.LATS[::2]:
                m, c2, _, _, _ = solve(geo, obs, ins, sig, lo, la, z)
                if m[0] > 0:
                    prof[k] = min(prof[k], c2)
                    if c2 < best[0]:
                        best = (c2, (lo, la, z))
    lo, la, z = best[1]
    m, c2, cg, ci, pred = solve(geo, obs, ins, sig, lo, la, z)
    ng, ni = 3 * len(obs), sum(len(d) for d in ins.values())
    npar = 3 + 3 + 3 * len(ins)  # lon, lat, depth + dV, te, tn + 3 per frame
    gonly, _, _ = G.fit_grid(obs, G.LONS, G.LATS, G.DEPTHS)
    dof = max(ng + ni - npar, 1)
    dlo, dhi = depth_interval(prof, max(c2 / dof, 1.0))
    row = dict(window=f"{win['interval']} {win['start']}..{win['end']}", gnss_start=str(t0.date()), gnss_end=str(t1.date()),
               n_gnss_sites=len(obs), n_insar_blocks=ni, lon=round(lo, 4), lat=round(la, 4), depth_km=z,
               dv_mm3=round(m[0], 3), depth_lo_km=dlo, depth_hi_km=dhi, chi2_dof_joint=round(c2 / dof, 2),
               chi2_gnss_part=round(cg, 1), chi2_insar_part=round(ci, 1),
               chi2_dof_insar_part=round(ci / max(ni, 1), 2), gnss_only_chi2=round(gonly["chi2"], 1),
               gnss_only_depth_km=gonly["depth_km"], gnss_only_dv_mm3=round(gonly["dv_mm3"], 3),
               gnss_only_chi2_dof=round(gonly["chi2_dof"], 2), gnss_only_lon=round(gonly["lon"], 4),
               gnss_only_lat=round(gonly["lat"], 4))
    return row, dict(frames=info, obs=obs, ins=ins, pred=pred, ng=ng)


def plot(name, res, row):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ins, pred, ng = res["ins"], res["pred"], res["ng"]
    fig, axes = plt.subplots(len(ins), 3, figsize=(12, 2.6 * len(ins)), squeeze=False)
    i0 = ng
    for k, (f, d) in enumerate(ins.items()):
        p = pred[i0:i0 + len(d)]
        i0 += len(d)
        v = np.nanpercentile(np.abs(d.los.values), 98)
        for j, (val, ttl) in enumerate(((d.los.values, "data"), (p, "model"), (d.los.values - p, "residual"))):
            ax = axes[k, j]
            sc = ax.scatter(d.lon, d.lat, c=val, s=30, cmap="RdBu_r", vmin=-v, vmax=v, marker="s", linewidths=0)
            ax.set_aspect(1 / np.cos(np.radians(63.9)))
            ax.plot(row["lon"], row["lat"], "k+", ms=10)
            ax.set_title(f"{f[:4]} {ttl} (mm LOS)", fontsize=9)
            ax.tick_params(labelsize=7)
        fig.colorbar(sc, ax=axes[k, :].tolist(), shrink=0.8)
    fig.suptitle(f"{row['window']}: Mogi {row['depth_km']} km, dV {row['dv_mm3']} Mm3, chi2/dof joint {row['chi2_dof_joint']}",
                 fontsize=10)
    fig.savefig(IMG / f"joint_fit_{name}.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def main():
    events = load_events()
    gn, _ = G.load_all(events)
    v, _, _ = linear_velocity(gn["SENG"], *G.PRE_UNREST)
    bg = {"SENG": float(v.u)}
    rows, chains = [], {}
    for w in CFG["windows"]:
        row, res = fit_window(w, gn, bg)
        if row is None:
            print(f"{w['interval']} {w['start']}: no usable InSAR chain")
            continue
        rows.append(row)
        chains[row["window"]] = res["frames"]
        plot(f"{w['interval']}_{w['start']}", res, row)
        print(row)
    pd.DataFrame(rows).to_csv(DATA / "joint_fits.csv", index=False)
    (DATA / "joint_chains.json").write_text(json.dumps(chains, indent=1))


if __name__ == "__main__":
    main()
