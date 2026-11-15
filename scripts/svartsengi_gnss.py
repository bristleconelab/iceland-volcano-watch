"""Pilot A: Svartsengi magma-domain inflation from NGL GNSS, Mogi point-source fits, volume per cycle.

Steps
 1. Load and clean daily GNSS (SENG, AUSV, GEVK, VOGC, NYLA); vertical-only equipment offsets are detected and removed.
 2. Fit a Mogi source (lon, lat, depth, dV) to the displacement over inflation-only intervals (no dikes), with a common
    horizontal translation as nuisance. Grid search over position and depth, weighted linear least squares for dV and
    translation. Reported per fit: chi2/dof, a depth interval from the chi2 profile with the errors rescaled by
    sqrt(chi2/dof), a leave-one-station-out jackknife, the parameter correlation matrix, and fits with the horizontal
    position held fixed. The Parks et al. (2025) source is fitted with the same weights.
 3. Volume series dV(t) from the SENG vertical with the source fixed, after removing SENG's pre-unrest vertical trend.
    Mogi dV is the volume change of the cavity (chamber), not the volume of magma that flowed in; with compressible
    magma the inflow is larger (e.g. Rivalta & Segall 2008), so volumes and rates here are lower bounds on inflow.
 4. Per event: volume lost = (median of the 5 days before onset) - (lowest 5-day running median in the 7 days after
    onset). Per cycle: recharge from that post-event level to the pre-onset level of the next event. Events whose dike
    came close to SENG are flagged by a rule fixed in config/svartsengi.json that uses seismicity only.
 5. Sensitivity: volumes over the depth range of the fits and the Parks source; pre-unrest background window.

Inputs: NGL .tenv3 (data/raw/ngl), IMO catalogues (data/raw/imo), config/svartsengi.json, docs/data/events.csv.
Outputs: docs/data/*.csv|json, docs/img/*.png, data/processed/gnss_clean/*.csv. Retrospective only, no forecasts.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ivwlib import ROOT, clean_series, linear_velocity, load_events, load_tenv3, mogi, to_local_km, window_mean  # noqa: E402,E501

CFG = json.loads((ROOT / "config" / "svartsengi.json").read_text())
SITES = ["SENG", "AUSV", "GEVK", "VOGC", "NYLA"]
IMG, DATA, PROC = ROOT / "docs" / "img", ROOT / "docs" / "data", ROOT / "data" / "processed" / "gnss_clean"
INTERVALS = [(i["name"], i["start"], i["end"]) for i in CFG["inflation_intervals"]]
SIG_H = CFG["gnss_errors_mm"]["horizontal"]
SIG_U = CFG["gnss_errors_mm"]["vertical"]
SIG_REL = CFG["gnss_errors_mm"]["relative_model_error"]
PRE_UNREST = (CFG["pre_unrest_window"]["start"], CFG["pre_unrest_window"]["end"])
PARKS = {k: CFG["parks2025_source"][k] for k in ("lon", "lat", "depth_km")}
LONS = np.arange(-22.56, -22.34, 0.004)
LATS = np.arange(63.82, 63.93, 0.002)
DEPTHS = np.arange(1.5, 8.01, 0.25)
SVARTSENGI_BOX = (-22.52, 63.80, -22.28, 63.93)


# ------------------------------------------------------------------------------------------------ data
def load_all(events):
    out, steps = {}, {}
    ev_dates = list(events["start"])  # dikes start abruptly; eruption ends cause no step
    for s in SITES:
        raw = load_tenv3(s)
        df, st = clean_series(raw, ev_dates)
        df.attrs.update(raw.attrs)
        out[s], steps[s] = df, st
        PROC.mkdir(parents=True, exist_ok=True)
        df.to_csv(PROC / f"{s}.csv", float_format="%.2f")
        print(f"{s}: {raw.index.min().date()}..{raw.index.max().date()} n={len(raw)} kept={len(df)} offsets={st}")
    return out, steps


def interval_obs(gn, t0, t1, bg):
    rows = []
    for s, df in gn.items():
        a, b = window_mean(df, t0), window_mean(df, t1)
        if a is None or b is None:
            continue
        d = b - a
        if s in bg:  # remove pre-unrest vertical trend (only SENG has one)
            d["u"] -= bg[s] * (pd.Timestamp(t1) - pd.Timestamp(t0)).days / 365.25
        rows.append(dict(site=s, lon=df.attrs["lon"], lat=df.attrs["lat"], de=d["e"], dn=d["n"], du=d["u"]))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------ Mogi fitting
def data_vector(obs):
    return np.r_[obs.de.values, obs.dn.values, obs.du.values]


def sigmas(d, n):
    return np.sqrt(np.r_[np.full(2 * n, SIG_H), np.full(n, SIG_U)] ** 2 + (SIG_REL * d) ** 2)


def design(obs, lon, lat, depth):
    x, y = to_local_km(obs.lon.values, obs.lat.values, lon, lat)
    gx, gy, gz = mogi(x, y, depth, 1.0)
    n = len(obs)
    G = np.zeros((3 * n, 3))
    G[:, 0] = np.r_[gx, gy, gz] * 1e3 * 1e6   # mm per 1e6 m3
    G[:n, 1] = 1.0  # common east translation
    G[n:2 * n, 2] = 1.0  # common north translation
    return G


def solve_linear(obs, lon, lat, depth, weighted=True):
    d = data_vector(obs)
    sig = sigmas(d, len(obs)) if weighted else np.ones_like(d)
    G = design(obs, lon, lat, depth)
    m, *_ = np.linalg.lstsq(G / sig[:, None], d / sig, rcond=None)
    pred = G @ m
    w = sigmas(d, len(obs))  # chi2 always with the same errors, so weighted and unweighted fits are comparable
    return m, pred, float(np.sum(((pred - d) / w) ** 2))


def fit_grid(obs, lons, lats, depths, fix_lonlat=None):
    """Grid search; returns best fit, prediction and the chi2 profile over depth (min over position)."""
    best, chi = (np.inf, None), np.full(len(depths), np.inf)
    lon_grid = [fix_lonlat[0]] if fix_lonlat else lons
    lat_grid = [fix_lonlat[1]] if fix_lonlat else lats
    for k, z in enumerate(depths):
        for lo in lon_grid:
            for la in lat_grid:
                m, _, c2 = solve_linear(obs, lo, la, z)
                if m[0] <= 0:
                    continue
                chi[k] = min(chi[k], c2)
                if c2 < best[0]:
                    best = (c2, dict(lon=lo, lat=la, depth_km=z, dv_mm3=m[0], te=m[1], tn=m[2]))
    c2, p = best
    m, pred, _ = solve_linear(obs, p["lon"], p["lat"], p["depth_km"])
    n = len(obs)
    npar = 4 if fix_lonlat else 6
    d = data_vector(obs)
    p.update(chi2=c2, dof=3 * n - npar, chi2_dof=c2 / max(3 * n - npar, 1),
             rms_mm=float(np.sqrt(np.mean((pred - d) ** 2))))
    return p, pred.reshape(3, n).T, chi


def depth_interval(chi, depths, chi2_dof):
    """Depths whose profile chi2 is within 1 of the minimum after rescaling errors by sqrt(max(chi2/dof, 1))."""
    scale = max(chi2_dof, 1.0)
    ok = depths[chi <= chi.min() + scale]
    return float(ok.min()), float(ok.max())


def correlation_matrix(obs, p):
    """Linearised parameter correlation at the best fit, parameters (x0 km E, y0 km N, depth km, dV, te, tn)."""
    d = data_vector(obs)
    sig = sigmas(d, len(obs))

    def pred(q):
        lon, lat = shift_lonlat(p["lon"], p["lat"], q[0], q[1])
        return design(obs, lon, lat, q[2]) @ q[3:]

    q0 = np.array([0.0, 0.0, p["depth_km"], p["dv_mm3"], p["te"], p["tn"]])
    h = np.array([0.02, 0.02, 0.02, 0.01 * max(p["dv_mm3"], 1e-3), 0.1, 0.1])
    J = np.column_stack([(pred(q0 + np.eye(6)[i] * h[i]) - pred(q0 - np.eye(6)[i] * h[i])) / (2 * h[i])
                         for i in range(6)]) / sig[:, None]
    cov = np.linalg.pinv(J.T @ J) * max(p["chi2_dof"], 1.0)
    sd = np.sqrt(np.diag(cov))
    return cov / np.outer(sd, sd), sd


def shift_lonlat(lon, lat, dx_km, dy_km):
    """Point dx_km east and dy_km north of (lon, lat) (small offsets, local ellipsoidal scale)."""
    x1, _ = to_local_km(lon + 0.01, lat, lon, lat)
    _, y1 = to_local_km(lon, lat + 0.01, lon, lat)
    return lon + dx_km / float(x1) * 0.01, lat + dy_km / float(y1) * 0.01


def fit_intervals(gn, bg):
    fits, resid_rows, corr_out, preds = [], [], {}, {}
    for name, t0, t1 in INTERVALS:
        obs = interval_obs(gn, t0, t1, bg)
        p, pred, chi = fit_grid(obs, LONS, LATS, DEPTHS)
        days = (pd.Timestamp(t1) - pd.Timestamp(t0)).days
        lo, hi = depth_interval(chi, DEPTHS, p["chi2_dof"])
        lo1, hi1 = depth_interval(chi, DEPTHS, 1.0)
        p.update(interval=name, start=t0, end=t1, days=days, n_sites=len(obs), sites="+".join(obs.site),
                 rate_m3s=p["dv_mm3"] * 1e6 / (days * 86400), depth_lo=lo, depth_hi=hi,
                 depth_lo_unscaled=lo1, depth_hi_unscaled=hi1)
        jk = []
        for s in obs.site:
            o2 = obs[obs.site != s].reset_index(drop=True)
            if len(o2) >= 3:
                q, _, _ = fit_grid(o2, LONS[::2], LATS[::2], DEPTHS[::2])
                jk.append((s, q["depth_km"], q["dv_mm3"], q["lon"], q["lat"]))
        p["jackknife"] = jk
        for lab, w in (("parks_src", True), ("parks_src_unweighted", False)):
            m, pr, c2 = solve_linear(obs, PARKS["lon"], PARKS["lat"], PARKS["depth_km"], weighted=w)
            p.update({f"dv_{lab}_mm3": m[0], f"rate_{lab}_m3s": m[0] * 1e6 / (days * 86400),
                      f"rms_{lab}_mm": float(np.sqrt(np.mean((pr - data_vector(obs)) ** 2))),
                      f"chi2_dof_{lab}": c2 / (3 * len(obs) - 3)})
        for lab, ll in (("fixed_parks_xy", (PARKS["lon"], PARKS["lat"])),):
            q, _, _ = fit_grid(obs, LONS, LATS, DEPTHS, fix_lonlat=ll)
            p.update({f"depth_{lab}": q["depth_km"], f"dv_{lab}_mm3": q["dv_mm3"], f"chi2_dof_{lab}": q["chi2_dof"]})
        corr, sd = correlation_matrix(obs, p)
        names = ["x0", "y0", "depth", "dV", "te", "tn"]
        corr_out[name] = dict(params=names, sd_scaled=[round(float(v), 3) for v in sd],
                              corr=[[round(float(v), 2) for v in row] for row in corr])
        fits.append(p)
        preds[name] = (obs, pred, p)
        for i, s in enumerate(obs.site):
            resid_rows.append(dict(interval=name, site=s, obs_e=obs.de[i], obs_n=obs.dn[i], obs_u=obs.du[i],
                                   mod_e=pred[i, 0], mod_n=pred[i, 1], mod_u=pred[i, 2]))
        print(f"{name} {t0}..{t1}: lon {p['lon']:.3f} lat {p['lat']:.3f} depth {p['depth_km']:.2f} "
              f"[{lo}-{hi}] km dV {p['dv_mm3']:.2f} Mm3 chi2/dof {p['chi2_dof']:.1f} rms {p['rms_mm']:.1f} mm | "
              f"Parks src dV {p['dv_parks_src_mm3']:.2f} (unweighted {p['dv_parks_src_unweighted_mm3']:.2f})")
    return pd.DataFrame(fits), pd.DataFrame(resid_rows), corr_out, preds


# ------------------------------------------------------------------------------------------------ volumes
def seng_gz(seng_lonlat, src):
    x, y = to_local_km(*seng_lonlat, src["lon"], src["lat"])
    return float(mogi(x, y, src["depth_km"], 1.0)[2]) * 1e9   # mm per Mm3


def volume_series(seng, bg_vu, sources):
    t_yr = (seng.index - pd.Timestamp(PRE_UNREST[0])).days.values / 365.25
    u_corr = seng["u"].values - bg_vu * t_yr
    vols = pd.DataFrame(index=seng.index)
    gz = {}
    for k, s in sources.items():
        gz[k] = seng_gz((seng.attrs["lon"], seng.attrs["lat"]), s)
        vols[k] = u_corr / gz[k]
    return vols, gz


def pre_level(series, onset, days=5):
    w = series[(series.index < onset) & (series.index >= onset - pd.Timedelta(days=days))]
    return float(w.median()) if len(w) else np.nan


def post_level(series, onset, days=7, win=5):
    """Lowest centred running median (window `win` days) using only days after the onset; returns (level, date)."""
    seg = series[(series.index > onset) & (series.index <= onset + pd.Timedelta(days=days + win // 2))]
    if len(seg) < 3:
        return np.nan, None
    rm = seg.rolling(win, center=True, min_periods=3).median()
    rm = rm[rm.index <= onset + pd.Timedelta(days=days)].dropna()
    if not len(rm):
        return np.nan, None
    return float(rm.min()), rm.idxmin()


def event_losses(vols, events):
    rows = []
    for _, e in events.iterrows():
        s = vols["this_work"]
        pre_old = s[(s.index < e.start) & (s.index >= e.start - pd.Timedelta(days=3))]
        post_old = s[(s.index > e.start) & (s.index <= e.start + pd.Timedelta(days=7))]
        if len(pre_old) == 0 or len(post_old) == 0:
            rows.append(dict(id=e.id, start=e.start.date(), note="no SENG data"))
            continue
        r = dict(id=e.id, start=e.start.date(), kind=e.kind, published_loss=e.published_loss_mm3,
                 published_sd=e.published_loss_sd_mm3,
                 loss_this_work_old_method=pre_old.iloc[-1] - post_old.min())
        for k in vols.columns:
            lvl, when = post_level(vols[k], e.start)
            r[f"loss_{k}"] = pre_level(vols[k], e.start) - lvl
            if k == "this_work":
                r["level_date"] = when.date() if when is not None else None
        rows.append(r)
    return pd.DataFrame(rows)


def cycle_recharge(vols, events, loss, seng_end):
    rec = []
    for i in range(len(events)):
        e = events.iloc[i]
        if pd.isna(loss.loc[i].get("level_date")):
            continue
        t_start = pd.Timestamp(loss.loc[i, "level_date"])
        last = i + 1 >= len(events)
        nxt = seng_end + pd.Timedelta(days=1) if last else events.iloc[i + 1].start
        r = dict(after_event=e.id, from_=t_start.date(), to=(nxt - pd.Timedelta(days=1)).date(),
                 days=(nxt - pd.Timedelta(days=1) - t_start).days,
                 next_event="none (ongoing)" if last else events.iloc[i + 1].id)
        for k in vols.columns:
            start_lvl, _ = post_level(vols[k], e.start)
            r[f"recharge_{k}"] = pre_level(vols[k], nxt) - start_lvl
        r["mean_rate_m3s_this_work"] = r["recharge_this_work"] * 1e6 / max(r["days"], 1) / 86400
        rec.append(r)
    return pd.DataFrame(rec)


def imo_comparison(vols):
    onset = pd.Timestamp("2025-07-16")
    out = []
    for item in CFG["imo_volumes"]:
        d = pd.Timestamp(item["date"])
        r = dict(date=item["date"], imo_mm3=item["volume_mm3"], imo_lo=item["lo"], imo_hi=item["hi"])
        for k in vols.columns:
            ref, when = post_level(vols[k], onset)
            w = vols[k][(vols.index >= d - pd.Timedelta(days=2)) & (vols.index <= d + pd.Timedelta(days=2))]
            r[k] = float(w.mean() - ref) if len(w) else np.nan
            if k == "this_work":
                r["ref_date"] = when.date()
        out.append(r)
    return pd.DataFrame(out)


def quarterly_inflow(vols, start="2025-07-17", end=None):
    """Mean inflow (m3/s, Mogi dV convention) in consecutive 3-month steps from `start`."""
    s = vols["this_work"]
    end = pd.Timestamp(end) if end else s.index.max()
    edges = pd.date_range(pd.Timestamp(start), end, freq=pd.DateOffset(months=3))
    if edges[-1] < end - pd.Timedelta(days=30):
        edges = edges.append(pd.DatetimeIndex([end]))
    rows = []
    for a, b in zip(edges[:-1], edges[1:]):
        va, vb = window_mean(s.to_frame("u").assign(e=0.0, n=0.0), a), window_mean(s.to_frame("u").assign(e=0.0, n=0.0), b)
        if va is None or vb is None:
            continue
        rows.append(dict(start=a.date(), end=b.date(), dv_mm3=round(vb["u"] - va["u"], 2),
                         rate_m3s=round((vb["u"] - va["u"]) * 1e6 / ((b - a).days * 86400), 2)))
    return pd.DataFrame(rows)


def seng_uplift_windows(seng, ends=("2026-09-15", "2026-09-30"), months=3):
    """Raw SENG uplift (cleaned, no background removed) over `months` before each end date, 5-day means."""
    out = []
    for e in ends:
        b = pd.Timestamp(e)
        a = b - pd.DateOffset(months=months)
        va, vb = window_mean(seng, a), window_mean(seng, b)
        out.append(dict(start=a.date(), end=b.date(), uplift_mm=round(float(vb["u"] - va["u"]), 1)))
    return pd.DataFrame(out)


def threshold_crossings(vols, rec, cycle_start):
    """First dates on which the ongoing cycle's recharge reaches the smallest and largest of the five previous
    complete-cycle recharges (retrospective description only)."""
    prev = rec[rec.next_event != "none (ongoing)"].tail(5)
    lo, hi = float(prev.recharge_this_work.min()), float(prev.recharge_this_work.max())
    s = vols["this_work"].rolling(5, center=True, min_periods=3).median()
    base, _ = post_level(vols["this_work"], pd.Timestamp(cycle_start))
    s = s[s.index > pd.Timestamp(cycle_start)] - base
    first = lambda v: str(s[s >= v].index.min().date()) if (s >= v).any() else None  # noqa: E731
    cur = float(rec.recharge_this_work.iloc[-1])
    return dict(previous_five=[round(float(v), 2) for v in prev.recharge_this_work], lo=round(lo, 2), hi=round(hi, 2),
                reached_lo=first(lo), reached_hi=first(hi), current=round(cur, 2),
                current_over_hi_pct=round(100 * (cur / hi - 1), 1))


# ------------------------------------------------------------------------------------------------ sensitivity
def volume_sensitivity(gn, fdf, vols_base, bg_vu):
    seng = gn["SENG"]
    ll = (seng.attrs["lon"], seng.attrs["lat"])
    med = dict(lon=float(np.median(fdf.lon)), lat=float(np.median(fdf.lat)))
    zs = np.arange(fdf.depth_km.min(), fdf.depth_km.max() + 1e-9, 0.25)
    rows = []
    for lab, src in [(f"this_work_xy_z{z:.2f}", dict(med, depth_km=float(z))) for z in zs] + [("parks2025", PARKS)]:
        g = seng_gz(ll, src)
        vols, _ = volume_series(seng, bg_vu, {"this_work": src})
        imo = imo_comparison(vols)
        rows.append(dict(source=lab, lon=src["lon"], lat=src["lat"], depth_km=src["depth_km"],
                         gz_mm_per_mm3=round(g, 2), volume_for_450mm=round(450 / g, 1),
                         **{f"vol_{d}": round(v, 2) for d, v in zip(imo.date, imo.this_work)}))
    return pd.DataFrame(rows)


def background_sensitivity(gn, srcs):
    seng = gn["SENG"]
    rows = []
    windows = [tuple(w) for w in CFG["background_windows_for_sensitivity"]]
    for w in windows:
        v, sv, n = linear_velocity(seng, *w)
        rows.append(dict(window=f"{w[0]}..{w[1]}", vu=round(float(v.u), 2), svu_white=round(float(sv.u), 2)))
    vn, _, _ = linear_velocity(gn["NYLA"], *PRE_UNREST)
    rows.append(dict(window="NYLA rate only (geothermal term set to 0)", vu=round(float(vn.u), 2), svu_white=np.nan))
    out = []
    for r in rows:
        vols, _ = volume_series(seng, r["vu"], {"this_work": srcs["this_work"]})
        imo = imo_comparison(vols)
        q = quarterly_inflow(vols)
        out.append(dict(**r, **{f"vol_{d}": round(v, 2) for d, v in zip(imo.date, imo.this_work)},
                        last_quarter_m3s=float(q.rate_m3s.iloc[-1]) if len(q) else np.nan))
    return pd.DataFrame(out)


def dike_distance(events, seng_lonlat):
    """5th percentile of the epicentral distance from SENG to reviewed IMO earthquakes in the first 24 h of each event
    (36 h from 00:00 when the onset time is not in events.csv), Svartsengi box."""
    cats = []
    for f in ("seiscomp_reykjanes_hengill_2023_2026.csv", "sil_reykjanes_hengill_2020_2026.csv"):
        p = ROOT / "data" / "raw" / "imo" / f
        if p.exists():
            d = pd.read_csv(p)
            d["time"] = pd.to_datetime(d["time"], utc=True, format="ISO8601").dt.tz_localize(None)
            d = d[(d.evaluation_mode == "manual") & (d["type"].isna() | (d["type"] == "earthquake"))]
            b = SVARTSENGI_BOX
            d = d[(d.longitude >= b[0]) & (d.longitude <= b[2]) & (d.latitude >= b[1]) & (d.latitude <= b[3])]
            cats.append((f.split("_")[0], d))
    rows = []
    for _, e in events.iterrows():
        if isinstance(e.start_time_utc, str) and e.start_time_utc:
            t0 = pd.Timestamp(f"{e.start.date()} {e.start_time_utc}")
            t1 = t0 + pd.Timedelta(hours=24)
        else:
            t0, t1 = e.start, e.start + pd.Timedelta(hours=36)
        r = dict(id=e.id, window_start=str(t0), window_end=str(t1))
        for name, d in cats:
            w = d[(d.time >= t0) & (d.time < t1)]
            if len(w) >= 10:
                x, y = to_local_km(w.longitude.values, w.latitude.values, *seng_lonlat)
                r.update(catalogue=name, n_events=len(w), d5_km=round(float(np.percentile(np.hypot(x, y), 5)), 2))
                break
        rows.append(r)
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------ figures
def plot_fits(preds):
    fig, axs = plt.subplots(1, len(preds), figsize=(4.6 * len(preds), 4.6))
    for ax, (name, (obs, pred, p)) in zip(axs, preds.items()):
        x, y = to_local_km(obs.lon.values, obs.lat.values, p["lon"], p["lat"])
        sc = 0.02
        ax.quiver(x, y, obs.de, obs.dn, color="k", angles="xy", scale_units="xy", scale=1 / sc, width=0.006,
                  label="observed")
        ax.quiver(x, y, pred[:, 0], pred[:, 1], color="tab:red", angles="xy", scale_units="xy", scale=1 / sc,
                  width=0.006, label="Mogi model")
        for xi, yi, s, uo, um in zip(x, y, obs.site, obs.du, pred[:, 2]):
            ax.annotate(f"{s}\nU {uo:.0f}/{um:.0f}", (xi, yi), fontsize=7, xytext=(3, -14), textcoords="offset points")
        px, py = to_local_km(PARKS["lon"], PARKS["lat"], p["lon"], p["lat"])
        ax.plot(0, 0, "r*", ms=12, label=f"fit, {p['depth_km']:.2f} km")
        ax.plot(px, py, "bx", ms=9, label="Parks et al. 2025")
        ax.set_title(f"{name}: {p['start']} to {p['end']}\ndV={p['dv_mm3']:.1f} Mm³, chi2/dof={p['chi2_dof']:.1f}",
                     fontsize=9)
        ax.set_xlim(min(x.min(), 0) - 3, max(x.max(), 0) + 4)
        ax.set_ylim(min(y.min(), 0) - 3, max(y.max(), 0) + 3)
        ax.set_aspect("equal")
        ax.set_xlabel("km E of fitted source")
        ax.grid(alpha=.3)
    axs[0].set_ylabel("km N")
    axs[0].legend(fontsize=7, loc="lower left")
    fig.suptitle("Svartsengi inflation: NGL GNSS horizontal displacement vs Mogi fit (U obs/model mm)", fontsize=11)
    fig.tight_layout()
    fig.savefig(IMG / "mogi_fits.png", dpi=130)


def plot_series(seng, vols, events):
    fig, ax = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
    raw = load_tenv3("SENG")
    ax[0].plot(raw.index, raw["u"], ".", ms=1.5, color="0.7", label="SENG up, NGL raw")
    ax[0].plot(seng.index, seng["u"], ".", ms=1.5, color="k", label="cleaned (vertical-only offsets removed)")
    ax[0].set_ylabel("Up (mm)")
    ax[0].legend(fontsize=8, loc="upper left")
    ax[1].plot(vols.index, vols["this_work"], color="tab:red", lw=1, label="this work source")
    ax[1].plot(vols.index, vols["parks2025"], color="tab:blue", lw=1, alpha=.7, label="Parks et al. 2025 source")
    ax[1].set_ylabel("Mogi dV from SENG up (Mm³)")
    ax[1].legend(fontsize=8, loc="upper left")
    for a in ax:
        for _, e in events.iterrows():
            a.axvline(e.start, color="r" if "eruption" in e.kind else "b", lw=.6,
                      ls="-" if "eruption" in e.kind else "--")
            if e.end > e.start:
                a.axvspan(e.start, e.end, color="pink", alpha=.4)
        a.grid(alpha=.3)
    ax[1].set_xlim(pd.Timestamp("2019-06-01"), seng.index.max() + pd.Timedelta(days=20))
    ax[0].set_title("Svartsengi: SENG vertical and Mogi volume change (red = dike+eruption, blue dashed = dike only)")
    fig.tight_layout()
    fig.savefig(IMG / "seng_volume_timeseries.png", dpi=130)


def plot_loss_recharge(loss, rec, events):
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    ok = loss.dropna(subset=["loss_this_work"])
    used = ok[ok.get("passes_rule", True).astype(bool)] if "passes_rule" in ok else ok
    ax[0].errorbar(ok.published_loss, ok.loss_this_work, xerr=ok.published_sd.fillna(0), fmt="o", color="0.6",
                   label="this work source, all events")
    ax[0].plot(used.published_loss, used.loss_this_work, "o", color="tab:red", label="events passing the dike rule")
    ax[0].plot(ok.published_loss, ok.loss_parks2025, "s", mfc="none", color="tab:blue", label="Parks source")
    for _, r in ok.iterrows():
        ax[0].annotate(str(r.id), (r.published_loss, r.loss_this_work), fontsize=8, xytext=(4, 2),
                       textcoords="offset points")
    lim = [0, max(ok.published_loss.max(), ok.loss_this_work.max()) * 1.1]
    ax[0].plot(lim, lim, "k:", lw=.8)
    ax[0].set_xlabel("published volume loss (Mm³; Parks 2025 Table 1, EGU26-6017)")
    ax[0].set_ylabel("this work, SENG-only Mogi (Mm³)")
    ax[0].legend(fontsize=7)
    ax[0].set_title("Volume lost per event (Mogi dV convention)")
    ax[0].grid(alpha=.3)
    if len(rec):
        xi = np.arange(len(rec))
        pub = events.set_index("id").published_loss_mm3
        ax[1].axhspan(17, 23, color="tab:green", alpha=.15,
                      label="17-23 Mm³: recharge before the last five events\n(Parks et al., EGU26-6017)")
        ax[1].bar(xi, rec.recharge_this_work, 0.6, color="tab:orange",
                  label="recharge from post-event level to next event (this work)")
        ax[1].plot(xi, [pub.get(a, np.nan) for a in rec.after_event], "kD", ms=5,
                   label="published loss in the event that opened the cycle")
        last = len(rec) - 1
        ax[1].annotate(f"ongoing, {rec.days.iloc[-1]} days", (xi[last], rec.recharge_this_work.iloc[-1]),
                       xytext=(-60, 8), textcoords="offset points", fontsize=8)
        ax[1].set_xticks(xi)
        ax[1].set_xticklabels([f"{a}→{b}" if b != "none (ongoing)" else f"{a}→now" for a, b in
                               zip(rec.after_event, rec.next_event)], fontsize=8)
        ax[1].set_xlabel("cycle (event id → next event id, see events.csv)")
        ax[1].set_ylabel("Mm³")
        ax[1].set_title("Recharge per cycle (retrospective, not a forecast)")
        ax[1].legend(fontsize=7, loc="upper left")
        ax[1].grid(alpha=.3, axis="y")
    fig.tight_layout()
    fig.savefig(IMG / "volume_loss_and_recharge.png", dpi=130)


# ------------------------------------------------------------------------------------------------ main
def main():
    events = load_events()
    gn, steps = load_all(events)
    (DATA / "gnss_offsets_removed.json").write_text(json.dumps(steps, indent=1))

    bg_rows = []
    for s in ("SENG", "NYLA"):
        v, sv, n = linear_velocity(gn[s], *PRE_UNREST)
        bg_rows.append(dict(site=s, start=PRE_UNREST[0], end=PRE_UNREST[1], n_days=n,
                            ve=v.e, vn=v.n, vu=v.u, sve=sv.e, svn=sv.n, svu=sv.u))
    bgdf = pd.DataFrame(bg_rows).round(2)
    bgdf.to_csv(DATA / "gnss_pre_unrest_velocity.csv", index=False)
    bg = {"SENG": float(bgdf.loc[bgdf.site == "SENG", "vu"].iloc[0])}

    fdf, resid, corr, preds = fit_intervals(gn, bg)
    fdf.drop(columns=["jackknife"]).round(4).to_csv(DATA / "mogi_fits.csv", index=False)
    resid.round(1).to_csv(DATA / "mogi_residuals.csv", index=False)
    (DATA / "mogi_jackknife.json").write_text(json.dumps(dict(zip(fdf.interval, fdf.jackknife)), indent=1,
                                                          default=float))
    (DATA / "mogi_correlation.json").write_text(json.dumps(corr, indent=1))
    plot_fits(preds)

    seng = gn["SENG"]
    srcs = {"this_work": dict(lon=float(np.median(fdf.lon)), lat=float(np.median(fdf.lat)),
                              depth_km=float(np.median(fdf.depth_km))), "parks2025": PARKS}
    vols, gz = volume_series(seng, bg["SENG"], srcs)
    print("SENG vertical Green's function (mm per Mm3):", {k: round(v, 2) for k, v in gz.items()}, srcs)
    vols.round(3).to_csv(DATA / "svartsengi_volume_seng.csv")

    loss = event_losses(vols, events)
    dd = dike_distance(events, (seng.attrs["lon"], seng.attrs["lat"]))
    rule = CFG["event_selection_rule"]
    loss = loss.merge(dd[["id", "catalogue", "n_events", "d5_km"]], on="id", how="left")
    loss["passes_rule"] = loss.d5_km >= 3.0
    loss["ratio"] = loss.loss_this_work / loss.published_loss
    rule_summary = {}
    for thr in rule["thresholds_reported_for_sensitivity_km"]:
        sel = loss[loss.d5_km >= thr]
        rule_summary[f"d5>={thr}"] = dict(events=[int(i) for i in sel.id],
                                          ratio_median=round(float(sel.ratio.median()), 2) if len(sel) else None,
                                          ratio_min=round(float(sel.ratio.min()), 2) if len(sel) else None,
                                          ratio_max=round(float(sel.ratio.max()), 2) if len(sel) else None)
    rule_summary["all"] = dict(ratio_median=round(float(loss.ratio.median()), 2))
    rec = cycle_recharge(vols, events, loss, seng.index.max())
    imo = imo_comparison(vols)
    q = quarterly_inflow(vols)
    upl = seng_uplift_windows(seng)
    thr = threshold_crossings(vols, rec, "2025-07-16")
    (DATA / "cycle_threshold_check.json").write_text(json.dumps(thr, indent=1))
    upl.to_csv(DATA / "seng_uplift_3months.csv", index=False)
    print("threshold:", thr, "\n", upl.to_string(index=False))
    vsens = volume_sensitivity(gn, fdf, vols, bg["SENG"])
    bsens = background_sensitivity(gn, srcs)
    loss.round(2).to_csv(DATA / "event_volume_loss.csv", index=False)
    rec.round(2).to_csv(DATA / "cycle_recharge.csv", index=False)
    imo.round(2).to_csv(DATA / "imo_volume_comparison.csv", index=False)
    q.to_csv(DATA / "quarterly_inflow.csv", index=False)
    vsens.to_csv(DATA / "volume_sensitivity_source.csv", index=False)
    bsens.to_csv(DATA / "volume_sensitivity_background.csv", index=False)
    (DATA / "event_selection_rule.json").write_text(json.dumps(dict(rule=rule, result=rule_summary), indent=1))
    for name, df in (("fits", fdf.drop(columns=["jackknife"])), ("loss", loss), ("recharge", rec), ("imo", imo),
                     ("quarterly", q), ("source sensitivity", vsens), ("background sensitivity", bsens)):
        print(f"--- {name}\n{df.round(2).to_string(index=False)}")
    print(json.dumps(rule_summary, indent=1))
    print(json.dumps(corr, indent=1))
    plot_series(seng, vols, events)
    plot_loss_recharge(loss, rec, events)


if __name__ == "__main__":
    main()
