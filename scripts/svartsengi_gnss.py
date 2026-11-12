"""Pilot A: Svartsengi magma-domain inflation from NGL GNSS, Mogi point-source fits, volume per cycle.

Steps
 1. Load and clean daily GNSS (SENG, AUSV, GEVK, VOGC, NYLA); vertical-only equipment offsets are detected and removed.
 2. Fit a Mogi source (lon, lat, depth, dV) to the displacement over inflation-only intervals (no dikes), with a common
    horizontal translation as nuisance (absorbs most of the plate motion of the near-field network). Grid search over
    position and depth, linear least squares for dV and translation. Leave-one-station-out jackknife for spread.
 3. Volume time series dV(t) from the SENG vertical with the source fixed, after removing SENG's pre-unrest vertical
    trend (2015.5-2019.9: geothermal-field subsidence plus regional motion).
 4. Per event: volume lost (from the day before onset to the minimum within 7 days) vs Parks et al. (2025) Table 1.
    Per cycle: volume recharged before each event vs the previous loss ("volume-threshold check", retrospective only).

Outputs: docs/data/*.csv, docs/img/*.png, data/processed/gnss_clean/*.csv
"""
from __future__ import annotations

import json
import sys

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from ivwlib import (NU, PARKS_SOURCE, ROOT, clean_series, linear_velocity, load_events, load_tenv3,  # noqa: E402
                    mogi, to_local_km, window_mean)

SITES = ["SENG", "AUSV", "GEVK", "VOGC", "NYLA"]
IMG, DATA, PROC = ROOT / "docs" / "img", ROOT / "docs" / "data", ROOT / "data" / "processed" / "gnss_clean"
# inflation-only intervals (no dike between start and end); first/last days chosen 4+ days from any event
INTERVALS = [
    ("I1", "2024-10-03", "2024-11-17"),   # VOGC starts 2024-09-30; AUSV has a gap Aug 2024 - Apr 2025
    ("I2", "2024-12-12", "2025-03-28"),
    ("I3", "2025-04-08", "2025-07-12"),
    ("I4", "2025-08-12", "2026-10-15"),   # current cycle, 14 months
]
SIG_H, SIG_U, SIG_REL = 3.0, 8.0, 0.05  # mm, mm, fraction of the signal (model error)
PRE_UNREST = ("2015-07-01", "2019-12-31")


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


def design(obs, lon, lat, depth):
    x, y = to_local_km(obs.lon.values, obs.lat.values, lon, lat)
    gx, gy, gz = mogi(x, y, depth, 1.0)
    n = len(obs)
    G = np.zeros((3 * n, 3))
    G[:, 0] = np.r_[gx, gy, gz] * 1e3 * 1e6   # mm per 1e6 m3
    G[:n, 1] = 1.0  # common east translation
    G[n:2 * n, 2] = 1.0  # common north translation
    return G


def fit_grid(obs, lons, lats, depths):
    d = np.r_[obs.de.values, obs.dn.values, obs.du.values]
    n = len(obs)
    sig = np.sqrt(np.r_[np.full(2 * n, SIG_H), np.full(n, SIG_U)] ** 2 + (SIG_REL * d) ** 2)
    best = (np.inf, None)
    chi = np.full((len(depths),), np.inf)
    for k, z in enumerate(depths):
        for lo in lons:
            for la in lats:
                G = design(obs, lo, la, z) / sig[:, None]
                m, *_ = np.linalg.lstsq(G, d / sig, rcond=None)
                if m[0] <= 0:
                    continue
                c2 = float(np.sum((G @ m - d / sig) ** 2))
                chi[k] = min(chi[k], c2)
                if c2 < best[0]:
                    best = (c2, dict(lon=lo, lat=la, depth_km=z, dv_mm3=m[0], te=m[1], tn=m[2]))
    c2, p = best
    G = design(obs, p["lon"], p["lat"], p["depth_km"])
    pred = G @ np.array([p["dv_mm3"], p["te"], p["tn"]])
    p.update(chi2=c2, dof=3 * n - 6, rms_mm=float(np.sqrt(np.mean((pred - d) ** 2))))
    return p, pred.reshape(3, n).T, chi


def main():
    events = load_events()
    gn, steps = load_all(events)
    (DATA / "gnss_offsets_removed.json").write_text(json.dumps(steps, indent=1))

    # ---- pre-unrest background (SENG, NYLA)
    bg_rows = []
    for s in ("SENG", "NYLA"):
        v, sv, n = linear_velocity(gn[s], *PRE_UNREST)
        bg_rows.append(dict(site=s, start=PRE_UNREST[0], end=PRE_UNREST[1], n_days=n,
                            ve=v.e, vn=v.n, vu=v.u, sve=sv.e, svn=sv.n, svu=sv.u))
    bgdf = pd.DataFrame(bg_rows).round(2)
    bgdf.to_csv(DATA / "gnss_pre_unrest_velocity.csv", index=False)
    print(bgdf.to_string(index=False))
    bg = {"SENG": float(bgdf.loc[bgdf.site == "SENG", "vu"].iloc[0])}

    # ---- Mogi geometry fits
    lons = np.arange(-22.56, -22.34, 0.004)
    lats = np.arange(63.82, 63.93, 0.002)
    depths = np.arange(1.5, 8.01, 0.25)
    fits, resid_rows = [], []
    fig, axs = plt.subplots(1, len(INTERVALS), figsize=(4.6 * len(INTERVALS), 4.6))
    for ax, (name, t0, t1) in zip(axs, INTERVALS):
        obs = interval_obs(gn, t0, t1, bg)
        p, pred, chi = fit_grid(obs, lons, lats, depths)
        days = (pd.Timestamp(t1) - pd.Timestamp(t0)).days
        p.update(interval=name, start=t0, end=t1, days=days, n_sites=len(obs), sites="+".join(obs.site),
                 rate_m3s=p["dv_mm3"] * 1e6 / (days * 86400))
        # depth range: depths whose best chi2 is within the 68 % band after scaling chi2 to dof
        scale = max(p["chi2"] / max(p["dof"], 1), 1.0)
        ok = depths[chi <= p["chi2"] + 1.0 * scale]
        p.update(depth_lo=float(ok.min()), depth_hi=float(ok.max()))
        # jackknife: drop one site
        jk = []
        for s in obs.site:
            o2 = obs[obs.site != s].reset_index(drop=True)
            if len(o2) < 3:
                continue
            q, _, _ = fit_grid(o2, lons[::2], lats[::2], depths[::2])
            jk.append((s, q["depth_km"], q["dv_mm3"], q["lon"], q["lat"]))
        p["jackknife"] = jk
        # fixed Parks source
        G = design(obs, PARKS_SOURCE["lon"], PARKS_SOURCE["lat"], PARKS_SOURCE["depth_km"])
        dvec = np.r_[obs.de.values, obs.dn.values, obs.du.values]
        m, *_ = np.linalg.lstsq(G, dvec, rcond=None)
        p.update(dv_parks_src_mm3=m[0], rms_parks_src_mm=float(np.sqrt(np.mean((G @ m - dvec) ** 2))),
                 rate_parks_src_m3s=m[0] * 1e6 / (days * 86400))
        fits.append(p)
        for i, s in enumerate(obs.site):
            resid_rows.append(dict(interval=name, site=s, obs_e=obs.de[i], obs_n=obs.dn[i], obs_u=obs.du[i],
                                   mod_e=pred[i, 0], mod_n=pred[i, 1], mod_u=pred[i, 2]))
        print(f"{name} {t0}..{t1}: lon {p['lon']:.3f} lat {p['lat']:.3f} depth {p['depth_km']:.2f} "
              f"[{p['depth_lo']}-{p['depth_hi']}] km dV {p['dv_mm3']:.2f} Mm3 rate {p['rate_m3s']:.2f} m3/s "
              f"rms {p['rms_mm']:.1f} mm | Parks src dV {p['dv_parks_src_mm3']:.2f} rms {p['rms_parks_src_mm']:.1f} | jk {jk}")
        # map panel
        x, y = to_local_km(obs.lon, obs.lat, p["lon"], p["lat"])
        sc = 0.02  # km per mm
        ax.quiver(x, y, obs.de, obs.dn, color="k", angles="xy", scale_units="xy", scale=1 / sc, width=0.006, label="observed")
        ax.quiver(x, y, pred[:, 0], pred[:, 1], color="tab:red", angles="xy", scale_units="xy", scale=1 / sc,
                  width=0.006, label="Mogi model")
        for xi, yi, s, uo, um in zip(x, y, obs.site, obs.du, pred[:, 2]):
            ax.annotate(f"{s}\nU {uo:.0f}/{um:.0f}", (xi, yi), fontsize=7, xytext=(3, -14), textcoords="offset points")
        px, py = to_local_km(PARKS_SOURCE["lon"], PARKS_SOURCE["lat"], p["lon"], p["lat"])
        ax.plot(0, 0, "r*", ms=12, label=f"fit, {p['depth_km']:.1f} km")
        ax.plot(px, py, "bx", ms=9, label="Parks et al. 2025")
        ax.set_title(f"{name}: {t0} to {t1}\ndV={p['dv_mm3']:.1f} Mm³, {p['rate_m3s']:.1f} m³/s", fontsize=9)
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
    fdf = pd.DataFrame(fits)
    fdf.drop(columns=["jackknife"]).round(4).to_csv(DATA / "mogi_fits.csv", index=False)
    pd.DataFrame(resid_rows).round(1).to_csv(DATA / "mogi_residuals.csv", index=False)
    (DATA / "mogi_jackknife.json").write_text(json.dumps({f["interval"]: f["jackknife"] for f in fits}, indent=1,
                                                          default=float))

    # ---- volume series from SENG vertical, source fixed
    seng = gn["SENG"]
    t_yr = (seng.index - pd.Timestamp(PRE_UNREST[0])).days.values / 365.25
    u_corr = seng["u"].values - bg["SENG"] * t_yr
    srcs = {
        "this_work": dict(lon=float(np.median(fdf.lon)), lat=float(np.median(fdf.lat)),
                          depth_km=float(np.median(fdf.depth_km))),
        "parks2025": PARKS_SOURCE,
    }
    vols = pd.DataFrame(index=seng.index)
    gz = {}
    for k, s in srcs.items():
        x, y = to_local_km(seng.attrs["lon"], seng.attrs["lat"], s["lon"], s["lat"])
        gz[k] = float(mogi(x, y, s["depth_km"], 1.0)[2])  # m per m3
        vols[k] = u_corr / 1e3 / gz[k] / 1e6  # Mm3
    print("SENG vertical Green's function (mm per Mm3):", {k: round(v * 1e9, 2) for k, v in gz.items()}, srcs)
    vols.round(3).to_csv(DATA / "svartsengi_volume_seng.csv")

    # per-event loss
    rows = []
    for _, e in events.iterrows():
        pre = vols[(vols.index < e.start) & (vols.index >= e.start - pd.Timedelta(days=3))]
        post = vols[(vols.index > e.start) & (vols.index <= e.start + pd.Timedelta(days=7))]
        if len(pre) == 0 or len(post) == 0:
            rows.append(dict(id=e.id, start=e.start.date(), note="no SENG data"))
            continue
        r = dict(id=e.id, start=e.start.date(), kind=e.kind, published_loss=e.published_loss_mm3,
                 published_sd=e.published_loss_sd_mm3)
        for k in srcs:
            r[f"loss_{k}"] = pre[k].iloc[-1] - post[k].min()
        r["min_date"] = post["this_work"].idxmin().date()
        rows.append(r)
    loss = pd.DataFrame(rows)
    # per-cycle recharge: from the post-event minimum to the day before the next event
    rec = []
    for i in range(len(events)):
        e = events.iloc[i]
        nxt = events.iloc[i + 1].start if i + 1 < len(events) else seng.index.max() + pd.Timedelta(days=1)
        mind = loss.loc[i, "min_date"] if "min_date" in loss and pd.notna(loss.loc[i, "min_date"]) else None
        if mind is None:
            continue
        seg = vols[(vols.index >= pd.Timestamp(mind)) & (vols.index < nxt)]
        if len(seg) < 2:
            continue
        r = dict(after_event=e.id, from_=mind, to=seg.index[-1].date(), days=(seg.index[-1] - pd.Timestamp(mind)).days,
                 next_event=events.iloc[i + 1].id if i + 1 < len(events) else "none (ongoing)")
        for k in srcs:
            r[f"recharge_{k}"] = seg[k].iloc[-1] - seg[k].iloc[0]
            r[f"prev_loss_{k}"] = loss.loc[i, f"loss_{k}"]
        r["mean_rate_m3s_this_work"] = r["recharge_this_work"] * 1e6 / max(r["days"], 1) / 86400
        rec.append(r)
    rec = pd.DataFrame(rec)
    loss.round(2).to_csv(DATA / "event_volume_loss.csv", index=False)
    rec.round(2).to_csv(DATA / "cycle_recharge.csv", index=False)
    print(loss.round(2).to_string(index=False))
    print(rec.round(2).to_string(index=False))

    # IMO comparison points for the current cycle (volume since 16/17 July 2025)
    # IMO counts accumulation from 16/17 July 2025, i.e. from the drawdown minimum at the start of the last eruption
    post = vols.loc[(vols.index >= "2025-07-16") & (vols.index <= "2025-07-23")]
    ref = post.loc[[post["this_work"].idxmin()]]
    imo = []
    if len(ref):
        for d, v, lo, hi in (("2026-07-21", 28.3, 27.5, 29.1), ("2026-08-18", 29.1, 28.3, 29.9)):
            w = vols.loc[(vols.index >= pd.Timestamp(d) - pd.Timedelta(days=2)) &
                         (vols.index <= pd.Timestamp(d) + pd.Timedelta(days=2))]
            if len(w):
                imo.append(dict(date=d, imo_mm3=v, imo_lo=lo, imo_hi=hi,
                                this_work=float(w["this_work"].mean() - ref["this_work"].iloc[0]),
                                parks_src=float(w["parks2025"].mean() - ref["parks2025"].iloc[0]),
                                ref_date=ref.index[0].date()))
    imodf = pd.DataFrame(imo).round(2)
    imodf.to_csv(DATA / "imo_volume_comparison.csv", index=False)
    print(imodf.to_string(index=False))

    # ---- figures: SENG vertical + volume series with events
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
            a.axvline(e.start, color="r" if "eruption" in e.kind else "b", lw=.6, ls="-" if "eruption" in e.kind else "--")
            if e.end > e.start:
                a.axvspan(e.start, e.end, color="pink", alpha=.4)
        a.grid(alpha=.3)
    ax[1].set_xlim(pd.Timestamp("2019-06-01"), seng.index.max() + pd.Timedelta(days=20))
    ax[0].set_title("Svartsengi: SENG vertical and Mogi volume change (red = dike+eruption, blue dashed = dike only)")
    fig.tight_layout()
    fig.savefig(IMG / "seng_volume_timeseries.png", dpi=130)

    # loss comparison + recharge vs previous loss
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    ok = loss.dropna(subset=["loss_this_work"])
    ax[0].errorbar(ok.published_loss, ok.loss_this_work, xerr=ok.published_sd.fillna(0), fmt="o", color="tab:red",
                   label="this work source")
    ax[0].plot(ok.published_loss, ok.loss_parks2025, "s", mfc="none", color="tab:blue", label="Parks source")
    for _, r in ok.iterrows():
        ax[0].annotate(str(r.id), (r.published_loss, r.loss_this_work), fontsize=8, xytext=(4, 2), textcoords="offset points")
    lim = [0, max(ok.published_loss.max(), ok.loss_this_work.max()) * 1.1]
    ax[0].plot(lim, lim, "k:", lw=.8)
    ax[0].set_xlabel("published volume loss (Mm³; Parks 2025 Table 1, EGU26-6017)")
    ax[0].set_ylabel("this work, SENG-only Mogi (Mm³)")
    ax[0].legend(fontsize=8)
    ax[0].set_title("Volume lost per event")
    ax[0].grid(alpha=.3)
    if len(rec):
        xi = np.arange(len(rec))
        pub = events.set_index("id").published_loss_mm3
        ax[1].axhspan(17, 23, color="tab:green", alpha=.15,
                      label="17-23 Mm³: recharge before the last five events\n(Parks et al., EGU26-6017)")
        ax[1].bar(xi, rec.recharge_this_work, 0.6, color="tab:orange",
                  label="recharge from post-event minimum to next event (this work)")
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
        ax[1].set_title("Volume-threshold check (retrospective, not a forecast)")
        ax[1].legend(fontsize=7, loc="upper left")
        ax[1].grid(alpha=.3, axis="y")
    fig.tight_layout()
    fig.savefig(IMG / "volume_loss_and_recharge.png", dpi=130)


if __name__ == "__main__":
    main()
