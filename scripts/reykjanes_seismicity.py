"""Pilot B: seismicity of the Svartsengi-Sundhnúkur sequence (IMO catalogue), 2020-2026.

 1. Dike of 10-11 Nov 2023: strike from the principal axis of the epicentres, along-strike distance vs time, length,
    depth range, speed of the seismic front.
 2. Daily rate in the Svartsengi box against the GNSS volume series (Pilot A) and the event list; rank correlation of
    the inter-event seismicity rate with the volume recharged since the previous event.
 3. Completeness (Mc: MAXC+0.2, goodness of fit, b-stability) and b-value (Aki-Utsu, bootstrap, b-positive) for
    the dike swarm, inter-event inflation periods, the current cycle and Hengill.

Input: CSV from scripts/fetch_imo_quakes.py. Only reviewed ('manual') earthquakes are used.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import seislib as sl  # noqa: E402
from ivwlib import PARKS_SOURCE, ROOT, load_events, to_local_km  # noqa: E402

IMG, DATA = ROOT / "docs" / "img", ROOT / "docs" / "data"
SVARTSENGI_BOX = (-22.52, 63.80, -22.28, 63.93)   # Svartsengi, Sundhnúkur crater row, Grindavík
HENGILL_BOX = (-21.50, 63.98, -21.15, 64.15)


def inbox(df, b):
    return df[(df.longitude >= b[0]) & (df.longitude <= b[2]) & (df.latitude >= b[1]) & (df.latitude <= b[3])]


def load(path):
    d = pd.read_csv(path)
    d["time"] = pd.to_datetime(d["time"], utc=True, format="ISO8601").dt.tz_localize(None)
    n0 = len(d)
    d = d[d.evaluation_mode == "manual"]
    if "type" in d:
        d = d[d["type"].isna() | (d["type"] == "earthquake")]
    d = d.dropna(subset=["magnitude"]).sort_values("time").reset_index(drop=True)
    print(f"{path}: {n0} rows, {len(d)} reviewed earthquakes, {d.time.min()} .. {d.time.max()}")
    return d


def stats_block(name, m_timeordered, dm=0.1, rng=0):
    m = np.asarray(m_timeordered)
    out = dict(set=name, n_all=int(m.size))
    if m.size < 50:
        return out
    maxc = sl.mc_maxc(m, dm, 0.2)
    gft, r = sl.mc_gft(m, dm)
    mbs = sl.mc_mbs(m, dm)
    mc = float(np.nanmax([maxc, gft, mbs]))  # conservative choice
    b, sb, n, a = sl.b_aki(m, mc, dm)
    bb = sl.b_bootstrap(m, mc, dm, 500, rng)
    bp, sbp, nbp = sl.b_positive(m, dm, 0.2)
    out.update(mc_maxc02=maxc, mc_gft=gft, gft_R=r, mc_mbs=mbs, mc_used=mc, n_ge_mc=int(n), b=b, b_shi_bolt=sb,
               b_boot_sd=bb, b_positive=bp, b_positive_sd=sbp, n_bpos=int(nbp), m_max=float(m.max()))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--label", default="IMO SIL")
    a = ap.parse_args()
    cat = load(a.catalog)
    ev = load_events()
    sv = inbox(cat, SVARTSENGI_BOX).reset_index(drop=True)

    # ------------------------------------------------------------------ 1. dike of 10 Nov 2023
    t0, t1 = pd.Timestamp("2023-11-10 12:00"), pd.Timestamp("2023-11-12 00:00")
    dk = inbox(cat, (-22.65, 63.75, -22.20, 63.95))
    dk = dk[(dk.time >= t0) & (dk.time < t1)].reset_index(drop=True)
    x, y = to_local_km(dk.longitude, dk.latitude, PARKS_SOURCE["lon"], PARKS_SOURCE["lat"])
    xy = np.c_[x, y]
    c = np.median(xy, axis=0)
    w, v = np.linalg.eigh(np.cov((xy - c).T))
    ax_ = v[:, np.argmax(w)]
    if ax_[1] < 0:
        ax_ = -ax_
    strike = float(np.degrees(np.arctan2(ax_[0], ax_[1])) % 180)
    s = (xy - c) @ ax_          # along strike, km (positive towards NE)
    nrm = (xy - c) @ np.array([ax_[1], -ax_[0]])
    dk["along_km"], dk["across_km"] = s, nrm
    main = dk[dk.time >= pd.Timestamp("2023-11-10 15:00")]
    lo, hi = np.percentile(main.along_km, [2.5, 97.5])
    # front: along-strike extent (2.5-97.5 %) in 1-hour bins
    hb = main.groupby(main.time.dt.floor("h"))
    front = pd.DataFrame({"n": hb.size(), "p05": hb.along_km.quantile(0.05), "p95": hb.along_km.quantile(0.95),
                          "depth_med": hb.depth.median()})
    front = front[front.n >= 20]
    dike = dict(window=f"{t0} .. {t1}", n_events=len(dk), n_after_15utc=len(main), strike_deg=round(strike, 1),
                length_km_2p5_97p5=round(hi - lo, 1), across_sd_km=round(float(np.std(main.across_km)), 2),
                depth_p10_p50_p90=[round(float(q), 1) for q in np.percentile(main.depth, [10, 50, 90])],
                m_max=float(dk.magnitude.max()), n_m_ge_3=int((dk.magnitude >= 3).sum()))
    # SW and NE extent growth between first hour with >=20 events and the hour of maximum extent
    if len(front) > 2:
        ext = front.p95 - front.p05
        dike["extent_first_hour_km"] = round(float(ext.iloc[0]), 1)
        dike["extent_max_km"] = round(float(ext.max()), 1)
        dike["hour_of_max_extent"] = str(ext.idxmax())
        hrs = (ext.idxmax() - ext.index[0]).total_seconds() / 3600
        dike["growth_km_per_h"] = round(float((ext.max() - ext.iloc[0]) / max(hrs, 1)), 2)
    print("dike:", dike)
    (DATA / "dike_2023-11-10.json").write_text(json.dumps(dike, indent=1, default=str))

    fig, axs = plt.subplots(1, 3, figsize=(16, 5.2), gridspec_kw=dict(width_ratios=[1, 1.4, 1]))
    hrs = (dk.time - pd.Timestamp("2023-11-10 00:00")).dt.total_seconds() / 3600
    sc = axs[0].scatter(dk.longitude, dk.latitude, c=hrs, s=2 + 3 * np.clip(dk.magnitude, 0, 5) ** 2, cmap="viridis")
    axs[0].plot(PARKS_SOURCE["lon"], PARKS_SOURCE["lat"], "r*", ms=12, label="Svartsengi Mogi source (Parks 2025)")
    axs[0].plot(-22.4337, 63.8822, "k^", label="SENG GNSS")
    axs[0].plot(-22.43, 63.843, "ks", mfc="none", label="Grindavík")
    axs[0].set_aspect(1 / np.cos(np.radians(63.87)))
    axs[0].legend(fontsize=7, loc="upper left")
    axs[0].set_title(f"{len(dk)} reviewed events, strike N{strike:.0f}E", fontsize=9)
    plt.colorbar(sc, ax=axs[0], label="hours after 10 Nov 00:00 UTC")
    axs[1].scatter(hrs, dk.along_km, s=2, c="k", alpha=.4)
    if len(front):
        fh = (front.index - pd.Timestamp("2023-11-10 00:00")).total_seconds() / 3600
        axs[1].plot(fh, front.p05, "b-", label="5 % along strike")
        axs[1].plot(fh, front.p95, "r-", label="95 % along strike")
    axs[1].set_xlabel("hours after 10 Nov 2023 00:00 UTC")
    axs[1].set_ylabel("km along strike (NE +)")
    axs[1].legend(fontsize=8)
    axs[1].set_title("Migration of the seismicity along the dike", fontsize=9)
    axs[1].grid(alpha=.3)
    axs[2].scatter(dk.along_km, dk.depth, s=2, c=hrs, cmap="viridis")
    axs[2].invert_yaxis()
    axs[2].set_xlabel("km along strike (NE +)")
    axs[2].set_ylabel("depth (km)")
    axs[2].set_title("Depth section", fontsize=9)
    axs[2].grid(alpha=.3)
    fig.suptitle(f"10-11 November 2023 Sundhnúkur dike, {a.label} reviewed catalogue", fontsize=11)
    fig.tight_layout()
    fig.savefig(IMG / "dike_2023-11-10.png", dpi=130)

    # ------------------------------------------------------------------ 2. rate vs inflation
    vol = pd.read_csv(DATA / "svartsengi_volume_seng.csv", parse_dates=["date"], index_col="date")["this_work"]
    # one completeness level for all inter-event periods: the larger Mc of the 2023-12..2025-07 inter-event set and of
    # the current cycle (the dike swarm itself is excluded; its short-term incompleteness would inflate Mc)
    st1 = stats_block("a", sv[(sv.time >= "2023-12-22") & (sv.time < "2025-07-16")].magnitude.values)
    st2 = stats_block("b", sv[sv.time >= "2025-08-06"].magnitude.values)
    mc_sv = float(max(st1.get("mc_used", 1.0), st2.get("mc_used", 1.0)))
    daily = sv[sv.magnitude >= mc_sv].set_index("time").resample("D").size()
    daily_all = sv.set_index("time").resample("D").size()
    rows = []
    for i in range(len(ev) - 1 + 1):
        e = ev.iloc[i]
        start = (e.end if e.end > e.start else e.start) + pd.Timedelta(days=2)
        end = ev.iloc[i + 1].start - pd.Timedelta(days=1) if i + 1 < len(ev) else daily.index.max()
        if end <= start:
            continue
        nd = (end - start).days
        q = sv[(sv.time >= start) & (sv.time < end) & (sv.magnitude >= mc_sv)]
        v0 = vol[(vol.index >= start - pd.Timedelta(days=1)) & (vol.index <= start + pd.Timedelta(days=1))]
        v1 = vol[(vol.index >= end - pd.Timedelta(days=2)) & (vol.index <= end)]
        if len(v0) == 0 or len(v1) == 0:
            continue
        # rate in the first and the last third of the interval
        thirds = np.array_split(np.arange(nd), 3)
        cnt = daily.reindex(pd.date_range(start.normalize(), periods=nd, freq="D"), fill_value=0).values
        last14 = cnt[-14:].mean() if nd >= 28 else np.nan
        first14 = cnt[:14].mean() if nd >= 28 else np.nan
        rows.append(dict(after_event=int(e.id), rate_first14=first14, rate_last14=last14, start=start.date(), end=end.date(), days=nd, n_ge_mc=len(q),
                         rate_per_day=len(q) / nd, rate_first_third=cnt[thirds[0]].mean(),
                         rate_last_third=cnt[thirds[-1]].mean(), dv_mm3=float(v1.mean() - v0.mean()),
                         inflow_m3s=float((v1.mean() - v0.mean()) * 1e6 / (nd * 86400))))
    rr = pd.DataFrame(rows)
    if len(rr) > 3:
        rho1, p1 = spearmanr(rr.inflow_m3s, rr.rate_per_day)
        rho2, p2 = spearmanr(rr.dv_mm3, rr.rate_last_third)
        print(f"Spearman rate vs inflow rate: rho={rho1:.2f} p={p1:.3f}; late rate vs recharged volume rho={rho2:.2f} p={p2:.3f}")
        rr.attrs = {}
        (DATA / "seismicity_vs_inflation_spearman.json").write_text(json.dumps(dict(
            mc=mc_sv, rate_vs_inflow_rho=rho1, rate_vs_inflow_p=p1, late_rate_vs_dv_rho=rho2, late_rate_vs_dv_p=p2,
            n_intervals=len(rr)), indent=1))
    rr.round(3).to_csv(DATA / "seismicity_vs_inflation.csv", index=False)
    print(rr.round(2).to_string(index=False))

    fig, ax = plt.subplots(2, 1, figsize=(13, 7.5), sharex=True)
    ax[0].bar(daily_all.index, daily_all.values, width=1, color="0.75", label="all reviewed events")
    ax[0].bar(daily.index, daily.values, width=1, color="k", label=f"M ≥ Mc = {mc_sv:.1f}")
    ax[0].set_yscale("symlog", linthresh=10)
    ax[0].set_ylabel("events per day (Svartsengi box)")
    ax[0].legend(fontsize=8, loc="upper right")
    ax[1].plot(vol.index, vol.values, color="tab:red")
    ax[1].set_ylabel("Mogi dV from SENG (Mm³)")
    for a_ in ax:
        for _, e in ev.iterrows():
            a_.axvline(e.start, color="r" if "eruption" in e.kind else "b", lw=.6,
                       ls="-" if "eruption" in e.kind else "--")
            if e.end > e.start:
                a_.axvspan(e.start, e.end, color="pink", alpha=.4)
        a_.grid(alpha=.3)
    ax[1].set_xlim(pd.Timestamp("2023-10-01"), daily.index.max())
    ax[0].set_title(f"Seismicity rate and magma-domain volume, Svartsengi ({a.label})")
    fig.tight_layout()
    fig.savefig(IMG / "seismicity_rate_vs_volume.png", dpi=130)

    # ------------------------------------------------------------------ 3. Mc / b-value table
    he = inbox(cat, HENGILL_BOX)
    sets = [
        ("Dike swarm 10-11 Nov 2023", dk.magnitude.values),
        ("Svartsengi box, inter-event 2023-12-22..2025-07-15",
         sv[(sv.time >= "2023-12-22") & (sv.time < "2025-07-16")].magnitude.values),
        ("Svartsengi box, current cycle 2025-08-06..", sv[sv.time >= "2025-08-06"].magnitude.values),
        ("Svartsengi box, 2020-01-01..2023-10-24 (pre-dike unrest)",
         sv[(sv.time >= "2020-01-01") & (sv.time < "2023-10-25")].magnitude.values),
        ("Hengill box, whole catalogue", he.magnitude.values),
    ]
    tab = pd.DataFrame([stats_block(n, m) for n, m in sets])
    tab.round(3).to_csv(DATA / "bvalue_table.csv", index=False)
    print(tab.round(2).to_string(index=False))

    fig, ax = plt.subplots(1, len(sets), figsize=(4 * len(sets), 4), sharey=True)
    for a_, (n, m), (_, r) in zip(ax, sets, tab.iterrows()):
        if len(m) < 50:
            continue
        bins, nc, cum = sl.fmd(m)
        a_.semilogy(bins, cum, "k.", ms=4, label="cumulative")
        a_.semilogy(bins, np.maximum(nc, 0.8), "^", color="0.6", ms=3, label="per bin")
        mm = np.linspace(r.mc_used, bins.max(), 20)
        a_.semilogy(mm, 10 ** (np.log10(r.n_ge_mc) - r.b * (mm - r.mc_used)), "r-",
                    label=f"b={r.b:.2f}±{r.b_boot_sd:.2f}\nMc={r.mc_used:.1f}, N={r.n_ge_mc}")
        a_.axvline(r.mc_used, color="r", ls=":")
        a_.set_title(n, fontsize=8)
        a_.legend(fontsize=7)
        a_.set_xlabel("magnitude")
        a_.grid(alpha=.3)
    ax[0].set_ylabel("N")
    fig.tight_layout()
    fig.savefig(IMG / "fmd_bvalues.png", dpi=130)

    # ------------------------------------------------------------------ map
    fig, ax = plt.subplots(figsize=(11, 6))
    rk = inbox(cat, (-22.9, 63.75, -21.0, 64.2))
    rk = rk[rk.time >= "2020-01-01"]
    yr = rk.time.dt.year + rk.time.dt.day_of_year / 366
    sc = ax.scatter(rk.longitude, rk.latitude, c=yr, s=0.3 + 0.6 * np.clip(rk.magnitude, 0, 5) ** 2, cmap="turbo",
                    rasterized=True)
    for b, nm in ((SVARTSENGI_BOX, "Svartsengi box"), (HENGILL_BOX, "Hengill box")):
        ax.plot([b[0], b[2], b[2], b[0], b[0]], [b[1], b[1], b[3], b[3], b[1]], "k-", lw=.8)
        ax.text(b[0], b[3] + 0.005, nm, fontsize=8)
    ax.set_aspect(1 / np.cos(np.radians(64)))
    plt.colorbar(sc, label="year")
    ax.set_title(f"Reviewed earthquakes 2020-2026, {a.label} ({len(rk)} events)")
    fig.tight_layout()
    fig.savefig(IMG / "seismicity_map.png", dpi=130)


if __name__ == "__main__":
    main()
