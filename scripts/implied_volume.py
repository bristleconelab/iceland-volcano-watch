"""What volume would a single station report, depending on where it stands? Radar "virtual stations".

For the windows with both tracks (I2, I3, I4 Aug-Nov 2025) the chained ascending and descending LOS are decomposed into
east and up (track_hv_check.decompose; north motion neglected, LOS weight ~0.1). A plane fitted to the up field more
than 10 km from the source is removed (datum and residual ramp). The up value of every 0.01 deg block is then divided
by the uplift that IMO's point source (Parks et al. 2025: 22.454 W, 63.869 N, 4 km) predicts there per Mm3. The result
is the cavity volume change a GNSS station at that block would report with IMO's model.

1. Ring 2.5-4.5 km from the source (where AUSV and GEVK are): implied volume against azimuth, fitted with
   mean + amplitude cos 2(azimuth - axis). A source elongated along the axis gives this smooth pattern; local
   subsidence (e.g. the Grindavik graben, south of the source) would show as a dip in one sector only.
2. Stations: the radar-implied volume at the station (median of the 3 nearest blocks) against the GNSS-implied
   volume over the same dates, with the GNSS positions at both ends from a one-sided local linear fit (the first
   15 days of the window at its start, the last 15 days at its end), so no data from before an event are used.

    python -I scripts/implied_volume.py -> docs/data/implied_volume_ring.csv, implied_volume_stations.csv,
                                           docs/img/paper/fig_implied_volume.{png,pdf}
Retrospective only; no forecasts.
"""
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams.update({"font.size": 7, "axes.linewidth": 0.6})

sys.path.insert(0, str(Path(__file__).resolve().parent))
import network_volume as N  # noqa: E402
import svartsengi_gnss as G  # noqa: E402
import svartsengi_joint as J  # noqa: E402
import track_hv_check as H  # noqa: E402
from ivwlib import mogi, to_local_km  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs" / "data"
IMG = ROOT / "docs" / "img" / "paper"
P = G.PARKS
RING = (2.5, 4.5)
SITES = ("SENG", "AUSV", "GEVK")


def local_position(df, t, side, days=15):
    """Position (e, n, u) at day t from a linear fit over the `days` after t (side=+1) or before t (side=-1),
    and its 1-sigma error (mm)."""
    t = pd.Timestamp(t)
    lo, hi = (t, t + pd.Timedelta(days=days)) if side > 0 else (t - pd.Timedelta(days=days), t)
    w = df[(df.index >= lo) & (df.index <= hi)]
    if len(w) < 8:
        return None
    x = (w.index - t).days.values.astype(float)
    A = np.c_[np.ones_like(x), x]
    out = {}
    for c in "enu":
        m, *_ = np.linalg.lstsq(A, w[c].values, rcond=None)
        r = w[c].values - A @ m
        out[c], out["s" + c] = float(m[0]), float(r.std(ddof=2) * np.sqrt(np.linalg.inv(A.T @ A)[0, 0]))
    return out


def gz_imo(x, y):
    return mogi(x, y, P["depth_km"], 1e6)[2] * 1e3          # mm per Mm3


def up_field(names):
    blk = {}
    for f, nm in names.items():
        los, tr = J.chained_los(f, nm)
        blk[f] = J.blocks(los, tr, f)
    df = H.decompose(*[blk[f] for f in J.CFG["frames"]])
    x, y = to_local_km(df.lon.values, df.lat.values, P["lon"], P["lat"])
    far = np.hypot(x, y) > 10
    A = np.c_[np.ones(far.sum()), x[far], y[far]]
    m, *_ = np.linalg.lstsq(A, df.u.values[far], rcond=None)
    u = df.u.values - (m[0] + m[1] * x + m[2] * y)
    return df.assign(x=x, y=y, u0=u, V=u / gz_imo(x, y), az=(np.degrees(np.arctan2(x, y)) + 360) % 360,
                     r=np.hypot(x, y))


def main():
    gc, bg = N.load_clean(), N.background()
    ring_rows, st_rows, fields = [], [], {}
    for w in J.CFG["windows"]:
        names = {f: J.chain(f, w["start"], w["end"]) for f in J.CFG["frames"]}
        if any(len(n) < 2 for n in names.values()):
            continue
        key = f"{w['interval']} {w['start']}"
        df = up_field(names)
        fields[key] = df
        ring = df[(df.r > RING[0]) & (df.r < RING[1])]
        th = np.radians(ring.az.values)
        A = np.c_[np.ones(len(ring)), np.cos(2 * th), np.sin(2 * th)]
        c, *_ = np.linalg.lstsq(A, ring.V.values, rcond=None)
        res = ring.V.values - A @ c
        axis = (np.degrees(np.arctan2(c[2], c[1])) / 2) % 180
        for a0 in range(0, 360, 30):
            s = (ring.az >= a0) & (ring.az < a0 + 30)
            mid = np.radians(a0 + 15)
            ring_rows.append(dict(window=key, sector=f"{a0}-{a0 + 30}", n_blocks=int(s.sum()),
                                  implied_median=round(float(np.median(ring.V[s])), 2) if s.any() else np.nan,
                                  cos2_model=round(float(c[0] + c[1] * np.cos(2 * mid) + c[2] * np.sin(2 * mid)), 2),
                                  ring_mean=round(float(c[0]), 2), ring_amplitude=round(float(np.hypot(c[1], c[2])), 2),
                                  axis_deg=round(float(axis), 1), residual_rms=round(float(res.std()), 2)))
        t0 = min(J.pdate(n[0].split("_")[0]) for n in names.values())
        t1 = max(J.pdate(n[-1].split("_")[1]) for n in names.values())
        for s in SITES:
            d = gc[s]
            a, b = local_position(d, t0, +1), local_position(d, t1, -1)
            sx, sy = to_local_km(d.attrs["lon"], d.attrs["lat"], P["lon"], P["lat"])
            k = np.argsort(np.hypot(df.x.values - sx, df.y.values - sy))[:3]
            row = dict(window=key, start=t0.date(), end=t1.date(), site=s,
                       azimuth=round(float((np.degrees(np.arctan2(sx, sy)) + 360) % 360), 1),
                       dist_km=round(float(np.hypot(sx, sy)), 2),
                       insar_implied=round(float(np.median(df.V.values[k])), 2))
            if a is not None and b is not None:
                vu = bg[s][2] if s in bg else 0.0
                du = b["u"] - a["u"] - vu * (t1 - t0).days / 365.25
                g = float(gz_imo(np.array([sx]), np.array([sy]))[0])
                row.update(gnss_implied=round(du / g, 2), gnss_sigma=round(np.hypot(a["su"], b["su"]) / g, 2))
                row["gnss_over_insar"] = round(row["gnss_implied"] / row["insar_implied"], 3)
            st_rows.append(row)
        print(key, f"ring mean {c[0]:.1f} amp {np.hypot(c[1], c[2]):.1f} axis N{axis:.0f}E rms {res.std():.1f}",
              flush=True)
    ring_df, st_df = pd.DataFrame(ring_rows), pd.DataFrame(st_rows)
    ring_df.to_csv(DATA / "implied_volume_ring.csv", index=False)
    st_df.to_csv(DATA / "implied_volume_stations.csv", index=False)
    print(st_df.to_string(index=False))
    plot(fields, ring_df, st_df, gc)


def plot(fields, ring_df, st_df, gc):
    key = "I4 2025-08-12"
    df = fields[key]
    mm = 1 / 25.4
    fig, axs = plt.subplots(1, 2, figsize=(175 * mm, 72 * mm), constrained_layout=True,
                            gridspec_kw=dict(width_ratios=[1.1, 1]))
    ax = axs[0]
    near = df[df.r < 6.0]
    ix = np.floor((near.lon.values + 22.75) / 0.01 + 1e-6).astype(int)
    iy = np.floor((near.lat.values - 63.80) / 0.01 + 1e-6).astype(int)
    grid = np.full((iy.max() + 1, ix.max() + 1), np.nan)
    grid[iy, ix] = near.V.values
    xe, ye = -22.75 + 0.01 * np.arange(ix.max() + 2), 63.80 + 0.01 * np.arange(iy.max() + 2)
    m0 = ring_df[ring_df.window == key].ring_mean.iloc[0]
    pc = ax.pcolormesh(xe, ye, grid, cmap="RdBu_r", vmin=m0 * 0.6, vmax=m0 * 1.4, shading="flat", rasterized=True)
    fig.colorbar(pc, ax=ax, shrink=0.85, label="Implied volume (Mm³)")
    for s in SITES:
        d = gc[s]
        ax.plot(d.attrs["lon"], d.attrs["lat"], "^", ms=5, mfc="white", mec="#222222", mew=0.8)
        ax.annotate(s, (d.attrs["lon"], d.attrs["lat"]), xytext=(3, 3), textcoords="offset points", fontsize=6)
    ax.plot(P["lon"], P["lat"], "*", ms=8, mfc="white", mec="#222222", mew=0.6)
    t = np.linspace(0, 2 * np.pi, 200)
    for rr in RING:
        ax.plot(P["lon"] + rr * np.sin(t) / (111.32 * np.cos(np.radians(P["lat"]))), P["lat"] + rr * np.cos(t) / 111.13,
                color="#555555", lw=0.6, ls=":")
    ax.set_aspect(1 / np.cos(np.radians(63.87)))
    ax.set_xlim(-22.60, -22.30)
    ax.set_ylim(63.81, 63.93)
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")
    ax.set_title("a  Aug-Nov 2025, IMO point source", loc="left", fontsize=7)
    ax = axs[1]
    for k, col in (("I2 2024-12-12", "#999999"), ("I3 2025-04-08", "#1f6fb4"), ("I4 2025-08-12", "#d0582a")):
        f = fields[k]
        ring = f[(f.r > RING[0]) & (f.r < RING[1])]
        mean = ring_df[ring_df.window == k].ring_mean.iloc[0]
        ax.plot(ring.az, ring.V / mean, ".", ms=2.5, color=col, alpha=0.6, label=k.split()[0] + " " + k.split()[1][:7])
        r = ring_df[ring_df.window == k]
        mids = np.arange(15, 360, 30)
        ax.plot(mids, r.cos2_model / mean, "-", color=col, lw=1.0)
    for s in ("AUSV", "GEVK"):
        az = st_df[st_df.site == s].azimuth.iloc[0]
        ax.axvline(az, color="#222222", lw=0.6, ls="--")
        ax.text(az + 3, 1.42, s, fontsize=6, rotation=90, va="top")
    ax.set_xlim(0, 360)
    ax.set_ylim(0.55, 1.45)
    ax.set_xticks(range(0, 361, 60))
    ax.set_xlabel("Azimuth from source (°)")
    ax.set_ylabel("Implied volume / ring mean")
    ax.legend(frameon=False, fontsize=6, loc="lower left", markerscale=3)
    ax.set_title("b  Ring 2.5-4.5 km; lines: cos 2θ fit", loc="left", fontsize=7)
    IMG.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMG / "fig_implied_volume.png", dpi=300)
    fig.savefig(IMG / "fig_implied_volume.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
