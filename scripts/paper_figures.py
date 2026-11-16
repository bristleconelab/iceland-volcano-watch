"""Figures for the station-dependence result (Volcanica widths: 85 mm single, 175 mm double column).

Fig. 1  Up field Apr-Jul 2025 from ascending + descending Sentinel-1 (track_hv_check.decompose) with the GNSS stations,
        the published point source (Parks et al. 2025), the joint point source and the best 4 km line; residuals of
        the up field for the point source and for the line (geometry fixed from the joint fits, dV + offset + plane
        fitted to the up field).
Fig. 2  Cavity volume from the 22 Jul 2025 post-eruption level to 21 Jul 2026 by station set (vertical motion):
        IMO's point source (any-four range for the five stations) against the 6 km line of point sources along N60E
        (range over the three two-track windows' joint-fit geometries), with IMO's 95 % interval.
Fig. 3  Weekly five-station volume since 22 Jul 2025 (vertical): IMO's point source and the 6 km line (I3 geometry),
        SENG alone, and IMO's two published values.

    python -I scripts/paper_figures.py -> docs/img/paper/fig{1,2,3}.{png,pdf}
Inputs are the CSVs written by network_volume.py and elongated_source.py and the cached LiCSAR chains.
"""
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import elongated_source as E  # noqa: E402
import network_volume as N  # noqa: E402
import svartsengi_gnss as G  # noqa: E402
import svartsengi_joint as J  # noqa: E402
import track_hv_check as H  # noqa: E402
from ivwlib import to_local_km  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs" / "data"
OUT = ROOT / "docs" / "img" / "paper"
MM = 1 / 25.4
TWO_TRACK = ("I2 2024-12-12", "I3 2025-04-08", "I4 2025-08-12")
INK, POINT, LINE, IMO_C = "#222222", "#1f6fb4", "#d0582a", "#6b6b6b"
plt.rcParams.update({"font.size": 7, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "savefig.dpi": 300, "font.family": "DejaVu Sans"})


def geom(window, L):
    e = pd.read_csv(DATA / "elongated_source.csv")
    r = e[(e.window == window) & (e.L_km == L)].iloc[0]
    return dict(lon=r.lon, lat=r.lat, depth_km=r.depth_km, L=r.L_km, strike=r.strike)


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.png")
    fig.savefig(OUT / f"{name}.pdf")
    plt.close(fig)


def up_field():
    w = next(x for x in J.CFG["windows"] if x["interval"] == "I3")
    blk = {}
    for f in J.CFG["frames"]:
        los, tr = J.chained_los(f, J.chain(f, w["start"], w["end"]))
        blk[f] = J.blocks(los, tr, f)
    return H.decompose(blk[J.CFG["frames"][0]], blk[J.CFG["frames"][1]])


def residual(df, p):
    x, y = to_local_km(df.lon.values, df.lat.values, p["lon"], p["lat"])
    gz = E.line_greens(x, y, p["depth_km"], p["L"], p["strike"])[2]
    A = np.c_[gz, np.ones_like(x), x, y]
    m, *_ = np.linalg.lstsq(A, df.u.values, rcond=None)
    return df.u.values - A @ m, float(np.sqrt(np.mean((df.u.values - A @ m) ** 2)))


def line_ends(p):
    a = np.radians(p["strike"])
    dlon = 1 / (111.32 * np.cos(np.radians(p["lat"])))
    dx, dy = p["L"] / 2 * np.sin(a), p["L"] / 2 * np.cos(a)
    return ([p["lon"] - dx * dlon, p["lon"] + dx * dlon], [p["lat"] - dy / 111.13, p["lat"] + dy / 111.13])


def fig1(gn):
    df = up_field()
    pt, ln = geom("I3 2025-04-08", 0.0), geom("I3 2025-04-08", 4.0)
    r_pt, rms_pt = residual(df, pt)
    r_ln, rms_ln = residual(df, ln)
    fig, axs = plt.subplots(1, 3, figsize=(175 * MM, 66 * MM), sharey=True, constrained_layout=True)
    box = (-22.62, -22.28, 63.80, 63.95)
    panels = ((df.u.values, "a  Up, 9 Apr-10 Jul 2025", "viridis", (0, np.nanpercentile(df.u, 99))),
              (r_pt, f"b  Point source, rms {rms_pt:.1f} mm", "RdBu_r", (-40, 40)),
              (r_ln, f"c  4 km line N60E, rms {rms_ln:.1f} mm", "RdBu_r", (-40, 40)))
    ix = np.floor((df.lon.values + 22.75) / 0.01 + 1e-6).astype(int)
    iy = np.floor((df.lat.values - 63.80) / 0.01 + 1e-6).astype(int)
    xe = -22.75 + 0.01 * np.arange(ix.max() + 2)
    ye = 63.80 + 0.01 * np.arange(iy.max() + 2)
    for ax, (v, title, cmap, lim) in zip(axs, panels):
        grid = np.full((iy.max() + 1, ix.max() + 1), np.nan)
        grid[iy, ix] = v
        sc = ax.pcolormesh(xe, ye, grid, cmap=cmap, vmin=lim[0], vmax=lim[1], shading="flat", rasterized=True)
        cb = fig.colorbar(sc, ax=ax, orientation="horizontal", shrink=0.8, pad=0.02, aspect=30)
        cb.set_label("mm", labelpad=1)
        for s, d in gn.items():
            if box[0] <= d.attrs["lon"] <= box[1] and box[2] <= d.attrs["lat"] <= box[3]:
                ax.plot(d.attrs["lon"], d.attrs["lat"], "^", ms=4, mfc="white", mec=INK, mew=0.7)
                ax.annotate(s, (d.attrs["lon"], d.attrs["lat"]), xytext=(3, 3), textcoords="offset points", fontsize=6)
        ax.plot(G.PARKS["lon"], G.PARKS["lat"], "*", ms=7, mfc="white", mec=INK, mew=0.6)
        ax.plot(pt["lon"], pt["lat"], "o", ms=3.5, mfc=POINT, mec="white", mew=0.4)
        ax.plot(*line_ends(ln), "-", color=LINE, lw=1.6)
        ax.set_xlim(box[0], box[1])
        ax.set_ylim(box[2], box[3])
        ax.set_aspect(1 / np.cos(np.radians(63.87)))
        ax.set_title(title, loc="left", fontsize=7)
        ax.set_xlabel("Longitude (°E)")
    axs[0].set_ylabel("Latitude (°N)")
    save(fig, "fig1_up_field_residuals")
    return rms_pt, rms_ln


def line_geom(window):
    r = pd.read_csv(DATA / "elongated_source.csv").query("window == @window and L_km == 6.0").iloc[0]
    return dict(lon=r.lon, lat=r.lat, depth_km=r.depth_km, L=6.0, strike=r.strike)


def fig2():
    """Volume 22 Jul 2025 - 21 Jul 2026 by station set: IMO's point source against the 6 km line (vertical)."""
    imo = pd.read_csv(DATA / "imo_volume_comparison.csv").iloc[0]
    st = pd.read_csv(DATA / "network_volume_station.csv").query("source == 'parks2025' and date == '2026-07-21'")
    nv = pd.read_csv(DATA / "network_volume_imo.csv").query("source == 'parks2025' and date == '2026-07-21' "
                                                             "and variant == 'vertical'")
    ln = pd.read_csv(DATA / "elongated_source_imo.csv").query("L_km == 6.0")
    ln = ln[ln.window.isin(TWO_TRACK)]
    rows = (("SENG alone", "SENG", "seng_alone"), ("GEVK alone", "GEVK", "gevk_alone"),
            ("AUSV alone", "AUSV", "ausv_alone"), ("Five stations", None, "all5_vertical"))
    fig, ax = plt.subplots(figsize=(85 * MM, 55 * MM), constrained_layout=True)
    ax.axvspan(imo.imo_lo, imo.imo_hi, color=IMO_C, alpha=0.25, lw=0)
    ax.axvline(imo.imo_mm3, color=IMO_C, lw=0.8)
    ax.text(imo.imo_mm3 + 0.3, len(rows) - 0.45, "IMO", color=IMO_C, fontsize=6)
    for i, (lab, site, col) in enumerate(rows):
        yv = len(rows) - 1 - i
        if site:
            ax.plot(float(st[st.site == site].implied_mm3.iloc[0]), yv + 0.12, "o", ms=3.5, color=POINT,
                    label="IMO point source" if i == 0 else None)
        else:
            full = float(nv[nv.left_out == "none"].dv_mm3.iloc[0])
            lo, hi = nv[nv.left_out != "none"].dv_mm3.min(), nv[nv.left_out != "none"].dv_mm3.max()
            ax.errorbar(full, yv + 0.12, xerr=[[full - lo], [hi - full]], fmt="o", ms=3.5, color=POINT,
                        elinewidth=0.9, capsize=1.5)
            ax.text(hi + 0.3, yv + 0.12, "any four", fontsize=5.5, color=POINT, va="center")
        v = ln[col]
        ax.errorbar(v.median(), yv - 0.12, xerr=[[v.median() - v.min()], [v.max() - v.median()]], fmt="s", ms=3.2,
                    color=LINE, elinewidth=0.9, capsize=1.5,
                    label="6 km source along N60°E" if i == 0 else None)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows][::-1])
    ax.set_xlabel("Cavity volume change, 22 Jul 2025 - 21 Jul 2026 (Mm³)")
    ax.set_xlim(15, 40)
    ax.legend(frameon=False, loc="center right", fontsize=6)
    save(fig, "fig2_volume_by_station_set")


def fig3():
    """Weekly five-station volume (vertical) with IMO's point source and the 6 km line (I3 geometry)."""
    s = pd.read_csv(DATA / "network_volume_series.csv", parse_dates=["date"])
    gclean, bg = N.load_clean(), N.background()
    ln = line_geom("I3 2025-04-08")
    pt, lv = [], []
    for t in s.date:
        obs = N.displacements(gclean, bg, t)
        pt.append(N.fit(obs, G.PARKS, "vertical")[0])
        lv.append(E.fit_network(obs, ln, "vertical")[0])
    imo = pd.read_csv(DATA / "imo_volume_comparison.csv", parse_dates=["date"])
    fig, ax = plt.subplots(figsize=(85 * MM, 58 * MM), constrained_layout=True)
    ax.plot(s.date, s.seng_only_mm3, color=INK, lw=1.0, ls="--", label="SENG alone, IMO point source")
    ax.plot(s.date, pt, color=POINT, lw=1.2, label="Five stations, IMO point source")
    ax.plot(s.date, lv, color=LINE, lw=1.2, label="Five stations, 6 km source (Apr-Jul 2025 fit)")
    ax.errorbar(imo.date, imo.imo_mm3, yerr=[imo.imo_mm3 - imo.imo_lo, imo.imo_hi - imo.imo_mm3], fmt="s", ms=3,
                color=IMO_C, capsize=1.5, elinewidth=0.9, label="IMO (95 %)")
    ax.set_ylabel("Volume since 22 Jul 2025 (Mm³)")
    ax.legend(frameon=False, fontsize=6, loc="upper left")
    ax.set_ylim(0, None)
    fig.autofmt_xdate()
    save(fig, "fig3_volume_series")
    return float(pt[-1]), float(lv[-1])


def main():
    gn = N.load_clean()
    rms = fig1(gn)
    fig2()
    last = fig3()
    print("fig1 rms point/line mm:", [round(x, 2) for x in rms], "| fig3 last point/line", [round(x, 2) for x in last])


if __name__ == "__main__":
    main()
