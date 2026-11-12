"""Pilot C: geothermal-field deformation at Hengill (Hellisheiði, Nesjavellir) and Svartsengi from LiCSAR + NGL GNSS.

InSAR: summer-to-summer (~1 yr) LiCSAR interferograms (scripts/licsar_annual.py), scaled to mm/yr. For each year the
median of up to three pairs is used; the multi-year rate is the median of the yearly maps. Field rates are measured
relative to a ring 4-7 km from the field centre, which removes the regional (plate-spreading) gradient to first order,
in the same spirit as the "after plate-motion correction" values in the literature.
Where both geometries exist, LOS rates are decomposed into east and vertical (north ignored: the LOS has a north
sensitivity of < 0.13).

GNSS: HUSM (Húsmúli re-injection area, 2012-2026) with offsets at its equipment changes (NGL steps file) and an annual
term estimated by least squares; vertical rate relative to REYK (Reykjavík).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import rasterio
from scipy.ndimage import uniform_filter

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ivwlib import ROOT, load_tenv3, to_local_km  # noqa: E402

LIC = ROOT / "data" / "raw" / "licsar"
IMG, DATA = ROOT / "docs" / "img", ROOT / "docs" / "data"
FRAMES = {"asc": "016A_02504_162118", "desc": "155D_02484_191814"}
FIELDS = {  # approximate centres of the production fields (power-plant / well-field area)
    "Hellisheiði": (-21.40, 64.04),
    "Nesjavellir": (-21.26, 64.11),
    "Svartsengi": (-22.43, 63.88),
    "Reykjanes": (-22.68, 63.83),
}
CC_MIN = 0.25
PRE_UNREST_YEARS = [2015, 2016, 2017, 2018]  # primaries; the 2018->2019 pair ends before the Jan 2020 unrest


def load_frame(frame):
    d = LIC / frame
    pairs = sorted(p.stem for p in (d / "los_mm").glob("*.tif"))
    stack, years, dts = [], [], []
    for p in pairs:
        with rasterio.open(d / "los_mm" / f"{p}.tif") as ds:
            a = ds.read(1)
            tr, shape = ds.transform, ds.shape
        ccf = d / "cc" / f"{p}.tif"
        if ccf.exists():
            with rasterio.open(ccf) as ds:
                cc = ds.read(1)
            a = np.where(cc >= CC_MIN, a, np.nan)
        t0, t1 = pd.Timestamp(p[:8]), pd.Timestamp(p[9:])
        dt = (t1 - t0).days
        stack.append(a * 365.25 / dt)
        years.append(t0.year)
        dts.append(dt)
    meta = {}
    for k in ("E", "N", "U"):
        with rasterio.open(d / "meta" / f"{k}.tif") as ds:
            m = ds.read(1)
            meta[k] = np.where(m == 0, np.nan, m)
    return pairs, np.array(stack), np.array(years), tr, shape, meta


def grid_lonlat(tr, shape):
    rows, cols = np.mgrid[0:shape[0], 0:shape[1]]
    lon = tr.c + (cols + 0.5) * tr.a
    lat = tr.f + (rows + 0.5) * tr.e
    return lon, lat


def field_rate(vmap, lon, lat, centre, r_in=1.5, ring=(4.0, 7.0)):
    """Peak (5th percentile of a 3x3-smoothed map within r_in km) minus ring median; also the centre mean."""
    x, y = to_local_km(lon, lat, *centre)
    r = np.hypot(x, y)
    sm = uniform_filter(np.nan_to_num(vmap), 3) / np.maximum(uniform_filter(np.isfinite(vmap).astype(float), 3), 1e-6)
    sm[~np.isfinite(vmap)] = np.nan
    ringv = vmap[(r >= ring[0]) & (r <= ring[1])]
    inner = sm[r <= r_in]
    if np.isfinite(ringv).sum() < 50 or np.isfinite(inner).sum() < 10:
        return np.nan, np.nan, np.nan
    ref = np.nanmedian(ringv)
    return float(np.nanpercentile(inner, 5) - ref), float(np.nanmedian(inner) - ref), float(np.isfinite(inner).mean())


def yearly_maps(stack, years):
    out = {}
    for y in sorted(set(years)):
        s = stack[years == y]
        out[y] = np.nanmedian(s, axis=0) if len(s) > 1 else s[0]
    return out


def husm_rate():
    rows = []
    ref = load_tenv3("REYK")
    d = load_tenv3("HUSM")
    j = d.join(ref[["e", "n", "u"]], rsuffix="_ref", how="inner")
    t = (j.index - pd.Timestamp("2012-01-01")).days.values / 365.25
    # NGL steps.txt equipment changes; keep only steps with data on both sides of the previous one
    steps = [pd.Timestamp("2016-07-26"), pd.Timestamp("2017-02-09")]
    steps = [s for i, s in enumerate(steps)
             if i + 1 == len(steps) or ((j.index >= s) & (j.index < steps[i + 1])).sum() > 30]
    for comp in ("e", "n", "u"):
        y = (j[comp] - j[comp + "_ref"]).values
        cols = [np.ones_like(t), t, np.sin(2 * np.pi * t), np.cos(2 * np.pi * t)]
        cols += [(j.index >= s).astype(float) for s in steps]
        A = np.c_[tuple(cols)]
        ok = np.isfinite(y)
        m, *_ = np.linalg.lstsq(A[ok], y[ok], rcond=None)
        r = y[ok] - A[ok] @ m
        cov = np.linalg.inv(A[ok].T @ A[ok]) * (r @ r) / (ok.sum() - A.shape[1])
        rows.append(dict(component=comp, rate_mm_yr=m[1], sigma_white=np.sqrt(cov[1, 1]),
                         offsets_mm=[round(float(v), 1) for v in m[4:]], rms_mm=float(np.sqrt(np.mean(r ** 2))),
                         n=int(ok.sum()), start=str(j.index.min().date()), end=str(j.index.max().date())))
    # period rates for the vertical (offset-corrected series)
    y = (j["u"] - j["u_ref"]).values.copy()
    um = next(r for r in rows if r["component"] == "u")
    for s, off in zip(steps, um["offsets_mm"]):
        y[j.index >= s] -= off
    per = []
    for a, b in (("2012-06-01", "2015-12-31"), ("2017-03-01", "2019-12-31"), ("2020-01-01", "2026-12-31")):
        w = (j.index >= a) & (j.index <= b)
        tt = t[w]
        if w.sum() > 100:
            A = np.c_[np.ones_like(tt), tt, np.sin(2 * np.pi * tt), np.cos(2 * np.pi * tt)]
            m, *_ = np.linalg.lstsq(A, y[w], rcond=None)
            per.append(dict(period=f"{a[:7]}..{b[:7]}", vu_rel_REYK=round(float(m[1]), 2), n=int(w.sum())))
    return rows, per, j.index, y


def main():
    res, maps = [], {}
    for geom, frame in FRAMES.items():
        if not (LIC / frame / "los_mm").exists():
            print("missing", frame)
            continue
        pairs, stack, years, tr, shape, meta = load_frame(frame)
        lon, lat = grid_lonlat(tr, shape)
        ym = yearly_maps(stack, years)
        print(f"{geom} {frame}: {len(pairs)} pairs, years {sorted(ym)}")
        for y, vm in ym.items():
            for f, c in FIELDS.items():
                pk, md, cov = field_rate(vm, lon, lat, c)
                res.append(dict(geom=geom, frame=frame, year=f"{y}-{y + 1}", field=f, peak_los_mm_yr=pk,
                                centre_los_mm_yr=md, coverage=cov, n_pairs=int((years == y).sum())))
        pre = [ym[y] for y in PRE_UNREST_YEARS if y in ym]
        allm = list(ym.values())
        maps[geom] = dict(pre=np.nanmedian(pre, axis=0), all=np.nanmedian(allm, axis=0), lon=lon, lat=lat,
                          tr=tr, meta=meta, years_pre=[y for y in PRE_UNREST_YEARS if y in ym], years_all=sorted(ym))
    df = pd.DataFrame(res)
    df.round(2).to_csv(DATA / "geothermal_insar_yearly.csv", index=False)

    # summary: medians of yearly field rates, pre-unrest and all years
    summ = []
    for (g, f), grp in df.groupby(["geom", "field"]):
        for lab, sel in (("pre-unrest 2015-2019", grp.year.str[:4].astype(int).isin(PRE_UNREST_YEARS)),
                         ("all years", np.ones(len(grp), bool))):
            q = grp[sel].dropna(subset=["peak_los_mm_yr"])
            if len(q):
                summ.append(dict(geom=g, field=f, period=lab, n_years=len(q),
                                 peak_los_median=q.peak_los_mm_yr.median(), peak_los_min=q.peak_los_mm_yr.min(),
                                 peak_los_max=q.peak_los_mm_yr.max(), centre_los_median=q.centre_los_mm_yr.median()))
    sdf = pd.DataFrame(summ).round(1)
    # vertical/east decomposition of the pre-unrest maps
    if "asc" in maps and "desc" in maps and maps["asc"]["pre"].shape == maps["desc"]["pre"].shape:
        A, D = maps["asc"], maps["desc"]
        ea, ua, ed, ud = A["meta"]["E"], A["meta"]["U"], D["meta"]["E"], D["meta"]["U"]
        det = ea * ud - ed * ua
        for key in ("pre", "all"):
            vu = (ea * D[key] - ed * A[key]) / det
            ve = (ud * A[key] - ua * D[key]) / det
            maps.setdefault("vert", {})[key] = vu
            maps.setdefault("east", {})[key] = ve
        for f, c in FIELDS.items():
            for key, lab in (("pre", "pre-unrest 2015-2019"), ("all", "all years")):
                pk, md, cov = field_rate(maps["vert"][key], A["lon"], A["lat"], c)
                sdf = pd.concat([sdf, pd.DataFrame([dict(geom="vertical (asc+desc)", field=f, period=lab,
                                                         peak_los_median=round(pk, 1), centre_los_median=round(md, 1))])])
    sdf.to_csv(DATA / "geothermal_insar_summary.csv", index=False)
    print(sdf.to_string(index=False))

    # GNSS HUSM
    rows, per, tidx, yu = husm_rate()
    (DATA / "husm_gnss.json").write_text(json.dumps(dict(full=rows, periods=per), indent=1, default=float))
    print(json.dumps(per), [(r["component"], round(r["rate_mm_yr"], 2), r["offsets_mm"]) for r in rows])

    # ---- figures: maps
    panels = [(g, k) for g in ("asc", "desc", "vert") if g in maps for k in ("pre",)]
    fig, axs = plt.subplots(len(panels), 1, figsize=(13, 3.6 * len(panels)))
    axs = np.atleast_1d(axs)
    ref = maps["asc"]
    ext = [ref["tr"].c, ref["tr"].c + ref["tr"].a * ref["lon"].shape[1], ref["tr"].f + ref["tr"].e * ref["lon"].shape[0],
           ref["tr"].f]
    for ax, (g, k) in zip(axs, panels):
        vm = maps[g][k] - np.nanmedian(maps[g][k])
        im = ax.imshow(vm, extent=ext, cmap="RdBu", vmin=-25, vmax=25, interpolation="nearest")
        for f, c in FIELDS.items():
            ax.plot(*c, "k+", ms=10)
            ax.annotate(f, c, fontsize=8, xytext=(4, 4), textcoords="offset points")
        ax.set_aspect(1 / np.cos(np.radians(64)))
        lab = {"asc": f"ascending {FRAMES['asc']} LOS", "desc": f"descending {FRAMES['desc']} LOS",
               "vert": "vertical (asc+desc decomposition)"}[g]
        yrs = maps["asc"]["years_pre"] if g == "vert" else maps[g]["years_pre"]
        ax.set_title(f"{lab}, median of summer-to-summer rates, primaries {yrs[0]}-{yrs[-1]} (pre-unrest)", fontsize=9)
        plt.colorbar(im, ax=ax, label="mm/yr (+ towards satellite / up)", shrink=.8)
    fig.tight_layout()
    fig.savefig(IMG / "geothermal_velocity_maps.png", dpi=130)

    # yearly field series
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
    for i, (geom, ls) in enumerate((("asc", "-o"), ("desc", "--s"))):
        q = df[df.geom == geom]
        for f, col in zip(FIELDS, ("tab:red", "tab:orange", "tab:blue", "tab:green")):
            s = q[q.field == f]
            ax[0].plot(s.year.str[:4].astype(int) + 0.5, s.peak_los_mm_yr, ls, color=col, ms=4,
                       label=f"{f} {geom}" if True else None)
    ax[0].axhline(-24, color="k", ls=":", lw=1)
    ax[0].text(2015.2, -23, "Hellisheiði ~24 mm/yr LOS (TSX, Juncu et al. 2015/2017)", fontsize=7)
    ax[0].set_ylabel("peak LOS rate rel. to 4-7 km ring (mm/yr)")
    ax[0].set_xlabel("year of the summer-to-summer interval (start + 0.5)")
    ax[0].legend(fontsize=6, ncol=2)
    ax[0].grid(alpha=.3)
    ax[0].set_ylim(-45, 45)
    ax[0].set_title("Geothermal fields: yearly LiCSAR rates (Svartsengi from 2020 = magmatic)", fontsize=9)
    ax[1].plot(tidx, yu, ".", ms=1, color="k")
    ax[1].set_title("HUSM (Húsmúli) vertical relative to REYK, equipment offsets removed", fontsize=9)
    ax[1].set_ylabel("mm")
    ax[1].grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(IMG / "geothermal_yearly_rates.png", dpi=130)


if __name__ == "__main__":
    main()
