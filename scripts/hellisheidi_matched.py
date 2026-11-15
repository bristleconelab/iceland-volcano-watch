"""Hellisheiði subsidence: recompute the LiCSAR rate with the definitions of the published 24-26 mm/yr.

The published numbers are TerraSAR-X track 41 ascending LOS rates (Juncu et al. 2017, JGR 122:692, 2012-2015,
"up to 26 mm/yr"; EGU2015-13455, from Oct 2009, "maximum of ~24 mm/yr"). Both are the *peak* of a StaMPS
small-baseline velocity field (3 m pixels, noise 4.9 mm/yr), *after subtracting a plate-motion model*, relative to a
reference point ~10 km SW of the plant near the GFEL GNSS station (Juncu et al. 2017, Fig. 7 hexagon). The LOS unit
vector is [-0.50, -0.12, 0.86] (E, N, U), incidence ~30.6 deg.

Pilot C measured something else: Sentinel-1 track 16 ascending (incidence 38 deg), 2015-2019 annual pairs, the 5th
percentile of a 3x3-smoothed map within 1.5 km of (-21.40 E, 64.04 N), minus the median of a 4-7 km ring. This script
moves those choices one at a time towards the published definition:

  D0  pilot C as published (centre -21.40/64.04, r 1.5 km, 5th pct, ring 4-7 km)
  D1  centre at the Hellisheiði source of Juncu et al. 2017, Table 3 (joint: 21.371 W, 64.033 N)
  D2  D1 but referenced to the GFEL hexagon area (median within 1 km of 63.963 N, 21.530 W) instead of the ring
  D3  D2 after removing a plane fitted outside the deforming areas (stand-in for the plate-motion model)
  D3min  D3 with the minimum (not the 5th percentile) of the 3x3-smoothed map: the "up to" statistic
  TSX  D3 decomposed (asc + desc) into E and U, re-projected onto the TerraSAR-X LOS (north term ignored)

It also checks the annual pairs against GNSS (HUSM minus REYK, projected into LOS, for the exact pair dates) and
writes yearly D0/D2/D3 series to see whether the rate fell after 2015.
Outputs: docs/data/hellisheidi_matched.json, docs/data/hellisheidi_matched_yearly.csv,
docs/img/hellisheidi_matched.png. Inputs: LiCSAR crops from scripts/licsar_annual.py, NGL HUSM and REYK.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geothermal_deformation import FRAMES, LIC, grid_lonlat, load_frame, yearly_maps  # noqa: E402
from ivwlib import ROOT, load_tenv3, to_local_km, window_mean  # noqa: E402

IMG, DATA = ROOT / "docs" / "img", ROOT / "docs" / "data"
PILOT_C_CENTRE = (-21.40, 64.04)
JUNCU_CENTRE = (-21.371, 64.033)          # Juncu et al. 2017 Table 3, joint inversion, Hellisheiði spheroid
GFEL_REF = (-21.530, 63.963)              # reference hexagon in Juncu et al. 2017 Figs. 7-8 (read off the map)
TSX_LOS = np.array([-0.50, -0.12, 0.86])  # Juncu et al. 2017, section 2.2 (E, N, U)
# areas left out of the plane fit: the two geothermal fields, the eastern-Hengill deep source (Juncu Table 3) and
# the 2017-18 Ölkelduháls uplift (Ducrocq et al. 2021), radius in km
PLANE_MASK = [(JUNCU_CENTRE, 5.0), ((-21.267, 64.107), 4.0), ((-21.25, 64.057), 7.0)]
PLANE_WINDOW = (-21.80, 63.90, -21.10, 64.16)  # lon0, lat0, lon1, lat1
PRE = [2015, 2016, 2017, 2018]
GFEL_R_KM, REF_NMIN = 2.0, 30
MIN_BOWL_COVERAGE = 0.3
DEFS = ("D0", "D1", "D2", "D3", "D3gfel", "D3min")
PUBLISHED = dict(juncu2017_tsx_los_peak=-26.0, egu2015_tsx_los_peak=-24.0, juncu2017_gnss_hh47_up=-23.0)


# --------------------------------------------------------------------------------------------- map statistics
def smooth3(vmap):
    """3x3 mean that ignores NaN; NaN where the input is NaN."""
    ok = np.isfinite(vmap)
    sm = uniform_filter(np.where(ok, vmap, 0.0), 3) / np.maximum(uniform_filter(ok.astype(float), 3), 1e-6)
    sm[~ok] = np.nan
    return sm


def disc(lon, lat, centre, r_km, r0_km=0.0):
    x, y = to_local_km(lon, lat, *centre)
    r = np.hypot(x, y)
    return (r >= r0_km) & (r <= r_km)


def peak(vmap, lon, lat, centre, r_km=1.5, stat="p5"):
    """Most negative value near `centre`: 5th percentile ('p5') or minimum ('min') of the 3x3-smoothed map."""
    v = smooth3(vmap)[disc(lon, lat, centre, r_km)]
    v = v[np.isfinite(v)]
    if v.size < 10:
        return np.nan
    return float(np.percentile(v, 5) if stat == "p5" else v.min())


def ref_value(vmap, lon, lat, centre, r_km, r0_km=0.0, nmin=20):
    v = vmap[disc(lon, lat, centre, r_km, r0_km)]
    v = v[np.isfinite(v)]
    return float(np.median(v)) if v.size >= nmin else np.nan


def remove_plane(vmap, lon, lat, mask=PLANE_MASK, window=PLANE_WINDOW):
    """Fit a + b*x + c*y to the pixels inside `window` and outside every masked disc; return map minus plane."""
    x, y = to_local_km(lon, lat, *JUNCU_CENTRE)
    use = (lon >= window[0]) & (lon <= window[2]) & (lat >= window[1]) & (lat <= window[3]) & np.isfinite(vmap)
    for c, r in mask:
        use &= ~disc(lon, lat, c, r)
    if use.sum() < 100:
        return vmap * np.nan, None
    A = np.c_[np.ones(use.sum()), x[use], y[use]]
    coef, *_ = np.linalg.lstsq(A, vmap[use], rcond=None)
    return vmap - (coef[0] + coef[1] * x + coef[2] * y), coef


def definitions(vmap, lon, lat):
    """All definitions for one LOS (or vertical) map; NaN where the reference or the bowl is poorly covered."""
    out = {}
    out["D0"] = peak(vmap, lon, lat, PILOT_C_CENTRE) - ref_value(vmap, lon, lat, PILOT_C_CENTRE, 7.0, 4.0, 50)
    ring1 = ref_value(vmap, lon, lat, JUNCU_CENTRE, 7.0, 4.0, 50)
    out["D1"] = peak(vmap, lon, lat, JUNCU_CENTRE) - ring1
    gref = ref_value(vmap, lon, lat, GFEL_REF, GFEL_R_KM, nmin=REF_NMIN)
    out["D2"] = peak(vmap, lon, lat, JUNCU_CENTRE) - gref
    dp, _ = remove_plane(vmap, lon, lat)
    out["D3"] = peak(dp, lon, lat, JUNCU_CENTRE)                       # far field (plane residual) = 0
    out["D3gfel"] = out["D3"] - ref_value(dp, lon, lat, GFEL_REF, GFEL_R_KM, nmin=REF_NMIN)
    out["D3min"] = peak(dp, lon, lat, JUNCU_CENTRE, stat="min")
    out["ring_minus_gfel"] = ring1 - gref
    m = disc(lon, lat, JUNCU_CENTRE, 1.5)
    out["bowl_coverage"] = float(np.isfinite(vmap[m]).mean())
    out["gfel_ref_pixels"] = int(np.isfinite(vmap[disc(lon, lat, GFEL_REF, GFEL_R_KM)]).sum())
    return out


def decompose(asc, desc, ga, gd):
    """LOS (asc, desc) -> (east, up), north ignored. ga/gd are dicts of E/U unit-vector maps."""
    det = ga["E"] * gd["U"] - gd["E"] * ga["U"]
    up = (ga["E"] * desc - gd["E"] * asc) / det
    east = (gd["U"] * asc - ga["U"] * desc) / det
    return east, up


def project(east, up, los=TSX_LOS):
    return los[0] * east + los[2] * up


# --------------------------------------------------------------------------------------------- GNSS check
def husm_reyk_offsets_removed():
    """HUSM minus REYK (mm), with the 2017-02-09 HUSM equipment offset (NGL steps file) removed per component."""
    h, r = load_tenv3("HUSM"), load_tenv3("REYK")
    j = h.join(r[["e", "n", "u"]], rsuffix="_r", how="inner")
    d = pd.DataFrame({c: j[c] - j[c + "_r"] for c in "enu"}, index=j.index)
    step = pd.Timestamp("2017-02-09")
    for c in "enu":
        a = d[c][(d.index >= step - pd.Timedelta(days=60)) & (d.index < step)]
        b = d[c][(d.index >= step) & (d.index < step + pd.Timedelta(days=60))]
        if len(a) > 10 and len(b) > 10:
            # offset = jump of local linear fits evaluated at the step date
            ta = (a.index - step).days.values
            tb = (b.index - step).days.values
            off = np.polyval(np.polyfit(tb, b.values, 1), 0) - np.polyval(np.polyfit(ta, a.values, 1), 0)
            d.loc[d.index >= step, c] -= off
    return d


def pixel(vmap, lon, lat, site_lonlat, r_km=1.0, nmin=20):
    """Median of coherent pixels within r_km of a site (NaN if fewer than nmin)."""
    return ref_value(vmap, lon, lat, site_lonlat, r_km, nmin=nmin)


def gnss_vs_insar(geom, frame):
    """Per annual pair: InSAR (HUSM area - REYK area) vs GNSS (HUSM - REYK) projected into this frame's LOS.

    HUSM sits on poorly coherent ground, so a 1 km disc (>= 20 coherent pixels) is used at both sites."""
    pairs, stack, years, tr, shape, meta = load_frame(frame)
    lon, lat = grid_lonlat(tr, shape)
    g = husm_reyk_offsets_removed()
    hs, rk = load_tenv3("HUSM"), load_tenv3("REYK")
    hl, rl = (hs.attrs["lon"], hs.attrs["lat"]), (rk.attrs["lon"], rk.attrs["lat"])
    e, n, u = (pixel(meta[k], lon, lat, hl, 0.3, 1) for k in "ENU")
    rows = []
    for p, rate in zip(pairs, stack):
        t0, t1 = pd.Timestamp(p[:8]), pd.Timestamp(p[9:])
        a, b = window_mean(g, t0, 3), window_mean(g, t1, 3)
        if a is None or b is None:
            continue
        dd = (b - a) / ((t1 - t0).days / 365.25)
        gl = float(e * dd["e"] + n * dd["n"] + u * dd["u"])
        il = pixel(rate, lon, lat, hl) - pixel(rate, lon, lat, rl)
        rows.append(dict(geom=geom, pair=p, gnss_los_mm_yr=round(gl, 2), insar_los_mm_yr=round(il, 2),
                         diff=round(il - gl, 2)))
    return rows


CLOSURE_DIR = ROOT / "data" / "raw" / "licsar_closure" / "016A_02504_162118"
CLOSURE_SETS = {  # long same-summer pair: chain of 6/12-day pairs over the same dates (public LiCSAR products)
    # 2015: 20150603_20150907 has an empty unw file and 20150615_20150627 too, so the chain starts on 27 Jun
    "20150627_20150907": ["20150627_20150709", "20150709_20150721", "20150721_20150802", "20150802_20150814",
                          "20150814_20150826", "20150826_20150907"],
    "20250604_20250926": ["20250604_20250610", "20250610_20250622", "20250622_20250704", "20250704_20250716",
                          "20250716_20250728", "20250728_20250809", "20250809_20250821", "20250821_20250902",
                          "20250902_20250914", "20250914_20250926"],
}


def read_masked(pair, d=CLOSURE_DIR):
    import rasterio
    with rasterio.open(d / "los_mm" / f"{pair}.tif") as ds:
        a, tr = ds.read(1), ds.transform
    with rasterio.open(d / "cc" / f"{pair}.tif") as ds:
        cc = ds.read(1)
    if not np.isfinite(cc).any():  # LiCSAR serves an empty coherence file for some pairs: keep unw unmasked
        return a, tr
    return np.where(cc >= 0.25, a, np.nan), tr


def _cc_empty(pair, d=CLOSURE_DIR):
    import rasterio
    with rasterio.open(d / "cc" / f"{pair}.tif") as ds:
        return not np.isfinite(ds.read(1)).any()


def closure_test():
    """Phase-bias (fading) check: sum of a chain of short pairs minus the direct long pair, same two dates.

    Values are displacements in mm (LOS, + towards the satellite). Only pixels coherent in every pair are used.
    Reported at the Hellisheiði bowl (median within 1.5 km of the Juncu centre) relative to the far field
    (median of the plane-fit area), together with the bowl signal of the long pair itself."""
    out = {}
    for long, chain in CLOSURE_SETS.items():
        if not (CLOSURE_DIR / "los_mm" / f"{long}.tif").exists():
            continue
        L, tr = read_masked(long)
        S = np.zeros_like(L)
        for p in chain:
            S = S + read_masked(p)[0]
        lon, lat = grid_lonlat(tr, L.shape)
        ok = np.isfinite(L) & np.isfinite(S)
        far = (lon >= PLANE_WINDOW[0]) & (lon <= PLANE_WINDOW[2]) & (lat >= PLANE_WINDOW[1]) & (lat <= PLANE_WINDOW[3])
        for c, r in PLANE_MASK:
            far &= ~disc(lon, lat, c, r)
        bowl = disc(lon, lat, JUNCU_CENTRE, 1.5)
        C = np.where(ok, S - L, np.nan)
        Lm = np.where(ok, L, np.nan)
        days = (pd.Timestamp(long[9:]) - pd.Timestamp(long[:8])).days
        cb, cf = np.nanmedian(C[bowl]), np.nanmedian(C[far])
        lb, lf = np.nanmedian(Lm[bowl]), np.nanmedian(Lm[far])
        out[long] = dict(days=days, n_chain=len(chain), chain_without_cc=[p for p in chain if _cc_empty(p)], bowl_pixels=int(np.isfinite(C[bowl]).sum()),
                         far_pixels=int(np.isfinite(C[far]).sum()),
                         closure_bowl_minus_far_mm=round(float(cb - cf), 2),
                         closure_far_median_mm=round(float(cf), 2),
                         long_pair_bowl_minus_far_mm=round(float(lb - lf), 2),
                         closure_iqr_far_mm=round(float(np.subtract(*np.nanpercentile(C[far], [75, 25]))), 2))
    return out


def ngl_vertical_rate(site, start, end):
    """Vertical rate (mm/yr) with an annual term, NGL IGS20."""
    d = load_tenv3(site)
    w = d[(d.index >= start) & (d.index <= end)]
    t = (w.index - pd.Timestamp(start)).days.values / 365.25
    A = np.c_[np.ones_like(t), t, np.sin(2 * np.pi * t), np.cos(2 * np.pi * t)]
    m, *_ = np.linalg.lstsq(A, w["u"].values, rcond=None)
    return float(m[1]), int(len(w))


def med(s):
    s = s.dropna()
    return round(float(s.median()), 1) if len(s) else None


# --------------------------------------------------------------------------------------------- main
def main():
    res = dict(published=PUBLISHED, centre_pilot_c=PILOT_C_CENTRE, centre_juncu=JUNCU_CENTRE, ref_gfel=GFEL_REF,
               gfel_ref_radius_km=GFEL_R_KM, min_bowl_coverage=MIN_BOWL_COVERAGE)
    rows, ym_all, geo = [], {}, {}
    for geom, frame in FRAMES.items():
        if not (LIC / frame / "los_mm").exists():
            raise SystemExit(f"missing {frame}: run scripts/licsar_annual.py first")
        pairs, stack, years, tr, shape, meta = load_frame(frame)
        lon, lat = grid_lonlat(tr, shape)
        ym_all[geom], geo[geom] = yearly_maps(stack, years), dict(meta=meta, lon=lon, lat=lat)
        for p, vm in zip(pairs, stack):
            d = definitions(vm, lon, lat)
            rows.append(dict(geom=geom, pair=p, year=int(p[:4]), **{k: round(v, 2) for k, v in d.items()}))
        res[f"los_unit_{geom}_ENU"] = [round(float(pixel(meta[c], lon, lat, JUNCU_CENTRE, 0.3, 1)), 3) for c in "ENU"]
    pdf = pd.DataFrame(rows)
    pdf["usable"] = pdf.bowl_coverage >= MIN_BOWL_COVERAGE
    pdf.to_csv(DATA / "hellisheidi_matched_pairs.csv", index=False)
    use = pdf[pdf.usable]

    # per-pair medians: pre-unrest primaries (2015-2018) and all years
    summ = {}
    for geom in ("asc", "desc"):
        q = use[use.geom == geom]
        summ[geom] = {}
        for lab, s in (("pre_2015_2019", q[q.year.isin(PRE)]), ("all_years", q)):
            summ[geom][lab] = {c: med(s[c]) for c in DEFS + ("ring_minus_gfel",)}
            summ[geom][lab] |= {f"n_{c}": int(s[c].notna().sum()) for c in DEFS}
            summ[geom][lab]["n_pairs_listed"] = int((pdf.geom == geom).sum() if lab == "all_years"
                                                    else ((pdf.geom == geom) & pdf.year.isin(PRE)).sum())
    res["los_pair_medians"] = summ

    # yearly medians of usable pairs (asc) for the time question
    yr = use.groupby(["geom", "year"])[list(DEFS)].median().round(1).reset_index()
    yr.to_csv(DATA / "hellisheidi_matched_yearly.csv", index=False)
    q = yr[yr.geom == "asc"].copy()
    q["t"] = q.year + 0.5
    trend = {}
    for c in ("D0", "D2", "D3"):
        ok = np.isfinite(q[c])
        b, _ = np.polyfit(q.t[ok], q[c][ok], 1)
        trend[c] = dict(slope_mm_yr_per_yr=round(float(b), 2), n_years=int(ok.sum()),
                        median_2015_2019=med(q[c][q.year < 2020]), median_2022_2026=med(q[c][q.year > 2020]))
    res["asc_time_change"] = trend

    # asc+desc -> east/up for the pre-unrest median maps, then onto the TSX LOS
    lon, lat = geo["asc"]["lon"], geo["asc"]["lat"]
    A = np.nanmedian([ym_all["asc"][y] for y in PRE if y in ym_all["asc"]], axis=0)
    D = np.nanmedian([ym_all["desc"][y] for y in PRE if y in ym_all["desc"]], axis=0)
    east, up = decompose(A, D, geo["asc"]["meta"], geo["desc"]["meta"])
    tsx = project(east, up)
    res["pre_unrest_median_maps"] = {name: {k: round(v, 1) for k, v in definitions(m, lon, lat).items()}
                                     for name, m in (("asc_los", A), ("desc_los", D), ("vertical", up),
                                                     ("east", east), ("tsx_los_equivalent", tsx))}
    tsx3, _ = remove_plane(tsx, lon, lat)
    sm = smooth3(tsx3)
    m4 = disc(lon, lat, JUNCU_CENTRE, 4.0) & np.isfinite(sm)
    i = int(np.nanargmin(np.where(m4, sm, np.nan)))
    res["tsx_equiv_minimum_location"] = dict(lon=round(float(lon.flat[i]), 3), lat=round(float(lat.flat[i]), 3),
                                             km_from_juncu_centre=round(float(np.hypot(*to_local_km(
                                                 lon.flat[i], lat.flat[i], *JUNCU_CENTRE))), 2))
    # absolute vertical, tied to REYK (NGL IGS20 vertical over the same years); no plane removed
    vr, nr = ngl_vertical_rate("REYK", "2015-06-01", "2019-09-30")
    rk = load_tenv3("REYK")
    up_reyk = pixel(up, lon, lat, (rk.attrs["lon"], rk.attrs["lat"]))
    res["vertical_tied_to_reyk"] = dict(
        reyk_ngl_up_mm_yr=round(vr, 2), reyk_days=nr, insar_up_at_reyk=round(up_reyk, 2),
        bowl_p5_abs=round(peak(up, lon, lat, JUNCU_CENTRE) - up_reyk + vr, 1),
        bowl_min_abs=round(peak(up, lon, lat, JUNCU_CENTRE, stat="min") - up_reyk + vr, 1))

    # GNSS check of the annual pairs
    g = []
    for geom, frame in FRAMES.items():
        g += gnss_vs_insar(geom, frame)
    gdf = pd.DataFrame(g)
    gdf.to_csv(DATA / "hellisheidi_gnss_check.csv", index=False)
    res["gnss_check_husm_minus_reyk"] = {}
    for geom, s in gdf.dropna(subset=["diff"]).groupby("geom"):
        res["gnss_check_husm_minus_reyk"][geom] = dict(
            n_pairs=int(len(s)), median_diff=round(float(s["diff"].median()), 2),
            mean_diff=round(float(s["diff"].mean()), 2),
            sem_diff=round(float(s["diff"].std(ddof=1) / np.sqrt(len(s))), 2) if len(s) > 1 else None,
            rms_diff=round(float(np.sqrt((s["diff"] ** 2).mean())), 2))
    res["closure_test"] = closure_test()
    (DATA / "hellisheidi_matched.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(json.dumps(res, indent=1, ensure_ascii=False))
    print(pdf[pdf.geom == "asc"].to_string(index=False))

    # figure: TSX-equivalent map (plane removed) and yearly series
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    win = (lon >= -21.75) & (lon <= -21.10) & (lat >= 63.90) & (lat <= 64.16)
    r_, c_ = np.where(win)
    sl = (slice(r_.min(), r_.max() + 1), slice(c_.min(), c_.max() + 1))
    ext = [lon[sl].min(), lon[sl].max(), lat[sl].min(), lat[sl].max()]
    im = ax[0].imshow(tsx3[sl], extent=ext, cmap="RdBu", vmin=-25, vmax=25, interpolation="nearest")
    ax[0].plot(*JUNCU_CENTRE, "k+", ms=12, label="Juncu et al. 2017 source")
    ax[0].plot(*PILOT_C_CENTRE, "kx", ms=8, label="pilot C centre")
    ax[0].plot(*GFEL_REF, "kh", ms=8, label="GFEL reference (Juncu Fig. 7)")
    th = np.linspace(0, 2 * np.pi, 100)
    for (x0, y0), r in PLANE_MASK:
        ax[0].plot(x0 + r * np.cos(th) / (111.32 * np.cos(np.radians(y0))), y0 + r * np.sin(th) / 110.57, "k:", lw=.6)
    ax[0].set_aspect(1 / np.cos(np.radians(64)))
    ax[0].set_title("S1 asc+desc 2015-2019 median, projected on the TerraSAR-X LOS,\n"
                    "plane removed outside the dotted circles (far field = 0)", fontsize=9)
    ax[0].legend(fontsize=7, loc="lower right")
    plt.colorbar(im, ax=ax[0], label="mm/yr (+ towards satellite)", shrink=.85)
    labels = {"D0": "D0 pilot C (ring)", "D2": "D2 GFEL reference", "D3": "D3 plane removed"}
    for c, st in (("D0", "-o"), ("D2", "-s"), ("D3", "-^")):
        ax[1].plot(q.t, q[c], st, ms=4, label=labels[c])
    ax[1].axhspan(-26, -24, color="grey", alpha=.3, label="published TSX LOS 24-26 (2009/2012-2015)")
    ax[1].set_ylabel("S1 ascending LOS rate at the bowl (mm/yr)")
    ax[1].set_xlabel("summer-to-summer interval (start + 0.5 yr); median of usable pairs")
    ax[1].grid(alpha=.3)
    ax[1].legend(fontsize=7)
    ax[1].set_title("Hellisheiði, yearly LiCSAR rate under three definitions", fontsize=9)
    fig.tight_layout()
    fig.savefig(IMG / "hellisheidi_matched.png", dpi=130)


if __name__ == "__main__":
    main()
