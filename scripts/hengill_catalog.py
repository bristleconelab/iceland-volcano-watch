"""Pilot D: the Zenodo "Hengill seismic catalog" (University of Campania, 2025) against IMO, and Húsmúli re-injection.

 1. Provenance. The Zenodo file hen_conv.dat (5 columns, no header beyond the record text) is compared row by row
    with the COSEISMIQ high-quality double-difference catalogue (Grigoli et al. 2022). The time encoding of the
    Zenodo file is checked against the COSEISMIQ origin times.
 2. Comparison with IMO in the common window (2018-12-01 to 2021-01-31) and box (HENGILL_BOX of Pilot B): event
    matching, magnitude offsets, Mc, b with 95 % bootstrap interval, b-positive, daily rate above Mc.
 3. Húsmúli re-injection area (3 km circle): IMO SIL yearly rates 2006-2025 at M >= 1.0 next to the HUSM GNSS
    vertical (relative to REYK); the COSEISMIQ rate inside and outside the circle in 2018-2021.

Inputs: scripts/fetch_hengill_catalogs.py (data/raw/hengill/), scripts/fetch_imo_quakes.py for
data/raw/imo/{sil,seiscomp}_hengill_1991_2026.csv, NGL HUSM and REYK (Pilot C).
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
import seislib as sl  # noqa: E402
from geothermal_deformation import husm_rate  # noqa: E402
from ivwlib import ROOT  # noqa: E402

IMG, DATA = ROOT / "docs" / "img", ROOT / "docs" / "data"
RAW = ROOT / "data" / "raw"
HENGILL_BOX = (-21.50, 63.98, -21.15, 64.15)   # same box as Pilot B
HUSMULI = (64.055, -21.395)                     # centre of the Húsmúli induced-seismicity cluster
HUSMULI_R_KM = 3.0
T0, T1 = pd.Timestamp("2018-12-01"), pd.Timestamp("2021-02-01")
INJECTION_START = pd.Timestamp("2011-09-01")
R: dict = {}


def inbox(df, b, lat="latitude", lon="longitude"):
    return df[(df[lon] >= b[0]) & (df[lon] <= b[2]) & (df[lat] >= b[1]) & (df[lat] <= b[3])]


def load_imo(name):
    d = pd.read_csv(RAW / "imo" / name)
    d["time"] = pd.to_datetime(d["time"], utc=True, format="ISO8601").dt.tz_localize(None)
    d = d[d["evaluation_mode"] == "manual"]
    if "type" in d:
        d = d[d["type"].isna() | (d["type"] == "earthquake")]
    if "magnitude_type" in d:
        d = d[d["magnitude_type"] == "ML_SIL"]
    return d.dropna(subset=["magnitude"]).sort_values("time").reset_index(drop=True)


# ----------------------------------------------------------------------------------------- 1. provenance
z = pd.DataFrame(np.loadtxt(RAW / "hengill" / "hen_conv.dat"), columns=["t", "lat", "lon", "c4", "c5"])
q = pd.read_csv(RAW / "hengill" / "catalogue_HQDD.csv")
q["time"] = pd.to_datetime(q["origtime"], utc=True, format="ISO8601").dt.tz_localize(None)
# pair each Zenodo row with the nearest COSEISMIQ epicentre (coordinates are unique in both files)
from scipy.spatial import cKDTree  # noqa: E402
# latitude, longitude and depth (depth scaled so 1 km weighs like 1e-3 degree)
dist, idx = cKDTree(np.c_[q.latitude, q.longitude, q.depth * 1e-3]).query(np.c_[z.lat, z.lon, z.c4 * 1e-3])
zs, qs = z, q.iloc[idx].reset_index(drop=True)
same_xy = pd.Series((np.abs(zs.lat - qs.latitude) < 5e-6) & (np.abs(zs.lon - qs.longitude) < 5e-6))
R["provenance_unique_partners"] = int(len(np.unique(idx)))
days = (qs.time - pd.Timestamp("2018-01-01")).dt.total_seconds() / 86400
err_s = (zs.t - (days + 1) / 365) * 365 * 86400   # best-fitting encoding: (days since 2018-01-01 + 1) / 365
ym = qs.time.dt.strftime("%Y-%m")
day_off = err_s.groupby(ym).median()
R["provenance"] = {
    "zenodo_rows": int(len(z)), "coseismiq_hqdd_rows": int(len(q)),
    "rows_with_same_lat_lon": int(same_xy.sum()),
    "col4_equals_hqdd_depth": int((np.abs(zs.c4 - qs.depth) < 1e-3).sum()),
    "col5_equals_hqdd_magnitude": int((np.abs(zs.c5 - qs.magnitude) < 1e-4).sum()),
    "col4_equals_hqdd_magnitude": int((np.abs(zs.c4 - qs.magnitude) < 1e-4).sum()),
    "col5_equals_abs_hqdd_magnitude": int((np.abs(zs.c5 - qs.magnitude.abs()) < 1e-4).sum()),
    "hqdd_negative_magnitudes": int((q.magnitude < 0).sum()),
    "zenodo_negative_magnitudes": int((z.c5 < 0).sum()),
    "time_encoding": "t = (days since 2018-01-01 + 1) / 365, i.e. decimal years from 2018 with 365-day years",
    "months_off_by_one_day": [k for k, v in day_off.items() if abs(v) > 43200],
    "rows_off_by_about_one_day": int(((err_s.abs() - 86400).abs() < 3600).sum()),
    "time_error_s_other_months_p5_p50_p95": [float(np.percentile(err_s[abs(err_s) < 43200], p)) for p in (5, 50, 95)],
    "time_resolution_s": {"t<1": 1e-7 * 365 * 86400, "t>=1": 1e-6 * 365 * 86400},
    "hqdd_time_range": [str(q.time.min()), str(q.time.max())],
    "hqdd_mag_range": [float(q.magnitude.min()), float(q.magnitude.max())],
}

# ----------------------------------------------------------------------------------------- 2. against IMO
cq = inbox(q, HENGILL_BOX).sort_values("time").reset_index(drop=True)
sc = load_imo("seiscomp_hengill_1991_2026.csv")
sil = load_imo("sil_hengill_1991_2026.csv")
win = lambda d: inbox(d[(d.time >= T0) & (d.time < T1)], HENGILL_BOX).reset_index(drop=True)  # noqa: E731
scw, silw = win(sc), win(sil)
ndays = (T1 - T0).days


def match(a, b, tol_s=3.0, tol_km=3.0):
    """For each event in a, the nearest-in-time event of b within tol_s and tol_km."""
    bb = b[["time", "latitude", "longitude", "magnitude"]].rename(
        columns={"latitude": "lat_b", "longitude": "lon_b", "magnitude": "m_b"}).assign(t_b=lambda d: d.time)
    m = pd.merge_asof(a.sort_values("time"), bb.sort_values("time"), on="time", direction="nearest",
                      tolerance=pd.Timedelta(seconds=tol_s)).dropna(subset=["m_b"])
    m = m[sl.haversine_km(m.latitude, m.longitude, m.lat_b, m.lon_b) <= tol_km]
    return m.sort_values("time").drop_duplicates("t_b")


def stats(name, m_time_ordered, rng=0):
    m = np.asarray(m_time_ordered, float)
    maxc = sl.mc_maxc(m, 0.1, 0.2)
    gft, rfit = sl.mc_gft(m)
    mc = float(np.nanmax([maxc, gft]))
    b, sig, n, _ = sl.b_aki(m, mc)
    rg = np.random.default_rng(rng)
    mm = sl.bin_mag(m, 0.1)
    mm = mm[mm >= mc - 0.05 + 1e-9]
    boots = np.array([sl.b_aki(rg.choice(mm, mm.size), mc)[0] for _ in range(1000)])
    bp, bps, nbp = sl.b_positive(m)
    return dict(set=name, n_all=int(m.size), mc_maxc02=maxc, mc_gft=gft, gft_fit_pct=rfit, mc=mc, n_ge_mc=int(n),
                b=float(b), b_shi_bolt=float(sig), b_ci95=[float(np.percentile(boots, 2.5)),
                                                         float(np.percentile(boots, 97.5))],
                b_positive=float(bp), b_positive_sigma=float(bps), rate_per_day_ge_mc=float(n / ndays))


R["window"] = {"start": str(T0.date()), "end_exclusive": str(T1.date()), "days": ndays, "box_WSEN": HENGILL_BOX}
R["counts_in_box"] = {"coseismiq_hqdd": int(len(cq)), "imo_seiscomp": int(len(scw)), "imo_sil": int(len(silw))}
mS = match(scw, cq)
mL = match(silw, cq)
for lab, m_ in (("seiscomp", mS), ("sil", mL)):
    d_ = m_["magnitude"] - m_["m_b"]           # IMO minus COSEISMIQ
    a_, b_ = np.polyfit(m_["m_b"], m_["magnitude"], 1)
    R[f"match_{lab}"] = {"imo_events": int(len(scw if lab == "seiscomp" else silw)),
                         "matched_to_coseismiq": int(len(m_)),
                         "imo_minus_coseismiq_median": float(d_.median()),
                         "iqr": [float(d_.quantile(.25)), float(d_.quantile(.75))],
                         "ols_imo_on_coseismiq": [float(b_), float(a_)]}
R["stats"] = [stats("COSEISMIQ HQDD (ML)", cq.magnitude), stats("IMO SeisComP (ML_SIL)", scw.magnitude),
              stats("IMO SIL", silw.magnitude)]
# share of IMO events above IMO's Mc that COSEISMIQ also has
for lab, w, m_ in (("seiscomp", scw, mS), ("sil", silw, mL)):
    mc = next(s["mc"] for s in R["stats"] if lab in s["set"].lower())
    R[f"match_{lab}"]["share_of_imo_ge_mc_found"] = float((m_["magnitude"] >= mc - 0.05).sum()
                                                          / (w["magnitude"] >= mc - 0.05).sum())

# ----------------------------------------------------------------------------------------- 3. Húsmúli
def in_circle(d, lat="latitude", lon="longitude"):
    return d[sl.haversine_km(d[lat], d[lon], *HUSMULI) <= HUSMULI_R_KM]


MC_LONG = 1.0
sh = in_circle(sil)
sh = sh[(sh.magnitude >= MC_LONG - 0.05) & (sh.time >= "2006-01-01") & (sh.time < "2026-01-01")]
yearly = sh.groupby(sh.time.dt.year).size().reindex(range(2006, 2026), fill_value=0)
R["husmuli_sil_yearly_M>=1"] = {int(k): int(v) for k, v in yearly.items()}
per = {"2006-01..2010-12 (before injection)": ("2006-01-01", "2011-01-01"),
       "2011-01..2011-08": ("2011-01-01", "2011-09-01"),
       "2011-09..2011-12": ("2011-09-01", "2012-01-01"),
       "2012-06..2015-12 (HUSM +2.1 mm/yr)": ("2012-06-01", "2016-01-01"),
       "2017-03..2019-12 (HUSM -5.1 mm/yr)": ("2017-03-01", "2020-01-01"),
       "2020-01..2025-12 (HUSM -9.9 mm/yr)": ("2020-01-01", "2026-01-01")}
R["husmuli_sil_rate_per_year_M>=1"] = {
    k: float(((sh.time >= a) & (sh.time < b)).sum() / ((pd.Timestamp(b) - pd.Timestamp(a)).days / 365.25))
    for k, (a, b) in per.items()}
R["husmuli_sil_max_mag_by_period"] = {k: float(sh[(sh.time >= a) & (sh.time < b)].magnitude.max())
                                      for k, (a, b) in per.items()}
cqh = in_circle(cq)
R["coseismiq_husmuli"] = {"events": int(len(cqh)), "share_of_box": float(len(cqh) / len(cq)),
                          "median_depth_km": float(cqh.depth.median()),
                          "median_depth_rest_of_box_km": float(cq.drop(cqh.index).depth.median())}
rows, husm_per, tidx, husm_u = husm_rate()
R["husm_vertical_rel_REYK"] = husm_per

# ----------------------------------------------------------------------------------------- figures
fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
ax = axs[0]
ax.plot(mS["m_b"], mS["magnitude"], ".", ms=2.5, alpha=.4, color="tab:blue",
        label=f"SeisComP ML_SIL (n = {len(mS)})")
ax.plot(mL["m_b"], mL["magnitude"], ".", ms=2.5, alpha=.4, color="tab:orange", label=f"SIL (n = {len(mL)})")
ax.plot([-1, 4.5], [-1, 4.5], "k:", lw=1)
ax.set_xlabel("COSEISMIQ HQDD magnitude (ML)")
ax.set_ylabel("IMO magnitude, same event")
ax.set_title("Matched events, 3 s and 3 km")
ax.legend(fontsize=8)
ax = axs[1]
for (nm, mags), col in zip((("COSEISMIQ HQDD", cq.magnitude), ("IMO SeisComP", scw.magnitude),
                            ("IMO SIL", silw.magnitude)), ("k", "tab:blue", "tab:orange")):
    bins, nc, cum = sl.fmd(np.asarray(mags), 0.1)
    ax.semilogy(bins, cum, "o", ms=3, color=col, label=f"{nm} (N = {len(mags)})")
    s_ = next(s for s in R["stats"] if s["set"].startswith(nm.split()[0]) and nm.split()[-1] in s["set"])
    ax.axvline(s_["mc"], color=col, ls="--", lw=.8)
ax.set_xlabel("magnitude (each catalogue in its own scale)")
ax.set_ylabel("cumulative number")
ax.set_title("Hengill box, Dec 2018 - Jan 2021")
ax.legend(fontsize=8)
ax = axs[2]
cd = cq.set_index("time").resample("7D").size()
hd = cqh.set_index("time").resample("7D").size()
ax.bar(cd.index, cd.values, width=6, color="0.7", label="COSEISMIQ, whole box")
ax.bar(hd.index, hd.values, width=6, color="k", label=f"COSEISMIQ, Húsmúli {HUSMULI_R_KM:.0f} km circle")
ax.set_ylabel("events per week (all magnitudes)")
ax.set_title("Weekly counts")
ax.legend(fontsize=8)
fig.autofmt_xdate()
fig.tight_layout()
fig.savefig(IMG / "hengill_catalog_comparison.png", dpi=130)
plt.close(fig)

fig, axs = plt.subplots(2, 1, figsize=(11, 6.5), sharex=True)
ax = axs[0]
mo = sh.set_index("time").resample("MS").size()
ax.bar(mo.index, mo.values, width=25, color="k")
ax.set_yscale("symlog", linthresh=10)
ax.set_ylabel(f"SIL events per month, M ≥ {MC_LONG}")
ax.set_title(f"Húsmúli ({HUSMULI_R_KM:.0f} km circle): IMO SIL seismicity and HUSM GNSS")
ax2 = axs[1]
ok = np.isfinite(husm_u)
ax2.plot(tidx[ok], husm_u[ok] - np.nanmedian(husm_u[ok][:60]), ".", ms=1.5, color="tab:green")
ax2.set_ylabel("HUSM vertical rel. REYK (mm)")
for p in husm_per:
    a, b = p["period"].split("..")
    ax2.annotate(f"{p['vu_rel_REYK']:+.1f} mm/yr", (pd.Timestamp(a + "-01") + (pd.Timestamp(b + "-01")
                 - pd.Timestamp(a + "-01")) / 2, ax2.get_ylim()[1] * 0.85 if ax2.get_ylim()[1] > 0 else 0),
                 ha="center", fontsize=8)
for a_ in axs:
    a_.axvline(INJECTION_START, color="tab:red", lw=1)
    a_.axvspan(T0, T1, color="tab:blue", alpha=.12)
axs[0].text(INJECTION_START, axs[0].get_ylim()[1] * 0.6, " re-injection starts (Sep 2011)", color="tab:red",
            fontsize=8)
axs[0].text(T0, axs[0].get_ylim()[1] * 0.6, " COSEISMIQ", color="tab:blue", fontsize=8)
ax2.set_xlim(pd.Timestamp("2006-01-01"), pd.Timestamp("2026-06-01"))
fig.tight_layout()
fig.savefig(IMG / "husmuli_seismicity_gnss.png", dpi=130)
plt.close(fig)

(DATA / "hengill_catalog.json").write_text(json.dumps(R, indent=1, default=float))
print(json.dumps(R, indent=1, default=float))
