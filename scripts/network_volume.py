"""Why is our SENG volume series ~1.2 times IMO's when both use the same Mogi source?

IMO's published volume series is an automated GNSS-only (24 h solutions) inversion with a Mogi source fixed at 4 km
depth at -22.454 E, 63.869 N (Parks et al. 2025, Sect. 3.2). Our headline series converts the SENG vertical alone into
volume. This script applies IMO's model to every open NGL station instead of SENG alone:

1. Reference = mean position over +/-2 days around the post-event level date of the July 2025 eruption (2025-07-22,
   the date used in docs/data/imo_volume_comparison.csv). Target = mean over +/-2 days around each date.
2. Pre-unrest velocities (2015-2019, all three components) are removed where a station has them (SENG, NYLA);
   AUSV, GEVK and VOGC started after the unrest began and have none.
3. Volume change of the fixed source by weighted linear least squares, with errors as in svartsengi_gnss.py, in
   three variants: (a) 3-D + common east/north translation, (b) 3-D without translation, (c) vertical only.
   The same fits are repeated leaving out one station at a time.
4. Per station, the volume implied by its vertical alone (u / G_z), to show which station drives the difference.

The joint GNSS + InSAR source (-22.456, 63.872, 4.0 km) is run as well. Outputs:
    docs/data/network_volume_imo.csv   values at the IMO comparison dates, every variant
    docs/data/network_volume_station.csv  per-station implied volumes at those dates
    docs/data/network_volume_series.csv   daily series, variant (a) and SENG-only, Parks source
Retrospective only; no forecasts.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import svartsengi_gnss as G  # noqa: E402
from ivwlib import mogi, to_local_km, window_mean  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs" / "data"
REF = pd.Timestamp("2025-07-22")
SOURCES = {"parks2025": G.PARKS, "joint": dict(lon=-22.456, lat=63.872, depth_km=4.0)}
VARIANTS = ("3d_translation", "3d", "vertical")


def load_clean():
    out = {}
    for s in G.SITES:
        df = pd.read_csv(G.PROC / f"{s}.csv", index_col=0, parse_dates=True)
        raw = pd.read_csv(ROOT / "data" / "raw" / "ngl" / f"{s}.tenv3", sep=r"\s+", nrows=50)
        lon = float(raw["_longitude(deg)"].median())
        df.attrs.update(site=s, lat=float(raw["_latitude(deg)"].median()), lon=lon - 360 if lon > 180 else lon)
        out[s] = df
    return out


def background():
    v = pd.read_csv(DATA / "gnss_pre_unrest_velocity.csv").set_index("site")
    return {s: v.loc[s, ["ve", "vn", "vu"]].to_numpy(float) for s in v.index}


def displacements(gn, bg, t):
    rows = []
    for s, df in gn.items():
        a, b = window_mean(df, REF), window_mean(df, t)
        if a is None or b is None:
            continue
        d = (b - a).to_numpy(float)
        if s in bg:
            d = d - bg[s] * (pd.Timestamp(t) - REF).days / 365.25
        rows.append(dict(site=s, lon=df.attrs["lon"], lat=df.attrs["lat"], de=d[0], dn=d[1], du=d[2]))
    return pd.DataFrame(rows)


def greens(obs, src):
    x, y = to_local_km(obs.lon.values, obs.lat.values, src["lon"], src["lat"])
    gx, gy, gz = mogi(x, y, src["depth_km"], 1.0)
    return gx * 1e9, gy * 1e9, gz * 1e9   # mm per Mm3


def fit(obs, src, variant):
    gx, gy, gz = greens(obs, src)
    n = len(obs)
    if variant == "vertical":
        d, g = obs.du.to_numpy(), gz
        sig = np.sqrt(G.SIG_U ** 2 + (G.SIG_REL * d) ** 2)
        A = g[:, None]
    else:
        d = np.r_[obs.de, obs.dn, obs.du]
        sig = G.sigmas(d, n)
        A = np.r_[gx, gy, gz][:, None]
        if variant == "3d_translation":
            T = np.zeros((3 * n, 2))
            T[:n, 0] = 1.0
            T[n:2 * n, 1] = 1.0
            A = np.c_[A, T]
    m, *_ = np.linalg.lstsq(A / sig[:, None], d / sig, rcond=None)
    r = (A @ m - d) / sig
    dof = max(len(d) - A.shape[1], 1)
    return float(m[0]), float(r @ r / dof)


def main():
    gn, bg = load_clean(), background()
    imo = pd.read_csv(DATA / "imo_volume_comparison.csv")
    rows, st_rows = [], []
    for _, it in imo.iterrows():
        obs = displacements(gn, bg, it.date)
        for sk, src in SOURCES.items():
            gz = greens(obs, src)[2]
            for _, o in obs.iterrows():
                k = obs.index[obs.site == o.site][0]
                st_rows.append(dict(date=it.date, source=sk, site=o.site, up_mm=round(o.du, 1),
                                    gz_mm_per_mm3=round(gz[k], 2), implied_mm3=round(o.du / gz[k], 2)))
            for v in VARIANTS:
                dv, c2 = fit(obs, src, v)
                rows.append(dict(date=it.date, imo_mm3=it.imo_mm3, imo_lo=it.imo_lo, imo_hi=it.imo_hi, source=sk,
                                 variant=v, left_out="none", n_sites=len(obs), dv_mm3=round(dv, 2),
                                 ratio_to_imo=round(dv / it.imo_mm3, 3), chi2_dof=round(c2, 2)))
                for s in obs.site:
                    sub = obs[obs.site != s].reset_index(drop=True)
                    dv_s, c2_s = fit(sub, src, v)
                    rows.append(dict(date=it.date, imo_mm3=it.imo_mm3, imo_lo=it.imo_lo, imo_hi=it.imo_hi,
                                     source=sk, variant=v, left_out=s, n_sites=len(sub), dv_mm3=round(dv_s, 2),
                                     ratio_to_imo=round(dv_s / it.imo_mm3, 3), chi2_dof=round(c2_s, 2)))
    res = pd.DataFrame(rows)
    res.to_csv(DATA / "network_volume_imo.csv", index=False)
    pd.DataFrame(st_rows).to_csv(DATA / "network_volume_station.csv", index=False)

    ser = []
    for t in pd.date_range(REF + pd.Timedelta(days=7), "2026-10-15", freq="7D"):
        obs = displacements(gn, bg, t)
        if len(obs) < 3 or "SENG" not in set(obs.site):
            continue
        dv, c2 = fit(obs, G.PARKS, "3d_translation")
        gz = greens(obs, G.PARKS)[2]
        k = obs.index[obs.site == "SENG"][0]
        ser.append(dict(date=t.date(), n_sites=len(obs), network_dv_mm3=round(dv, 2), chi2_dof=round(c2, 2),
                        seng_only_mm3=round(obs.du[k] / gz[k], 2)))
    pd.DataFrame(ser).to_csv(DATA / "network_volume_series.csv", index=False)

    main_rows = res[res.left_out == "none"]
    print(main_rows[["date", "imo_mm3", "source", "variant", "n_sites", "dv_mm3", "ratio_to_imo", "chi2_dof"]]
          .to_string(index=False))
    print(pd.DataFrame(st_rows).query("source == 'parks2025'").to_string(index=False))


if __name__ == "__main__":
    main()
