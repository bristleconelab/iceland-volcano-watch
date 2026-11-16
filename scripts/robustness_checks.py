"""Small robustness numbers quoted in the README and the manuscript, written to CSV so they can be checked.

1. Endpoint sensitivity: single-station volumes from 22 Jul 2025 to 21 Jul 2026 (IMO's point source, vertical) with
   5-day means (network_volume.py) against one-sided 15-day linear fits (implied_volume.local_position).
2. Residual uplift lobe west-south-west of the source in the Apr-Jul 2025 up field (paper_figures.py Fig. 1): maximum
   and 90th percentile of the residual in the box 22.58-22.46 W, 63.84-63.875 N, and the residual at SENG's block,
   for the point source, the 4 and 6 km lines and the pressurised prolate CDM (each at its joint-fit geometry, volume refitted
   to the up field with offset and plane).
    python -I scripts/robustness_checks.py -> docs/data/robustness_endpoints.csv, docs/data/robustness_residual_lobe.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cdm_source as C  # noqa: E402
import elongated_source as E  # noqa: E402
import implied_volume as IV  # noqa: E402
import network_volume as N  # noqa: E402
import paper_figures as PF  # noqa: E402
import svartsengi_gnss as G  # noqa: E402
from ivwlib import mogi, to_local_km, window_mean  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs" / "data"


def endpoints():
    gc, bg = N.load_clean(), N.background()
    t0, t1 = pd.Timestamp("2025-07-22"), pd.Timestamp("2026-07-21")
    rows = []
    for s in ("SENG", "AUSV", "GEVK"):
        d = gc[s]
        sx, sy = to_local_km(d.attrs["lon"], d.attrs["lat"], G.PARKS["lon"], G.PARKS["lat"])
        g = float(mogi(np.atleast_1d(sx), np.atleast_1d(sy), 4.0, 1e6)[2][0] * 1e3)
        vu = bg[s][2] if s in bg else 0.0
        yrs = (t1 - t0).days / 365.25
        du5 = window_mean(d, t1).u - window_mean(d, t0).u - vu * yrs
        a, b = IV.local_position(d, t0, +1), IV.local_position(d, t1, -1)
        dul = b["u"] - a["u"] - vu * yrs
        rows.append(dict(site=s, gz_mm_per_mm3=round(g, 2), vol_5day_mean=round(du5 / g, 2),
                         vol_one_sided_fit=round(dul / g, 2), fit_sigma=round(np.hypot(a["su"], b["su"]) / g, 2),
                         difference=round((dul - du5) / g, 2)))
    pd.DataFrame(rows).to_csv(DATA / "robustness_endpoints.csv", index=False)
    return pd.DataFrame(rows)


def lobe():
    df = PF.up_field()
    w = ((df.lon < -22.46) & (df.lon > -22.58) & (df.lat > 63.84) & (df.lat < 63.875)).values
    s = ((abs(df.lon + 22.4337) < 0.006) & (abs(df.lat - 63.8822) < 0.006)).values
    cd = pd.read_csv(DATA / "cdm_source.csv").set_index("window").loc["I3 2025-04-08"]
    rows = []
    for label, p in (("point (joint)", PF.geom("I3 2025-04-08", 0.0)), ("line L=4", PF.geom("I3 2025-04-08", 4.0)),
                     ("line L=6", PF.geom("I3 2025-04-08", 6.0)), ("pressurised prolate CDM", None)):
        if p is None:
            q = dict(lon=cd.lon, lat=cd.lat, depth_km=cd.depth_km, ax=cd.ax_km, ay=cd.ay_az_km, strike=cd.strike)
            x, y = to_local_km(df.lon.values, df.lat.values, q["lon"], q["lat"])
            gz = C.cdm_greens(np.asarray(x), np.asarray(y), q)[2]
            A = np.c_[gz, np.ones_like(gz), x, y]
            m, *_ = np.linalg.lstsq(A, df.u.values, rcond=None)
            r = df.u.values - A @ m
        else:
            r, _ = PF.residual(df, p)
        rows.append(dict(model=label, lobe_max_mm=round(float(r[w].max()), 1),
                         lobe_p90_mm=round(float(np.percentile(r[w], 90)), 1),
                         residual_at_seng_mm=round(float(r[s].mean()), 1), rms_mm=round(float(np.sqrt(np.mean(r ** 2))), 2),
                         peak_up_mm=round(float(np.nanpercentile(df.u, 99.5)), 0)))
    pd.DataFrame(rows).to_csv(DATA / "robustness_residual_lobe.csv", index=False)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    print(endpoints().to_string(index=False))
    print(lobe().to_string(index=False))
