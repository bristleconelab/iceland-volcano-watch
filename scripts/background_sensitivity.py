"""Sensitivity of the station volumes to an unknown background trend at AUSV and GEVK (second review, item 8).

AUSV and GEVK began recording after the unrest started, so no pre-unrest vertical trend can be removed. NYLA's
pre-unrest rate is -5.5 mm/yr and geothermal subsidence occurs on the peninsula. Here a vertical background of
-10, -5, 0, +5 mm/yr is assumed at both stations, and the volume from 22 Jul 2025 to 21 Jul 2026 is recomputed with
IMO's point source and with a 6 km line of point sources (joint-fit I3 geometry): each station alone, all five.
    python -I scripts/background_sensitivity.py -> docs/data/background_sensitivity.csv
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import elongated_source as E  # noqa: E402
import network_volume as N  # noqa: E402
import svartsengi_gnss as G  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def main():
    obs0 = N.displacements(N.load_clean(), N.background(), "2026-07-21")
    yrs = (pd.Timestamp("2026-07-21") - N.REF).days / 365.25
    es = pd.read_csv(ROOT / "docs" / "data" / "elongated_source.csv")
    r6 = es[(es.window == "I3 2025-04-08") & (es.L_km == 6.0)].iloc[0]
    sources = {"IMO point source": dict(G.PARKS, L=0.0, strike=0),
               "line L=6 N60E (I3)": dict(lon=r6.lon, lat=r6.lat, depth_km=r6.depth_km, L=6.0, strike=r6.strike)}
    rows = []
    for v in (-10.0, -5.0, 0.0, 5.0):
        obs = obs0.copy()
        m = obs.site.isin(["AUSV", "GEVK"])
        obs.loc[m, "du"] = obs.loc[m, "du"] - v * yrs      # remove an assumed background of v mm/yr
        for name, p in sources.items():
            r = dict(background_mm_per_yr=v, source=name)
            for s in ("SENG", "AUSV", "GEVK"):
                o = obs[obs.site == s].iloc[0]
                r[f"{s.lower()}_alone"] = round(float(o.du / E.gz_at((o.lon, o.lat), p)), 2)
            r["all5_vertical"], r["chi2_dof"] = [round(x, 2) for x in E.fit_network(obs, p, "vertical")]
            r["seng_over_ausv"] = round(r["seng_alone"] / r["ausv_alone"], 3)
            rows.append(r)
    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "docs" / "data" / "background_sensitivity.csv", index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
