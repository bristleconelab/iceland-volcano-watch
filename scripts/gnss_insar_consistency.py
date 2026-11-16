"""Do GNSS and InSAR agree at the stations? A check that the station dependence of the volume is not a GNSS artefact.

For every radar window and track, the GNSS displacement over the exact chain dates (1-day means at both ends) is
projected onto the line of sight at the station, and compared with the median chained LOS of the pixels within 350 m.
The InSAR datum is arbitrary per chain, so the comparison uses differences between stations (SENG - GEVK,
SENG - AUSV), where the datum cancels. Absolute values are written too, for reference only.
    python -I scripts/gnss_insar_consistency.py -> docs/data/gnss_insar_consistency.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

sys.path.insert(0, str(Path(__file__).resolve().parent))
import network_volume as N  # noqa: E402
import svartsengi_joint as J  # noqa: E402
from ivwlib import window_mean  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SITES = ("SENG", "AUSV", "GEVK")


def near(arr, tr, lon, lat, r_m=350.0):
    rows, cols = np.mgrid[0:arr.shape[0], 0:arr.shape[1]]
    xs, ys = rasterio.transform.xy(tr, rows.ravel(), cols.ravel())
    xs, ys = np.asarray(xs).reshape(arr.shape), np.asarray(ys).reshape(arr.shape)
    d = np.hypot((xs - lon) * 111320 * np.cos(np.radians(lat)), (ys - lat) * 111130)
    v = arr[(d <= r_m) & np.isfinite(arr)]
    return (float(np.median(v)), int(v.size)) if v.size else (np.nan, 0)


def main():
    gc = N.load_clean()
    rows = []
    for w in J.CFG["windows"]:
        for f in J.CFG["frames"]:
            names = J.chain(f, w["start"], w["end"])
            if len(names) < 2:
                continue
            los, tr = J.chained_los(f, names)
            unit = {k: J.read(J.RAW / f / "meta" / f"{k}.tif")[0] for k in "ENU"}
            t0, t1 = J.pdate(names[0].split("_")[0]), J.pdate(names[-1].split("_")[1])
            v = {}
            for s in SITES:
                a, b = window_mean(gc[s], t0, 1), window_mean(gc[s], t1, 1)
                if a is None or b is None:
                    continue
                lo, la = gc[s].attrs["lon"], gc[s].attrs["lat"]
                u = {k: near(unit[k], tr, lo, la)[0] for k in "ENU"}
                g = b - a
                ins, npx = near(los, tr, lo, la)
                v[s] = (g.e * u["E"] + g.n * u["N"] + g.u * u["U"], ins, npx)
            for o in ("GEVK", "AUSV"):
                if "SENG" in v and o in v:
                    gd, idf = v["SENG"][0] - v[o][0], v["SENG"][1] - v[o][1]
                    rows.append(dict(window=f"{w['interval']} {w['start']}", frame=f, start=t0.date(), end=t1.date(),
                                     pair=f"SENG-{o}", gnss_los_mm=round(gd, 1), insar_los_mm=round(idf, 1),
                                     ratio_insar_over_gnss=round(idf / gd, 3), seng_gnss_abs=round(v["SENG"][0], 1),
                                     seng_insar_abs=round(v["SENG"][1], 1), n_pixels=min(v["SENG"][2], v[o][2])))
    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "docs" / "data" / "gnss_insar_consistency.csv", index=False)
    print(df.to_string(index=False))
    print("median ratio", df.ratio_insar_over_gnss.median(), "range", df.ratio_insar_over_gnss.min(),
          df.ratio_insar_over_gnss.max())


if __name__ == "__main__":
    main()
