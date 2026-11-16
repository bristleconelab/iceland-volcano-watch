"""Ring harmonics of the radar 'virtual station' volumes, observed and predicted (second-round review, M2).

implied_volume.py fits only mean + cos 2(theta - axis) on the 2.5-4.5 km ring. A mislocated or asymmetric source adds
a first harmonic (cos theta). Here, for the windows with both tracks, the ring values are fitted with
mean + A1 cos(theta - t1) + A2 cos 2(theta - t2), for the data and for synthetic up fields of candidate sources at
their joint-fit geometry and volume, processed the same way (plane beyond 10 km removed, divided by IMO's point-source
response). Model up fields are taken directly (no east/up decomposition), which is adequate for the harmonic shape.
Sources: point (joint fit), lines of point sources L = 4 and 6 km, pressurised prolate CDM (cdm_source.csv).
Rings 2.0-4.0, 2.5-4.5 (main) and 3.0-5.0 km; observed amplitudes get a 1-sigma spread from 500 bootstrap resamples
of the ring blocks (blocks are spatially correlated, so this is a lower bound).
    python -I scripts/implied_volume_harmonics.py -> docs/data/implied_volume_harmonics.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cdm_source as C  # noqa: E402
import elongated_source as E  # noqa: E402
import implied_volume as IV  # noqa: E402
import svartsengi_joint as J  # noqa: E402
from ivwlib import to_local_km  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs" / "data"


def harmonics(az, v):
    th = np.radians(az)
    A = np.c_[np.ones_like(th), np.cos(th), np.sin(th), np.cos(2 * th), np.sin(2 * th)]
    c, *_ = np.linalg.lstsq(A, v, rcond=None)
    rms = float((v - A @ c).std())
    a1, a2 = np.hypot(c[1], c[2]) / c[0], np.hypot(c[3], c[4]) / c[0]
    t1 = np.degrees(np.arctan2(c[2], c[1])) % 360
    t2 = (np.degrees(np.arctan2(c[4], c[3])) / 2) % 180
    return dict(mean=round(float(c[0]), 2), A1_pct=round(100 * a1, 1), A1_toward=round(float(t1), 0),
                A2_pct=round(100 * a2, 1), A2_axis=round(float(t2), 0), rms=round(rms, 2))


def boot(az, v, n=500, seed=3):
    rng = np.random.default_rng(seed)
    a1, a2 = [], []
    for _ in range(n):
        k = rng.integers(0, len(v), len(v))
        h = harmonics(az[k], v[k])
        a1.append(h["A1_pct"])
        a2.append(h["A2_pct"])
    return round(float(np.std(a1)), 1), round(float(np.std(a2)), 1)


def process(df, u):
    far = df.r.values > 10
    A = np.c_[np.ones(far.sum()), df.x.values[far], df.y.values[far]]
    m, *_ = np.linalg.lstsq(A, u[far], rcond=None)
    uu = u - (m[0] + m[1] * df.x.values + m[2] * df.y.values)
    return uu / IV.gz_imo(df.x.values, df.y.values)


def main():
    es = pd.read_csv(DATA / "elongated_source.csv")
    cd = pd.read_csv(DATA / "cdm_source.csv").set_index("window")
    rows = []
    for w in J.CFG["windows"]:
        names = {f: J.chain(f, w["start"], w["end"]) for f in J.CFG["frames"]}
        if any(len(n) < 2 for n in names.values()):
            continue
        key = f"{w['interval']} {w['start']}"
        df = IV.up_field(names)
        models = []
        for L in (0.0, 4.0, 6.0):
            r = es[(es.window == key) & (es.L_km == L)].iloc[0]
            models.append(("point (joint)" if L == 0 else f"line L={L:g}", "line",
                           dict(lon=r.lon, lat=r.lat, depth_km=r.depth_km, L=L, strike=r.strike), r.dv_mm3))
        c = cd.loc[key]
        models.append(("pressurised prolate CDM", "cdm", dict(lon=c.lon, lat=c.lat, depth_km=c.depth_km, ax=c.ax_km,
                                                              ay=c.ay_az_km, strike=c.strike), c.dv_mm3))
        Vm = {}
        for label, kind, p, dv in models:
            x, y = to_local_km(df.lon.values, df.lat.values, p["lon"], p["lat"])
            x, y = np.asarray(x), np.asarray(y)
            gz = (E.line_greens(x, y, p["depth_km"], p["L"], p["strike"])[2] if kind == "line"
                  else C.cdm_greens(x, y, p)[2])
            Vm[label] = process(df, gz * dv)
        for ring_km in ((2.0, 4.0), (2.5, 4.5), (3.0, 5.0)):
            ring = ((df.r > ring_km[0]) & (df.r < ring_km[1])).values
            az = df.az.values[ring]
            h = harmonics(az, df.V.values[ring])
            s1, s2 = boot(az, df.V.values[ring])
            rows.append(dict(window=key, ring_km=f"{ring_km[0]:g}-{ring_km[1]:g}", field="observed", n_blocks=int(ring.sum()),
                             **h, A1_sigma=s1, A2_sigma=s2))
            for label, V in Vm.items():
                rows.append(dict(window=key, ring_km=f"{ring_km[0]:g}-{ring_km[1]:g}", field=f"model: {label}",
                                 n_blocks=int(ring.sum()), **harmonics(az, V[ring])))
        print(pd.DataFrame(rows).tail(5 * 3).to_string(index=False), flush=True)
    pd.DataFrame(rows).to_csv(DATA / "implied_volume_harmonics.csv", index=False)


if __name__ == "__main__":
    main()
