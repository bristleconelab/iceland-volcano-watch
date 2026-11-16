"""Volume lost from the Svartsengi source in each dike event, from all GNSS stations with a dike (Okada) term.

The single-station (SENG) losses cannot separate the deflating source from the opening dike next to it. Here, per
event:
1. Dike geometry from the reviewed IMO earthquakes of the first 24 h (36 h from 00:00 when the onset time is unknown),
   Svartsengi box: centre and strike from the principal axis of the epicentres, length = 2.5-97.5 % extent along
   strike, vertical, top at 0 km for eruptions and 1 km for dike-only events, bottom at 4 km (sensitivity 3 and 5 km).
2. GNSS displacement across the event at every station with data: mean over the 2 days centred 2 days before onset to
   the mean over 2 days centred on the post-event level date of docs/data/event_volume_loss.csv.
3. Linear least squares for: source volume change (point source at the joint-fit position −22.456, 63.872, 4.0 km),
   uniform dike opening, common east/north translation. GNSS errors as in svartsengi_gnss.py.
Loss = −(source volume change). Compared with Parks et al. 2025 Table 1 / EGU26-6017.
    python -I scripts/event_dike_losses.py -> docs/data/event_dike_losses.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import svartsengi_gnss as G  # noqa: E402
from ivwlib import load_events, mogi, to_local_km, window_mean  # noqa: E402
from okada import rect_source  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SRC = dict(lon=-22.456, lat=63.872, depth_km=4.0)


def seismicity():
    d = pd.read_csv(ROOT / "data" / "raw" / "imo" / "seiscomp_reykjanes_hengill_2023_2026.csv")
    d["time"] = pd.to_datetime(d["time"], utc=True, format="ISO8601").dt.tz_localize(None)
    d = d[(d.evaluation_mode == "manual") & (d["type"].isna() | (d["type"] == "earthquake"))]
    b = G.SVARTSENGI_BOX
    return d[(d.longitude >= b[0]) & (d.longitude <= b[2]) & (d.latitude >= b[1]) & (d.latitude <= b[3])]


def dike_geometry(w):
    x, y = to_local_km(w.longitude.values, w.latitude.values, SRC["lon"], SRC["lat"])
    P = np.c_[x, y]
    c = P.mean(axis=0)
    _, _, vt = np.linalg.svd(P - c, full_matrices=False)
    ax = vt[0] if vt[0][1] >= 0 else -vt[0]          # point the axis north-ish
    s = (P - c) @ ax
    lo, hi = np.percentile(s, [2.5, 97.5])
    mid = c + ax * (lo + hi) / 2
    strike = float(np.degrees(np.arctan2(ax[0], ax[1])) % 180)
    return dict(e=float(mid[0]), n=float(mid[1]), strike=strike, length=float(hi - lo), n_eq=len(w))


def displacement(gn, t_pre, t_post):
    rows = []
    for s, df in gn.items():
        a, b = window_mean(df, t_pre, 1), window_mean(df, t_post, 1)
        if a is None or b is None:
            continue
        d = b - a
        rows.append(dict(site=s, lon=df.attrs["lon"], lat=df.attrs["lat"], de=d["e"], dn=d["n"], du=d["u"]))
    return pd.DataFrame(rows)


def solve(obs, geom, top, bottom):
    n = len(obs)
    x, y = to_local_km(obs.lon.values, obs.lat.values, SRC["lon"], SRC["lat"])
    me, mn, mu = mogi(x, y, SRC["depth_km"], 1e6)
    H = bottom - top
    de, dn, du = rect_source(x, y, geom["e"], geom["n"], top + H / 2, geom["length"], H, geom["strike"], 90.0,
                             opening=1.0)  # 1 m opening -> metres
    A = np.zeros((3 * n, 4))
    A[:, 0] = np.r_[me, mn, mu] * 1e3                    # mm per Mm3
    A[:, 1] = np.r_[de, dn, du] * 1e3                    # mm per m of opening
    A[:n, 2], A[n:2 * n, 3] = 1.0, 1.0
    d = G.data_vector(obs)
    sig = G.sigmas(d, n)
    m, *_ = np.linalg.lstsq(A / sig[:, None], d / sig, rcond=None)
    cov = np.linalg.inv((A / sig[:, None]).T @ (A / sig[:, None]))
    r = (A @ m - d) / sig
    dof = max(3 * n - 4, 1)
    scale = max(float(r @ r) / dof, 1.0)
    return dict(loss=-m[0], loss_sd=float(np.sqrt(cov[0, 0] * scale)), opening_m=m[1],
                dike_mm3=m[1] * geom["length"] * H, chi2_dof=float(r @ r) / dof)


def main():
    events = load_events()
    gn, _ = G.load_all(events)
    ev = pd.read_csv(ROOT / "docs" / "data" / "event_volume_loss.csv")
    meta = pd.read_csv(ROOT / "docs" / "data" / "events.csv")
    eq = seismicity()
    rows = []
    for _, e in ev.iterrows():
        m = meta[meta.id == e.id].iloc[0]
        if isinstance(m.start_time_utc, str) and m.start_time_utc:
            t0 = pd.Timestamp(f"{m.start} {m.start_time_utc}")
            t1 = t0 + pd.Timedelta(hours=24)
        else:
            t0 = pd.Timestamp(m.start)
            t1 = t0 + pd.Timedelta(hours=36)
        w = eq[(eq.time >= t0) & (eq.time < t1)]
        if len(w) < 10:
            continue
        geom = dike_geometry(w)
        obs = displacement(gn, (pd.Timestamp(e.start) - pd.Timedelta(days=2)).strftime("%Y-%m-%d"), e.level_date)
        if len(obs) < 3:
            continue
        top = 0.0 if "eruption" in e.kind else 1.0
        res = {b: solve(obs, geom, top, b) for b in (3.0, 4.0, 5.0)}
        r4 = res[4.0]
        rows.append(dict(id=e.id, start=e.start, kind=e.kind, n_sites=len(obs), sites="+".join(obs.site),
                         n_eq=geom["n_eq"], dike_strike=round(geom["strike"], 1), dike_length_km=round(geom["length"], 1),
                         loss_mm3=round(r4["loss"], 2), loss_sd=round(r4["loss_sd"], 2),
                         loss_bottom3=round(res[3.0]["loss"], 2), loss_bottom5=round(res[5.0]["loss"], 2),
                         opening_m=round(r4["opening_m"], 2), dike_mm3=round(r4["dike_mm3"], 1),
                         chi2_dof=round(r4["chi2_dof"], 2), loss_seng_only=e.loss_this_work,
                         published_loss=e.published_loss, published_dike=m.published_dike_mm3))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(ROOT / "docs" / "data" / "event_dike_losses.csv", index=False)


if __name__ == "__main__":
    main()
