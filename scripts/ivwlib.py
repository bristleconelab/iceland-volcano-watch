"""Shared helpers for iceland-volcano-watch: GNSS loading and cleaning, Mogi point source, event list.

GNSS: Nevada Geodetic Laboratory (NGL) daily IGS20 .tenv3 files, read as text. Displacements in mm relative to the
first epoch. NGL does not remove offsets that are missing from its steps file; `clean_series` detects vertical-only
offsets (antenna-height style jumps with no horizontal motion and no volcanic event) and gross outliers.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
NGL = ROOT / "data" / "raw" / "ngl"
EVENTS_CSV = ROOT / "docs" / "data" / "events.csv"

# Mogi point source reported by Parks et al. (2025, EPSL 658:119324, Fig. 6): fixed for their GNSS volume series.
PARKS_SOURCE = dict(lon=-22.454, lat=63.869, depth_km=4.0)
NU = 0.25  # Poisson's ratio used by Parks et al. and most Mogi studies


# ------------------------------------------------------------------------------------------------ events
def load_events() -> pd.DataFrame:
    ev = pd.read_csv(EVENTS_CSV, parse_dates=["start", "end"])
    return ev


# ------------------------------------------------------------------------------------------------ geometry
def to_local_km(lon, lat, lon0: float, lat0: float):
    """Equirectangular projection (km) about (lon0, lat0); error < 0.1 % over the ~30 km used here."""
    lon, lat = np.asarray(lon, float), np.asarray(lat, float)
    x = (lon - lon0) * np.cos(np.radians(lat0)) * 111.32
    y = (lat - lat0) * 110.57
    return x, y


def mogi(x_km, y_km, depth_km: float, dv_m3: float, nu: float = NU):
    """Surface displacement (m) of a Mogi point source at (0, 0, -depth) in an elastic half-space.

    u = (1 - nu) * dV / pi * (x, y, d) / R^3   (e.g. Segall 2010, eq. 7.14 with dV the source volume change)."""
    x, y, d = np.asarray(x_km) * 1e3, np.asarray(y_km) * 1e3, depth_km * 1e3
    r3 = (x ** 2 + y ** 2 + d ** 2) ** 1.5
    c = (1 - nu) * dv_m3 / np.pi
    return c * x / r3, c * y / r3, c * d / r3


# ------------------------------------------------------------------------------------------------ GNSS
def load_tenv3(site: str, ngl_dir: Path = NGL) -> pd.DataFrame:
    df = pd.read_csv(ngl_dir / f"{site}.tenv3", sep=r"\s+")
    t = pd.to_datetime(df["YYMMMDD"], format="%y%b%d")
    e = (df["_e0(m)"] + df["__east(m)"]) * 1e3
    n = (df["____n0(m)"] + df["_north(m)"]) * 1e3
    u = (df["u0(m)"] + df["____up(m)"]) * 1e3
    out = pd.DataFrame({
        "date": t, "e": e - e.iloc[0], "n": n - n.iloc[0], "u": u - u.iloc[0],
        "se": df["sig_e(m)"] * 1e3, "sn": df["sig_n(m)"] * 1e3, "su": df["sig_u(m)"] * 1e3,
    })
    out.attrs["lat"] = float(df["_latitude(deg)"].median())
    lon = float(df["_longitude(deg)"].median())
    out.attrs["lon"] = lon - 360 if lon > 180 else lon
    out.attrs["site"] = site
    return out.set_index("date")


def clean_series(df: pd.DataFrame, event_dates, min_step_mm: float = 150.0, max_h_mm: float = 25.0,
                 win: int = 5, guard_days: int = 3) -> tuple[pd.DataFrame, list[dict]]:
    """Remove vertical-only offsets and gross outliers.

    A vertical offset is accepted when the median of the next `win` days minus the median of the previous `win` days
    exceeds `min_step_mm` in Up, the same horizontal difference is below `max_h_mm`, and no volcanic event lies within
    `guard_days`. Volcanic events (dikes) move stations horizontally by decimetres to metres, so a large vertical jump
    with no horizontal motion and no event is treated as an equipment/metadata offset.
    Returns the cleaned frame and the list of removed offsets."""
    d = df.copy()
    ev = pd.to_datetime(pd.Series(list(event_dates)))
    steps = []
    i = win
    idx = d.index
    while i < len(d) - win:
        before = d.iloc[i - win:i]
        after = d.iloc[i:i + win]
        du = after["u"].median() - before["u"].median()
        dh = np.hypot(after["e"].median() - before["e"].median(), after["n"].median() - before["n"].median())
        near_event = ((ev - idx[i]).abs() <= pd.Timedelta(days=guard_days)).any()
        if abs(du) > min_step_mm and dh < max_h_mm and not near_event:
            # locate the jump precisely: the largest single-day change in the rolling-median series
            seg = d["u"].iloc[i - win:i + win].rolling(3, center=True, min_periods=1).median()
            j = i - win + int(np.argmax(np.abs(np.diff(seg.values)))) + 1
            off = d["u"].iloc[j:j + win].median() - d["u"].iloc[j - win:j].median()
            d.loc[d.index >= idx[j], "u"] -= off
            steps.append({"date": idx[j].date().isoformat(), "offset_up_mm": round(float(off), 1)})
            i = j + win
            continue
        i += 1
    # gross outliers: > 5 robust sigma (and > 30 mm) from an 11-day running median, per component
    keep = np.ones(len(d), bool)
    for c, floor in (("e", 15.0), ("n", 15.0), ("u", 30.0)):
        med = d[c].rolling(11, center=True, min_periods=3).median()
        res = d[c] - med
        mad = 1.4826 * res.abs().rolling(61, center=True, min_periods=10).median()
        keep &= ~((res.abs() > np.maximum(5 * mad, floor)).values)
    return d[keep], steps


def window_mean(df: pd.DataFrame, center, half_days: int = 2) -> pd.Series | None:
    """Mean position over center +/- half_days (None if no data)."""
    c = pd.Timestamp(center)
    w = df[(df.index >= c - pd.Timedelta(days=half_days)) & (df.index <= c + pd.Timedelta(days=half_days))]
    if len(w) == 0:
        return None
    return w[["e", "n", "u"]].mean()


def linear_velocity(df: pd.DataFrame, start, end) -> tuple[pd.Series, pd.Series, int]:
    """Least-squares linear velocity (mm/yr) and 1-sigma formal error (scaled by residual scatter)."""
    w = df[(df.index >= pd.Timestamp(start)) & (df.index <= pd.Timestamp(end))]
    t = (w.index - w.index[0]).days.values / 365.25
    A = np.c_[np.ones_like(t), t]
    v, s = {}, {}
    for c in ("e", "n", "u"):
        coef, res, *_ = np.linalg.lstsq(A, w[c].values, rcond=None)
        r = w[c].values - A @ coef
        cov = np.linalg.inv(A.T @ A) * (r @ r) / max(len(t) - 2, 1)
        v[c], s[c] = coef[1], np.sqrt(cov[1, 1])
    return pd.Series(v), pd.Series(s), len(w)
