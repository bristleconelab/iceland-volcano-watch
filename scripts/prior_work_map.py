"""Prior-work and coverage map: study areas of published work, open data used here, and what this project measures.

Background: LiCSAR (Copernicus DEM-based) heights cropped by licsar_annual.py, shown as hillshade; sea is masked.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np
import rasterio

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ivwlib import PARKS_SOURCE, ROOT, load_tenv3  # noqa: E402

HGT = ROOT / "data" / "raw" / "licsar" / "016A_02504_162118" / "meta" / "hgt.tif"

STUDIES = [  # (label, bbox W S E N, colour) - approximate footprints of the published analyses
    ("Parks et al. 2025 EPSL: GNSS+InSAR Mogi/sill, 9 dikes 2023-24", (-22.62, 63.78, -22.25, 63.97), "tab:red"),
    ("Flóvenz et al. 2022 Nat. Geosci.: Svartsengi 2020-21 cycles, poroelastic", (-22.55, 63.83, -22.35, 63.93), "tab:purple"),
    ("Juncu et al. 2017/2018: Hengill TSX + GNSS, 2012-2015", (-21.50, 63.98, -21.15, 64.15), "tab:green"),
    ("Reykjanes field S1 2015-2017 (~16 mm/yr LOS)", (-22.74, 63.80, -22.62, 63.86), "tab:olive"),
]


def main():
    with rasterio.open(HGT) as ds:
        h = ds.read(1)
        b = ds.bounds
    h = np.where(h <= 3, np.nan, h)  # sea and coastal noise
    gy, gx = np.gradient(np.nan_to_num(h, nan=0.0))
    shade = np.clip(0.6 + 0.02 * (gx - gy), 0, 1)
    shade[~np.isfinite(h)] = np.nan
    fig, ax = plt.subplots(figsize=(13, 6.5))
    ext = [b.left, b.right, b.bottom, b.top]
    ax.set_facecolor("#cfe3f3")
    ax.imshow(shade, extent=ext, cmap="gray", vmin=0, vmax=1, interpolation="bilinear")
    for lab, bb, col in STUDIES:
        ax.add_patch(Rectangle((bb[0], bb[1]), bb[2] - bb[0], bb[3] - bb[1], fill=False, ec=col, lw=1.8, label=lab))
    for s in ["SENG", "AUSV", "GEVK", "VOGC", "NYLA", "HUSM", "REYK"]:
        d = load_tenv3(s)
        ax.plot(d.attrs["lon"], d.attrs["lat"], "k^", ms=7)
        ax.annotate(f"{s} {d.index.min().year}-{d.index.max().year}", (d.attrs["lon"], d.attrs["lat"]), fontsize=7,
                    xytext=(4, 3), textcoords="offset points")
    ax.plot(PARKS_SOURCE["lon"], PARKS_SOURCE["lat"], "r*", ms=14, label="Svartsengi Mogi source (Parks 2025; 4 km)")
    ax.plot(-22.446, 63.872, "y*", mec="k", ms=12, label="Mogi source, this work (GNSS, 3.75 km)")
    for nm, (lo, la) in {"Hellisheiði": (-21.40, 64.04), "Nesjavellir": (-21.26, 64.11), "Svartsengi PP": (-22.43, 63.88),
                         "Reykjanes PP": (-22.68, 63.83)}.items():
        ax.plot(lo, la, "s", color="orange", mec="k", ms=7)
        ax.annotate(nm, (lo, la), fontsize=8, xytext=(4, -10), textcoords="offset points", color="darkorange")
    ax.text(-22.74, 63.805, "LiCSAR 016A_02504_162118 (asc) and 155D_02484_191814 (desc) both cover this whole window",
            fontsize=8)
    ax.set_xlim(b.left, b.right)
    ax.set_ylim(b.bottom, b.top)
    ax.set_aspect(1 / np.cos(np.radians(64)))
    ax.legend(fontsize=7, loc="lower right", framealpha=.9)
    ax.set_title("Reykjanes Peninsula: prior work footprints, open data used here (NGL GNSS ▲, LiCSAR frames), "
                 "geothermal power plants ■", fontsize=10)
    fig.tight_layout()
    fig.savefig(ROOT / "docs" / "img" / "prior_work_map.png", dpi=130)


if __name__ == "__main__":
    main()
