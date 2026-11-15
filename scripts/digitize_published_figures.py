"""Read the Hellisheiði subsidence bowl off two published InSAR velocity maps, with the same statistics as our LiCSAR maps.

Why: the published "up to 24-26 mm/yr LOS" (TerraSAR-X, 2009/2012-2015) is a single extreme value after a
plate-motion correction. The papers do not publish the velocity grids, so the maps are digitised from the figures:
  * Juncu et al. 2017, JGR 122:692-709, doi:10.1002/2016JB013626, Fig. 7 (p. 700): TerraSAR-X T41 ascending mean LOS
    velocity 2012-2015, left panel as observed, right panel after subtracting the plate model; reference = hexagon.
  * Ducrocq et al. 2021, Front. Earth Sci. 9:725109, doi:10.3389/feart.2021.725109, Fig. 5 (p. 7): mean LOS
    velocity 2015-2017 for TerraSAR-X T41 and Sentinel-1 T16 ascending (the LiCSAR track we use); reference = the
    black box near GFEL. No plate correction.

Method: each panel is rendered at 300 dpi with poppler's pdftoppm, colours are inverted to mm/yr by hue against the
figure's own colour bar (hillshading changes brightness, not hue), line art is suppressed with a 3x3 median, and the
panel is georeferenced with an affine transform fitted to three GNSS station markers (HUSM, OLKE, HVER; NGL positions).
Statistics: p5 and minimum within 1.5 km of the Juncu et al. Hellisheiði source, relative to (a) the paper's own
reference and (b) the median of a 4-7 km ring (pilot C's reference).
Limits: hue inversion saturates where the palette darkens (Ducrocq bar below about -19 mm/yr; Juncu bar below about
-24 mm/yr), so minima beyond those values are lower bounds in magnitude. Output: docs/data/published_figures_digitized.json.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from matplotlib.colors import rgb_to_hsv
from PIL import Image
from scipy.ndimage import generic_filter

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ivwlib import ROOT, to_local_km  # noqa: E402

PAPERS = ROOT / "data" / "raw" / "papers"
OUT = ROOT / "docs" / "data" / "published_figures_digitized.json"
JUNCU_CENTRE = (-21.371, 64.033)
STATIONS = {"HUSM": (-21.417253, 64.066649), "OLKE": (-21.219893, 64.063126), "HVER": (-21.184810, 64.017153)}

# colour bars: crop box (300 dpi page pixels), then pixel -> value along the bar
BARS = {
    "juncu": dict(pdf="juncu2017/juncu2017.pdf", page=9, box=(272, 2768, 545, 2850), axis="x", rows=(52, 60),
                  span=(16, 257), value=lambda x: (x - 135.5) / 5.0, vmax=8.0),
    "ducrocq": dict(pdf="ducrocq2021/ducrocq2021.pdf", page=7, box=(2104, 649, 2242, 1494), axis="y",
                    cols=(20, 30), span=(23, 829), value=lambda y: 23 - (y - 21) / ((830 - 21) / 46), vmax=23.0),
}
# panels: crop box (300 dpi), station marker centres and reference (in crop pixels)
PANELS = {
    "juncu2017_fig7_tsx_2012_2015_observed": dict(
        bar="juncu", box=(259, 2230, 900, 2809), markers={"HUSM": (270.5, 220.0), "OLKE": (444.5, 227.5),
                                                          "HVER": (476.0, 320.0)},
        ref=("ring_px", (178, 438), 10, 25)),
    "juncu2017_fig7_tsx_2012_2015_plate_corrected": dict(
        bar="juncu", box=(1663, 2230, 2318, 2809), markers={"HUSM": (287.0, 220.0), "OLKE": (461.0, 227.5),
                                                            "HVER": (492.5, 320.0)},
        ref=("ring_px", (195, 437), 10, 25)),
    "ducrocq2021_fig5_tsx_2015_2017": dict(
        bar="ducrocq", box=(300, 289, 739, 666), markers={"HUSM": (110.8, 134.4), "OLKE": (235.4, 139.5),
                                                          "HVER": (258.7, 206.3)},
        ref=("box_px", (34, 249, 76, 288))),
    "ducrocq2021_fig5_s1_t16_2015_2017": dict(
        bar="ducrocq", box=(300, 682, 739, 1059), markers={"HUSM": (111.6, 135.1), "OLKE": (237.3, 141.2),
                                                           "HVER": (259.6, 207.3)},
        ref=("box_px", (34, 249, 76, 288))),
}


def render(pdf: Path, page: int, outdir: Path) -> Image.Image:
    exe = shutil.which("pdftoppm")
    if exe is None:
        raise SystemExit("pdftoppm (poppler) is needed")
    stem = outdir / f"{pdf.stem}_{page}"
    subprocess.run([exe, "-r", "300", "-f", str(page), "-l", str(page), "-png", str(pdf), str(stem)], check=True)
    return Image.open(next(outdir.glob(f"{stem.name}*.png"))).convert("RGB")


def bar_lut(img: Image.Image, spec: dict):
    """Hue and value along a colour bar; drop unsaturated/dark entries and values above vmax."""
    a = np.asarray(img.crop(spec["box"])).astype(float) / 255
    s0, s1 = spec["span"]
    pos = np.arange(s0, s1)
    if spec["axis"] == "x":
        rgb = a[spec["rows"][0]:spec["rows"][1], s0:s1].mean(0)
    else:
        rgb = a[s0:s1, spec["cols"][0]:spec["cols"][1]].mean(1)
    hsv = rgb_to_hsv(rgb)
    val = np.array([spec["value"](p) for p in pos])
    keep = (hsv[:, 1] > 0.35) & (hsv[:, 2] > 0.3) & (val <= spec["vmax"])
    return hsv[keep, 0], val[keep]


def hue_to_value(rgb, lut_h, lut_v, smin=0.35, vmin=0.25):
    """Nearest-hue lookup (hue is circular); NaN for grey, dark or white pixels."""
    hsv = rgb_to_hsv(rgb)
    d = np.abs(hsv[..., 0][..., None] - lut_h)
    d = np.minimum(d, 1 - d)
    v = lut_v[np.argmin(d, axis=-1)]
    return np.where((hsv[..., 1] > smin) & (hsv[..., 2] > vmin), v, np.nan)


def fit_affine(px_xy: dict, lonlat: dict):
    """Affine map from local km (about JUNCU_CENTRE) to pixel x, y, fitted exactly to the station markers."""
    names = sorted(px_xy)
    km = np.array([to_local_km(*lonlat[n], *JUNCU_CENTRE) for n in names])
    A = np.c_[km, np.ones(len(names))]
    px = np.array([px_xy[n] for n in names])
    cx, *_ = np.linalg.lstsq(A, px[:, 0], rcond=None)
    cy, *_ = np.linalg.lstsq(A, px[:, 1], rcond=None)
    M = np.array([cx[:2], cy[:2]])           # px per km
    return M, np.array([cx[2], cy[2]])       # pixel of the Juncu centre


def pixel_km(shape, M, c0):
    """km coordinates (x east, y north) of every pixel, about JUNCU_CENTRE."""
    yy, xx = np.mgrid[:shape[0], :shape[1]]
    inv = np.linalg.inv(M)
    d = np.stack([xx - c0[0], yy - c0[1]], -1)
    km = d @ inv.T
    return km[..., 0], km[..., 1]


def panel_stats(vmap, x_km, y_km, ref):
    r = np.hypot(x_km, y_km)
    inner = vmap[r <= 1.5]
    inner = inner[np.isfinite(inner)]
    ring = np.nanmedian(vmap[(r >= 4) & (r <= 7)])
    yy, xx = np.mgrid[:vmap.shape[0], :vmap.shape[1]]
    if ref[0] == "ring_px":
        (rx, ry), r0, r1 = ref[1], ref[2], ref[3]
        rr = np.hypot(xx - rx, yy - ry)
        refv = np.nanmedian(vmap[(rr > r0) & (rr < r1)])
    else:
        x0, y0, x1, y1 = ref[1]
        refv = np.nanmedian(vmap[y0:y1, x0:x1])
    p5, mn = float(np.percentile(inner, 5)), float(inner.min())
    return dict(n_inner=int(inner.size), ref_reading=round(float(refv), 1), ring_4_7_reading=round(float(ring), 1),
                p5_rel_paper_ref=round(p5 - refv, 1), min_rel_paper_ref=round(mn - refv, 1),
                p5_rel_ring=round(p5 - ring, 1), min_rel_ring=round(mn - ring, 1))


def main():
    out = {"note": "digitised from published figures; approximate, see script docstring", "panels": {}}
    with tempfile.TemporaryDirectory() as td:
        pages, luts = {}, {}
        for name, spec in BARS.items():
            key = (spec["pdf"], spec["page"])
            if key not in pages:
                pages[key] = render(PAPERS / spec["pdf"], spec["page"], Path(td))
            luts[name] = bar_lut(pages[key], spec)
            h, v = luts[name]
            out[f"bar_{name}"] = dict(n=int(len(v)), vmin=round(float(v.min()), 1), vmax=round(float(v.max()), 1))
        for pname, p in PANELS.items():
            bar = BARS[p["bar"]]
            img = pages[(bar["pdf"], bar["page"])]
            rgb = np.asarray(img.crop(p["box"])).astype(float) / 255
            vmap = generic_filter(hue_to_value(rgb, *luts[p["bar"]]), np.nanmedian, size=3, mode="nearest")
            M, c0 = fit_affine(p["markers"], STATIONS)
            x_km, y_km = pixel_km(vmap.shape, M, c0)
            st = panel_stats(vmap, x_km, y_km, p["ref"])
            st["px_per_km"] = round(float(np.sqrt(abs(np.linalg.det(M)))), 2)
            out["panels"][pname] = st
    OUT.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore", message="All-NaN slice")
    main()
