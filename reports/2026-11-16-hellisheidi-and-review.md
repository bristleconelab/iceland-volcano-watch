# Report, 2026-11-16: the Hellisheiði rate, and an external review of Pilot A

Numbers come from `docs/data/hellisheidi_matched.json`, `published_figures_digitized.json`, `mogi_fits.csv`,
`mogi_correlation.json`, `volume_sensitivity_*.csv`, `event_volume_loss.csv` and `event_selection_rule.json`.
Errors found in earlier versions are listed in [docs/corrections.md](../docs/corrections.md).

## 1. Hellisheiði: why −14.6 is not "60 % of 24-26"

### What the published numbers are

| Source | Quantity | Period | Reference |
|---|---|---|---|
| Juncu et al. 2017, JGR 122:692-709, p. 704 | "up to 26 mm/yr LOS": the most negative value of a TerraSAR-X track 41 ascending StaMPS velocity field, 3 m pixels, after a plate-motion model is subtracted (Fig. 7, p. 700). LOS unit vector [−0.50, −0.12, 0.86] (E, N, U), incidence about 30.6° | 2012-2015 | a point near GFEL, about 10 km SW of the plant |
| EGU2015-13455 | "maximum of ~24 mm/yr" LOS, same sensor, after plate-motion correction | from Oct 2009 | not stated |
| Juncu et al. 2017, Fig. 3 (p. 696) | campaign GNSS HH47, −23 mm/yr vertical, relative to stable North America (ITRF08) | 2012-2015 | |

### What pilot C measured

- LiCSAR frame 016A_02504_162118, which is Sentinel-1 track 16 ascending. At the bowl the LOS unit vector is
  [−0.604, −0.124, 0.787], incidence 38°.
- Summer-to-summer pairs with 2015-2018 primaries.
- The 5th percentile of a 3×3-smoothed map within 1.5 km of (−21.40 E, 64.04 N), minus the median of a 4-7 km ring.
- No plate model; the ring removes a planar gradient.

So the quantity, the statistic and the reference all differ from the published ones. Only the sensor geometry was
discussed before.

### Matched definitions

Each column is the same statistic, the 5th percentile (p5) within 1.5 km of the Juncu et al. source (21.371 W,
64.033 N), applied to a different map. Published maps were digitised by colour (hue) against their own colour bars and
georeferenced on the HUSM, OLKE and HVER markers.

| Map | p5 minus ring median (4-7 km) | p5 relative to the paper's own reference |
|---|---|---|
| TerraSAR-X T41 2012-2015, plate-corrected (Juncu Fig. 7 right) | **−17.6** | −19.8 (min −21.8) |
| TerraSAR-X T41 2012-2015, as observed (Juncu Fig. 7 left) | −16.2 | −15.0 |
| TerraSAR-X T41 2015-2017 (Ducrocq et al. 2021 Fig. 5, p. 7) | −15.4 | −17.3 |
| Sentinel-1 T16 2015-2017 (Ducrocq et al. 2021 Fig. 5) | **−15.3** | −15.6 |
| **This work**, LiCSAR T16 2015-2019, median of 7 usable annual pairs | **−13.1** | −16.8 (GFEL area, 2 km) |
| This work, pilot C definition (centre −21.40/64.04) | −13.9 (pairs); −14.6 (yearly maps, README) | |

What each choice is worth, in mm/yr at the bowl:

- **Plate-motion correction: about −5.** Juncu's own panels, relative to their reference, read −15.0 observed and
  −19.8 corrected.
- **"Up to" vs p5: about −2 at the resolution we can read.** On the corrected panel the minimum is −21.8 and p5 is
  −19.8. The single most extreme 3 m pixel, 26, is beyond what a 3×3 median of the printed figure keeps.
- **Viewing geometry: below 1.**
  - Our asc+desc east and up maps re-projected onto the TerraSAR-X LOS give −15.8, against −15.5 for the
    Sentinel-1 LOS (plane removed, far field = 0).
  - The larger vertical weight of TerraSAR-X is offset by its smaller east weight.
- **Phase bias (pair length): below 1 per summer.** We summed a chain of 6/12-day pairs and subtracted the direct
  long pair over the same dates, bowl minus far field:
  - +0.22 mm over 72 days (2015, 6 pairs);
  - +0.55 mm over 114 days (2025, 10 pairs).

  The far-field spread of that closure is about 1 mm. Our rates come from single annual pairs, which carry no
  chain-sum bias. Their own fading bias was not measured directly. 2019 and 2023 short pairs return HTTP 403 from
  LiCSAR.
- **InSAR vs GNSS.** HUSM minus REYK projected onto the LOS, for the exact pair dates:
  - ascending: InSAR − GNSS = +0.16 ± 2.2 mm/yr (mean ± s.e., 8 pairs);
  - descending: −5.2 ± 3.7 (3 pairs).

  No significant bias in the ascending annual pairs. HUSM is on poorly coherent ground, so only 8 of 12 ascending
  pairs could be used.
- **Absolute vertical.** Our 2015-2019 vertical, tied to REYK's NGL rate (−2.66 mm/yr), gives p5 = −19.2 mm/yr. The
  2012-2015 campaign point HH47 gives −23. The two differ in period, point vs 1.5 km statistic, and tie point.

### Did the rate fall after 2015?

- Production at Hellisheiði was about 40 Mt/yr in 2012-2014 and was cut to about 32 Mt/yr in 2015 (Juncu et al. 2017,
  p. 696). The plant was expanded in 2016 with production wells in Hverahlíð (Ducrocq et al. 2021, p. 3).
- The same-sensor TerraSAR-X maps, read the same way, give −16.2 (2012-2015) and −15.4 (2015-2017). That difference is
  inside the digitising error: about 1 mm/yr from the colour bar, plus about 0.2 km of georeferencing.
- Our LiCSAR yearly series shows no trend from 2015 to 2026: slope +0.01 mm/yr per year for D0.
- So a production-driven drop is **not needed and not seen**. Open production and pressure series after 2015 were
  not found: the annual-report pages of the operator did not resolve, so this is unverified.

### Verdict: explained difference

The 24-26 mm/yr is the most extreme pixel of a plate-corrected, point-referenced TerraSAR-X field. With the same
statistic, the published maps give −15 to −18 mm/yr, and the published Sentinel-1 track-16 map gives −15.3. Our LiCSAR
value is −13.1 to −14.6.

A residual of about 2 mm/yr (15 %) remains against the published Sentinel-1 map and is not resolved further. Possible
causes, none tested:
- different periods (2015-17 vs 2015-19);
- our coherence mask (cc ≥ 0.25 keeps 38-89 % of the bowl pixels in the 7 ascending 2015-2018 pairs);
- SBAS stack vs median of annual pairs.

## 2. External review of Pilot A: point by point

| # | Point | Status | Evidence |
|---|---|---|---|
| 1 | χ²/dof ≫ 1, depth intervals not meaningful | **Fixed (reporting); partly disputed** | χ²/dof = 9.3, 24.5, 19.0, 30.3 (`mogi_fits.csv`). The depth intervals already rescaled the errors by sqrt(χ²/dof) (Δχ² ≤ χ²/dof): 3.0-4.5, 3.0-5.5, 2.75-4.5, 4.0-5.75 km. Without rescaling they are 3.5-3.75, 3.75-4.25, 3.25-3.5, 4.75 km. Rescaling assumes a correct model; the misfit is systematic, so the intervals are now labelled as not formal. |
| 2 | 29.2 vs 29.1 near-coincidental | **Confirmed, fixed** | SENG G_z = 18.14 / 14.35 / 9.51 mm per Mm³ at 3.25 / 3.75 / 4.75 km, so 450 mm = 24.8 / 31.4 / 47.3 Mm³ (reviewer: 24.7 / 31.3 / 47.2). Current-cycle volume at 21 Jul / 16 Aug 2026: 29.1 / 30.2 Mm³ for the median source; 24.6-41.6 / 25.5-43.2 over depth 3.5-4.75 km and the Parks source. IMO: 28.3 (27.5-29.1) / 29.1 (28.3-29.9). Wording is now "consistent within the source-geometry uncertainty". |
| 3 | Mogi ΔV is the chamber volume change | **Fixed (text)** | Method, README and outline now say that rates are lower bounds on magma inflow when magma is compressible (Rivalta & Segall 2008, GRL 35(4), doi:10.1029/2007GL032521). |
| 4 | Background held constant | **Tested** | Alternative windows give SENG −15.76 to −16.82 mm/yr. The 16 Aug 2026 volume changes by ≤ 0.05 Mm³ and the last-quarter inflow by ≤ 0.01 m³/s. Extreme case, geothermal term zero after 2019 (NYLA rate −5.51 only): volume −0.85 Mm³ (30.22 → 29.37), last quarter 0.47 → 0.45 m³/s (`volume_sensitivity_background.csv`). |
| 5 | Weak geometry | **Quantified** | Correlation depth-ΔV 0.87-0.95, east position-east translation 0.45-0.63. Scaled 1σ: depth 0.74-1.02 km, ΔV 35-42 %, position 0.3-0.6 km. With the horizontal position fixed at Parks': depth 4.0-4.25 km, χ²/dof 8.2-30.4 (`mogi_correlation.json`, `mogi_fits.csv`). |
| 6 | 110.57 km/deg latitude | **Confirmed, fixed** | pyproj geodesic (WGS84); meridian degree at 63.87° N = 111.48 km (test). Fits and volumes changed (corrections.md). |
| 7 | Post-hoc event selection | **Fixed; result worse** | Rule written before use, seismicity only (d5 ≥ 3 km). No event passes: d5 = 1.65-2.36 km. At 2 km: events 3, 4, 8, 9, 10, ratio −0.26 to 1.22, median 0.65. The single-station loss is not validated by a blind rule; "agrees" is withdrawn. |
| 8 | Minimum-picking bias | **Fixed** | 5-day medians before onset and the lowest 5-day running median after. Losses for events 4, 6, 7, 8, 9, 11 fall by 1.87, 1.26, 1.18, 0.61, 0.85, 0.81 Mm³ (old-method column kept in `event_volume_loss.csv`). |
| 9 | Unequal weighting | **Fixed** | Parks source weighted: 9.65 / 15.67 / 15.25 / 23.63 Mm³ vs our 9.33 / 16.70 / 13.86 / 31.50 Mm³ (I1-I4). Unweighted values kept. |
| 10 | No tests / CI | **Fixed** | `tests/test_core.py`, 22 synthetic tests (Mogi vs Segall 2010 eq. 7.14 and the 2(1−ν)ΔV uplift-volume identity, projection, rate fit, cleaning, LOS decomposition and sign, plane and peak, grid fit, levels, digitiser, pair selection, vendored seislib hash). `.github/workflows/tests.yml` runs them. CI has not run because nothing was pushed. |
| 11 | Long main(), constants in code | **Fixed** | `svartsengi_gnss.py` split into fitting, volume, sensitivity and plotting functions. Intervals, errors, the Parks source, the IMO values (with URL and caption dates) and the selection rule moved to `config/svartsengi.json`. |
| 12 | Vendored seislib | **Fixed** | Header names seismo-lab commit 5229c22 and the upstream sha256; a test checks the body hash. |
