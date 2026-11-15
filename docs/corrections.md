# Corrections

Errors found in earlier versions of this repository, with the evidence and what changed. Newest first.
Numbers come from the files named; rerun the scripts listed in [method.md](method.md) to reproduce them.

## 2026-11-16

### Hellisheiði "mismatch": a definition difference, not a measurement difference

**Was:** "−14.6 mm/yr LOS vs published 24-26 mm/yr: lower, ~60 % of the published peak … not separated" (README §C,
pilot report).

**Evidence:**
- Juncu et al. 2017 (JGR 122:692-709, p. 704) give 26 mm/yr as "up to" the largest value of the TerraSAR-X track 41
  velocity field. That field is plate-motion-corrected and referenced to a point near GFEL (Fig. 7, p. 700). The
  24 mm/yr of EGU2015-13455 is the same kind of number.
- We read their Fig. 7 with our own statistic (`scripts/digitize_published_figures.py`, results in
  `docs/data/published_figures_digitized.json`): 5th percentile within 1.5 km of their source, minus the median of a
  4-7 km ring.
  - Plate-corrected panel: −17.6 mm/yr.
  - Observed panel: −16.2 mm/yr.
  - Both relative to their reference point: −19.8 and −15.0. So the plate correction alone adds about −5 mm/yr at
    the bowl.
- Ducrocq et al. 2021 (Front. Earth Sci. 9:725109, Fig. 5, p. 7) show Sentinel-1 track 16 for 2015-2017. That is the
  same track as our LiCSAR frame 016A. Read the same way it gives −15.3 mm/yr, against our −14.6 (pilot C) and −13.1
  (same centre as Juncu). TerraSAR-X for 2015-2017 in the same figure gives −15.4.
- The pair-length (phase-bias) check is small. Sum of short pairs minus the direct pair, bowl minus far field:
  +0.22 mm over 72 days and +0.55 mm over 114 days (`hellisheidi_matched.json`, `closure_test`).

**Now:** an explained difference. See README §C and `reports/2026-11-16-hellisheidi-and-review.md`.

### Local coordinates: latitude scale was the equatorial value

**Was:** `ivwlib.to_local_km` used 110.57 km per degree of latitude and 111.32·cos(lat0) per degree of longitude. Its
comment claimed "error < 0.1 %".

**Evidence:**
- At 63.87° N the meridian degree is 111.48 km (111132.92 − 559.82 cos 2φ + 1.175 cos 4φ m). North distances were
  therefore 0.8 % short and east distances 0.27 % short.
- Test: `tests/test_core.py::TestProjection`.

**Now:** geodesic azimuthal-equidistant coordinates from `pyproj.Geod` (WGS84). Every script that uses it was rerun.
Changes:

| Quantity | Before | After |
|---|---|---|
| I1 fit | 3.50 km, 8.41 Mm³ | 3.75 km, 9.33 Mm³ |
| I3 fit | 12.17 Mm³ | 13.86 Mm³ |
| Median source depth | 3.75 km | 3.875 km |
| SENG G_z | 14.4 mm/Mm³ | 13.58 mm/Mm³ |
| Volumes (all) | | about +6 % |
| Hellisheiði and other InSAR field rates | | change by ≤ 0.2 mm/yr |
| Dike length | 15.3 km | 15.4 km |

### Volume per event: minimum-picking bias; post-hoc event selection

**Was:**
- The loss was the last value before onset minus the lowest daily value within 7 days.
- Events 4, 6, 7, 8, 9, 11 were called "agrees, ratio 0.87-1.21". That set was chosen after the comparison with
  Parks et al. had been seen.

**Now:**
- The loss is the median of the 5 days before onset minus the lowest 5-day running median after onset. It is
  0.6-1.9 Mm³ smaller for those six events (`event_volume_loss.csv`, both columns kept).
- A selection rule that uses only seismicity was fixed in `config/svartsengi.json` before it was applied: 5th
  percentile of the distance from SENG to the dike earthquakes in the first 24 h ≥ 3 km.
- **No event passes it.** The d5 values are 1.65-2.36 km.
- At 2 km, events 3, 4, 8, 9, 10 pass. Their loss ratios are −0.26 to 1.22, median 0.65.
- The single-station loss method is therefore **not validated** for any event by a blind rule. The earlier
  "agrees" is withdrawn.

### "Matches IMO" → "consistent with IMO within the source-geometry uncertainty"

**Was:** 28.0 / 29.2 Mm³ vs IMO 28.3 / 29.1 Mm³, "agrees".

**Now:**
- 29.1 / 30.2 Mm³ (median source, 3.875 km). The 18 Aug 2026 IMO value covers data to 16 Aug 2026 (figure caption),
  so we compare at 16 Aug.
- Over the depth range of our fits (3.50-4.75 km) and with the Parks et al. source, the same SENG uplift gives
  24.6-41.6 and 25.5-43.2 Mm³ (`volume_sensitivity_source.csv`).
- The near-equality reported before was partly coincidence. SENG's sensitivity is 18.1 / 14.4 / 9.5 mm per Mm³ at
  3.25 / 3.75 / 4.75 km.

### Numbers in the README that no script produced

| Quantity | Was | Reproduced | Output file |
|---|---|---|---|
| Inflow per quarter | 1.41 → 0.79 → 0.60 → 0.53 → 0.44 m³/s | 1.23 → 0.97 → 0.56 → 0.51 → 0.47 m³/s | `quarterly_inflow.csv`, 3-month steps from 17 Jul 2025 |
| SENG uplift, 3 months to 15 / 30 Sep 2026 | 54.6 / 46.7 mm | 49.2 / 41.1 mm | `seng_uplift_3months.csv` |
| Recharge range entered / top passed | ~15 Dec 2025 / ~24 Feb 2026 | 6 Dec 2025 / 16 Feb 2026 | `cycle_threshold_check.json` |

### Parks-source comparison used a different weighting

**Was:** our source was fitted by weighted least squares, the Parks et al. source by unweighted least squares.

**Now:** both are weighted (`mogi_fits.csv`: `dv_parks_src_mm3`; the unweighted value is kept as
`dv_parks_src_unweighted_mm3`).
- I1-I4 with the Parks source: 9.65 / 15.67 / 15.25 / 23.63 Mm³ (was 11.25 / 19.49 / 19.33 / 28.65 unweighted).
- With our source: 9.33 / 16.70 / 13.86 / 31.50 Mm³.

### Citations and licences

- **Juncu et al. "2018, JVGR".** The article is Juncu, D. et al. 2020, *JVGR* 391:106337,
  doi:10.1016/j.jvolgeores.2018.03.019 (online 2018; Crossref).
  - It reports "up to 2 cm of surface displacements during 2011-2012" at the start of injection (abstract). That is
    not "~2 cm uplift since 2011".
  - HUSM GNSS starts in May 2012, so our "+2.1 mm/yr 2012-15 matches" compared different things. It is withdrawn.
- **"Later HUSM subsidence not in the literature".** Ducrocq et al. 2021 report localised subsidence at Húsmúli in
  2017-2018: diameter 6 km, up to −8 mm vertical (p. 10). That fits the sign of our HUSM −5.1 mm/yr in 2017-19.
  2020-2026 is still not covered by a paper we found.
- **LiCSAR licence.** LiCSAR products are under the UK Open Government Licence v3 (CEDA record "LiCSAR
  interferometry products", Lazecký & Maghsoudi 2022), with a set acknowledgement. They are not CC BY 4.0 as assumed
  in the task brief. See README "Data and licences".

### LiCSAR pair counts

**Was:** "27 ascending + 9 descending pairs".

**Now:** 4 ascending and 3 descending of those files are empty (unwrapped phase all zero, or no coherence). So 23 + 6
pairs carry data. The rates were already computed only from finite pixels, so no rate changes. `hellisheidi_matched.py`
now also gates on coverage.

### Smaller wording fixes

- method.md said "SENG sits 1.8 km from the source". That is the distance to Parks et al.'s source. Ours is 1.3 km
  away.
- method.md said "±26 % between our source and Parks et al.'s". After the projection fix it is +19 %. Over our own
  depth range it is −16 % to +43 %.
- README §C said the first-cycle magmatic uplift rate was "~3,300 mm/yr". It is ~3,000 mm/yr: 7.0 m³/s × 13.58 mm
  per Mm³.
