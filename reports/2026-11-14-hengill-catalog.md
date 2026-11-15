# Report, 2026-11-14: the Zenodo "Hengill seismic catalog"

> **Note (2026-11-16):** some numbers and verdicts in this report were superseded (coordinate fix, loss method,
> event selection, Hellisheiði definitions, citations). See [docs/corrections.md](../docs/corrections.md); the
> current values are in the README.

Scope: check the catalogue (Zenodo 10.5281/zenodo.16311162), compare it with IMO, relate Húsmúli seismicity to
re-injection and the HUSM GNSS series, and see whether it helps with the Hellisheiði InSAR mismatch.
Script `scripts/hengill_catalog.py`; numbers in `docs/data/hengill_catalog.json`; method in docs/method.md, Pilot D.

## Licence, method, magnitude type

- Zenodo record: CC BY 4.0, creator "University of Campania", published 2025-07-22. The description lists five
  columns and nothing else: no method, no network, no magnitude type.
- The rows are the COSEISMIQ high-quality double-difference catalogue (Grigoli et al. 2022, *Scientific Data*
  9:220; figshare, CC0): 8,691 rows in both, epicentres equal within single-precision rounding (8,685 rows within
  5×10⁻⁶ degrees). Method: automatic picking on a dense temporary network, then double-difference relocation.
  Magnitude: local ML (Wood-Anderson simulation), not tied to IMO's scale.

## Defects in the Zenodo copy

| Defect | Evidence |
|---|---|
| Columns 4 and 5 described the wrong way round | column 4 = HQDD depth (8,653 rows within 1 m), column 5 = HQDD magnitude |
| Minus sign of magnitudes lost | 2,054 negative magnitudes in HQDD, 0 in the Zenodo file; column 5 = abs(ML) for 8,686 rows |
| "time (sec)" is not seconds | t = (days since 2018-01-01 + 1)/365, 7 significant digits: 3 s resolution before 2019, 32 s after |
| One-day time errors | 1,517 rows (Dec 2018, Dec 2019-Feb 2020, Dec 2020) are one day off |

Because of these, the analysis uses the original HQDD file.

## Results

- **(a) vs IMO, Hengill box, Dec 2018-Jan 2021.** COSEISMIQ 7,969 events, IMO SeisComP 1,928, SIL 2,017.
  Mc 0.2 / 0.4 / 0.8 in each catalogue's own scale. Matched-event offsets (IMO − COSEISMIQ): SeisComP −0.51,
  SIL −0.15, so COSEISMIQ is complete about 0.7 units lower. b (95 % CI): 0.88 (0.86-0.90), 0.87 (0.82-0.94),
  0.91 (0.85-0.98): no difference. Rate above Mc: 6.2, 1.0 and 0.94 events per day. 88 % of IMO events above IMO's
  Mc are in COSEISMIQ.
- **(b) Húsmúli.** 47 % of COSEISMIQ events in the box are within 3 km of Húsmúli, median depth 2.0 km. IMO SIL at
  M ≥ 1.0: 29 per year in 2006-2010, 2,751 per year in Sep-Dec 2011 (largest 4.04), 158 per year in 2012-15
  (HUSM +2.1 mm/yr), 35 in 2017-19 (−5.1 mm/yr), 83 in 2020-25 (−9.9 mm/yr).
- **(c)** The catalogue does not explain the Hellisheiði InSAR mismatch: it measures earthquakes, not subsidence,
  and covers 2018-2021 only. It does not explain the HUSM change either: the Húsmúli rate drops when HUSM turns to
  subsidence and rises again while the subsidence speeds up.

## Side finding

IMO SeisComP has reviewed events in the Hengill box from 2015-06-20 (28,016 to 2026-11-13). Earlier docs said
"back to Oct 2023"; that was only where our Reykjanes fetch started. The docs are corrected.

## Not verified

- Injection and production volumes for Húsmúli and Hellisheiði (not open), so the cause of the HUSM sign change
  stays untested.
- The University of Campania study that produced `hen_conv.dat` was not found; why the signs and times changed
  is unknown. The authors were not contacted.
- COSEISMIQ's ML was not calibrated to Mw or to IMO beyond the median offsets above.
- The COSEISMIQ HQ (not DD) and lower-quality catalogues were not analysed.
- The 0.36 SeisComP-SIL offset implied here was not reconciled with the 0.62 measured on Reykjanes (Pilot B).
- There is no test suite in this repository; the scripts were run end to end.
