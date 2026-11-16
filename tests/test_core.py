"""Unit tests on synthetic cases with known answers. Run: .venv/bin/python -I -m unittest discover -s tests

No downloaded data are needed; every input is generated here."""
from __future__ import annotations

import hashlib
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import ivwlib  # noqa: E402
from ivwlib import clean_series, linear_velocity, mogi, to_local_km  # noqa: E402


class TestProjection(unittest.TestCase):
    def test_origin_is_zero(self):
        x, y = to_local_km(-22.44, 63.87, -22.44, 63.87)
        self.assertAlmostEqual(float(x), 0.0, places=9)
        self.assertAlmostEqual(float(y), 0.0, places=9)

    def test_meridian_degree_at_64N(self):
        # meridian arc 63.37 -> 64.37 deg on WGS84; series 111132.92 - 559.82 cos 2phi + 1.175 cos 4phi (m per deg)
        phi = np.radians(63.87)
        ref = (111132.92 - 559.82 * np.cos(2 * phi) + 1.175 * np.cos(4 * phi)) / 1e3
        _, y = to_local_km(-22.44, 64.37, -22.44, 63.37)
        self.assertAlmostEqual(float(y), ref, delta=0.01)   # 10 m over 111 km
        self.assertGreater(float(y), 111.4)                  # the old formula used 110.57

    def test_parallel_degree_at_64N(self):
        # length of 0.01 deg of longitude at 63.87 N: (pi/180) N cos(phi) / 100, N = prime-vertical radius
        a, f = 6378137.0, 1 / 298.257223563
        e2 = f * (2 - f)
        phi = np.radians(63.87)
        n = a / np.sqrt(1 - e2 * np.sin(phi) ** 2)
        ref = np.radians(0.01) * n * np.cos(phi) / 1e3
        x, _ = to_local_km(-22.43, 63.87, -22.44, 63.87)
        self.assertAlmostEqual(float(x), ref, delta=1e-5)

    def test_distance_is_geodesic(self):
        lon, lat = np.array([-22.0, -21.4, -22.7]), np.array([64.1, 63.9, 63.8])
        x, y = to_local_km(lon, lat, -22.44, 63.87)
        _, _, d = ivwlib.GEOD.inv(np.full(3, -22.44), np.full(3, 63.87), lon, lat)
        np.testing.assert_allclose(np.hypot(x, y), d / 1e3, rtol=1e-12)

    def test_shapes_broadcast(self):
        lon, lat = np.meshgrid(np.linspace(-22.6, -22.3, 4), np.linspace(63.8, 63.9, 3))
        x, y = to_local_km(lon, lat, -22.44, 63.87)
        self.assertEqual(x.shape, (3, 4))
        self.assertTrue(np.all(np.diff(x, axis=1) > 0) and np.all(np.diff(y, axis=0) > 0))


class TestMogi(unittest.TestCase):
    def test_on_axis_uplift(self):
        # u_z(0) = (1 - nu) dV / (pi d^2)  (Segall 2010, eq. 7.14); nu = 0.25, dV = 1e6 m3, d = 4 km -> 14.92 mm
        _, _, uz = mogi(0.0, 0.0, 4.0, 1e6, nu=0.25)
        self.assertAlmostEqual(float(uz) * 1e3, 0.75e6 / (np.pi * 16e6) * 1e3, places=9)
        self.assertAlmostEqual(float(uz) * 1e3, 14.92, places=2)

    def test_horizontal_to_vertical_ratio(self):
        ux, uy, uz = mogi(3.0, 4.0, 2.5, 1e6)
        self.assertAlmostEqual(float(np.hypot(ux, uy) / uz), 5.0 / 2.5, places=12)
        self.assertAlmostEqual(float(ux / uy), 3.0 / 4.0, places=12)

    def test_uplift_volume_integral(self):
        # integral of u_z over the free surface = 2 (1 - nu) dV  (1.5 dV for nu = 0.25)
        d, dv = 3.0, 2e6
        r = np.linspace(0, 400 * d, 400001)
        uz = mogi(r, np.zeros_like(r), d, dv)[2]
        vol = np.trapezoid(uz * 2 * np.pi * r * 1e3, r * 1e3)
        tail = 2 * np.pi * (1 - 0.25) * dv / np.pi * (d * 1e3) / (r[-1] * 1e3)   # analytic tail beyond r_max
        self.assertAlmostEqual((vol + tail) / dv, 1.5, places=4)

    def test_linear_in_volume_and_sign(self):
        a = np.array(mogi(1.0, 2.0, 4.0, 1e6))
        b = np.array(mogi(1.0, 2.0, 4.0, -3e6))
        np.testing.assert_allclose(b, -3 * a)


class TestOkada(unittest.TestCase):
    """Okada (1985, BSSA 75:1135, doi:10.1785/BSSA0750041135) Table 2 check list, p. 1149, lambda = mu (nu = 0.25),
    L = 3, W = 2, unit dislocation. Values as printed (4 significant digits); cross-checked against the check list in
    F. Beauducel's okada85_checklist.m (github.com/IPGP/deformation-lib)."""
    TABLE2 = {  # case: (x, y, d, dip), {type: (ux, uy, uz)}
        2: ((2, 3, 4, 70), {"strike": (-8.689e-3, -4.298e-3, -2.747e-3), "dip": (-4.682e-3, -3.527e-2, -3.564e-2),
                            "tensile": (-2.660e-4, +1.056e-2, +3.214e-3)}),
        3: ((0, 0, 4, 90), {"strike": (0.0, +5.253e-3, 0.0), "dip": (0.0, 0.0, 0.0),
                            "tensile": (+1.223e-2, 0.0, -1.606e-2)}),
        4: ((0, 0, 4, -90), {"strike": (0.0, -1.303e-3, 0.0), "dip": (0.0, 0.0, 0.0),
                             "tensile": (+3.507e-3, 0.0, -7.740e-3)}),
    }
    UNIT = {"strike": (1, 0, 0), "dip": (0, 1, 0), "tensile": (0, 0, 1)}

    def test_table2_check_list(self):
        from okada import okada85
        for case, ((x, y, d, dip), rows) in self.TABLE2.items():
            for kind, ref in rows.items():
                got = okada85(x, y, d, dip, 3.0, 2.0, *self.UNIT[kind], nu=0.25)
                for g, r in zip(got, ref):
                    # printed to 4 significant digits: half a unit in the last digit
                    tol = 0.5 * 10 ** (np.floor(np.log10(abs(r))) - 3) if r else 1e-9
                    self.assertAlmostEqual(float(g), r, delta=tol, msg=f"case {case} {kind}")

    def test_geographic_wrapper_reproduces_case2(self):
        from okada import rect_source
        L, W, dip, d = 3.0, 2.0, 70.0, 4.0
        zc = d - W / 2 * np.sin(np.radians(dip))
        ref = self.TABLE2[2][1]
        for kind, (ux, uy, uz) in ref.items():
            s, ds, op = self.UNIT[kind]
            # strike 90 (x = east, Okada +y = north): centroid at (L/2, W/2 cos dip)
            e, n, z = rect_source(2.0, 3.0, L / 2, W / 2 * np.cos(np.radians(dip)), zc, L, W, 90.0, dip, s, ds, op)
            np.testing.assert_allclose([float(e), float(n), float(z)], [ux, uy, uz], rtol=5e-4, atol=6e-7)  # 4 significant digits
            # strike 0 (x = north, Okada +y = west): observation at N = 2, E = -3
            e, n, z = rect_source(-3.0, 2.0, -W / 2 * np.cos(np.radians(dip)), L / 2, zc, L, W, 0.0, dip, s, ds, op)
            np.testing.assert_allclose([float(n), float(-e), float(z)], [ux, uy, uz], rtol=5e-4, atol=6e-7)  # 4 significant digits

    def test_small_sill_is_a_horizontal_point_crack(self):
        # horizontal tensile point source: u = 3 dV d^2 (x, y, d) / (2 pi R^5)  (Okada 1985 eq. 4-6 with dip = 0)
        from okada import rect_source
        a, op, d = 0.05, 1.0, 4.0
        dv = a * a * op
        e, n = np.array([0.0, 1.5, 3.0, -2.0]), np.array([0.0, 2.0, -1.0, 4.0])
        ue, un, uz = rect_source(e, n, 0.0, 0.0, d, a, a, 30.0, 0.0, opening=op)
        R5 = (e ** 2 + n ** 2 + d ** 2) ** 2.5
        np.testing.assert_allclose(uz, 3 * dv * d ** 3 / (2 * np.pi * R5), rtol=1e-3)
        np.testing.assert_allclose(ue, 3 * dv * d ** 2 * e / (2 * np.pi * R5), rtol=1e-3, atol=1e-12)
        np.testing.assert_allclose(un, 3 * dv * d ** 2 * n / (2 * np.pi * R5), rtol=1e-3, atol=1e-12)

    def test_sill_uplift_volume_equals_opening_volume(self):
        # for a horizontal tensile crack the surface uplift integrates to the crack volume, for any nu
        from okada import rect_source
        g = np.linspace(-60, 60, 1201)
        e, n = np.meshgrid(g, g)
        _, _, uz = rect_source(e, n, 0.0, 0.0, 3.0, 4.0, 2.0, 0.0, 0.0, opening=1.0, nu=0.3)
        vol = uz.sum() * (g[1] - g[0]) ** 2
        tail = 8.0 * 3.0 / 60.0 * 0.5   # crude bound of what lies beyond 60 km (point-crack tail ~ 8*d/(r_max) / 2)
        self.assertAlmostEqual(vol / 8.0, 1.0, delta=tail / 8.0 + 0.01)

    def test_vertical_dike_symmetry(self):
        from okada import rect_source
        e = np.array([1.0, -1.0, 2.5, -2.5])
        n = np.array([0.7, 0.7, -1.2, -1.2])
        ue, un, uz = rect_source(e, n, 0.0, 0.0, 2.0, 5.0, 3.0, 0.0, 90.0, opening=1.0)
        np.testing.assert_allclose(uz[0], uz[1], rtol=1e-12)          # uplift symmetric about the dike plane
        np.testing.assert_allclose(ue[0], -ue[1], rtol=1e-12)         # opening pushes the two sides apart
        self.assertGreater(float(ue[0]), 0)
        np.testing.assert_allclose(un[2], un[3], rtol=1e-12)


class TestJoint(unittest.TestCase):
    def test_geometry_shift_matches_direct_projection(self):
        # the joint fit places a trial source by shifting about one origin; within ~30 km this must agree with the
        # exact local projection about the source to < 0.1 % of the distance
        import pandas as pd
        from ivwlib import to_local_km
        from svartsengi_joint import Geometry
        lon = np.array([-22.70, -22.30, -22.45, -22.20]); lat = np.array([63.82, 63.95, 63.88, 63.99])
        geo = Geometry(pd.DataFrame(dict(lon=lon, lat=lat)), {}, -22.454, 63.869)
        for slon, slat in ((-22.50, 63.85), (-22.40, 63.90)):
            sx, sy = geo.src(slon, slat)
            ex, ey = to_local_km(lon, lat, slon, slat)
            d = np.hypot(ex, ey)
            err = np.hypot(geo.xg - sx - ex, geo.yg - sy - ey)
            self.assertTrue(np.all(err / d < 1e-3), msg=f"{err / d}")

    def test_chain_picks_longest_span_with_fewest_pairs(self):
        import tempfile
        from pathlib import Path
        import svartsengi_joint as J
        with tempfile.TemporaryDirectory() as t:
            d = Path(t) / "F" / "los_mm"; d.mkdir(parents=True)
            for n in ("20250101_20250113", "20250113_20250125", "20250101_20250107", "20250107_20250113",
                      "20250125_20250131", "20250201_20250213"):
                (d / f"{n}.tif").touch()
            old = J.RAW; J.RAW = Path(t)
            try:
                c = J.chain("F", "2025-01-01", "2025-01-31")
            finally:
                J.RAW = old
            self.assertEqual(c, ["20250101_20250113", "20250113_20250125", "20250125_20250131"])


class TestCovariance(unittest.TestCase):
    def test_exponential_covariance_recovered(self):
        # synthetic field with known exponential covariance (s = 3 mm, lam = 2.5 km) plus white noise (1 mm)
        from joint_covariance import fit_covariance
        rng = np.random.default_rng(3)
        g = np.arange(0, 20, 0.5)
        x, y = [a.ravel() for a in np.meshgrid(g, g)]
        r = np.hypot(x[:, None] - x[None, :], y[:, None] - y[None, :])
        C = 9.0 * np.exp(-r / 2.5) + 1.0 * np.eye(len(x))
        v = np.linalg.cholesky(C) @ rng.standard_normal(len(x))
        c = fit_covariance(x, y, v)
        self.assertAlmostEqual(c["lam_km"], 2.5, delta=1.0)
        self.assertAlmostEqual(np.sqrt(c["s2"]), 3.0, delta=1.0)

    def test_whitened_solve_equals_weighted_lstsq_for_diagonal_covariance(self):
        import pandas as pd
        import svartsengi_gnss as G
        from joint_covariance import whitened_solve
        obs = pd.DataFrame(dict(site=["A", "B"], lon=[0, 0], lat=[0, 0], de=[1.0, 2.0], dn=[0.5, -1.0], du=[3.0, 4.0]))
        ins = {"F": pd.DataFrame(dict(los=np.array([1.0, 2.0, 0.5])))}
        rng = np.random.default_rng(0)
        A = rng.standard_normal((9, 3))
        d = np.r_[G.data_vector(obs), ins["F"].los.values]
        sig_i = 2.0
        m, c2, _ = whitened_solve(A, d, obs, ins, {"F": sig_i * np.eye(3)})
        s = np.r_[G.sigmas(d[:6], 2), np.full(3, sig_i)]
        m2, *_ = np.linalg.lstsq(A / s[:, None], d / s, rcond=None)
        np.testing.assert_allclose(m, m2, rtol=1e-10)
        self.assertAlmostEqual(c2, float(np.sum(((A @ m2 - d) / s) ** 2)), places=8)


class TestRateFitting(unittest.TestCase):
    def series(self, rate_u=-16.4, n_days=1500, noise=3.0, seed=1):
        rng = np.random.default_rng(seed)
        t = pd.date_range("2015-07-01", periods=n_days, freq="D")
        yr = np.arange(n_days) / 365.25
        return pd.DataFrame({"e": 10 * yr + rng.normal(0, noise, n_days), "n": -2 * yr + rng.normal(0, noise, n_days),
                             "u": rate_u * yr + rng.normal(0, 3 * noise, n_days)}, index=t)

    def test_linear_velocity_recovers_rates(self):
        v, s, n = linear_velocity(self.series(), "2015-07-01", "2019-12-31")
        self.assertEqual(n, 1500)
        self.assertAlmostEqual(v.e, 10.0, delta=5 * s.e)
        self.assertAlmostEqual(v.u, -16.4, delta=5 * s.u)
        self.assertLess(s.u, 0.5)

    def test_clean_series_removes_vertical_only_step(self):
        df = self.series(noise=1.0)
        df.loc[df.index >= "2017-03-01", "u"] += 470.0
        out, steps = clean_series(df, event_dates=[])
        self.assertEqual(len(steps), 1)
        self.assertEqual(steps[0]["date"], "2017-03-01")
        self.assertAlmostEqual(steps[0]["offset_up_mm"], 470.0, delta=5)
        v, _, _ = linear_velocity(out, "2015-07-01", "2019-12-31")
        self.assertAlmostEqual(v.u, -16.4, delta=1.0)

    def test_clean_series_keeps_step_with_horizontal_motion(self):
        df = self.series(noise=1.0)
        m = df.index >= "2017-03-01"
        df.loc[m, "u"] += 470.0
        df.loc[m, "e"] += 300.0   # a dike-like step moves the station sideways: must not be "cleaned"
        _, steps = clean_series(df, event_dates=[])
        self.assertEqual(steps, [])

    def test_clean_series_respects_event_guard(self):
        df = self.series(noise=1.0)
        df.loc[df.index >= "2017-03-01", "u"] += 470.0
        _, steps = clean_series(df, event_dates=[pd.Timestamp("2017-03-02")])
        self.assertEqual(steps, [])


class TestMogiFit(unittest.TestCase):
    def test_grid_fit_recovers_source(self):
        import svartsengi_gnss as sv
        src = dict(lon=-22.448, lat=63.872, depth_km=4.0, dv=12.0)   # Mm3
        sites = pd.DataFrame(dict(site=list("ABCDE"), lon=[-22.434, -22.417, -22.472, -22.331, -22.738],
                                  lat=[63.882, 63.844, 63.835, 63.979, 63.974]))
        x, y = to_local_km(sites.lon.values, sites.lat.values, src["lon"], src["lat"])
        ux, uy, uz = mogi(x, y, src["depth_km"], src["dv"] * 1e6)
        obs = sites.assign(de=ux * 1e3 + 4.0, dn=uy * 1e3 - 2.0, du=uz * 1e3)
        lons = np.round(np.arange(-22.46, -22.435, 0.004), 3)
        lats = np.round(np.arange(63.866, 63.879, 0.002), 3)
        p, pred, chi = sv.fit_grid(obs, lons, lats, np.arange(3.0, 5.01, 0.25))
        self.assertAlmostEqual(p["lon"], src["lon"], places=3)
        self.assertAlmostEqual(p["lat"], src["lat"], places=3)
        self.assertAlmostEqual(p["depth_km"], 4.0)
        self.assertAlmostEqual(p["dv_mm3"], 12.0, places=6)
        self.assertAlmostEqual(p["te"], 4.0, places=6)
        self.assertLess(p["chi2"], 1e-12)

    def test_depth_interval_rescaling(self):
        import svartsengi_gnss as sv
        depths = np.arange(3.0, 5.01, 0.25)
        chi = 10 + 8 * (depths - 4.0) ** 2           # minimum 10 at 4 km
        self.assertEqual(sv.depth_interval(chi, depths, 1.0), (3.75, 4.25))      # delta chi2 <= 1
        self.assertEqual(sv.depth_interval(chi, depths, 9.0), (3.0, 5.0))        # errors x3 -> delta chi2 <= 9


class TestLevels(unittest.TestCase):
    def test_pre_and_post_level_on_a_noisy_step(self):
        import svartsengi_gnss as sv
        rng = np.random.default_rng(3)
        t = pd.date_range("2024-01-01", periods=40, freq="D")
        s = pd.Series(np.where(t < pd.Timestamp("2024-01-20"), 20.0, 8.0) + rng.normal(0, 0.8, 40), index=t)
        onset = pd.Timestamp("2024-01-20")
        loss_new = sv.pre_level(s, onset) - sv.post_level(s, onset)[0]
        loss_old = s[(s.index < onset)].iloc[-1] - s[(s.index > onset) & (s.index <= onset + pd.Timedelta(days=7))].min()
        self.assertAlmostEqual(loss_new, 12.0, delta=1.0)
        # the raw minimum is biased high in magnitude on noise; the running median much less
        self.assertLess(abs(loss_new - 12.0), abs(loss_old - 12.0) + 0.5)


class TestLosAndMaps(unittest.TestCase):
    def test_decompose_and_project_round_trip(self):
        import hellisheidi_matched as hm
        ga = {"E": np.full((2, 2), -0.604), "U": np.full((2, 2), 0.787)}
        gd = {"E": np.full((2, 2), 0.553), "U": np.full((2, 2), 0.824)}
        east, up = np.array([[-3.0, 1.0], [0.0, 2.0]]), np.array([[-15.0, 0.0], [-5.0, 4.0]])
        asc, desc = ga["E"] * east + ga["U"] * up, gd["E"] * east + gd["U"] * up
        e2, u2 = hm.decompose(asc, desc, ga, gd)
        np.testing.assert_allclose(e2, east, atol=1e-12)
        np.testing.assert_allclose(u2, up, atol=1e-12)
        np.testing.assert_allclose(hm.project(east, up), -0.50 * east + 0.86 * up)

    def test_los_sign_convention(self):
        # LiCSBAS convention used in licsar_annual.py: d_LOS = -unw * lambda / (4 pi); positive = towards satellite.
        # Uplift (+U) with an upward-pointing LOS unit vector (U component > 0) must give positive LOS.
        import hellisheidi_matched as hm
        self.assertGreater(float(hm.project(np.array(0.0), np.array(10.0))), 0)

    def test_plane_removal_and_peak(self):
        import hellisheidi_matched as hm
        lon, lat = np.meshgrid(np.arange(-21.80, -21.10, 0.002), np.arange(63.90, 64.16, 0.001))
        x, y = to_local_km(lon, lat, *hm.JUNCU_CENTRE)
        plane = 3.0 + 0.4 * x - 0.2 * y
        bowl = -20.0 * np.exp(-(x ** 2 + y ** 2) / (2 * 1.0 ** 2))
        out, coef = hm.remove_plane(plane + bowl, lon, lat)
        self.assertAlmostEqual(coef[1], 0.4, delta=0.01)
        pk = hm.peak(out, lon, lat, hm.JUNCU_CENTRE, stat="min")
        self.assertAlmostEqual(pk, -20.0, delta=0.6)   # 3x3 smoothing of a 1 km Gaussian costs < 3 %
        ring = hm.ref_value(plane + bowl, lon, lat, hm.JUNCU_CENTRE, 7.0, 4.0)
        self.assertAlmostEqual(ring, 3.0, delta=0.3)   # a planar ramp averages out over a full ring


class TestDigitize(unittest.TestCase):
    def test_hue_inversion_and_affine(self):
        import matplotlib
        import digitize_published_figures as dg
        cmap = matplotlib.colormaps["rainbow"]
        vals = np.linspace(-20, 20, 401)
        rgb = cmap((vals + 20) / 40)[:, :3]
        from matplotlib.colors import rgb_to_hsv
        h = rgb_to_hsv(rgb)[:, 0]
        test = np.array([-17.3, -5.0, 0.0, 12.4])
        got = dg.hue_to_value(cmap((test + 20) / 40)[:, :3][None] * 0.8, h, vals)[0]   # darkened like hillshade
        np.testing.assert_allclose(got, test, atol=0.15)
        lonlat = dict(dg.STATIONS)
        M_true, c_true = np.array([[18.0, 0.5], [-0.3, -18.2]]), np.array([300.0, 290.0])
        px = {}
        for n, ll in lonlat.items():
            k = np.array([float(v) for v in to_local_km(*ll, *dg.JUNCU_CENTRE)])
            px[n] = tuple(M_true @ k + c_true)
        M, c0 = dg.fit_affine(px, lonlat)
        np.testing.assert_allclose(M, M_true, atol=1e-9)
        np.testing.assert_allclose(c0, c_true, atol=1e-9)


class TestLicsarSelection(unittest.TestCase):
    def test_summer_pairs(self):
        from licsar_annual import summer_pairs
        pairs = ["20150603_20160609", "20150615_20160609", "20151010_20161005", "20150603_20150615",
                 "20160609_20170604", "20160609_20180604"]
        sel = summer_pairs(pairs, range(2015, 2017), per_year=3)
        # winter pair, 12-day pair and 2-year pair are rejected; ordered by closeness to 365 days (360 d before 372 d)
        self.assertEqual(sel, ["20150615_20160609", "20150603_20160609", "20160609_20170604"])


class TestVendor(unittest.TestCase):
    def test_seislib_matches_upstream(self):
        lines = (ROOT / "scripts" / "seislib.py").read_bytes().split(b"\n", 2)
        self.assertIn(b"855571eb03b542cff41bd95ecf54316917fefd9425150c63a37e9caab2e3be36", lines[1])
        self.assertEqual(hashlib.sha256(lines[2]).hexdigest(),
                         "855571eb03b542cff41bd95ecf54316917fefd9425150c63a37e9caab2e3be36")


if __name__ == "__main__":
    unittest.main()
