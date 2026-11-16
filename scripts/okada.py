"""Okada (1985) surface displacements of a finite rectangular dislocation in an elastic half-space.

Reference: Okada, Y. 1985. Surface deformation due to shear and tensile faults in a half-space. Bulletin of the
Seismological Society of America 75(4): 1135-1154. doi:10.1785/BSSA0750041135. Equations (25)-(30); the check list in
Table 2 (p. 1149) is reproduced by tests/test_core.py::TestOkada.

Native frame (`okada85`): x along strike, z up, the medium in z <= 0. The fault's reference point is (0, 0, -d) on its
DEEP edge; the fault spans 0 <= xi' <= L along strike and 0 <= eta' <= W up-dip, so its points are
(xi', eta' cos(dip), -d + eta' sin(dip)). With x along strike and the fault rising towards +y, it dips towards -y, which
is to the right of the strike direction (Aki-Richards convention). U1, U2, U3 are strike slip, dip slip and opening
of the hanging wall relative to the footwall. Displacements come out in the units of U; lengths only enter as ratios.

Singular points follow Okada (1992, BSSA 82:1018) practice: atan(xi eta / (q R)) = 0 when q = 0, I5 = 0 when xi = 0,
and ln(R + eta) -> -ln(R - eta), 1/(R + eta) -> 0 where R + eta vanishes (and the same for xi).

`rect_source` wraps the native routine in local east/north/up coordinates with the source given by its centroid,
strike (degrees clockwise from north) and dip.
"""
from __future__ import annotations

import numpy as np

EPS = 1e-12


def _corner(xi, eta, q, sd, cd, mu_lm):
    """Chinnery term f(xi, eta) for the three dislocation types; returns 9 arrays (ss x,y,z; ds x,y,z; ts x,y,z)."""
    yt = eta * cd + q * sd
    dt = eta * sd - q * cd
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    X = np.sqrt(xi ** 2 + q ** 2)

    r_eta = R + eta
    r_xi = R + xi
    bad_eta = r_eta < EPS
    bad_xi = r_xi < EPS
    ln_r_eta = np.where(bad_eta, -np.log(np.where(bad_eta, R - eta, 1.0)), np.log(np.where(bad_eta, 1.0, r_eta)))
    inv_r_reta = np.where(bad_eta, 0.0, 1.0 / (R * np.where(bad_eta, 1.0, r_eta)))
    inv_reta = np.where(bad_eta, 0.0, 1.0 / np.where(bad_eta, 1.0, r_eta))
    inv_r_rxi = np.where(bad_xi, 0.0, 1.0 / (R * np.where(bad_xi, 1.0, r_xi)))

    qzero = np.abs(q) < EPS
    atn = np.where(qzero, 0.0, np.arctan(xi * eta / np.where(qzero, 1.0, q * R)))

    if abs(cd) > EPS:
        xzero = np.abs(xi) < EPS
        num = eta * (X + q * cd) + X * (R + X) * sd
        den = xi * (R + X) * cd
        i5 = np.where(xzero, 0.0, mu_lm * 2.0 / cd * np.arctan(num / np.where(xzero, 1.0, den)))
        i4 = mu_lm / cd * (np.log(R + dt) - sd * ln_r_eta)
        i3 = mu_lm * (yt / (cd * (R + dt)) - ln_r_eta) + sd / cd * i4
        i1 = mu_lm * (-xi / (cd * (R + dt))) - sd / cd * i5
    else:
        rd = R + dt
        i1 = -mu_lm / 2.0 * xi * q / rd ** 2
        i3 = mu_lm / 2.0 * (eta / rd + yt * q / rd ** 2 - ln_r_eta)
        i4 = -mu_lm * q / rd
        i5 = -mu_lm * xi * sd / rd
    i2 = mu_lm * (-ln_r_eta) - i3

    c = -1.0 / (2 * np.pi)
    ss = (c * (xi * q * inv_r_reta + atn + i1 * sd),
          c * (yt * q * inv_r_reta + q * cd * inv_reta + i2 * sd),
          c * (dt * q * inv_r_reta + q * sd * inv_reta + i4 * sd))
    ds = (c * (q / R - i3 * sd * cd),
          c * (yt * q * inv_r_rxi + cd * atn - i1 * sd * cd),
          c * (dt * q * inv_r_rxi + sd * atn - i5 * sd * cd))
    t = 1.0 / (2 * np.pi)
    ts = (t * (q ** 2 * inv_r_reta - i3 * sd ** 2),
          t * (-dt * q * inv_r_rxi - sd * (xi * q * inv_r_reta - atn) - i1 * sd ** 2),
          t * (yt * q * inv_r_rxi + cd * (xi * q * inv_r_reta - atn) - i5 * sd ** 2))
    return ss + ds + ts


def okada85(x, y, d, dip_deg, L, W, U1=0.0, U2=0.0, U3=0.0, nu=0.25):
    """Surface displacement (ux, uy, uz) in Okada's native frame (see module docstring). d = depth of the deep edge."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    dip = np.radians(dip_deg)
    sd, cd = np.sin(dip), np.cos(dip)
    if abs(cd) < 1e-10:
        cd = 0.0
        sd = float(np.sign(sd))
    mu_lm = 1.0 - 2.0 * nu   # mu / (lambda + mu)
    p = y * cd + d * sd
    q = y * sd - d * cd
    out = np.zeros((9,) + np.broadcast(x, y).shape)
    for xi, eta, sgn in ((x, p, 1), (x, p - W, -1), (x - L, p, -1), (x - L, p - W, 1)):
        out += sgn * np.array(np.broadcast_arrays(*_corner(xi, eta, q, sd, cd, mu_lm)))
    ux = U1 * out[0] + U2 * out[3] + U3 * out[6]
    uy = U1 * out[1] + U2 * out[4] + U3 * out[7]
    uz = U1 * out[2] + U2 * out[5] + U3 * out[8]
    return ux, uy, uz


def rect_source(e, n, ec, nc, zc, length, width, strike_deg, dip_deg, strike_slip=0.0, dip_slip=0.0,
                opening=0.0, nu=0.25):
    """Surface displacement (uE, uN, uZ) at local points (e, n) of a rectangle with centroid (ec, nc) at depth zc > 0.

    Lengths in one unit (km here), displacement in the unit of the slip/opening. Strike clockwise from north; the
    plane dips to the right of strike. For a horizontal sill use dip_deg = 0 (strike then orients its length)."""
    phi, dip = np.radians(strike_deg), np.radians(dip_deg)
    s = np.array([np.sin(phi), np.cos(phi)])        # unit vector along strike (E, N)
    yh = np.array([-np.cos(phi), np.sin(phi)])      # Okada +y: up-dip side, to the left of strike
    d = zc + width / 2 * np.sin(dip)                 # deep-edge depth
    o = np.array([ec, nc]) - length / 2 * s - width / 2 * np.cos(dip) * yh
    de, dn = np.asarray(e, float) - o[0], np.asarray(n, float) - o[1]
    x = de * s[0] + dn * s[1]
    y = de * yh[0] + dn * yh[1]
    ux, uy, uz = okada85(x, y, d, dip_deg, length, width, strike_slip, dip_slip, opening, nu)
    return ux * s[0] + uy * yh[0], ux * s[1] + uy * yh[1], uz
