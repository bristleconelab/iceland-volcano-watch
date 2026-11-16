"""Compound dislocation model (CDM) and point ellipsoidal cavity ratios, ported from Nikkhoo et al. (2017).

Python port of CDM.m (with RDdispSurf, AngSetupFSC, AngDisDispSurf) and of the Eshelby shape tensor in pECM.m,
version 2022-10-24 from https://volcanodeformation.com/software. A CDM is three mutually orthogonal rectangular
tensile dislocations with one uniform opening; with suitable semi-axes it represents a pressurised ellipsoidal cavity
of any orientation (Nikkhoo, Walter, Lundgren & Prats-Iraola 2017, GJI 208:877-894, doi:10.1093/gji/ggw427). The
angular-dislocation solutions use the artefact-free configuration of Nikkhoo & Walter (2015, GJI 201:1119-1141).

Original code: Copyright (c) 2016 Mehdi Nikkhoo. Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software
without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions: The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE
AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Units: coordinates, depth and semi-axes in one length unit; displacements in the unit of the opening.
Checked in tests/test_core.py against the repository's Okada (1985) implementation (a CDM with one zero semi-axis is
a single rectangular dislocation) and against the Mogi solution (a small CDM with equal semi-axes far away).
"""
from __future__ import annotations

import numpy as np


def _rot(omega_x, omega_y, omega_z):
    cx, sx = np.cos(np.radians(omega_x)), np.sin(np.radians(omega_x))
    cy, sy = np.cos(np.radians(omega_y)), np.sin(np.radians(omega_y))
    cz, sz = np.cos(np.radians(omega_z)), np.sin(np.radians(omega_z))
    Rx = np.array([[1, 0, 0], [0, cx, sx], [0, -sx, cx]])
    Ry = np.array([[cy, 0, -sy], [0, 1, 0], [sy, 0, cy]])
    Rz = np.array([[cz, sz, 0], [-sz, cz, 0], [0, 0, 1]])
    return Rz @ Ry @ Rx


def _ang_dis_disp_surf(y1, y2, beta, b1, b2, b3, nu, a):
    sinB, cosB = np.sin(beta), np.cos(beta)
    cotB = 1.0 / np.tan(beta)
    z1 = y1 * cosB + a * sinB
    z3 = y1 * sinB - a * cosB
    r = np.sqrt(y1 ** 2 + y2 ** 2 + a ** 2)
    Fi = 2 * np.arctan2(y2, (r + a) / np.tan(beta / 2) - y1)
    n2 = 1 - 2 * nu
    v1b1 = b1 / 2 / np.pi * ((1 - n2 * cotB ** 2) * Fi + y2 / (r + a) * (n2 * (cotB + y1 / 2 / (r + a)) - y1 / r)
                             - y2 * (r * sinB - y1) * cosB / r / (r - z3))
    v2b1 = b1 / 2 / np.pi * (n2 * ((.5 + cotB ** 2) * np.log(r + a) - cotB / sinB * np.log(r - z3))
                             - 1 / (r + a) * (n2 * (y1 * cotB - a / 2 - y2 ** 2 / 2 / (r + a)) + y2 ** 2 / r)
                             + y2 ** 2 * cosB / r / (r - z3))
    v3b1 = b1 / 2 / np.pi * (n2 * Fi * cotB + y2 / (r + a) * (2 * nu + a / r)
                             - y2 * cosB / (r - z3) * (cosB + a / r))
    v1b2 = b2 / 2 / np.pi * (-n2 * ((.5 - cotB ** 2) * np.log(r + a) + cotB ** 2 * cosB * np.log(r - z3))
                             - 1 / (r + a) * (n2 * (y1 * cotB + .5 * a + y1 ** 2 / 2 / (r + a)) - y1 ** 2 / r)
                             + z1 * (r * sinB - y1) / r / (r - z3))
    v2b2 = b2 / 2 / np.pi * ((1 + n2 * cotB ** 2) * Fi - y2 / (r + a) * (n2 * (cotB + y1 / 2 / (r + a)) - y1 / r)
                             - y2 * z1 / r / (r - z3))
    v3b2 = b2 / 2 / np.pi * (-n2 * cotB * (np.log(r + a) - cosB * np.log(r - z3))
                             - y1 / (r + a) * (2 * nu + a / r) + z1 / (r - z3) * (cosB + a / r))
    v1b3 = b3 / 2 / np.pi * (y2 * (r * sinB - y1) * sinB / r / (r - z3))
    v2b3 = b3 / 2 / np.pi * (-y2 ** 2 * sinB / r / (r - z3))
    v3b3 = b3 / 2 / np.pi * (Fi + y2 * (r * cosB + a) * sinB / r / (r - z3))
    return v1b1 + v1b2 + v1b3, v2b1 + v2b2 + v2b3, v3b1 + v3b2 + v3b3


def _ang_setup_fsc(X, Y, b, PA, PB, nu):
    side = PB - PA
    eZ = np.array([0.0, 0.0, 1.0])
    beta = np.arccos(-side @ eZ / np.linalg.norm(side))
    eps = np.finfo(float).eps
    if abs(beta) < eps or abs(np.pi - beta) < eps:
        z = np.zeros_like(X)
        return z, z.copy(), z.copy()
    ey1 = np.array([side[0], side[1], 0.0])
    ey1 /= np.linalg.norm(ey1)
    ey3 = -eZ
    ey2 = np.cross(ey3, ey1)
    A = np.column_stack([ey1, ey2, ey3])          # transformation matrix (columns)
    pts = np.vstack([X - PA[0], Y - PA[1], np.full_like(X, -PA[2])])
    y1A, y2A, _ = A @ pts
    y1AB, y2AB, _ = A @ side
    y1B, y2B = y1A - y1AB, y2A - y2AB
    b1, b2, b3 = A @ b
    I = (beta * y1A) >= 0
    v1A, v2A, v3A = (np.empty_like(X) for _ in range(3))
    v1B, v2B, v3B = (np.empty_like(X) for _ in range(3))
    for mask, bb in ((I, -np.pi + beta), (~I, beta)):
        if mask.any():
            v1A[mask], v2A[mask], v3A[mask] = _ang_dis_disp_surf(y1A[mask], y2A[mask], bb, b1, b2, b3, nu, -PA[2])
            v1B[mask], v2B[mask], v3B[mask] = _ang_dis_disp_surf(y1B[mask], y2B[mask], bb, b1, b2, b3, nu, -PB[2])
    v = np.vstack([v1B - v1A, v2B - v2A, v3B - v3A])
    ue, un, uv = A.T @ v
    return ue, un, uv


def rd_disp_surf(X, Y, P1, P2, P3, P4, opening, nu):
    """Surface displacement of a rectangular tensile dislocation with vertices P1..P4 (z negative below surface)."""
    vn = np.cross(P2 - P1, P4 - P1)
    b = opening * vn / np.linalg.norm(vn)
    out = np.zeros((3, X.size))
    for PA, PB in ((P1, P2), (P2, P3), (P3, P4), (P4, P1)):
        out += np.vstack(_ang_setup_fsc(X, Y, b, PA, PB, nu))
    return out[0], out[1], out[2]


def cdm(X, Y, X0, Y0, depth, omega_x, omega_y, omega_z, ax, ay, az, opening, nu=0.25):
    """Surface displacements (ue, un, uv) and potency DV of a CDM (Nikkhoo et al. 2017), as CDM.m."""
    X = np.asarray(X, float).ravel()
    Y = np.asarray(Y, float).ravel()
    ax, ay, az = 2 * ax, 2 * ay, 2 * az
    R = _rot(omega_x, omega_y, omega_z)
    P0 = np.array([X0, Y0, -depth], float)
    P1 = P0 + ay * R[:, 1] / 2 + az * R[:, 2] / 2
    P2 = P1 - ay * R[:, 1]
    P3 = P2 - az * R[:, 2]
    P4 = P1 - az * R[:, 2]
    Q1 = P0 - ax * R[:, 0] / 2 + az * R[:, 2] / 2
    Q2 = Q1 + ax * R[:, 0]
    Q3 = Q2 - az * R[:, 2]
    Q4 = Q1 - az * R[:, 2]
    R1 = P0 + ax * R[:, 0] / 2 + ay * R[:, 1] / 2
    R2 = R1 - ax * R[:, 0]
    R3 = R2 - ay * R[:, 1]
    R4 = R1 - ay * R[:, 1]
    if any(v[2] > 0 for v in (P1, P2, P3, P4, Q1, Q2, Q3, Q4, R1, R2, R3, R4)):
        raise ValueError("Half-space solution: the CDM must be under the free surface")
    ue, un, uv = np.zeros(X.size), np.zeros(X.size), np.zeros(X.size)
    for size_ok, verts in ((ay != 0 and az != 0, (P1, P2, P3, P4)), (ax != 0 and az != 0, (Q1, Q2, Q3, Q4)),
                           (ax != 0 and ay != 0, (R1, R2, R3, R4))):
        if size_ok:
            e, n, v = rd_disp_surf(X, Y, *verts, opening, nu)
            ue, un, uv = ue + e, un + n, uv + v
    DV = (ax * ay + ax * az + ay * az) * opening
    return ue, un, uv, DV


# ------------------------------------------------------------------ pressurised ellipsoid: volume change / potency
def _rf(x, y, z, r=1e-16):
    xm, ym, zm = x, y, z
    A0 = (x + y + z) / 3
    Q = max(abs(A0 - x), abs(A0 - y), abs(A0 - z)) / (3 * r) ** (1 / 6)
    n, Am = 0, A0
    while abs(Am) <= Q / 4 ** n:
        lam = np.sqrt(xm * ym) + np.sqrt(xm * zm) + np.sqrt(ym * zm)
        Am, xm, ym, zm = (Am + lam) / 4, (xm + lam) / 4, (ym + lam) / 4, (zm + lam) / 4
        n += 1
    X = (A0 - x) / 4 ** n / Am
    Y = (A0 - y) / 4 ** n / Am
    Z = -X - Y
    E2, E3 = X * Y - Z ** 2, X * Y * Z
    return (1 - E2 / 10 + E3 / 14 + E2 ** 2 / 24 - 3 * E2 * E3 / 44) / np.sqrt(Am)


def _rd(x, y, z, r=1e-16):
    xm, ym, zm = x, y, z
    A0 = (x + y + 3 * z) / 5
    Q = max(abs(A0 - x), abs(A0 - y), abs(A0 - z)) / (r / 4) ** (1 / 6)
    n, Am, S = 0, A0, 0.0
    while abs(Am) <= Q / 4 ** n:
        lam = np.sqrt(xm * ym) + np.sqrt(xm * zm) + np.sqrt(ym * zm)
        S += (1 / 4 ** n) / np.sqrt(zm) / (zm + lam)
        Am, xm, ym, zm = (Am + lam) / 4, (xm + lam) / 4, (ym + lam) / 4, (zm + lam) / 4
        n += 1
    X = (A0 - x) / 4 ** n / Am
    Y = (A0 - y) / 4 ** n / Am
    Z = -(X + Y) / 3
    E2 = X * Y - 6 * Z ** 2
    E3 = (3 * X * Y - 8 * Z ** 2) * Z
    E4 = 3 * (X * Y - Z ** 2) * Z ** 2
    E5 = X * Y * Z ** 3
    return ((1 - 3 * E2 / 14 + E3 / 6 + 9 * E2 ** 2 / 88 - 3 * E4 / 22 - 9 * E2 * E3 / 52 + 3 * E5 / 26)
            / 4 ** n / Am ** 1.5 + 3 * S)


def shape_tensor(a1, a2, a3, nu):
    """Eshelby shape tensor components S_iijj (9 values, row order 11,22,33) for semi-axes a1 >= a2 >= a3 > 0."""
    if a1 == a2 == a3:
        d, o = (7 - 5 * nu) / 15 / (1 - nu), (5 * nu - 1) / 15 / (1 - nu)
        return np.array([d, o, o, o, d, o, o, o, d])
    pi = np.pi
    if a1 > a2 > a3 > 0:
        sin_t = np.sqrt(1 - a3 ** 2 / a1 ** 2)
        k = np.sqrt((a1 ** 2 - a2 ** 2) / (a1 ** 2 - a3 ** 2))
        c = 1 / sin_t ** 2
        F = _rf(c - 1, c - k ** 2, c)
        E = F - k ** 2 / 3 * _rd(c - 1, c - k ** 2, c)
        I1 = 4 * pi * a1 * a2 * a3 / (a1 ** 2 - a2 ** 2) / np.sqrt(a1 ** 2 - a3 ** 2) * (F - E)
        I3 = 4 * pi * a1 * a2 * a3 / (a2 ** 2 - a3 ** 2) / np.sqrt(a1 ** 2 - a3 ** 2) * (
            a2 * np.sqrt(a1 ** 2 - a3 ** 2) / a1 / a3 - E)
        I2 = 4 * pi - I1 - I3
        I12 = (I2 - I1) / (a1 ** 2 - a2 ** 2)
        I13 = (I3 - I1) / (a1 ** 2 - a3 ** 2)
        I11 = (4 * pi / a1 ** 2 - I12 - I13) / 3
        I23 = (I3 - I2) / (a2 ** 2 - a3 ** 2)
        I21, I22 = I12, (4 * pi / a2 ** 2 - I23 - I12) / 3
        I31, I32, I33 = I13, I23, (4 * pi / a3 ** 2 - I13 - I23) / 3
    elif a1 == a2 and a2 > a3 > 0:
        I1 = 2 * pi * a1 * a2 * a3 / (a1 ** 2 - a3 ** 2) ** 1.5 * (np.arccos(a3 / a1) - a3 / a1 * np.sqrt(1 - a3 ** 2 / a1 ** 2))
        I2, I3 = I1, 4 * pi - 2 * I1
        I13 = (I3 - I1) / (a1 ** 2 - a3 ** 2)
        I11 = pi / a1 ** 2 - I13 / 4
        I12, I23 = I11, I13
        I22 = pi / a2 ** 2 - I23 / 4
        I21, I31, I32 = I12, I13, I23
        I33 = (4 * pi / a3 ** 2 - 2 * I31) / 3
    elif a1 > a2 and a2 == a3 > 0:
        I2 = 2 * pi * a1 * a2 * a3 / (a1 ** 2 - a3 ** 2) ** 1.5 * (a1 / a3 * np.sqrt(a1 ** 2 / a3 ** 2 - 1) - np.arccosh(a1 / a3))
        I3, I1 = I2, 4 * pi - 2 * I2
        I12 = (I2 - I1) / (a1 ** 2 - a2 ** 2)
        I13 = I12
        I11 = (4 * pi / a1 ** 2 - 2 * I12) / 3
        I21 = I12
        I22 = pi / a2 ** 2 - I21 / 4
        I23, I32, I31 = I22, I22, I13
        I33 = (4 * pi / a3 ** 2 - I31 - I32) / 3
    else:
        raise ValueError("semi-axes must be sorted a1 >= a2 >= a3 > 0")
    k = 1 / 8 / pi / (1 - nu)
    m = (1 - 2 * nu) / 8 / pi / (1 - nu)
    return np.array([3 * k * a1 ** 2 * I11 + m * I1, k * a2 ** 2 * I12 - m * I1, k * a3 ** 2 * I13 - m * I1,
                     k * a1 ** 2 * I21 - m * I2, 3 * k * a2 ** 2 * I22 + m * I2, k * a3 ** 2 * I23 - m * I2,
                     k * a1 ** 2 * I31 - m * I3, k * a2 ** 2 * I32 - m * I3, 3 * k * a3 ** 2 * I33 + m * I3])


def pressurised_ellipsoid(ax, ay, az, nu=0.25):
    """Volume change dV and potencies (DVx, DVy, DVz) per unit p/mu of a uniformly pressurised ellipsoidal cavity
    (pECM.m, Nikkhoo et al. 2017). Returns dV, DV = sum of potencies, and the potencies along x, y, z."""
    mu = 1.0
    lam = 2 * mu * nu / (1 - 2 * nu)
    K = lam + 2 * mu / 3
    p = 1.0
    axes = np.array([ax, ay, az], float)
    order = np.argsort(-axes, kind="stable")
    ai = axes[order]
    S = shape_tensor(*ai, nu)
    Sm = np.array([[S[0] - 1, S[1], S[2]], [S[3], S[4] - 1, S[5]], [S[6], S[7], S[8] - 1]])
    eT = -np.linalg.solve(Sm, np.full(3, p / 3 / K))
    eT[np.sign(eT) != np.sign(p)] = 0
    V = 4 / 3 * np.pi * ax * ay * az
    dV = (eT.sum() - p / K) * V
    DV = dV + p * V / K
    pot = np.zeros(3)
    pot[order] = V * eT
    return dV, DV, pot
