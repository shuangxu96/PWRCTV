"""
Python version of PWRCTV (TGRS 2024).
"""

import time

import numpy as np
from numba import njit

from .operators import _EPS, diff_x, diff_y, psf2otf, imcorrfilter, svdsecon

MN_DTYPE = np.int64


def _da_gather_indices(H, W):
    MN = H * W
    x = np.arange(H)
    y = np.arange(W)
    # element (a, b) of the 2D grid holds the F-flat index of the neighbour;
    # ravel(order='F') puts that value at position a + H*b (F-flat of (a,b)),
    # so idxp[k] = F-flat of the neighbour of pixel k.
    idxp_x = (((x + 1) % H)[:, None] + H * y[None, :]).ravel(order="F")
    idxp_y = (x[:, None] + H * ((y + 1) % W)[None, :]).ravel(order="F")
    idxm_x = (((x - 1) % H)[:, None] + H * y[None, :]).ravel(order="F")
    idxm_y = (x[:, None] + H * ((y - 1) % W)[None, :]).ravel(order="F")
    return (idxp_x.astype(MN_DTYPE), idxp_y.astype(MN_DTYPE),
            idxm_x.astype(MN_DTYPE), idxm_y.astype(MN_DTYPE))


# --------------------------------------------------------------------------
# PWRCTV kernels (TGRS 2024) -- two-stage rho update, faithfully preserved
# --------------------------------------------------------------------------
# Layout contract: every (MN, B) matrix is passed as a FLAT C-contiguous view
# (reshape(-1)); every (MN, r) auxiliary array likewise.  U2d stays the
# F-contiguous 2D view of the U tensor (gather indices + FFT view depend on
# it); U_c is its C-contiguous copy used only for BLAS matmuls.  All loops
# are fused in-place, so no 33MB temporaries are allocated per iteration.

@njit(cache=True, fastmath=False)
def _pwr_a(U2d, V2d, Ef, Sf, M3f, M1f, M2f, mu, Df,
           W1, W2, rho1, rho2,
           idxp_x, idxp_y, idxm_x, idxm_y,
           F1, F2, A2f, numer, U_c, MN, r, B, tau0, tau1):
    """PWRCTV pre-FFT: W*=rho, weighted soft thresholds, adjoint diff, numer."""
    inv_mu = 1.0 / mu
    nMB = MN * B
    for i in range(MN):
        for j in range(r):
            W1[i, j] = W1[i, j] * rho1[i, j]
            W2[i, j] = W2[i, j] * rho2[i, j]
    for i in range(nMB):
        A2f[i] = Df[i] - Ef[i] - Sf[i] + M3f[i] * inv_mu
    for i in range(MN):
        for j in range(r):
            U_c[i, j] = U2d[i, j]
    for i in range(MN):
        ipx = idxp_x[i]
        ipy = idxp_y[i]
        for j in range(r):
            k = i * r + j
            dx = U2d[ipx, j] - U2d[i, j]
            dy = U2d[ipy, j] - U2d[i, j]
            z1 = dx + M1f[k] * inv_mu
            z2 = dy + M2f[k] * inv_mu
            a1 = np.abs(z1) - tau0 * inv_mu * W1[i, j]
            a2 = np.abs(z2) - tau1 * inv_mu * W2[i, j]
            if np.isnan(a1):
                a1 = -np.inf
            if np.isnan(a2):
                a2 = -np.inf
            F1[k] = np.sign(z1) * np.maximum(a1, 0.0)
            F2[k] = np.sign(z2) * np.maximum(a2, 0.0)
    A2_2d = A2f.reshape(MN, B)
    temp = A2_2d @ V2d
    for i in range(MN):
        imx = idxm_x[i]
        imy = idxm_y[i]
        for j in range(r):
            k = i * r + j
            kmx = imx * r + j
            kmy = imy * r + j
            g1 = F1[k] - M1f[k] * inv_mu
            g2 = F2[k] - M2f[k] * inv_mu
            rhs = (F1[kmx] - M1f[kmx] * inv_mu - g1
                   + F2[kmy] - M2f[kmy] * inv_mu - g2)
            numer[k] = rhs + temp[i, j]


@njit(cache=True, fastmath=False)
def _pwr_b(U2d, A2f, Ef, Sf, M3f, mu, Df,
           F1, F2,
           idxp_x, idxp_y,
           leq1, leq2, leq3f, Bmat, Cmat, V2d,
           beta, lam, normD,
           U_c, V_T, UVt, MN, r, B):
    """PWRCTV post-FFT: V via eig-based svdecon, E/S closed forms, residuals.

    Multiplier/mu update is NOT performed here; the caller applies it only in
    the else branch of the stopping rule (MATLAB semantics -- the iteration
    that flips update_rho skips the multiplier update).
    """
    inv_mu = 1.0 / mu
    nMB = MN * B
    for i in range(MN):
        for j in range(r):
            U_c[i, j] = U2d[i, j]
    A2_2d = A2f.reshape(MN, B)
    Bmat[:, :] = (A2_2d.T) @ U_c
    # svdecon(Bmat) with m > n (eig-based branch)
    Cmat[:, :] = Bmat.T @ Bmat
    Dw, Vw = np.linalg.eig(Cmat)
    dw = np.abs(np.real(Dw))
    ix = np.argsort(-dw)
    Vw = np.ascontiguousarray(np.real(Vw[:, ix]))
    sv = np.sqrt(dw[ix])
    uu = Bmat @ Vw
    for j in range(r):
        for i in range(B):
            uu[i, j] = uu[i, j] / sv[j]
    V2d[:, :] = uu @ Vw.T
    # V_T must be copied from the UPDATED V (numpy: UVt = U @ V_new.T)
    for i in range(r):
        for j in range(B):
            V_T[i, j] = V2d[j, i]
    UVt[:, :] = U_c @ V_T
    # fused E / S / leq3
    inv_den = 1.0 / (2.0 * beta + mu)
    scale = mu * inv_den
    if lam > 100:
        for i in range(MN):
            base = i * B
            for j in range(B):
                k = base + j
                e_new = (A2f[k] + Ef[k] - UVt[i, j]) * scale
                Ef[k] = e_new
                Sf[k] = 0.0
                leq3f[k] = Df[k] - UVt[i, j] - e_new
    else:
        for i in range(MN):
            base = i * B
            for j in range(B):
                k = base + j
                t = A2f[k] + Ef[k] - UVt[i, j]
                e_new = t * scale
                Ef[k] = e_new
                t2 = t + Sf[k] - e_new
                arg = np.abs(t2) - lam * inv_mu
                if np.isnan(arg):
                    arg = -np.inf
                s_new = np.sign(t2) * np.maximum(arg, 0.0)
                Sf[k] = s_new
                leq3f[k] = Df[k] - UVt[i, j] - e_new - s_new
    s1 = 0.0
    s2 = 0.0
    for i in range(MN):
        ipx = idxp_x[i]
        ipy = idxp_y[i]
        for j in range(r):
            k = i * r + j
            l1 = (U2d[ipx, j] - U2d[i, j]) - F1[k]
            l2 = (U2d[ipy, j] - U2d[i, j]) - F2[k]
            leq1[k] = l1
            leq2[k] = l2
            s1 += l1 * l1
            s2 += l2 * l2
    s3 = 0.0
    for i in range(nMB):
        s3 += leq3f[i] * leq3f[i]
    stopC1 = np.sqrt(s1) / normD
    stopC2 = np.sqrt(s2) / normD
    stopC3 = np.sqrt(s3) / normD
    return stopC1, stopC2, stopC3


@njit(cache=True, fastmath=False)
def _pwr_mult(M1f, M2f, M3f, mu, leq1, leq2, leq3f, rho, max_mu, MN, r, B):
    """Multiplier and mu update (else branch of the PWRCTV stopping rule)."""
    for i in range(MN):
        for j in range(r):
            k = i * r + j
            M1f[k] += mu * leq1[k]
            M2f[k] += mu * leq2[k]
    nMB = MN * B
    for i in range(nMB):
        M3f[i] += mu * leq3f[i]
    return min(max_mu, mu * rho)


def pwrctv(Nhsi, Pan, beta=0.5, lam=1.0, tau=(0.01, 0.01), r=6, q=10):
    """Numba backend of PWRCTV.  Returns (output_image, U, V)."""
    Nhsi = np.asarray(Nhsi, dtype=float)
    Pan = np.asarray(Pan, dtype=float)
    tau = np.asarray(tau, dtype=float).ravel()

    tol = 1e-5
    tol_rho = 100.0
    max_iter = 200
    rho = 1.5
    mu0 = 1.0
    max_mu = 1e6
    eps = 1e-4

    H, W, nband = Nhsi.shape
    sizeU = (H, W, r)
    D = Nhsi.reshape(H * W, nband, order="F")
    Dc = np.ascontiguousarray(D)
    MN = H * W

    norm_two = svdsecon(D, 1)[0]
    mu = mu0 / max(norm_two, _EPS)
    normD = np.linalg.norm(D, "fro")

    Eny_x = np.abs(psf2otf(np.array([[1.0], [-1.0]]), sizeU)) ** 2
    Eny_y = np.abs(psf2otf(np.array([[1.0, -1.0]]), sizeU)) ** 2
    determ = Eny_x + Eny_y

    u, s, vh = np.linalg.svd(D, full_matrices=False)
    U = u[:, :r] @ np.diag(s[:r])
    V = vh[:r, :].T

    S = np.zeros((MN, nband))
    E = np.zeros((MN, nband))
    M1 = np.zeros(MN * r)
    M2 = np.zeros(MN * r)
    M3 = np.zeros((MN, nband))
    W3 = 1.0

    G1 = diff_x(Pan, None)
    G2 = diff_y(Pan, None)

    W1 = np.repeat(np.abs(G1)[:, :, None], r, axis=2)
    W2 = np.repeat(np.abs(G2)[:, :, None], r, axis=2)
    W1 = (1.0 - W1) ** q
    W2 = (1.0 - W2) ** q
    # NOTE: numpy PWRCTV reads the weights as ravel(order='F').reshape(MN, r,
    # order='F'); element (i, j) with i = x + y*H equals the (x, y, j) entry.
    W1 = np.ascontiguousarray(W1.ravel(order="F").reshape(MN, r, order="F"))
    W2 = np.ascontiguousarray(W2.ravel(order="F").reshape(MN, r, order="F"))

    rho_XY1 = np.ones((H, W, r), order="F")
    rho_XY2 = np.ones((H, W, r), order="F")
    rho1 = rho_XY1.reshape(MN, r, order="F")
    rho2 = rho_XY2.reshape(MN, r, order="F")
    update_rho = False

    # ---- buffers ----
    U3 = np.zeros((H, W, r), order="F")
    U2d = U3.reshape(MN, r, order="F")
    U2d[:] = U
    U_c = np.zeros((MN, r))
    Vk = np.ascontiguousarray(V)
    V_T = np.zeros((r, nband))
    UVt = np.zeros((MN, nband))
    E_k = np.zeros((MN, nband))
    S_k = np.zeros((MN, nband))
    M3_k = np.zeros((MN, nband))
    M1_k = np.zeros(MN * r)
    M2_k = np.zeros(MN * r)
    A2 = np.zeros((MN, nband))
    leq3 = np.zeros((MN, nband))
    F1 = np.zeros((MN, r))
    F2 = np.zeros((MN, r))
    leq1 = np.zeros((MN, r))
    leq2 = np.zeros((MN, r))
    Bmat = np.zeros((nband, r))
    Cmat = np.zeros((r, r))
    numer2d = np.zeros((MN, r))
    numer3 = numer2d.reshape(H, W, r, order="F")
    idxp_x, idxp_y, idxm_x, idxm_y = _da_gather_indices(H, W)
    Df = Dc.reshape(-1)
    Ef = E_k.reshape(-1)
    Sf = S_k.reshape(-1)
    M3f = M3_k.reshape(-1)
    A2f = A2.reshape(-1)
    leq3f = leq3.reshape(-1)
    F1f = F1.reshape(-1)
    F2f = F2.reshape(-1)
    leq1f = leq1.reshape(-1)
    leq2f = leq2.reshape(-1)
    numerf = numer2d.reshape(-1)

    it = 0
    for it in range(1, max_iter + 1):
        _pwr_a(U2d, Vk, Ef, Sf, M3f, M1_k, M2_k, mu, Df,
               W1, W2, rho1, rho2,
               idxp_x, idxp_y, idxm_x, idxm_y,
               F1f, F2f, A2f, numerf, U_c, MN, r, nband, tau[0], tau[1])
        U3[:] = np.real(np.fft.ifftn(np.fft.fftn(numer3) / (determ + 1.0 + eps)))
        if update_rho:
            # MATLAB bug faithfully preserved: only band r is written, and
            # the OTHER bands are zeroed (numpy PWRCTV: rho_XY1 = zeros(...)).
            UG1 = diff_x(U3, None)
            UG2 = diff_y(U3, None)
            rho_XY1[:] = 0.0
            rho_XY2[:] = 0.0
            rho_XY1[:, :, r - 1] = np.abs(imcorrfilter(UG1[:, :, r - 1], G1, 2, 2))
            rho_XY2[:, :, r - 1] = np.abs(imcorrfilter(UG2[:, :, r - 1], G2, 2, 2))
        c1, c2, c3 = _pwr_b(
            U2d, A2f, Ef, Sf, M3f, mu, Df,
            F1f, F2f,
            idxp_x, idxp_y,
            leq1f, leq2f, leq3f, Bmat, Cmat, Vk,
            beta, lam, normD,
            U_c, V_T, UVt, MN, r, nband)

        if c1 < tol and c2 < tol and c3 < tol and update_rho:
            break
        elif c1 < tol_rho * tol and c2 < tol_rho * tol and c3 < tol_rho * tol \
                and not update_rho:
            update_rho = True
        else:
            mu = _pwr_mult(M1_k, M2_k, M3f, mu, leq1f, leq2f, leq3f,
                           rho, max_mu, MN, r, nband)

    output = (U2d @ Vk.T).reshape(H, W, nband, order="F")
    return output, U2d.copy(), Vk.copy()


def warm_up(mn=65536, r=4, nband=63, verbose=True):
    """Pre-compile the PWRCTV JIT kernels so the first real case runs warm.

    Uses production-size dummy buffers (same dtypes / layouts as the real
    calls) so numba compiles exactly the signatures used at runtime.  No
    numerical result is used.  Also writes the numba disk cache
    (cache=True), so later processes skip compilation too.
    """
    import time as _t
    t0 = _t.perf_counter()
    rng = np.random.default_rng(20261008)
    B = nband
    MN = mn
    M = N = int(round(np.sqrt(MN)))

    U2d = np.asfortranarray(rng.standard_normal((MN, r)))
    Vk = np.ascontiguousarray(rng.standard_normal((B, r)))
    Df = rng.standard_normal((MN, B)).reshape(-1)
    Ef = np.zeros(MN * B); Sf = np.zeros(MN * B); M3f = np.zeros(MN * B)
    A2f = np.zeros(MN * B); leq3f = np.zeros(MN * B)
    M1f = np.zeros(MN * r); M2f = np.zeros(MN * r)
    leq1f = np.zeros(MN * r); leq2f = np.zeros(MN * r)
    numerf = np.zeros(MN * r)
    U_c = np.zeros((MN, r)); V_T = np.zeros((r, B)); UVt = np.zeros((MN, B))
    Bmat = np.zeros((B, r)); Cmat = np.zeros((r, r))
    F1f = np.zeros(MN * r); F2f = np.zeros(MN * r)
    idxp_x, idxp_y, idxm_x, idxm_y = _da_gather_indices(M, N)
    mu = 1e-3
    # NOTE: real fast_pwrctv passes rho1/rho2 as F-order (rho_XY1.reshape
    # (MN, r, order='F')); layout is part of the numba signature, so the
    # warm-up buffers must match it or the first real call recompiles.
    W1 = np.ones((MN, r)); W2 = np.ones((MN, r))
    rho1 = np.asfortranarray(np.ones((MN, r)))
    rho2 = np.asfortranarray(np.ones((MN, r)))

    _pwr_a(U2d, Vk, Ef, Sf, M3f, M1f, M2f, mu, Df,
           W1, W2, rho1, rho2, idxp_x, idxp_y, idxm_x, idxm_y,
           F1f, F2f, A2f, numerf, U_c, MN, r, B, 0.4, 0.4)
    _pwr_b(U2d, A2f, Ef, Sf, M3f, mu, Df, F1f, F2f,
           idxp_x, idxp_y, leq1f, leq2f, leq3f, Bmat, Cmat, Vk,
           100.0, 1.0, 1.0, U_c, V_T, UVt, MN, r, B)
    _pwr_mult(M1f, M2f, M3f, mu, leq1f, leq2f, leq3f, 1.5, 1e6, MN, r, B)

    dt = _t.perf_counter() - t0
    if verbose:
        print(f"[warm_up] PWRCTV JIT kernels compiled in {dt:.2f}s "
              f"(MN={MN}, r={r}, B={B})")
    return dt
