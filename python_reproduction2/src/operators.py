"""
helper functions for PWRCTV
"""

from __future__ import annotations

import numpy as np

_EPS = np.finfo(float).eps


def _as_3d(X, sizeX):
    if X.ndim in (1, 2):
        return X.reshape(tuple(sizeX), order="F")
    return X


def diff_x(X, sizeX=None):
    """Forward difference along dim 1 (rows), circular boundary.

    MATLAB: d(1:end-1) = X(2:end) - X(1:end-1); d(end) = X(1) - X(end).
    Returns a column vector when sizeX is given, else a full-size array.
    """
    if sizeX is not None:
        X = _as_3d(X, sizeX)
        d = np.empty_like(X)
        d[:-1] = X[1:] - X[:-1]
        d[-1] = X[0] - X[-1]
        return d.reshape(-1, order="F")
    d = np.empty_like(X)
    d[:-1] = X[1:] - X[:-1]
    d[-1] = X[0] - X[-1]
    return d


def diff_y(X, sizeX=None):
    """Forward difference along dim 2 (columns), circular boundary."""
    if sizeX is not None:
        X = _as_3d(X, sizeX)
        d = np.empty_like(X)
        d[:, :-1] = X[:, 1:] - X[:, :-1]
        d[:, -1] = X[:, 0] - X[:, -1]
        return d.reshape(-1, order="F")
    d = np.empty_like(X)
    d[:, :-1] = X[:, 1:] - X[:, :-1]
    d[:, -1] = X[:, 0] - X[:, -1]
    return d


def psf2otf(psf, outSize):
    """PSF -> OTF by zero-padding to outSize[:2] and repmat over dim 3."""
    pad = np.zeros((int(outSize[0]), int(outSize[1])), dtype=float)
    pad[: psf.shape[0], : psf.shape[1]] = np.asarray(psf, dtype=float)
    otf2 = np.fft.fft2(pad)
    return np.repeat(otf2[:, :, np.newaxis], int(outSize[2]), axis=2)

def svdsecon(X, k):
    """Port of PRTV_plus_code/utils/svdsecon.m: largest k singular values."""
    X = np.asarray(X, dtype=float)
    m, n = X.shape
    assert k <= m and k <= n, "k needs to be smaller than size(X,1) and size(X,2)"
    if m <= n:
        C = X @ X.T
        d = np.linalg.eigvalsh(C)
        s = np.sqrt(np.abs(d[-k:][::-1]))
    else:
        C = X.T @ X
        d = np.linalg.eigvalsh(C)
        s = np.sqrt(np.abs(d[-k:][::-1]))
    return s


def boxfilter(imSrc, r1, r2):
    """Port of the cumulative-sum box filter inside imcorrfilter.m."""
    imSrc = np.asarray(imSrc, dtype=float)
    hei, wid = imSrc.shape[:2]
    imDst = np.zeros_like(imSrc)

    imCum = np.cumsum(imSrc, axis=0)
    imDst[0 : r1 + 1, :] = imCum[r1 : 2 * r1 + 1, :]
    imDst[r1 + 1 : hei - r1, :] = imCum[2 * r1 + 1 : hei, :] - imCum[0 : hei - 2 * r1 - 1, :]
    imDst[hei - r1 : hei, :] = np.tile(imCum[hei - 1, :], (r1, 1)) - imCum[hei - 2 * r1 - 1 : hei - r1 - 1, :]

    imCum = np.cumsum(imDst, axis=1)
    imDst[:, 0 : r2 + 1] = imCum[:, r2 : 2 * r2 + 1]
    imDst[:, r2 + 1 : wid - r2] = imCum[:, 2 * r2 + 1 : wid] - imCum[:, 0 : wid - 2 * r2 - 1]
    imDst[:, wid - r2 : wid] = np.tile(imCum[:, wid - 1, None], (1, r2)) - imCum[:, wid - 2 * r2 - 1 : wid - r2 - 1]

    imDst = imDst / ((2 * r1 + 1) * (2 * r2 + 1))
    return imDst


def imcorrfilter(X, Y, w1, w2):
    """Port of PRTV_plus_code/utils/imcorrfilter.m (local correlation)."""
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    mean_X = boxfilter(X, w1, w2)
    mean_Y = boxfilter(Y, w1, w2)
    var_X = boxfilter((X - mean_X) ** 2, w1, w2)
    var_Y = boxfilter((Y - mean_Y) ** 2, w1, w2)
    cov_XY = boxfilter((X - mean_X) * (Y - mean_Y), w1, w2)
    rho_XY = cov_XY / (np.sqrt(var_X) * np.sqrt(var_Y))
    return rho_XY
