import numpy as np
import warnings
from numba import guvectorize, float64


# msd_fft-------------
def autocorrFFT(X):
    M = X.shape[0]
    Fx = np.fft.rfft(X.T[0], 2 * M)
    Fy = np.fft.rfft(X.T[1], 2 * M)
    Fz = np.fft.rfft(X.T[2], 2 * M)
    corr = abs(Fx) ** 2 + abs(Fy) ** 2 + abs(Fz) ** 2
    res = np.fft.irfft(corr)
    res = (res[:M]).real
    return res


def msd_fft(X):
    M = X.shape[0]
    D = np.square(X).sum(axis=1)
    D = np.append(D, 0)
    S2 = autocorrFFT(X)
    S2 = S2 / np.arange(M, 0, -1)
    Q = 2 * D.sum()
    S1 = np.zeros(M)
    for m in range(M):
        Q = Q - D[m - 1] - D[M - m]
        S1[m] = Q / (M - m)
    return S1 - 2 * S2


# normal p modes---------------
@guvectorize([(float64[:, :], float64[:, :], float64[:, :])], '(n,p),(p,m)->(n,m)', target='parallel')  # target='cpu','gpu'
def batch_dot(a, b, ret):  # much more faster than np.tensordot or np.einsum
    r"""Vectorized universal function.
    :param a: np.ndarray, (...,N,P)
    :param b: np.ndarray, (...,P,M)
    axes will be assigned automatically to last 2 axes due to the signatures.
    this functions is actual np.einsum('...mp,....pn->...mn', a, b), or
    np.matmul(a, b). But this is much faster.
    :param ret: np.ndarray, generated automatically by guvectorize
    :return: np.ndarray, results. (...,N,M)
    """
    for i in range(ret.shape[0]):
        for j in range(ret.shape[1]):
            tmp = 0.
            for k in range(a.shape[1]):
                tmp += a[i, k] * b[k, j]
            ret[i, j] = tmp


def normal_modes(pos: np.ndarray, modes=None) -> np.ndarray:
    chain_length = pos.shape[-2]
    # modes = np.atleast_1d(np.asarray(modes)) - 1 / 2 if modes is not None else np.arange(1, chain_length + 1) - 1 / 2
    modes = np.atleast_1d(np.asarray(modes)) if modes is not None else np.arange(1, chain_length + 1)
    if 0 in modes:
        warnings.warn("Please use UNWRAPPED coordinates to calculate the 0th mode!")
    s = np.expand_dims(np.arange(1, chain_length + 1) - 1 / 2, 0) * np.pi / chain_length
    factors = np.sqrt(2 / chain_length) * np.cos(np.expand_dims(modes, -1) * s)
    # s = np.expand_dims(np.arange(1, chain_length + 1), 0) * np.pi / chain_length
    # factors = 1 / chain_length * np.cos(np.expand_dims(modes, -1) * s)
    """
    :param pos: np.ndarray, positions in (n_frames (optional), n_chains (optional), chain_length, n_dimensions)
    :param modes: iterable, modes to calculate. mode 1 ~ chain_length are calculated by default.
    :return: np.ndarray, normal modes (..., n_modes, n_dimensions)
    """
    return batch_dot(factors, pos)


# acf_fft-------------
def mat_ac(x, axes=None):
    r"""Matrix autocorrelation function.
    :param x: np.ndarray -> (n_frames, ...) of input
    :param axes: tuple or int, axes that summing up.
    :return: np.ndarray -> (n_frames, ...) of output
    :raises: ValueError, if axes contains 0.
    """
    fft = np.fft.rfft
    ifft = np.fft.irfft
    if np.issubdtype(x.dtype, np.dtype(complex)):
        fft = np.fft.fft
        ifft = np.fft.ifft
    n = x.shape[0]
    s = 2 * n  # 2 * n - 1 is fine.
    norm = np.arange(n, 0, -1).reshape(n, *[1] * (x.ndim - 1))
    if axes is None:
        return ifft(abs(fft(x, axis=0, n=s)) ** 2, axis=0, n=s)[:n].real / norm
    else:
        axes = np.atleast_1d(np.asarray(axes, dtype=int))
        if 0 in axes:
            raise ValueError("The 1st axis should be time axis!")
        norm = norm.reshape(n, *[1] * (x.ndim - 1 - axes.size))
        return ifft(np.sum(abs(fft(x, axis=0, n=s)) ** 2, axis=tuple(axes)), axis=0, n=s)[:n].real / norm
