#!/share/simulation/Anaconda3-2022.05/bin/python
import warnings
import sys
import os
from simpletraj.dcd.dcd import DCDReader
import numpy as np
from numba import guvectorize
from numba import float64
import matplotlib.pyplot as plt


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
    # chain_length = 50
    # modes = np.atleast_1d(np.asarray(modes)) - 1 / 2 if modes is not None else np.arange(1, chain_length + 1) - 1 / 2
    modes = np.atleast_1d(np.asarray(modes)) if modes is not None else np.arange(1, chain_length + 1)
    if 0 in modes:
        warnings.warn("Please use UNWRAPPED coordinates to calculate the 0th mode!")
    # s = np.expand_dims(np.arange(1, chain_length + 1), 0) * np.pi / chain_length
    s = np.expand_dims(np.arange(1, chain_length + 1) - 1 / 2, 0) * np.pi / chain_length
    factors = 1 / chain_length * np.cos(np.expand_dims(modes, -1) * s)
    """
    :param pos: np.ndarray, positions in (n_frames (optional), n_chains (optional), chain_length, n_dimensions)
    :param modes: iterable, modes to calculate. mode 1 ~ chain_length are calculated by default.
    :return: np.ndarray, normal modes (..., n_modes, n_dimensions)
    """
    return batch_dot(factors, pos)


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

# path_raw = os.getcwd()
# all_items = os.listdir(path_raw)
# path = [os.path.join(path_raw, _) for _ in all_items if os.path.isdir(os.path.join(path_raw, _)) and (_.startswith("chain") or _.startswith("pure"))]
# path = path[::-1]

dtau = float(input("input dtau or (default:0.005):") or "0.005")
N = int(input("input N or (default:50):") or "50")
filenames = [i for i in os.listdir(os.getcwd()) if i.startswith("particle") and i.endswith("dcd")]

for i in filenames:
    traj_raw = np.array([_.copy() for _ in DCDReader(i)], dtype=np.float32)
    period = int(i.split("_")[-1].split("Period")[0])

    print("Processing " + i + " ...")
    print("Period: " + str(period))
    print("")

    if "_" in os.getcwd():
        chain_ring = int(os.getcwd().split("_")[1]) + 7 * int(os.getcwd().split("_")[-1])
        chain_N = int(os.getcwd().split("_")[1])

        print("chain_ring: " + str(chain_ring))
        print("chain_N: " + str(chain_N))
        print("")

        traj = traj_raw.reshape(traj_raw.shape[0], -1, chain_ring, 3) # n_frames, n_chains, chain&ring_length, n_dimensions
        pos = traj[:, :, chain_N:] # n_frames, n_chains, chain_length, n_dimensions
    else:
        print("Processing pure polymer...")
        print(f"chain_N: {N}")
        pos = traj_raw.reshape(traj_raw.shape[0], -1, N, 3) # n_frames, n_chains, chain_length, n_dimensions
    mode = normal_modes(pos).mean(axis=1) # n_frames, all mode
    acf = mat_ac(mode, -1) #n_frames, all mode
    # print(acf.shape)
    
    t = period * dtau
    time = np.arange(mode.shape[0]) * t
    np.savetxt(f'acf{period}.dat', np.c_[time, acf], fmt='%.8f')
