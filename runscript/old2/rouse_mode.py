#!/share/simulation/Anaconda3-2022.05/bin/python
import warnings
import sys
import os
from simpletraj.dcd.dcd import DCDReader
import numpy as np
from numba import guvectorize
from numba import float64
import matplotlib.pyplot as plt

fontsize = 16
linesize = 2
def make_figure():
    plt.gcf().set_size_inches(8, 6)
    plt.tick_params(labelsize=fontsize, width=linesize)
    ax = plt.gca()
    ax.spines['bottom'].set_linewidth(linesize)
    ax.spines['left'].set_linewidth(linesize)
    ax.spines['top'].set_linewidth(linesize)
    ax.spines['right'].set_linewidth(linesize)
    plt.xlabel("r", fontdict={'size': fontsize})
    plt.ylabel("g(r)", fontdict={'size': fontsize})
    #plt.tick_params(direction='in')
    plt.tick_params(axis='both', which='both', direction='in')
    ls = plt.legend(fontsize=fontsize, loc='lower left')
    ls.get_frame().set_linewidth(linesize)
    ls.get_frame().set_edgecolor('black')


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
    modes = np.atleast_1d(np.asarray(modes)) - 1 / 2 if modes is not None else np.arange(1, chain_length + 1) - 1 / 2
    if 0 in modes:
        warnings.warn("Please use UNWRAPPED coordinates to calculate the 0th mode!")
    s = np.expand_dims(np.arange(1, chain_length + 1), 0) * np.pi / chain_length
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

path_raw = os.getcwd()
all_items = os.listdir(path_raw)
path = [os.path.join(path_raw, _) for _ in all_items if os.path.isdir(os.path.join(path_raw, _)) and (_.startswith("chain") or _.startswith("pure"))]
path = path[::-1]

dtau = float(input("input dtau or (default:0.005):") or "0.005")
N = int(input("input N or (default:50):") or "50")
for j in path:
    os.chdir(j)
    filenames = [i for i in os.listdir(os.getcwd()) if i.startswith("particle") and i.endswith("dcd")]
    for i in filenames:
        if "particle_100000Period" in i or "particle_50000Period" in i:
        # if "particle_20000Period" in i or "particle_50000Period" in i:
            period = int(i.split("_")[-1].split("Period")[0])
            traj_raw = np.array([_.copy() for _ in DCDReader(i)], dtype=np.float32)
            print("")
            print("Located at: " + j)
            print("Processing " + i + " ...")
            print("Period: " + str(period))
            if "_" in os.getcwd():
                chain_ring = int(os.getcwd().split("_")[1]) + 7 * int(os.getcwd().split("_")[-1])
                chain_N = int(os.getcwd().split("_")[1])

                print("chain_ring: " + str(chain_ring))
                print("chain_N: " + str(chain_N))
                print("")

                traj = traj_raw.reshape(traj_raw.shape[0], -1, chain_ring, 3) # n_frames, n_chains, chain&ring_length, n_dimensions
                pos = traj[:, :, :chain_N] # n_frames, n_chains, chain_length, n_dimensions
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

            for i in np.array([1,2,3,10])-1:
                acf_final = acf[:, i] / acf[0, i]
                if acf_final.min() > 0:
                    idx = acf_final.shape[0]
                else:
                    idx = np.argwhere(acf_final < 0)[0, 0]
                plt.plot(time[:idx], acf_final[:idx], label=f"p={i+1}")
            make_figure()
            plt.xlabel("Time")
            plt.ylabel(r"$\langle X_p(t)X_p(0) \rangle / \langle X_{p}^2 \rangle$")
            plt.xscale('log')
            # plt.xlim(10**2,10**7.5)
            plt.title(rf"{os.getcwd().split('/')[-1]}", fontdict={'size':fontsize})
            plt.savefig(rf"acf_{os.getcwd().split('/')[-1]}.png", dpi=600, bbox_inches='tight')
            # os.system("mv *.png ../")
            plt.show()
