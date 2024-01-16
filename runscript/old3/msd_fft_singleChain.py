#!/share/simulation/Anaconda3-2022.05/bin/python
import numpy as np
import sys
import os
from simpletraj.dcd.dcd import DCDReader


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


filenames = [i for i in os.listdir(os.getcwd()) if i.startswith("particle") and i.endswith("dcd")]

for i in filenames:
    period = int(i.split("_")[-1].split("Period")[0])
    print("Processing " + i + " ...")
    print("Period: " + str(period))

    # chain_ring = int(os.getcwd().split("_")[1]) + 7 * int(os.getcwd().split("_")[-1])
    # chain_N = int(os.getcwd().split("_")[1])

    # print("chain_ring: " + str(chain_ring))
    # print("chain_N: " + str(chain_N))
    print("")

    traj_raw = np.array([_.copy() for _ in DCDReader(i)], dtype=np.float32)
    traj = traj_raw.reshape(traj_raw.shape[0], -1, 50, 3)
    traj_final = traj[:, :, :].reshape(traj.shape[0], -1, 3)
    u = traj_final.swapaxes(0, 1)

    msd = np.asarray([msd_fft(_) for _ in u]).mean(axis=0)

    dt = 0.001 * period
    time = dt * np.arange(msd.shape[0])

    np.savetxt(f'msd{period}.dat', np.c_[time, msd], fmt='%.8f')
