#!/share/simulation/Anaconda3-2022.05/bin/python
from simpletraj.dcd.dcd import DCDReader
import numpy as np
import warnings
import os
# from numba import guvectorize, float64
import matplotlib.pyplot as plt


filenames = [i for i in os.listdir(os.getcwd()) if i.startswith("particle") and i.endswith("dcd")]
for i in filenames:
    if "particle_100000Period.dcd" in i:
        period = int(i.split("_")[-1].split("Period")[0])
        print("Processing " + i + " ...")
        print("Period: " + str(period))

        chain_ring = int(os.getcwd().split("_")[1]) + 7 * int(os.getcwd().split("_")[-1])
        chain_N = int(os.getcwd().split("_")[1])

        print("chain_ring: " + str(chain_ring))
        print("chain_N: " + str(chain_N))
        print("")

        traj_raw = np.array([_.copy() for _ in DCDReader(i)], dtype=np.float32)
        traj = traj_raw.reshape(traj_raw.shape[0], -1, chain_ring, 3)
        pos = traj[:, :, :chain_N].reshape(traj.shape[0], -1, 3)
        pos = pos.reshape(-1, 1, pos.shape[-2], 3)[:, :, :]
        mode = normal_modes(pos, 1).reshape(-1, 3)
        print(mode.shape)

        # 计算归一化的ACF
        acf = np.zeros(mode.shape[0])
        times = np.ones(mode.shape[0])
        for j in range(mode.shape[0]):
            for jj in range(j + 1, mode.shape[0]):
                acf[jj-j] += np.dot(mode[jj], mode[j])
                times[jj-j] += 1
        print(acf.shape)
        acf = acf / times

        t = period * 0.005
        time = np.arange(mode.shape[0]) * t
        plt.plot(time[1:], acf[1:])
        plt.xlabel("Time")
        plt.ylabel("ACF")
        plt.xscale('log')
        plt.show()
