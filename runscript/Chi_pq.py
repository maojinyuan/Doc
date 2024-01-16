#!/usr/bin/env python
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os
from simpletraj.dcd.dcd import DCDReader
import numpy as np
import warnings
from numba import guvectorize
from numba import float64
import msd_acf
from xml_io import uxml
from list_dir import get_current_folder_info

# get_current_folder_info() #return current_folder, father_folder, dcd_file, gala_file
cf, ff, dcdfile, galafile = get_current_folder_info()[:-1]

# dtau = [line.split()[-1] for line in open(galafile[-1], "r") if "app.setDt" in line and "#" not in line][-1]
# dtau = float(dtau.split("(")[-1].split(")")[0])
# print("dtau: " + str(dtau))

xml = uxml.read("out.xml")
N = xml.natoms
# print("N: " + str(N))
# print("----------------------")

# for i in dcdfile[2:3]:
for i in dcdfile:
    print("Processing " + i + " ...")
    slice_dcd = DCDReader(i)
    dtau = 0.01
    print(f"dtau = {dtau}")
    traj_raw = np.array([_.copy() for _ in slice_dcd], dtype=np.float32)
    period = slice_dcd.skip_timestep
    natoms = slice_dcd.numatoms
    if "chain" in cf and "ring" in cf:
        ring_N = os.path.basename(cf)
        chain_ring = int(ring_N.split("_")[1]) + 7 * int(ring_N.split("_")[-1])
        chain_N = int(ring_N.split("_")[1])
        print(f"N + ring = {chain_ring}, N = {chain_N}")

        traj = traj_raw.reshape(traj_raw.shape[0], -1, chain_ring, 3)  # n_frames, n_chains, chain&ring_length, n_dimensions
        pos = traj[:, :, :chain_N]  # n_frames, n_chains, chain_length, n_dimensions

        # mode = msd_acf.normal_modes(pos).mean(axis=1)  # n_frames, all mode, Xp (here, we do the n_chains average)
        # acf = msd_acf.mat_ac(mode, -1)  # n_frames, all mode (here, we sum up the Xp axis)

    # elif "jlx" in cf:
    #     types = xml_parser("out.xml").nodes["type"]
    #     crosslinker = np.where(types == "A")[0].shape[0]
    #     print("crosslinker: " + str(crosslinker))
    #     traj = traj_raw[:, -crosslinker:]  # n_frames, n_chains, crosslinker, n_dimensions
    #     pos = traj.reshape(traj.shape[0], -1, crosslinker, 3) # n_frames, n_chains, crosslinker, n_dimensions

    #     mode = msd_acf.normal_modes(pos).mean(axis=1) # n_frames, all mode, Xp (here, we do the n_chains average)
    #     acf = msd_acf.mat_ac(mode, -1) #n_frames, all mode (here, we sum up the Xp axis)

    elif "diff" in cf:
        traj = traj_raw.reshape(traj_raw.shape[0], -1, 400, 3)  # n_frames, n_chains, chain&ring_length, n_dimensions
        pos = traj[:, :, :50]  # n_frames, n_chains, chain_length, n_dimensions
    else:
        print(f"N = {N}")
        pos = traj_raw.reshape(traj_raw.shape[0], -1, N, 3)  # n_frames, n_chains, chain_length, n_dimensions

    mode = msd_acf.normal_modes(pos).mean(axis=1)  # n_frames, all mode, Xp (here, we do the n_chains average)
    xp = mode.mean(axis=0)  #all mode, Xp (here, we do the time average)

    # def chi_pq(X_p, X_q):
        # numerator = np.abs(np.dot(X_p, X_q))
        # denominator = np.sqrt(np.dot(X_p, X_p) * np.dot(X_q, X_q))
        # result = numerator / denominator
        # return result

    # result_matrix = np.empty((xp.shape[0], xp.shape[0]))
    # for i in range(xp.shape[0]):
        # for j in range(xp.shape[0]):
            # result_matrix[i, j] = chi_pq(xp[i], xp[j])


    # plt.imshow(result_matrix, cmap='hot', interpolation='nearest')
    # plt.colorbar()
    # plt.show()

    data_numerator = np.abs((xp[None, :, :] * xp[:, None, :]).sum(axis=-1))
    data_denominator = np.linalg.norm(xp[None, :, :], axis=-1) * np.linalg.norm(xp[:, None, :], axis=-1)
    data = data_numerator / data_denominator

    plt.imshow(data, cmap='hot', interpolation='nearest')
    plt.colorbar()
    plt.show()
