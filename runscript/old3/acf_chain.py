#!/usr/bin/env python
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

dtau = [line.split()[-1] for line in open(galafile[-1], "r") if "app.setDt" in line and "#" not in line][-1]
dtau = float(dtau.split("(")[-1].split(")")[0])
print("dtau: " + str(dtau))

xml = uxml.read("out.xml")
N = xml.natoms
print("N: " + str(N))
print("----------------------")

for i in dcdfile:
    print("Processing " + i + " ...")
    traj_raw = np.array([_.copy() for _ in DCDReader(i)], dtype=np.float32)
    period = int(i.split("_")[-1].split("Period")[0])
    if "chain" in cf and "ring" in cf:
        chain_ring = int(cf.split("_")[1]) + 7 * int(cf.split("_")[-1])
        chain_N = int(cf.split("_")[1])

        traj = traj_raw.reshape(traj_raw.shape[0], -1, chain_ring, 3)  # n_frames, n_chains, chain&ring_length, n_dimensions
        pos = traj[:, :, :chain_N]  # n_frames, n_chains, chain_length, n_dimensions

        mode = msd_acf.normal_modes(pos).mean(axis=1)  # n_frames, all mode, Xp (here, we do the n_chains average)
        acf = msd_acf.mat_ac(mode, -1)  # n_frames, all mode (here, we sum up the Xp axis)

    # elif "jlx" in cf:
        # types = xml_parser("out.xml").nodes["type"]
        # crosslinker = np.where(types == "A")[0].shape[0]
        # print("crosslinker: " + str(crosslinker))
        # traj = traj_raw[:, -crosslinker:]  # n_frames, n_chains, crosslinker, n_dimensions
        # pos = traj.reshape(traj.shape[0], -1, crosslinker, 3) # n_frames, n_chains, crosslinker, n_dimensions

        # mode = msd_acf.normal_modes(pos).mean(axis=1) # n_frames, all mode, Xp (here, we do the n_chains average)
        # acf = msd_acf.mat_ac(mode, -1) #n_frames, all mode (here, we sum up the Xp axis)

    else:
        pos = traj_raw.reshape(traj_raw.shape[0], -1, N, 3)  # n_frames, n_chains, chain_length, n_dimensions

        mode = msd_acf.normal_modes(pos).mean(axis=1)  # n_frames, all mode, Xp (here, we do the n_chains average)
        acf = msd_acf.mat_ac(mode, -1)  # n_frames, all mode (here, we sum up the Xp axis)

    t = period * dtau
    time = np.arange(mode.shape[0]) * t
    np.savetxt(f'acf{period}_chain.dat', np.c_[time, acf], fmt='%10.3f')
