#!/usr/bin/env python3
import numpy as np
import sys
import os
from simpletraj.dcd.dcd import DCDReader
import msd_acf
from list_dir import get_current_folder_info
from xml_io import uxml

# get_current_folder_info() #return current_folder, father_folder, dcd_file, gala_file
cf, ff, dcdfile, galafile = get_current_folder_info()[:-1]

# dtau = [line.split()[-1] for line in open(galafile[-1], "r") if "app.setDt" in line and "#" not in line][-1]
# dtau = float(dtau.split("(")[-1].split(")")[0])
# print("dtau: " + str(dtau))
# print("------------------------")
# dtau = [line.split()[-1] for line in open(galafile[-1], "r") if "app.setDt" in line and "#" not in line][-1]
# dtau = float(dtau.split("(")[-1].split(")")[0])
# print("dtau: " + str(dtau))

xml = uxml.read("out.xml")
N = xml.natoms
print("N: " + str(N))
print("----------------------")

for i in dcdfile:
    print("Processing " + i + " ...")
    slice_dcd = DCDReader(i)
    dtau = 0.01
    traj_raw = np.array([_.copy() for _ in slice_dcd], dtype=np.float32)
    period = slice_dcd.skip_timestep
    natoms = slice_dcd.numatoms
    # period = int(i.split("_")[-1].split("Period")[0])

    if "chain" in cf and "ring" in cf:
        # ring_N = os.path.basename(cf)
        # chain_ring = int(ring_N.split("_")[1]) + 7 * int(ring_N.split("_")[-1])
        # chain_N = int(ring_N.split("_")[1])

        chain_ring = int(os.path.basename(cf).split("_")[1]) + 7 * int(os.path.basename(cf).split("_")[-1])
        chain_N = int(os.path.basename(cf).split("_")[1])
        print(f"N + ring = {chain_ring}, N = {chain_N}")

        traj = traj_raw.reshape(traj_raw.shape[0], -1, chain_ring, 3)
        traj_final = traj[:, :, :chain_N].reshape(traj.shape[0], -1, 3) # n_frames, n_length, ...

        #ring and ringCM
        traj_ring_final = traj[:, :, chain_N:].reshape(traj.shape[0], -1, 7, 3)
        traj_ringCM_final = traj_ring_final.mean(axis=-2) # n_frames, n_ringCM, ...

    # if "jlx" in cf:
        # types = xml_parser("out.xml").nodes["type"]
        # crosslinker = np.where(types == "A")[0].shape[0]
        # print("crosslinker: " + str(crosslinker))
        # traj_final = traj_raw[:, -crosslinker:]

    # if "purepolymer" in cf:
        # traj_final = traj_raw
    
    # if "purepolymer_copied" in cf:
        # traj_final = traj_raw

    u = traj_final.swapaxes(0, 1)
    msd = np.asarray([msd_acf.msd_fft(_) for _ in u]).mean(axis=0)

    u_ring = traj[:, :, chain_N:].reshape(traj.shape[0], -1, 3).swapaxes(0, 1)
    msd_ring = np.asarray([msd_acf.msd_fft(_) for _ in u_ring]).mean(axis=0)

    u_ringCM = traj_ringCM_final.swapaxes(0, 1)
    msd_ringCM = np.asarray([msd_acf.msd_fft(_) for _ in u_ringCM]).mean(axis=0)

    dt = dtau * period
    time = dt * np.arange(msd.shape[0])

    np.savetxt(f'msd{period}.dat', np.c_[time, msd, msd_ring, msd_ringCM], fmt='%.8f')
    # np.savetxt(f'msd{period}.dat', np.c_[time, msd], fmt='%6.8f')
