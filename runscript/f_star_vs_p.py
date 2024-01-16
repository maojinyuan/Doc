#!/usr/bin/env python
import os
import sys
import numpy as np
from list_dir import get_current_folder_info as gc
from simpletraj.dcd.dcd import DCDReader
# from msd_acf import mat_ac
# from scipy.optimize import curve_fit
from list_dir import *
# import MDAnalysis as mda
# from MDAnalysis.analysis import polymer
# from matplotlib import pyplot as plt

path = sorted(list_folders(), key=lambda s: (int(s.split('b')[1]), int(s.split('p')[1].split('b')[0])))[:7]

nframes = 500
N = int(sys.argv[-1])
ree_f = []
for i in path:
    print(f"Processing {i}")
    num_p = int(i.split("p")[-1].split("b")[0])
    os.chdir(i)
    force_folders = list_folders()
    for single_ff in force_folders:
    # for single_ff in force_folders[6:]:
        if "force1.0_chain_p0b10" not in single_ff:
            print(f"Processing {single_ff}")
            dcdfile = f"{single_ff}/particle.dcd"
            traj_raw = np.array([_.copy() for _ in DCDReader(dcdfile)], dtype=np.float32)

            # ree
            traj_terminal = traj_raw[:, N + num_p * 7:(N + (num_p + 1) * 7), :].mean(axis=1) - traj_raw[:, N - 1, :]  # nframes, ...
            traj_terminal = traj_terminal[-nframes:, -1]
            ree_z_mean = np.abs(traj_terminal).mean()

            # force
            force = float(single_ff.split("force")[-1].split("_")[0])

            arr = np.column_stack((ree_z_mean, force))
            ree_f.append(arr)

    os.chdir("..")
        
ree_f = np.vstack(ree_f)
np.savetxt(rf"ree_f.dat", ree_f, fmt="%10.6f")
