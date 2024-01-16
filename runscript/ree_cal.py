#!/usr/bin/env python
import os
import numpy as np
from list_dir import get_current_folder_info as gc
from simpletraj.dcd.dcd import DCDReader
from msd_acf import mat_ac
from scipy.optimize import curve_fit
from list_dir import *
import MDAnalysis as mda
from MDAnalysis.analysis import polymer
from matplotlib import pyplot as plt

# def calculate_persistent_length(coordinates):
    # bond_vectors = coordinates[:, 1:, :] - coordinates[:, :-1, :]
    # ave_bond = np.linalg.norm(bond_vectors, axis=-1).mean(axis=-1)

    # cos_theta = np.sum(bond_vectors[:, :-1, :] * bond_vectors[:, 1:, :], axis=-1) / (np.linalg.norm(bond_vectors[:, :-1, :], axis=-1) * np.linalg.norm(bond_vectors[:, 1:, :], axis=-1))

    # # persistent_lengths = -1 / np.log(np.mean(cos_theta, axis=1))
    # persistent_lengths = -1 / np.log(np.mean(cos_theta, axis=1))
    # return bond_vectors, ave_bond, cos_theta, persistent_lengths

def bondVec_numpy(_a):
    def e_fit(x, _alpha):
        return np.exp(-x/_alpha)

    bond_vec_raw = _a[:,1:,:] - _a[:,:-1,:]
    bond_vec = bond_vec_raw / np.linalg.norm(bond_vec_raw, axis=-1)[:, :, None]

    lb = np.linalg.norm(bond_vec_raw, axis=-1).mean()

    ret = []
    for frame in bond_vec:
        ret.append(np.asarray([(frame[inv:] * frame[:-inv]).sum(axis=-1).mean() if inv != 0 else 1.0 for inv in range(frame.shape[0])]))
    ret = np.asarray(ret).mean(axis=0)

    a_fit = curve_fit(e_fit, np.arange(bond_vec.shape[1]), ret)[0][0]
    lp = a_fit * lb

    return np.arange(bond_vec.shape[1]), ret, lp

path = sorted(list_folders(), key=lambda s: (int(s.split('b')[1]), int(s.split('p')[1].split('b')[0])))[:7]

nframes = 500
N = 400
ree_f_n_lp_pring_brings = []
bond_acfs = []
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
            # print(ree_mean)

            # force
            force = float(single_ff.split("force")[-1].split("_")[0])
            # force_arr = np.ones(ree.shape[0]) * force

            # calculate min_distance_idx
            vec_CM_polymer = traj_raw[:, N + num_p * 7:(N + (num_p + 1) * 7), :].mean(axis=1)[:, None, :] - traj_raw[:, :N, :]
            distance_CM_polymer = np.linalg.norm(vec_CM_polymer[-nframes:], axis=-1)
            min_distance_idx = int(np.argmin(distance_CM_polymer, axis=1).mean())

            # n
            n_segment = N - min_distance_idx

            # # lp
            # u = mda.Universe(f"{single_ff}/xml/particle.0100000000.xml", f"{single_ff}/particle.dcd")
            # # total_chain = [u.select_atoms("type A T")]
            # idx = f"{min_distance_idx}:{N - 1}"
            # total_chain = [u.select_atoms(f"type A T and index {idx}")]
            # # print(total_chain)
            # persistence_length = polymer.PersistenceLength(total_chain)
            # persistence_length = persistence_length.run(start=-nframes, stop=-1, step=1)
            # # print(persistence_length.results.lp)
            # # lp_arr = np.ones(ree.shape[0]) * persistence_length.results.lp

            #lp method2
            # traj_raw_lpp = traj_raw[-nframes:, :min_distance_idx, :]
            traj_raw_lpb = traj_raw[-nframes:, min_distance_idx:N, :]
            # lpp = bondVec_numpy(traj_raw_lpp)[-1]
            lpb = bondVec_numpy(traj_raw_lpb)[-1]
            # y, lpb, x = bondVec_numpy(traj_raw_lpb)
            # print(min_distance_idx)
            print(lpb)
            print("")
            # lpp_arr = np.ones(ree.shape[0]) * lpp
            # lpb_arr = np.ones(ree.shape[0]) * lpb

            #bond_acf
            # bond_acfp = bondVec_numpy(traj_raw_lpp)[0]
            # bond_acfb = bondVec_numpy(traj_raw_lpb)[0]

            # xp = bondVec_numpy(traj_raw_lpp)[-1]
            # xb = bondVec_numpy(traj_raw_lpb)[-1]

            # print(bond_acfb)


            # lp_re = calculate_persistent_length(traj_raw[-nframes:, :min_distance_idx, :])
            # print(lp_re[-1].mean())


            # pring, bring
            pring = int(single_ff.split('p')[-1].split('b')[0])
            bring = int(single_ff.split('b')[-1])
            # pring_arr = np.ones(ree.shape[0]) * pring
            # bring_arr = np.ones(ree.shape[0]) * bring

            # time axis
            # dtau, period = 0.005, 100000
            # time = np.arange(ree.shape[0]) * dtau * period

            # ree_f_lp_pring_bring = np.column_stack((ree, force_arr, lpp_arr, lpb_arr, pring_arr, bring_arr))
            # ree_f_n_lp_pring_bring = np.column_stack((ree_z_mean, force, n_segment, lpp, lpb, pring, bring))
            # bond_acf = np.column_stack((bond_acfb))
            # bond_acfs.append(bond_acfp)
            # bond_acfs.append(bond_acfb)
            # ree_f_n_lp_pring_brings.append(ree_f_n_lp_pring_bring)
            # plt.plot(xp, bond_acfp, label=f"{single_ff}")
            # plt.plot(xb, bond_acfb, label=f"{single_ff}")
            # plt.scatter(x[:-10], y[:-10], label=f"{single_ff}")
    # plt.show()
    os.chdir("..")
        
# plt.legend()
# plt.xscale('log')
# plt.show()
# ree_f_n_lp_pring_brings = np.vstack(ree_f_n_lp_pring_brings)
# np.savetxt(rf"ree_f_n_lpp_lpb_pring_brings.dat", ree_f_n_lp_pring_brings, fmt="%10.6f")
# np.savetxt(rf"bond_acf_pb.dat", bond_acfs, fmt="%10.6f")
