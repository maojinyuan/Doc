#!/share/simulation/Anaconda3-2022.05/bin/python
import warnings
import sys
import os
from simpletraj.dcd.dcd import DCDReader
import numpy as np
from numba import guvectorize
from numba import float64
import matplotlib.pyplot as plt
import re

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
    # ls = plt.legend(fontsize=fontsize, loc='upper right')
    ls = plt.legend(fontsize=fontsize, loc='lower left')
    ls.get_frame().set_linewidth(linesize)
    ls.get_frame().set_edgecolor('black')



path_raw = os.getcwd()
all_items = os.listdir(path_raw)
path = [os.path.join(path_raw, _) for _ in all_items if os.path.isdir(os.path.join(path_raw, _)) and (_.startswith("chain") or _.startswith("pure"))]
path = path[::-1]

p_mode = 1 # should >= 0
fig, ax = plt.subplots()
for i in path:
    print(f"Processing {i}")
    os.chdir(i)
    # get path/*.dat
    dat_file = np.asarray([file for file in os.listdir() if file.startswith("acf") and file.endswith("dat")])
    digits = [int(re.search(r'\d+', s).group()) for s in dat_file]
    dat_file = dat_file[np.argsort(digits)]

    msd_u1 = np.loadtxt(dat_file[0])
    msd_u2 = np.loadtxt(dat_file[1])
    msd_u3 = np.loadtxt(dat_file[2])
    msd_u4 = np.loadtxt(dat_file[3])

    # print(msd_u1.shape)
    # print(msd_u2.shape)
    # print(msd_u3.shape)
    # print(msd_u4.shape)
    # exit()

    log_x = np.logspace(0, 7, 20)

    idx1 = np.unique(np.searchsorted(msd_u1[:, p_mode], log_x))
    idx2 = np.unique(np.searchsorted(msd_u2[:, p_mode], log_x))
    idx3 = np.unique(np.searchsorted(msd_u3[:, p_mode], log_x))
    idx4 = np.unique(np.searchsorted(msd_u4[:, p_mode], log_x))

    msd_u1_new = msd_u1[idx1[:-2], :]
    msd_u2_new = msd_u2[idx2[:-2], :]
    msd_u3_new = msd_u3[idx3[:-2], :]
    msd_u4_new = msd_u4[idx4[:-2], :]

    msd_u = np.concatenate((msd_u1_new, msd_u2_new, msd_u3_new, msd_u4_new), axis=0)
    # msd_u = np.insert(msd_u, 0, msd_u1[0, :], axis=0)

    idx = np.unique(np.searchsorted(msd_u[:, 0], log_x))
    msd = msd_u[idx[:-2], :]
    # slope = np.log(np.diff(msd[2:, 1])) / np.log(np.diff(msd[2:, 0]))
    
    #slope
    # slope = np.log(np.diff(msd[2:, 1])) / np.log(np.diff(msd[2:, 0]))
    # np.savetxt(f'../slope{i.split("/")[-1]}.dat', np.c_[slope], fmt='%.8f')

    if i.split("/")[-1] == "purepolymer":
        ax.loglog(msd[:, 0], msd[:, 1], color="black", label=f"{i.split('/')[-1]}")
    else:
        ax.loglog(msd[:, 0], msd[:, 1], label=f"{i.split('/')[-1]}")
    make_figure()
    plt.xlabel("t")
    plt.ylabel("$g_1(t)$")

os.chdir(path_raw)

plt.savefig(f"acf_{path_raw.split('/')[-1]}.png", dpi=600, bbox_inches='tight')
plt.show()
