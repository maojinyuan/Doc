#!/usr/bin/env python
import numpy as np
import sys
import os
import re
import matplotlib.pyplot as plt
from figure_para import *
import msd_acf
from list_dir import get_current_folder_info

fontsize=12
linesize=12
# get_current_folder_info() #return current_folder, father_folder, dcd_file, gala_file
cf, ff, dcdfile, galafile = get_current_folder_info()[:-1]

# get abs path
cal_dir_prefix = f"chain"
path_all = sorted([os.path.join(cf, _) for _ in os.listdir(cf) if os.path.isdir(os.path.join(cf, _)) and (_.startswith(cal_dir_prefix))])
path = sorted(path_all, key=lambda x: int(x.split("ring_")[-1]))[1:]
# all_items = os.listdir(cf)
# path = [os.path.join(cf, _) for _ in all_items if os.path.isdir(os.path.join(cf, _)) and (_.startswith("chain") or _.startswith("pure") or _.startswith("gel"))]
# path = path[::-1]

fig, ax = plt.subplots()
for i in path:
    print(f"Processing {i}")
    os.chdir(i)
    # get path/*.dat
    dat_file = np.asarray([file for file in os.listdir() if file.startswith("msd") and file.endswith("dat")])
    digits = [int(re.search(r'\d+', s).group()) for s in dat_file]
    dat_file = dat_file[np.argsort(digits)]

    msd_u1 = np.loadtxt(dat_file[0])
    msd_u2 = np.loadtxt(dat_file[1])
    msd_u3 = np.loadtxt(dat_file[2])
    msd_u4 = np.loadtxt(dat_file[3])

    log_x = np.logspace(0, 6, 30)

    idx1 = np.unique(np.searchsorted(msd_u1[:, 0], log_x))
    idx2 = np.unique(np.searchsorted(msd_u2[:, 0], log_x))
    idx3 = np.unique(np.searchsorted(msd_u3[:, 0], log_x))
    idx4 = np.unique(np.searchsorted(msd_u4[:, 0], log_x))

    msd_u1_new = msd_u1[idx1[:-2], :]
    msd_u2_new = msd_u2[idx2[:-2], :]
    msd_u3_new = msd_u3[idx3[:-2], :]
    msd_u4_new = msd_u4[idx4[:-2], :]

    msd_u = np.concatenate((msd_u1_new, msd_u2_new, msd_u3_new, msd_u4_new), axis=0)

    idx = np.unique(np.searchsorted(msd_u[:, 0], log_x))
    msd = msd_u[idx[:-2], :]

    #slope
    # slope = np.log(np.diff(msd[2:, 1])) / np.log(np.diff(msd[2:, 0]))
    # np.savetxt(f'../slope{i.split("/")[-1]}.dat', np.c_[slope], fmt='%.8f')

    if "purepolymer_copied" in i:
        ax.loglog(msd[:, 0], msd[:, 2], color="black", label=f"{i.split('/')[-1]}")
    else:
        # ax.loglog(msd[:, 0], msd[:, 1], label=f"{i.split('/')[-1].split('50_')[-1]}_chain", color=color_list[path.index(i)])
        ax.loglog(msd[:, 0], msd[:, 2], label=f"Ring={os.path.basename(i).split('_')[-1]}", color=color_list[path.index(i)])
        # ax.loglog(msd[:, 0], msd[:, 2], label=f"{i.split('/')[-1].split('50_')[-1]}_ringCM", color=color_list[path.index(i)], linestyle="--")
make_figure(20,2)
os.chdir(cf)
ax.set_xlabel("t")
# ax.set_ylabel("$g_1(t)$")
# ax.set_ylabel("$g_1^{ringCM}(t)$")
ax.set_ylabel("$g_1^{ring}(t)$")
# ax.set_title(f"{cf.split("")}")
# ax.set_title(f"Chain=100")
legend = ax.legend(fontsize=18, loc="lower right")
legend.get_frame().set_linewidth(2)
legend.get_frame().set_edgecolor('black')
plt.savefig(f"msd_ring.png", dpi=600, bbox_inches='tight')
# plt.savefig(f"msd_ringCM.png", dpi=600, bbox_inches='tight')
# plt.savefig(f"msd_{cf.split('/')[-1]}.png", dpi=600, bbox_inches='tight')
plt.show()
