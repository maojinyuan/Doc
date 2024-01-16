#!/usr/bin/env python3
import numpy as np
from scipy.optimize import curve_fit
import re
from figure_para import *
from list_dir import get_current_folder_info


def func(t, a, b):
    return np.exp(-(t / a) ** b)
    # return np.exp(-a * t ** b)


# get_current_folder_info() #return current_folder, father_folder, dcd_file, gala_file
cf, ff, dcdfile, galafile = get_current_folder_info()[:-1]

# get abs path
all_items = os.listdir(cf)
path = [os.path.join(cf, _) for _ in all_items if os.path.isdir(os.path.join(cf, _)) and (_.startswith("chain") or _.startswith("pure") or _.startswith("gel"))]
path = path[::-1]
last_str = path.pop()
path.insert(0, last_str)

# p_mode = int(input("input pMode(>=1, default 1): ") or 1) # should >= 1
pmodes = np.array([1, 2, 3, 4, 5])
# pmodes = np.array([1, 2, 5, 10])
fig, ax = plt.subplots()
for i in path:
    for p_mode in pmodes:
        print(f"Processing {os.path.basename(i)} pmode={p_mode}...")
        os.chdir(i)
        # get path/*.dat
        dat_file = np.asarray([file for file in os.listdir() if file.startswith("acf") and file.endswith("dat")])
        digits = [int(re.search(r'\d+', s).group()) for s in dat_file]
        dat_file = dat_file[np.argsort(digits)]

        acf_u1 = np.loadtxt(dat_file[0])
        acf_u2 = np.loadtxt(dat_file[1])
        acf_u3 = np.loadtxt(dat_file[2])
        acf_u4 = np.loadtxt(dat_file[3])

        log_x = np.logspace(0, 8, 100)

        idx1 = np.unique(np.searchsorted(acf_u1[:, 0], log_x))
        idx2 = np.unique(np.searchsorted(acf_u2[:, 0], log_x))
        idx3 = np.unique(np.searchsorted(acf_u3[:, 0], log_x))
        idx4 = np.unique(np.searchsorted(acf_u4[:, 0], log_x))

        acf_u1_new = acf_u1[idx1[:-2], :]
        acf_u2_new = acf_u2[idx2[:-2], :]
        acf_u3_new = acf_u3[idx3[:-2], :]
        acf_u4_new = acf_u4[idx4[:-2], :]

        acf_u = np.concatenate((acf_u1_new, acf_u2_new, acf_u3_new, acf_u4_new), axis=0)

        idx = np.unique(np.searchsorted(acf_u[:, 0], log_x))
        acf = acf_u[idx[:-2], :]

        time = acf[:, 0]
        acf_final = acf[:, p_mode] / acf[0, p_mode]

        # fitting
        popt, pcov = curve_fit(func, time, acf_final)
        print(popt)
        # fitting_x = np.logspace(0, 6.5, 200)
        # fitting_y = [func(i, popt[0], popt[1]) for i in fitting_x]

        # plot
        if "purepolymer" in i:
            plt.scatter(time, acf_final, label=f"p={p_mode}")
            # plt.plot(fitting_x, fitting_y, label=f"pmode={p_mode} fitting", linestyle=':')
        else:
            plt.scatter(time, acf_final, label=f"p={p_mode}")
            # plt.plot(fitting_x, fitting_y, label=f"pmode={p_mode} fitting", linestyle=':')

    os.chdir(cf)
    make_figure()
    plt.xlabel("t")
    plt.ylabel(r"$\langle X_p(t)X_p(0) \rangle / \langle X_{p}^2 \rangle$")
    plt.xscale('log')
    # plt.yscale('log')
    plt.title(rf"{os.path.basename(i)} (p={pmodes})")
    # figure_show()
    # plt.savefig(f"{i}_acf_p{pmodes}.png", dpi=600, bbox_inches='tight')
    plt.show()
    
    # exit()
