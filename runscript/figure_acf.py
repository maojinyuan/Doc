#!/usr/bin/env python3
import numpy as np
from scipy.optimize import curve_fit
import re
from figure_para import *
from list_dir import get_current_folder_info


def func(t, a, b):
    return np.exp(-(t / a) ** b)

# get_current_folder_info() #return current_folder, father_folder, dcd_file, gala_file
cf, ff, dcdfile, galafile = get_current_folder_info()[:-1]

# get abs path
all_items = os.listdir(cf)
path = [os.path.join(cf, _) for _ in all_items if os.path.isdir(os.path.join(cf, _)) and (_.startswith("chain") or _.startswith("pure") or _.startswith("gel"))]
path = path[::-1]
last_str = path.pop()
path.insert(0, last_str)

p_mode = int(input("input pMode(>=1, default 1): ") or 1) # should >= 1
fig, ax = plt.subplots()
for i in path:
    print(f"Processing {i}")
    os.chdir(i)
    # get path/*.dat
    dat_file = np.asarray([file for file in os.listdir() if file.startswith("acf") and file.endswith("dat")])
    digits = [int(re.search(r'\d+', s).group()) for s in dat_file]
    dat_file = dat_file[np.argsort(digits)]

    acf_u1 = np.loadtxt(dat_file[0])
    acf_u2 = np.loadtxt(dat_file[1])
    acf_u3 = np.loadtxt(dat_file[2])
    acf_u4 = np.loadtxt(dat_file[3])

    log_x = np.logspace(0, 7, 50)

    idx1 = np.unique(np.searchsorted(acf_u1[:, 0], log_x))
    idx2 = np.unique(np.searchsorted(acf_u2[:, 0], log_x))
    idx3 = np.unique(np.searchsorted(acf_u3[:, 0], log_x))
    idx4 = np.unique(np.searchsorted(acf_u4[:, 0], log_x))

    acf_u1_new = acf_u1[idx1[:-2], :]
    acf_u2_new = acf_u2[idx2[:-2], :]
    acf_u3_new = acf_u3[idx3[:-2], :]
    acf_u4_new = acf_u4[idx4[:-2], :]

    acf_u = np.concatenate((acf_u1_new, acf_u2_new, acf_u3_new, acf_u4_new), axis=0)
    # if "polymer" in i:
        # acf_u = np.concatenate((acf_u1_new, acf_u2_new, acf_u3_new, acf_u4_new), axis=0)
    # if "chain_50_ring_25" in i:
        # acf_u = np.concatenate((acf_u2_new, acf_u4_new), axis=0)
    # if "chain_50_ring_10" in i:
        # acf_u = np.concatenate((acf_u2_new, acf_u3_new), axis=0)
    # if "chain_50_ring_5" in i:
        # acf_u = np.concatenate((acf_u1_new, acf_u3_new), axis=0)

    idx = np.unique(np.searchsorted(acf_u[:, 0], log_x))
    acf = acf_u[idx[:-2], :]


    # acf = np.loadtxt(dat_file[-1])
    #slope
    # slope = np.log(np.diff(acf[2:, 1])) / np.log(np.diff(acf[2:, 0]))
    # np.savetxt(f'../slope{i.split("/")[-1]}.dat', np.c_[slope], fmt='%.8f')

    time = acf[:, 0]
    acf_final = acf[:, p_mode] / acf[0, p_mode]
    popt, pcov = curve_fit(func, time, acf_final)
    print(popt)
    fitting_x = np.logspace(0, 6.5, 200)
    fitting_y = [func(i, popt[0],popt[1]) for i in fitting_x]

    # idx = -1

    # if acf_final.min() > 0:
        # idx = acf_final.shape[0]
    # else:
        # idx = np.argwhere(acf_final < 0)[0, 0]

    # if "purepolymer" in i:
        # plt.scatter(time[:idx], acf_final[:idx], color="black", label=f"{i.split('/')[-1]}")
    # else:
        # plt.scatter(time[:idx], acf_final[:idx], label=f"{i.split('/')[-1]}")
    if "purepolymer" in i:
        plt.scatter(time, acf_final, color="black", label=f"{i.split('/')[-1]}")
        plt.plot(fitting_x, fitting_y, color="black", linestyle=":")
    else:
        plt.scatter(time, acf_final, label=f"{i.split('/')[-1]}")
        plt.plot(fitting_x, fitting_y, linestyle=":")
os.chdir(cf)
make_figure()
plt.xlabel("t")
plt.ylabel(r"$\langle X_p(t)X_p(0) \rangle / \langle X_{p}^2 \rangle$")
plt.xscale('log')
plt.title(rf"{cf.split('/')[-1]} (pmode={p_mode})", fontdict={'size':20})
# plt.savefig(f"acf_{cf.split('/')[-1]}_pmode{p_mode}.png", dpi=600, bbox_inches='tight')
plt.show()
