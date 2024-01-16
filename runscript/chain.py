#!/usr/bin/env python3
import math
import sys
import numpy as np
from scipy.optimize import curve_fit
import re
from figure_para import *
from list_dir import *

# color_list = ['#63b2ee', '#76da91', '#f8cb7f', '#f89588', '#7cd6cf', '#9192ab', '#7898e1', '#efa666', '#eddd86','#9987ce', '#63b2ee', '#76da91']
color_list = ['#63b2ee', '#76da91', '#f89588', '#7cd6cf', '#9192ab', '#7898e1', '#efa666', '#eddd86', '#9987ce', '#63b2ee', '#76da91', '#f89588']
markers = ['o', 'v', '<', 's', 'D', '^', '>', 'p', 'H', '*', 'x']


def func(x, tau_p, beta_):
    return -(x / tau_p) ** beta_


# get_current_folder_info() #return current_folder, father_folder, dcd_file, gala_file
cf, ff, dcdfile, galafile = get_current_folder_info()[:-1]

# get abs path
chain = int(sys.argv[1])
cal_dir_prefix = f"chain_{chain}"
path_all = sorted([os.path.join(cf, _) for _ in os.listdir(cf) if os.path.isdir(os.path.join(cf, _)) and (_.startswith(cal_dir_prefix))])
path = sorted(path_all, key=lambda x: int(x.split("_")[-1]))

# p_modes = np.array([1, 2, 3, 5, 10, 15, 25, 40])
p_modes = np.array([1])
print(f"p mode = {p_modes}")

# acf_pp = np.loadtxt(f"{cf}/chain_50_ring_0/acf10.dat")
# acf_final_pp = acf_pp[:, :-1] / acf_pp[0, :-1]
# time_pp = acf_pp[:, 0]

tau_p_beta_p = []
fig, ax = plt.subplots()
# for i in path[::-1]:
for i in path:
    print(f"Processing {i}")
    os.chdir(i)
    ring = int(os.path.basename(i).split('_')[-1])
    phi = np.around(np.multiply(ring / chain, 100)).astype(int)
    for p_mode in p_modes.tolist():
        # get path/*.dat
        dat_file = np.asarray([file for file in os.listdir() if file.startswith("acf") and file.endswith("dat")])
        digits = [int(re.search(r'\d+', s).group()) for s in dat_file]
        dat_file = dat_file[np.argsort(digits)]

        acf_u1 = np.loadtxt(dat_file[0])
        acf_u2 = np.loadtxt(dat_file[1])
        acf_u3 = np.loadtxt(dat_file[2])
        acf_u4 = np.loadtxt(dat_file[3])

        # cut_num = -1
        # acf_u1_new = acf_u1[:cut_num, 1:] / acf_u1[0, 1:]
        # acf_u2_new = acf_u2[:cut_num, 1:] / acf_u2[0, 1:] 
        # acf_u3_new = acf_u3[:cut_num, 1:] / acf_u3[0, 1:]
        # acf_u4_new = acf_u4[:cut_num, 1:] / acf_u4[0, 1:]

        # acf_u = np.concatenate((acf_u1_new, acf_u2_new, acf_u3_new, acf_u4_new), axis=0)
        # time = np.concatenate((acf_u1[:cut_num, 0], acf_u2[:cut_num, 0], acf_u3[:cut_num, 0], acf_u4[:cut_num, 0])).reshape((-1,1))

        # acf = np.concatenate((time, acf_u), axis=-1)
        # unique_indices = np.unique(acf[:, 0], return_index=True)[1]
        # acf = acf[unique_indices]
        # # acf = acf[np.linspace(0, acf.shape[0]-1, num=100, dtype=int)]

        # time = acf[:, 0]
        # acf_final = acf[:, p_mode]

        if p_mode <= 3:
            acf = acf_u3[:, :-1]
        elif p_mode <= 14:
            acf = acf_u2[:, :-1]
        else:
            acf = acf_u1[:, :-1]

        acf_final = acf[:, p_mode] / acf[0, p_mode]
        time = acf[:, 0]
        # 找到acf_final中第一个小于0.2的值对应的time的index
        cut_idx1 = np.where(acf_final < 0.4)[0][0]
        # cut_idx2 = np.where(acf_final < 0)[0][0]
        if "ring_100" in i:
            cut_idx2 = acf_final.shape[0]
        else:
            # cut_idx2 = np.where(acf_final < 0.3)[0][0]
            cut_idx2 = np.where(acf_final < 0.01)[0][0]
        # 用这个index截取time和acf_final
        time_fit = time[:cut_idx1]
        acf_final_fit = acf_final[:cut_idx1]
        # 做拟合
        popt, pcov = curve_fit(func, time_fit, np.log(acf_final_fit))
        tau_p_beta_p.append(popt)
        tau_p, beta_p = popt
        print(popt)

        log_x = np.logspace(-2, 5, 200)
        fitting_x = log_x
        fitting_y = np.exp([func(i, popt[0], popt[1]) for i in fitting_x])
        plt.scatter(time[:cut_idx2], acf_final[:cut_idx2], label=f"pmode={p_mode}", marker=markers[np.where(p_modes == p_mode)[0][0]], edgecolor=color_list[np.where(p_modes == p_mode)[0][0] + 3], facecolor='none')
        plt.plot(fitting_x, fitting_y, color='black', linestyle="--")
    # plt.plot(time_pp, acf_final_pp[:, ring], label=f"pmode={ring}", color='black')
    make_figure()
    fontsize = 16
    legend = plt.legend(fontsize=fontsize, loc='upper right')
    legend.get_frame().set_linewidth(linesize)
    legend.get_frame().set_edgecolor('black')
    plt.xlim(0.04, 0.9 * 10 ** 5)
    plt.ylim(-0.1, 1.1)
    plt.xlabel("t$/\\tau$")
    plt.ylabel(r"$\langle X_p(t)X_p(0) \rangle / \langle X_{p}^2 \rangle$")
    plt.xscale('log')
    plt.title(rf"Ring={ring},Chain={chain} ($\phi={phi}$%)", fontdict={'size': 20})

    img_dir = "img_pmode"
    if not os.path.exists(os.path.join(cf, img_dir)):
        os.makedirs(os.path.join(cf, img_dir))

    # plt.savefig(f"{cf}/img_pmode/phi={phi}%.png", dpi=600, bbox_inches='tight')
    plt.show()
    plt.close()
tau_p_beta_p = np.asarray(tau_p_beta_p).reshape((len(path), -1, 2))
# np.save(f"{cf}/tau_p_beta_p.npy", tau_p_beta_p)
