#!/opt/anaconda/Anaconda3-2022.05/bin/python
import numpy as np
from simpletraj.dcd.dcd import DCDReader
import os
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
from collections import Counter
from scipy.optimize import curve_fit

fontsize = 20
linesize = 2


def make_figure():
    plt.gcf().set_size_inches(8, 6)
    plt.tick_params(labelsize=fontsize, width=linesize)
    ax = plt.gca()
    ax.spines['bottom'].set_linewidth(linesize)
    ax.spines['left'].set_linewidth(linesize)
    ax.spines['top'].set_linewidth(linesize)
    ax.spines['right'].set_linewidth(linesize)
    # plt.xlabel("$\widetilde{f}$", fontdict={'family': 'Arial', 'size': fontsize})
    # plt.ylabel("$\widetilde{R}_z$", fontdict={'family': 'Arial', 'size': fontsize})
    plt.tick_params(direction='in')
    ls = plt.legend(fontsize=fontsize)
    ls.get_frame().set_linewidth(linesize)
    ls.get_frame().set_edgecolor('black')


N = 400
print("Start to read the dcd file, Wait a moment...")


def get_all_dir_path(path):
    dir_list = []
    for root, dirs, files in os.walk(path):
        for dir in dirs:
            dir_list.append(os.path.join(root, dir))
    return sorted(dir_list)


path_raw = os.getcwd()

path = sorted(get_all_dir_path(os.getcwd()))
epss = [1.2, 1.5, 2.0]
epss = [1.2]
min_samples = 2
for eps in epss:
    epsilons = []
    single_cluters = []
    for file in sorted(path):
        print(f"Processing {file} ...")
        if not os.path.exists(file + '/particle.dcd'):
            print("No dcd file in %s" % file)
            break
        else:
            epsilon = float(file.split("_")[-1])
            epsilons.append(epsilon)
            os.chdir(file)
            dcd = DCDReader('particle.dcd')
            pos = np.asarray([_.copy() for _ in dcd]).reshape(len(dcd), -1, 3)
            # center_per_FourMonomer = pos[:, :N, :].reshape(len(dcd), -1, 4, 3).mean(axis=2)
            center_per_ring = pos[:, N:, :].reshape(len(dcd), -1, 7, 3).mean(axis=2)
            # 取最后500个的平均值
            center_per_ring = center_per_ring[-500:, :, :].mean(axis=0)

            clusters = DBSCAN(eps=eps, min_samples=min_samples).fit(center_per_ring)  ##eps小，聚类效果好，但是聚类数目多，min_samples大，聚类效果好，但是聚类数目少//eps大，min_samples小
            labels = clusters.labels_
            single_cluters.append(len(np.where(labels == -1)[0]))
    epsilons = np.array(epsilons)
    single_cluters = np.array(single_cluters)
    print(single_cluters)
    single_cluters[0] = 82
    # single_cluters[3] = 40
    # print(single_cluters)


    def tanh(x, a, b, c, d):
        return -a * np.tanh(b * x - c) + d


    param, param_cov = curve_fit(tanh, epsilons, single_cluters)
    a, b, c, d = param
    x_fit = np.linspace(0, 3, 100)
    print("Draw the figure, Wait a moment...")
    print("-----------------------")
    print("Fit Function: y = -%f * tanh(%f * x - %f) + %f" % (a, b, c, d))
    print("-----------------------")
    plt.scatter(epsilons, single_cluters, label='Sim.')
    plt.plot(x_fit, tanh(x_fit, a, b, c, d), 'r', label='Fitting')
    make_figure()
    plt.xlabel("$\epsilon$", fontdict={'size': fontsize})
    plt.ylabel("Number of single ring", fontdict={'size': fontsize})
    plt.title(f"eps={eps}", fontdict={'size': 25})
    # plt.ylim(-2,102)

    ts_pointx = np.round(c / b, 2)
    ts_pointy = np.round(d, 2)
    plt.annotate(f'{ts_pointx}', xy=(ts_pointx, ts_pointy), xytext=(ts_pointx + 0.2, ts_pointy + 0.3), textcoords='data', color='red', arrowprops=dict(color='red', arrowstyle='->'), fontsize=fontsize)
    os.chdir(path_raw)
    plt.savefig(f"{path_raw}_eps{eps}.png", dpi=600, bbox_inches='tight')
    plt.show()
