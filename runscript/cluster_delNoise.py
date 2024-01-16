#!/opt/anaconda/Anaconda3-2022.05/bin/python
import numpy as np
from simpletraj.dcd.dcd import DCDReader
import os
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
from collections import Counter
from scipy.optimize import curve_fit
import seaborn as sns
import pandas as pd

fontsize = 13
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
    ls = plt.legend(fontsize=fontsize, loc='upper left')
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


path = sorted(get_all_dir_path(os.getcwd()))
# path = ["/home/mjy/rigid_flexible_chain/1rigid_chain/LJ/LJ_1.0"]

epss = [1.5, 2.0]
min_samples = 2

for eps in epss:
    epsilons = []
    for file in sorted(path):
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
            center_per_ring = center_per_ring[-2000:, :, :].mean(axis=0)

            clusters = DBSCAN(eps=eps, min_samples=min_samples).fit(center_per_ring)  ##eps小，聚类效果好，但是聚类数目多，min_samples大，聚类效果好，但是聚类数目少//eps大，min_samples小
            labels = clusters.labels_
            labels_delNoise = np.delete(labels, np.where(labels == -1))
            # delNoise
            cluster_size_num = np.array(list(Counter(labels_delNoise).values()))
            hist, bin = np.histogram(cluster_size_num, bins=range(1, 100 + 1))
            # 用pandas转换成dataframe，其中行首为epsilon，列首为bin[:-1],值为hist

            hist_epsilon = pd.DataFrame(hist, index=(bin[:-1] + bin[1:]) / 2, columns=[epsilon]).T
            # 循环将hist_epsilon添加到hist_all中
            if epsilon == epsilons[0]:
                hist_all = hist_epsilon
            else:
                hist_all = pd.concat([hist_all, hist_epsilon], axis=0)
    # print(hist_all)
    sns.heatmap(hist_all[::-1], cmap='Blues', cbar_kws={'label': 'Cluster Number'})
    #设置x轴的标签呈斜向
    plt.xticks(rotation=45, fontsize=10)
    #设置y轴的标签呈斜向
    plt.yticks(rotation=30, fontsize=10)
    plt.xlabel("Cluster Size", fontdict={'size': 15})
    plt.ylabel("epsilon", fontdict={'size': 15})
    plt.title(f"eps={eps}", fontdict={'size': 25})
    plt.savefig(rf"/home/mjy/rigid_flexible_chain/1rigid_chain/LJ/clusterDelNoise_eps={eps}.png", dpi=600, bbox_inches='tight')
    # plt.show()
