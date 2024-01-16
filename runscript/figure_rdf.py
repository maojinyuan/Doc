#!/opt/Anaconda3-2021.11/bin/python

import numpy as np
import matplotlib.pyplot as plt
import os
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
    ls = plt.legend(fontsize=fontsize)
    ls.get_frame().set_linewidth(linesize)
    ls.get_frame().set_edgecolor('black')

path_raw = os.getcwd()
all_items = os.listdir(path_raw)
path = [os.path.join(path_raw, _) for _ in all_items if os.path.isdir(os.path.join(path_raw, _)) and (_.startswith("chain") or _.endswith("polymer"))]
path = sorted(path, key=lambda x: int(re.findall(r'\d+', x)[-1]))
print(path)
for i in path:
    os.chdir(i)
    data = np.load("rdf.npy")
    if "polymer" in i:
        plt.plot(data[:, 0], data[:, 1], color="black", label=i.split("/")[-1])
    else:
        plt.plot(data[:, 0], data[:, 1], label=i.split("/")[-1])

os.chdir(path_raw)
make_figure()
plt.savefig(f"rdf_chain50_ring10.png", dpi=600, bbox_inches='tight')
plt.show()
