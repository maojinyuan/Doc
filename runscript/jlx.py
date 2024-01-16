#!/usr/bin/env python3
import numpy as np
from xml_io import xml_parser
import matplotlib.pyplot as plt
import os
from list_dir import get_current_folder_info

def pbc(a, box):
    return a - box * np.round(a / box)

# get_current_folder_info() #return current_folder, father_folder, dcd_file, gala_file
cf, ff, dcdfile, galafile = get_current_folder_info()

dt = [line.split()[-1] for line in open(galafile[-1], "r") if "app.setDt" in line and "#" not in line][-1]
dt = float(dt.split("(")[-1].split(")")[0])
interval = [line for line in open(galafile[-1], "r") if "xml.setPeriod" in line and "#" not in line[0]][-1]
interval = int(interval.split("(")[-1].split(")")[0])
print("dt: " + str(dt))
print("interval: " + str(interval))

fnames = [i for i in sorted(os.listdir(f"{os.getcwd()}/xml/")) if i.startswith("particle") and i.endswith("xml")][-200:]
# fnames = [i for i in sorted(os.listdir(os.getcwd())) if i.startswith("particle") and i.endswith("xml")][-20:]
print(fnames)

os.chdir(f"{os.getcwd()}/xml/")

pos_As = []
for fname in fnames:
    xml = xml_parser(fname)
    pos = xml.nodes['position']
    pos = pbc(pos, xml.box)
    bond = xml.nodes['bond']

    bond_TA = np.unique(bond[np.where(bond[:, 0] == 'T-A')][:, -1])
    pos_A = pos[bond_TA.astype(int)]
    pos_As.append(pos_A)
pos_As = np.asarray(pos_As, dtype=object)
natom = pos_As.shape[1]

time = dt * interval
t = np.arange(len(fnames)) * time

particle_idx = 100
position_idx = 0
# pos_1particle_x_mean = pos_As[:, particle_idx, 0].mean(axis=0)
pos_1particle_x_mean = pos_As[:, :, position_idx].mean(axis=0)
# pos_1particle_x_mean = np.ones(nframe) * pos_1particle_x_mean
# pos_1particle_x_mean = np.ones(pos_1particle_x_mean.shape[0]) * pos_1particle_x_mean
# distance = pos_As[:, particle_idx, 0]-pos_1particle_x_mean
distance = pos_As[:, :, position_idx]-pos_1particle_x_mean[None,:] #nframe, natom
# distance = distance.T #natom, nframe
nk = distance.shape[0] * 15
h = np.zeros(nk, dtype=np.int_)
for i in range(natom):
    hist, edges = np.histogram(distance[:, i], density=True, bins=nk, range=(distance.min(), distance.max()))
    h = h + hist
    # pos_1particle = pos_As[:, 10, :]
    # pos_1particle_x = pos_1particle[:,0]
h = h / natom
r_mid = np.linspace(distance.min(), distance.max(), nk)
plt.plot(r_mid, h)
    # plt.plot(t, pos_1particle_x)
    # plt.xlabel("Time")
    # plt.ylabel("X axis position")
# plt.yscale("log")
plt.xlabel("Distence from center particle")
plt.ylabel("$P(r)$")
plt.title(f"{os.getcwd().split('/')[-2]}_Xdirection")
os.chdir("..")
plt.savefig(f"{os.getcwd().split('/')[-1]}_Xdirection.png", dpi=600, bbox_inches='tight')
plt.show()
