#!/usr/bin/env python3
import numpy as np
from xml_io import xml_parser
import matplotlib.pyplot as plt
import os

def pbc(a, box):
    return a - box * np.round(a / box)

nframe = 100 + 1

os.chdir(f"{os.getcwd()}/xml")
fnames = np.sort(os.listdir())[-nframe:-1]

pos_As = []
for fname in fnames:
    xml = xml_parser(fname)
    pos = xml.nodes['position']
    pos = pbc(pos, xml.box)
    bond = xml.nodes['bond']

    bond_TA = np.unique(bond[np.where(bond[:, 0] == 'T-A')][:, -1])
    pos_A = pos[bond_TA.astype(int)]
    pos_As.append(pos_A)
pos_As = np.asarray(pos_As)

dt = 0.01
interval = 100000
time = dt * interval
t = np.arange(len(fnames)) * time
particle_idx = 54
pos_1particle_x = pos_As[:, particle_idx, 0]
pos_1particle_y = pos_As[:, particle_idx, 1]
pos_1particle_z = pos_As[:, particle_idx, 2]
plt.plot(t, pos_1particle_x, label="position x")
plt.plot(t, pos_1particle_y, label="position y")
plt.plot(t, pos_1particle_z, label="position z")
plt.xlabel("Time")
plt.ylabel("Axis position")
plt.legend()
plt.show()
