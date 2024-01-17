#!/usr/bin/env python
import numpy as np
import sys
from xml_io import uxml
import os

def Rotate(p, theta, p1, p2):
    u = (p2 - p1) / np.linalg.norm(p2 - p1)
    p = p - p1
    c, s = np.cos(theta), np.sin(theta)
    x = ((c + (1 - c) * u[0] * u[0]) * p[0] + \
         ((1 - c) * u[0] * u[1] - u[2] * s) * p[1] + \
         ((1 - c) * u[0] * u[2] + u[1] * s) * p[2])
    y = (((1 - c) * u[0] * u[1] + u[2] * s) * p[0] + \
         (c + (1 - c) * u[1] * u[1]) * p[1] + \
         ((1 - c) * u[1] * u[2] - u[0] * s) * p[2])
    z = (((1 - c) * u[0] * u[2] - u[1] * s) * p[0] + \
         ((1 - c) * u[1] * u[2] + u[0] * s) * p[1] + \
         (c + (1 - c) * u[2] * u[2]) * p[2])
    return p1 + np.array([x, y, z])


class system(object):
    def __init__(self):
        pass


def generate_ring():
    n = 7
    itheta = 2 * np.pi / n
    bondlen = 0.97
    scale = (bondlen / 2) / np.sin(itheta / 2)
    startpoint = np.array([1, 0, 0]) * scale

    p1 = np.array([0, 0, 1])
    p2 = np.array([0, 0, 0])

    onering = np.array([Rotate(startpoint, itheta * _, p1, p2) for _ in range(n)])
    S.onering = onering


chain_N, ring_m, ring_n = map(int, sys.argv[1:])
print(f"chain: {chain_N}")
print(f"ring: {ring_m}")


def generate_PR():
    chainlen_axis = chain_N
    nr = ring_m
    # nr = int(chainlen_axis * 0.25)
    bondlen_axis = 0.97

    # lx, ly, lz = 100, 100, (chainlen_axis * bondlen_axis) * 1.2
    lx, ly, lz = 80, 80, 107

    center_of_axis = np.arange(chainlen_axis) * bondlen_axis
    dshift_termial = 0.5
    center_of_axis[0] = center_of_axis[0] - dshift_termial
    center_of_axis[-1] = center_of_axis[-1] + dshift_termial

    center_of_ring = np.array([np.array([0, 0, _]) for _ in np.linspace(bondlen_axis + 3.0, (chainlen_axis - 1) * bondlen_axis - 3.0, nr)])
    position_of_axis = np.array([np.array([0, 0, _]) for _ in center_of_axis])
    position_of_ring = np.array([S.onering.copy() + _ for _ in center_of_ring]).reshape(-1, 3)

    S.position = np.r_[position_of_axis, position_of_ring]
    shift = chainlen_axis * bondlen_axis / 2.0
    S.position[:, 2] = S.position[:, 2] - shift

    bond_axis = np.array([(_, _ + 1) for _ in np.arange(chainlen_axis - 1)])
    bond_onering = np.array([[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 0]])
    bond_ring = np.array([bond_onering + _ * 7 for _ in np.arange(nr)]) + bond_axis.max() + 1

    # general info
    S.box = np.array([lx, ly, lz])

    # add types
    # S.type = np.array(['T'] + ['A'] * (chain_N - 2) + ['T'] + ['R'] * left_ring_m * ring_n + ['M'] * ring_n + ['R'] * int(nr - left_ring_m - 1), dtype=object)
    S.type = np.array(['T'] + ['A'] * (chain_N - 2) + ['T'] + ['R'] * ring_n * ring_m, dtype=object)

    # add bonds
    S.bond = np.r_[bond_axis.reshape(-1, 2), bond_ring.reshape(-1, 2)]

    # add angles
    angle_axis = np.array([(_, _ + 1, _ + 2) for _ in np.arange(chainlen_axis - 2)])
    angle_onering = np.array([(np.arange(3) + i) % ring_n for i in range(ring_n)])
    # angle_onering = np.array([[0, 1, 2], [1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 6], [5, 6, 0], [6, 0, 1]])
    angle_ring = np.array([angle_onering + _ * ring_n for _ in np.arange(nr)]) + chainlen_axis
    S.angle = np.r_[angle_axis.reshape(-1, 3), angle_ring.reshape(-1, 3)]

    # add dihedral
    dihedral_axis = np.array([(_, _ + 1, _ + 2, _ + 3) for _ in np.arange(chainlen_axis - 3)])
    dihedral_onering = np.array([(np.arange(4) + i) % ring_n for i in range(ring_n)])
    dihedral_ring = np.array([dihedral_onering + _ * ring_n for _ in np.arange(nr)]) + chainlen_axis
    S.dihedral = np.r_[dihedral_axis.reshape(-1, 4), dihedral_ring.reshape(-1, 4)]

    bondname_mapping = {
        'A-A': 'polymer',
        'T-A': 'terminal',
        'A-T': 'terminal',
        'R-R': 'R-R',
        # 添加其他可能的映射
    }
    bondnames = []
    for bond in S.bond:
        typea = S.type[bond[0]]
        typeb = S.type[bond[1]]

        bondname = typea + '-' + typeb

        if bondname in bondname_mapping:
            bondname = bondname_mapping[bondname]

        bondnames.append(bondname)

    anglename_mapping = {
        'C-R-R': 'R-R-R',
        'R-R-C': 'R-R-R',
        'R-R-R': 'R-R-R',
        'M-M-M': 'R-R-R',
        'R-C-R': 'R-C-R',
        'T-A-A': 'A-A-A',
        'A-A-A': 'A-A-A',
        'A-A-T': 'A-A-A',
        # 添加其他可能的映射
    }
    anglenames = []
    for angle in S.angle:
        typea = S.type[angle[0]]
        typeb = S.type[angle[1]]
        typec = S.type[angle[2]]
        anglename = typea + '-' + typeb + '-' + typec

        if anglename in anglename_mapping:
            anglename = anglename_mapping[anglename]

        anglenames.append(anglename)

    dihedralname_mapping = {
        'T-A-A-A': 'chain_dihedral',
        'A-A-A-A': 'chain_dihedral',
        'A-A-A-T': 'chain_dihedral',
        # 添加其他可能的映射
    }
    dihedralnames = []
    for dihedral in S.dihedral:
        typea = S.type[dihedral[0]]
        typeb = S.type[dihedral[1]]
        typec = S.type[dihedral[2]]
        typed = S.type[dihedral[3]]
        dihedralname = typea + '-' + typeb + '-' + typec + '-' + typed

        if dihedralname in dihedralname_mapping:
            dihedralname = dihedralname_mapping[dihedralname]

        dihedralnames.append(dihedralname)

    S.bond = np.c_[np.array(bondnames, dtype=object), S.bond]
    S.angle = np.c_[np.array(anglenames, dtype=object), S.angle]
    S.dihedral = np.c_[np.array(dihedralnames, dtype=object), S.dihedral]

    # write to xml
    writer = uxml.dump()
    writer.box = S.box
    writer.nodes['position'] = S.position
    writer.nodes['bond'] = S.bond
    writer.nodes['type'] = S.type
    writer.nodes['angle'] = S.angle
    writer.nodes['dihedral'] = S.dihedral

    # add diameter
    DiameterTerminalBead = 3.0
    writer.nodes['diameter'] = np.array([1.0 if _ != 'T' else DiameterTerminalBead for _ in writer.nodes['type']])

    writer.write('out.xml', xml_type='galamost')


# main
global S
S = system()
generate_ring()
generate_PR()
