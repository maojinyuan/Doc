#!/usr/bin/env python
import re
import numpy as np
import os
from list_dir import list_folders
from simpletraj.dcd.dcd import DCDReader
# from equilibrium_Zdirection_pb import equilibrium_Zdirection_pb
# from eq_0andFixRing import eq_0andFixRing
# from average_n import average_n

N = 400
xml_file_number = 3000
target_folder_names = ['p0b0', 'p10b0', 'p20b0', 'p50b0', 'p10b10', 'p20b20', 'p50b50', 'p10b20', 'p50b20', 'p20b10', 'p20b50', 'p0b10', 'p0b20', 'p0b50']

def get_all_dir_path(path):
    dir_list = []
    for root, dirs, files in os.walk(path):
        for dir in dirs:
            dir_list.append(os.path.join(root, dir))
    return sorted(dir_list)

def xr_Rl_n():
    target_folder_names_final = []
    xr_list = []
    Rl_list = []
    idx_list = []
    std_xr = []
    std_Rl = []
    std_idx = []
    for file in list_folders()[:5]:
        xr_mean_list = []
        Rl_mean_list = []
        idx_mean_list = []
        du = get_all_dir_path(file)
        for force in du:
            if force.split("/")[1] in target_folder_names:
                p = int(re.search(r'p(\d+)b', force).group(1))
                try:
                    dcd = DCDReader(f'{force}/particle.dcd')
                except FileNotFoundError:
                    continue
                pos = np.asarray([_.copy() for _ in dcd], dtype=np.float32).reshape(len(dcd), -1, 3)
                com = pos[0, N+7*p:N+7*p+7, :].mean(axis=0)

                v_xr = pos[:, N - 1, :] - com[None, :]
                xr_mean = np.abs(v_xr[-xml_file_number:, ].mean(axis=0))
                xr_mean_list.append(xr_mean[-1])

                v_Rl = pos[:, 0, :] - com[None, :]
                Rl_mean = np.linalg.norm(v_Rl[-xml_file_number:], axis=-1).mean()
                Rl_mean_list.append(Rl_mean)

                v_n = pos[-xml_file_number:, :N] - com[None, :]
                idx = np.linalg.norm(v_n, axis=-1).argmin(axis=-1)
                idx_mean = N - np.mean(idx)
                idx_mean_list.append(idx_mean)

                target_folder_names_final.append(force.split("/")[1])


        xr_mean_list = np.asarray(xr_mean_list).reshape((-1, 3, 12))
        std_xr.extend(np.std(xr_mean_list, axis=1))
        xr_list.extend(xr_mean_list.mean(axis=1))

        Rl_mean_list = np.asarray(Rl_mean_list).reshape((-1, 3, 12))
        std_Rl.extend(np.std(Rl_mean_list, axis=1))
        Rl_list.extend(Rl_mean_list.mean(axis=1))

        idx_mean_list = np.asarray(idx_mean_list).reshape((-1, 3, 12))
        std_idx.extend(np.std(idx_mean_list, axis=1))
        idx_list.extend(idx_mean_list.mean(axis=1))

    xr_list = np.array(xr_list)
    Rl_list = np.array(Rl_list)
    idx_list = np.array(idx_list)

    std_xr = np.array(std_xr)
    std_Rl = np.array(std_Rl)
    std_idx = np.array(std_idx)

    unique_elements, unique_indices = np.unique(target_folder_names_final, return_index=True)
    target_folder_names_final = unique_elements[np.argsort(unique_indices)]

    return xr_list, Rl_list, idx_list, std_xr, std_Rl, std_idx, target_folder_names_final

xr, Rl, idx, std_xr, std_Rl, std_idx, target_folder_names_final = xr_Rl_n()
print(xr.shape)
print(Rl.shape)
print(idx.shape)
print(f"xml_file_number: {xml_file_number}")
print(f"N: {N}")
print(f"target_folder_names: {target_folder_names_final}")

np.savetxt("xr.txt", xr, fmt="%10.6f")
np.savetxt("Rl.txt", Rl, fmt="%10.6f")
np.savetxt("idx.txt", idx, fmt="%10.6f")

np.savetxt("xr_std.txt", std_xr, fmt="%10.6f")
np.savetxt("Rl_std.txt", std_Rl, fmt="%10.6f")
np.savetxt("idx_std.txt", std_idx, fmt="%10.6f")
