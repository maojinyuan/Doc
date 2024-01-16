#!/usr/bin/env python3
import re
import numpy as np
from simpletraj.dcd.dcd import DCDReader
import os
from xml_io import xml_parser
import pandas as pd

def average_n():
    # N = int(input("N = (default 400):") or "400")
    N = 400
    p = int(re.search(r'p(\d+)b', os.path.basename(os.getcwd())).group(1))
    # xml_file_number = int(input("xml_file_number = (default 100):") or "100")
    xml_file_number = 2000


    def get_all_dir_path(path):
        dir_list = []
        for root, dirs, files in os.walk(path):
            for dir in dirs:
                dir_list.append(os.path.join(root, dir))
        return sorted(dir_list)


    idx_mean_list = []
    for file in sorted(get_all_dir_path(os.getcwd())):
        print(file)
        # 如果报错FileNotFoundError: [Errno 2] No such file or directory: 'particle.dcd',或者rr在的话，就break
        # if not os.path.exists(file + '/particle.dcd'):
            # pass
        # elif "nointeraction" in file:
            # print(file)
        #     pass
        # else:
            # os.chdir(file)
        current_directory = file
        contents = os.listdir(file)
        folders = [item for item in contents if os.path.isdir(os.path.join(current_directory, item))]
        for folder in folders:
            dcd = DCDReader(f'{file}/{folder}/particle.dcd')
         
            pos = np.asarray([_.copy() for _ in dcd]).reshape(len(dcd), -1, 3)
            dcd.close()
            com = pos[0, N+7*p:N+7*p+7, :].mean(axis=0)
            v = pos[-xml_file_number:, :N] - com[None, :]
            idx = np.linalg.norm(v, axis=-1).argmin(axis=-1)
            # print(idx)
            idx_mean = N - np.mean(idx)
            idx_mean_list.append(idx_mean)
            # print(idx_mean)
    idx = np.array(idx_mean_list).reshape((3, -1)).T
    return idx
    # df = pd.DataFrame(idx).to_string(index=False)
    # print(os.path.basename(os.path.abspath(__file__)))
    # print(f"N: {N}")
    # print(f"p: {p}")
    # print(f"xml_file_number: {xml_file_number}")
    # print(df)
    # print("")
