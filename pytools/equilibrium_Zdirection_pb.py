#!/share/simulation/Anaconda3-2022.05/bin/python
import re
import numpy as np
from simpletraj.dcd.dcd import DCDReader
import os
from xml_io import xml_parser
import pandas as pd


def equilibrium_Zdirection_pb():

    #N = int(input("N = (default 400):") or "400")
    N = 400
    # p = int(input("p_ring = (default 20):") or "20")
    p = int(re.search(r'p(\d+)b', os.path.basename(os.getcwd())).group(1))
    #xml_file_number = int(input("xml_file_number = (default 100):") or "100")
    xml_file_number = 2000

    def get_all_dir_path(path):
        dir_list = []
        for root, dirs, files in os.walk(path):
            for dir in dirs:
                dir_list.append(os.path.join(root, dir))
        return sorted(dir_list)

    z_direction_mean_list = []
    # i = 0
    for file in sorted(get_all_dir_path(os.getcwd())):
        print(file)
        # if "0relax" in file or "du" in file.split("/")[-1]:
        # i += 1
            # print(f"Processing {file}")
            # pass 
            #os.chdir(file)
            #fname = [xml_file for xl_file in sorted(os.listdir(os.getcwd())) if xml_file.startswith("particle") and xml_file.endswith(".xml")][:]
            #positions = np.reshape([xml_parser(xml).nodes['position'] for xml in fname], (len(fname), -1, 3))
            #com = positions[0, (N + p * 7):(N + (p + 1) * 7)].mean(axis=0)
            #v = positions[:, N - 1, :] - com[None, :]
            ##v = positions[:, 0, :] - com[None, :]
            #z_direction_mean = np.abs(v[-xml_file_number:].mean(axis=0))
            #print(z_direction_mean[-1])
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
            v = pos[:, N - 1, :] - com[None, :]
            #v = pos[:, 0, :] - com[None, :]
            z_direction_mean = np.abs(v[-xml_file_number:, ].mean(axis=0))
            # print(z_direction_mean)
            z_direction_mean_list.append(z_direction_mean[-1])
            #print(file)
    z_direction = np.array(z_direction_mean_list).reshape((3, -1)).T
    # df = pd.DataFrame(z_direction).to_string(index=False)
    # print(os.path.basename(os.path.abspath(__file__)))
    # print(f"N: {N}")
    # print(f"p: {p}")
    # print(f"xml_file_number: {xml_file_number}")
    # print(df)
    # print("")
    # result = df.to_numpy()
    return z_direction
