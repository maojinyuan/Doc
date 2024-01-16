import os

def list_folders(path=None):
    """
    List the subfolders in the specified path.

    Args:
        path (str): The path to list the subfolders from. Defaults to the current working directory.

    Returns:
        list: A sorted list of subfolder names.
    """
    folder_names = []
    if path is None:
        path = os.getcwd()
    for entry in os.scandir(path):
        if entry.is_dir():
            folder_names.append(entry.name)
    folder_names.sort()
    return folder_names

def find_files(file_names):
    matching_files = []
    for file in os.listdir('.'):
        if os.path.isfile(file) and any(name in file for name in file_names):
            matching_files.append(os.path.abspath(file))
    matching_files.sort()
    return matching_files

def find_specific_files(folder=None, prefix='acf', suffix='dat'):
    if folder is None:
        folder = '.'

    matching_files = []
    for file in os.listdir(folder):
        full_path = os.path.join(folder, file)
        if os.path.isfile(full_path) and file.startswith(prefix) and file.endswith(suffix):
            matching_files.append(os.path.abspath(full_path))
    matching_files.sort()
    return matching_files

def get_current_folder_info():
    current_folder = os.getcwd()
    father_folder = os.path.abspath('..')
    dcd_files = [f for f in os.listdir(current_folder) if f.endswith('.dcd')]
    dcd_file_paths = [os.path.join(current_folder, f) for f in dcd_files]
    gala_files = [f for f in os.listdir(current_folder) if f.endswith('.gala')]
    dat_files = [f for f in os.listdir(current_folder) if f.endswith('.dat')]
    return current_folder, father_folder, dcd_files, gala_files, dat_files
