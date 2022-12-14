import os, shutil

thuc_dir = '/aims-nas/luc/thuc_runs_40k_SM/'


for root, dirs, files in os.walk(thuc_dir):

    if (not dirs) and ('QualityControl' not in root) and ('LOGS' not in root) and \
        ('AnnAGNPS_Cell_Data_Section.csv' not in files) or \
        ('UPAREA.ASC' not in files):
    # If there are no subdirectories and a cell data section file is not present then it means that it failed
        # print(f'Deleting {root}')
        # shutil.rmtree(root)
        # Renaming the directory to failed_thuc
        id = len(thuc_dir)
        new_dir_name = root[id::].replace('thuc_', 'failed_')
        print(f'Renaming {root[id::]} to {new_dir_name}')
        os.rename(root[id::], new_dir_name)