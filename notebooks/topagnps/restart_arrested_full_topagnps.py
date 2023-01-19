"""
Given a master directory containing subdirectories of THUC runs,
this script will find the ones that
have a TOPAGNPS.XML control file with ROW/COL information
will delete all temporary files and restart the delineation
"""

import sys, os, shutil
from glob import glob
import socket
sys.path.append('/home/luc/projects/pyagnps/')
import geopandas as gpd
import pandas as pd

from src.pyagnps import topagnps
from src.pyagnps.utils import log_to_file, get_current_time, remove_all_files_from_dir_except_from_list, move_files_from_dir_to_dir

import time
import json

# Files we ultimately want to keep
keep_files =['AgFlow_LS_Factor.asc',
             'AgFlow_LS_Factor.prj',
             'AnnAGNPS_Cell_IDs.asc',
             'AnnAGNPS_Cell_IDs.PRJ',
             'AnnAGNPS_Reach_IDs.asc',
             'AnnAGNPS_Reach_IDs.PRJ',
             'FLOVEC.ASC',
             'FLOVEC.PRJ',
             'NETFUL.ASC',
             'NETFUL.PRJ',
             'UPAREA.ASC',
             'UPAREA.PRJ',
             'AgFlow_Reach_Data.csv',
             'AnnAGNPS_Cell_Data_Section.csv',
             'AnnAGNPS_Reach_Data_Section.csv',
             'command_line_output.txt',
             'TopAGNPS_log.CSV',
             'TopAGNPS_status.CSV',
             'TopAGNPS_wrn.CSV',
             'TOPAGNPS.XML']

bigbang = time.time()

nodename = socket.gethostname()

path_to_TOPAGNPS_bin = '/aims-nas/luc/bins/TopAGNPS_v6.00.a.018_release_64-bit_Linux' # absolute or with respect to a sub directory in path_to_dir
path_to_thucs = '/aims-nas/luc/data/tophuc_S_M_40000_closed_holes_with_container_thuc_merged_bbox_area_first_kept.gpkg'
root_dir = '/aims-nas/luc/thuc_runs_40k_SM/'
dir_runs_name = '40000_SM_res_10_buff_500' # Directory to place simulations if they were succesfully restarted and finished to completion

run_dir = '/home/luc/tmp/' # Directory where thuc runs have potentially failed and need to be restarted
if not(os.path.exists(run_dir) and os.path.isdir(run_dir)):
    os.makedirs(run_dir)

path_to_log_dir = f'{root_dir}/LOGS/'
if not(os.path.exists(path_to_log_dir) and os.path.isdir(path_to_log_dir)):
    os.makedirs(path_to_log_dir)

path_to_qc_dir = f'{root_dir}/QualityControl/'
if not(os.path.exists(path_to_qc_dir) and os.path.isdir(path_to_qc_dir)):
    os.makedirs(path_to_qc_dir)

path_to_time_log = f'{path_to_log_dir}/{nodename}_batch_time_log.txt'
path_to_general_log = f'{path_to_log_dir}/{nodename}_batch_general_log.txt'

# path_to_thuc_runlist = f'{path_to_log_dir}/lmrb.csv'
path_to_thuc_faillist = f'{path_to_log_dir}/{nodename}_fail_list.csv'

thucs = gpd.read_file(path_to_thucs) # GeoDataFrame containing the thucs and their geometry

if not(os.path.exists(path_to_time_log) and os.path.isfile(path_to_time_log)):
    log_to_file(path_to_time_log, 'thuc,time_s') # Initialize completion time log for thucs

# Get list of candidate thucs to run inside run_dir using glob
runlist = [(os.path.basename(path).split('_')[1], path) for path in glob(f'{run_dir}/thuc_*')]

files_to_keep_in_thuc_dir = ['TOPAGNPS.XML']


for thuc_id, path_to_run_dir in runlist:

    everythingwentwell = False # Initialize the variable to know if everything went well

    start = time.time()

    thucid_dir_name = os.path.basename(path_to_run_dir)
    path_to_dir = f'{root_dir}/{dir_runs_name}/{thucid_dir_name}'

    thuc_select = thucs[thucs['tophucid']==thuc_id]

    contained_files = [os.path.basename(path) for path in glob(f'{path_to_run_dir}/*')]

    # Read topagnps control file
    topagnpsXML = topagnps.read_topagnps_xml_control_file(path_to_run_dir+'/TOPAGNPS.XML')

    # Find dem_filename
    dem_filename = os.path.basename(glob(f'{path_to_run_dir}/thuc*.asc')[0])

    if 'OUTROW' and 'OUTCOL' not in topagnpsXML.keys():
        # delete path_to_run_dir and its contents and continue
        shutil.rmtree(path_to_run_dir)
        continue

    elif all([file in contained_files for file in ['UPAREA.ASC', 'UPAREA.OUT', 'RELIEF.ASC', 'RELIEF.OUT']]):
        # topagnps can be run again with READOUT option 1
        remove_all_files_from_dir_except_from_list(path_to_run_dir, ['UPAREA.OUT', 'RELIEF.OUT', dem_filename, dem_filename.replace('.asc', '.prj')])
        # no need to change topagnpsXML

    else: # some files are missing but we know the OUTROW and OUTCOL so we do a full processing from the beginning
        topagnpsXML['READOUT'] = 0
        remove_all_files_from_dir_except_from_list(path_to_run_dir, [dem_filename, dem_filename.replace('.asc', '.prj')])
        

    now = get_current_time()
    log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: [RETRY] Creating TopAGNPS control file')

    topagnps.create_topagnps_xml_control_file(topagnpsXML, path_to_run_dir+'/TOPAGNPS.XML')

    try:
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: [RETRY] TopAGNPS full processing')
        topagnps.run_topagnps(path_to_run_dir, path_to_TOPAGNPS_bin)
    except:
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: [RETRY] Error! Failed to run TopAGNPS to completion')
        log_to_file(path_to_time_log, f'{thuc_id},0')
        log_to_file(path_to_thuc_faillist, f'{thuc_id}')
        continue

    try:
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Computing quality control for THUC')
        path_to_cell_IDs_asc = f'{path_to_run_dir}/AnnAGNPS_Cell_IDs.asc'
        path_to_topagnps_wrn = f'{path_to_run_dir}/TopAGNPS_wrn.CSV'

        quality = topagnps.quality_control_areas_vs_boundary(path_to_cell_IDs_asc, thuc_select)
        quality['thuc'] = thuc_id

        touch_edges = topagnps.check_topagnps_wrn_log(path_to_topagnps_wrn)
        quality['touching_edges'] = touch_edges

        log_to_file(f'{path_to_qc_dir}/{thuc_id}.json', json.dumps(quality, indent=2))

        everythingwentwell = True

    except:
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: [RETRY] Error! Failed to compute quality control')
        log_to_file(path_to_time_log, f'{thuc_id},0')
        log_to_file(path_to_thuc_faillist, f'{thuc_id}')
        continue

    if everythingwentwell:

        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Deleting unnecessary files...')

        keep_files.append(f'{dem_filename}')
        keep_files.append(f'{dem_filename}'.replace('.asc','.prj'))
        file_del_errors = remove_all_files_from_dir_except_from_list(path_to_run_dir, keep_files)

        now = get_current_time()
        end = time.time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Finished normally in {end-start} seconds')
        log_to_file(path_to_time_log, f'{thuc_id},{end-start}')

        # Move files from run directory to output directory
        if path_to_dir != path_to_run_dir:
            now = get_current_time()
            log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Moving files to output directory...')
            move_files_erros = move_files_from_dir_to_dir(path_to_run_dir, path_to_dir)

            if len(move_files_erros) == 0:
                now = get_current_time()
                log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Files moved successfully!')
                # Remove run directory
                log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Removing run directory...')
                shutil.rmtree(path_to_run_dir)
            else:
                now = get_current_time()
                log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Error! Failed to move files from {path_to_run_dir} to {path_to_dir}')

    else:

        now = get_current_time()
        end = time.time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: [RETRIES] Finished with errors in {end-start} seconds')
        log_to_file(path_to_time_log, f'{thuc_id},{end-start}')

end = time.time()

print(f'Finished batch RETRIES! Overall process took {(end-bigbang)/3600} hours')