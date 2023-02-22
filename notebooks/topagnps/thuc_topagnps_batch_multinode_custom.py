import sys, os, shutil
import socket
sys.path.append('/home/luc/projects/pyagnps/')
import geopandas as gpd
import pandas as pd

from src.pyagnps import topagnps
from src.pyagnps.utils import log_to_file, get_current_time, remove_all_files_from_dir_except_from_list, move_files_from_dir_to_dir, copy_files_from_dir_to_dir

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
             'valgrind_topagnps_preprocessing.txt',
             'valgrind_topagnps_full_processing.txt',
             'TopAGNPS_log.CSV',
             'TopAGNPS_status.CSV',
             'TopAGNPS_wrn.CSV',
             'TOPAGNPS.XML']

# Keep log files in case of failure
keep_files_failure = ['command_line_output.txt',
                      'valgrind_topagnps_preprocessing.txt',
                      'valgrind_topagnps_full_processing.txt',
                      'RELIEF.ASC',
                      'RELIEF.PRJ',
                      'TopAGNPS_log.CSV',
                      'TopAGNPS_status.CSV',
                      'TopAGNPS_wrn.CSV'
                      'TopAGNPS_err.CSV',
                      'TOPAGNPS.XML']


bigbang = time.time()

nodename = socket.gethostname()

CSA = 10 # Critical Source Area in ha
MSCL = 250 # Minimum Stream Channel Length in m
RESOLUTION = 30 # Resolution of the DEM in m
BUFFER = 500 # Buffer around the thuc in m

path_to_TOPAGNPS_bin = '/aims-nas/luc/bins/TopAGNPS_v6.00.b.023_release_64-bit_Linux' # absolute or with respect to a sub directory in path_to_dir
path_to_thucs = '/aims-nas/luc/data/tophuc_S_M_40000_closed_holes_with_container_thuc_merged_bbox_area_first_kept.gpkg'
root_dir = '/aims-nas/luc/thuc_runs_40k_SM/'
dir_runs_name = f'40000_SM_res_{RESOLUTION}_csa_{CSA}_mscl_{MSCL}_buff_{BUFFER}' # Name of the directory where the runs will be stored

run_dir = '/home/luc/tmp/' # Directory where the computation will be carried out
if not(os.path.exists(run_dir) and os.path.isdir(run_dir)):
    os.makedirs(run_dir)

path_to_log_dir = f'{root_dir}/LOGS/{dir_runs_name}'
if not(os.path.exists(path_to_log_dir) and os.path.isdir(path_to_log_dir)):
    os.makedirs(path_to_log_dir)

path_to_qc_dir = f'{root_dir}/QualityControl/{dir_runs_name}'
if not(os.path.exists(path_to_qc_dir) and os.path.isdir(path_to_qc_dir)):
    os.makedirs(path_to_qc_dir)

path_to_time_log = f'{path_to_log_dir}/{nodename}_batch_time_log.txt'
path_to_general_log = f'{path_to_log_dir}/{nodename}_batch_general_log.txt'

path_to_thuc_runlist = f'{path_to_log_dir}/lmrb.csv'
path_to_thuc_faillist = f'{path_to_log_dir}/{nodename}_fail_list.csv'

thucs = gpd.read_file(path_to_thucs) # GeoDataFrame containing the thucs and their geometry
thucs = thucs.sort_values(by=['bbox_area_sqkm'], ascending=True)

# runlist = thucs['tophucid'].to_list()

runlist = ['1060', '1148']

# runlist = pd.read_csv(path_to_thuc_runlist, dtype=object)
# runlist = runlist.iloc[:,0].to_list() # Get the list of thucs that need to be 

if not(os.path.exists(path_to_time_log) and os.path.isfile(path_to_time_log)):
    log_to_file(path_to_time_log, 'thuc,time_s') # Initialize completion time log for thucs


for _, tuc in thucs.iterrows():

    everythingwentwell = False # Initialize the variable to know if everything went well

    thuc_id = tuc['tophucid']

    if thuc_id not in runlist:
        continue

    start = time.time()

    # Make sure the path exists
    thucid_dir_name = f'thuc_{thuc_id}_{dir_runs_name}'
    path_to_dir = f'{root_dir}/{dir_runs_name}/{thucid_dir_name}'
    path_to_run_dir = f'{run_dir}/{thucid_dir_name}'

    if not(os.path.exists(f'{root_dir}/{dir_runs_name}') and os.path.isdir(f'{root_dir}/{dir_runs_name}')):
        os.makedirs(f'{root_dir}/{dir_runs_name}')

    if os.path.exists(path_to_dir) and os.path.isdir(path_to_dir):
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: THUC previously computed: SKIPPING')
        continue
    else:
        topagnps.create_topagnps_directory(f'{root_dir}/{dir_runs_name}', thucid_dir_name)
        topagnps.create_topagnps_directory(run_dir, thucid_dir_name)

    thuc_select = thucs[thucs['tophucid']==thuc_id]

    try:
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Downloading DEM')
        dem, path_to_asc = topagnps.download_dem(thuc_select, path_to_run_dir, name=thucid_dir_name, resolution_m=RESOLUTION, buffer_m=BUFFER)

        dem_filename = path_to_asc.rsplit('/',1)[-1] # Part of the string after the last / = "thuc_1173_rest_10_m.asc"
        print(dem_filename)

        topagnpsXML = {'DEMPROC': 2,
                    'FORMAT': 0,
                    'CSA': CSA,
                    'MSCL': MSCL,
                    'KEEPFILES': 1,
                    'FILENAME': dem_filename}

        print(topagnpsXML)

        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Creating TopAGNPS control file')

        topagnps.create_topagnps_xml_control_file(topagnpsXML, path_to_run_dir+'/TOPAGNPS.XML')

    except:
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Failed to download DEM')
        # Delete path_to_dir so that it can be reattempted some other time
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Deleting run and output directory to reattempt later')
        shutil.rmtree(path_to_dir)
        shutil.rmtree(path_to_run_dir)
        # log_to_file(path_to_thuc_faillist, f'{thuc_id}')
        #continue

    try:
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: TopAGNPS pre-processing step')

        # topagnps.run_topagnps(path_to_run_dir, path_to_TOPAGNPS_bin)
        topagnps.run_topagnps(path_to_run_dir, path_to_TOPAGNPS_bin, memtrack=True, output_memtrack_filename="valgrind_topagnps_preprocessing.txt")
    except:
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Error! Failed to run pre-processing')
        log_to_file(path_to_time_log, f'{thuc_id},0')
        log_to_file(path_to_thuc_faillist, f'{thuc_id}')
        #continue

    try:
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Finding outlet...')
        xout, yout, rowout, colout, raster_crs = topagnps.find_outlet_uparea_shape_intersection(path_to_run_dir+'/UPAREA.ASC', thuc_select)

        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Outlet found! x = {xout}, y = {yout} ({raster_crs}) / row = {rowout}, col = {colout}')
    except:
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Error! Failed to find outlet')
        log_to_file(path_to_time_log, f'{thuc_id},0')
        log_to_file(path_to_thuc_faillist, f'{thuc_id}')
        #continue

    try:
        
        # Renaming RELIEF(.ASC|.PRJ) to be used as input DEM
        new_dem_name = f'RELIEF_{dem_filename}'
        new_prj_name = new_dem_name.replace('.asc','.prj')
      
        os.rename(f'{path_to_run_dir}/RELIEF.ASC', f'{path_to_run_dir}/{new_dem_name}')
        os.rename(f'{path_to_run_dir}/RELIEF.PRJ', f'{path_to_run_dir}/{new_prj_name}')

        topagnpsXML = {'DEMPROC': 0,
                       'FORMAT': 0,
                       'CSA': CSA,
                       'MSCL': MSCL,
                       'KEEPFILES': 1,
                       'OUTROW': rowout,
                       'OUTCOL': colout,
                       'READOUT': 1,
                       'FILENAME': new_dem_name}

        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Updating control file')

        topagnps.create_topagnps_xml_control_file(topagnpsXML, path_to_run_dir+'/TOPAGNPS.XML')

        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Finishing TopAGNPS processing with outlet information')

        # topagnps.run_topagnps(path_to_run_dir, path_to_TOPAGNPS_bin)
        topagnps.run_topagnps(path_to_run_dir, path_to_TOPAGNPS_bin, memtrack=True, output_memtrack_filename="valgrind_topagnps_full_processing.txt")
    except:
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Error! Failed to finish TopAGNPS processing')
        log_to_file(path_to_time_log, f'{thuc_id},0')
        log_to_file(path_to_thuc_faillist, f'{thuc_id}')
        #continue

    try:
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Computing quality control for THUC')
        path_to_cell_IDs_asc = f'{path_to_run_dir}/AnnAGNPS_Cell_IDs.asc'
        path_to_topagnps_wrn = f'{path_to_run_dir}/TopAGNPS_wrn.CSV'

        quality = topagnps.quality_control_areas_vs_boundary(path_to_cell_IDs_asc, thuc_select)
        quality['thuc'] = thuc_id

        touch_edges = topagnps.check_topagnps_wrn_log(path_to_topagnps_wrn)
        quality['touching_edges'] = touch_edges

        log_to_file(f'{path_to_qc_dir}/{thuc_id}_res_{RESOLUTION}_csa_{CSA}_mscl_{MSCL}.json', json.dumps(quality, indent=2))

        everythingwentwell = True

    except:
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Error! Failed to compute quality control')
        log_to_file(path_to_time_log, f'{thuc_id},0')
        log_to_file(path_to_thuc_faillist, f'{thuc_id}')
        #continue

    if everythingwentwell:

        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Deleting unnecessary files...')

        keep_files_tmp = keep_files.copy()

        # Original DEM
        keep_files_tmp.append(f'{dem_filename}')
        keep_files_tmp.append(f'{dem_filename}'.replace('.asc','.prj'))

        # RELIEF map
        keep_files_tmp.append(new_dem_name)
        keep_files_tmp.append(new_prj_name)

        file_del_errors = remove_all_files_from_dir_except_from_list(path_to_run_dir, keep_files_tmp)

        now = get_current_time()
        end = time.time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Finished normally in {end-start} seconds')
        log_to_file(path_to_time_log, f'{thuc_id},{end-start}')

        # Move files from run directory to output directory
        if path_to_dir != path_to_run_dir:
            now = get_current_time()
            log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Moving files to output directory...')
            move_files_errors = move_files_from_dir_to_dir(path_to_run_dir, path_to_dir)

            if len(move_files_errors) == 0:
                now = get_current_time()
                log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Files moved successfully!')
                # Remove run directory
                log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Removing run directory...')
                shutil.rmtree(path_to_run_dir)
            else:
                now = get_current_time()
                log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Error! Failed to move files from {path_to_run_dir} to {path_to_dir}')

    else:            

        # test if path_to_dir exists if it does not it was removed before
        if not os.path.exists(path_to_dir):
            continue

        if not os.path.exists(path_to_run_dir):
            continue

        # Removing useless files from run directory
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Deleting unnecessary files...')

        keep_files_failure_tmp = keep_files.copy()

        # Original DEM
        keep_files_failure_tmp.append(f'{dem_filename}')
        keep_files_failure_tmp.append(f'{dem_filename}'.replace('.asc','.prj'))

        # RELIEF map
        keep_files_failure_tmp.append(new_dem_name)
        keep_files_failure_tmp.append(new_prj_name)
        
        file_del_errors = remove_all_files_from_dir_except_from_list(path_to_run_dir, keep_files_failure_tmp)

        # Move log files from run directory to output directory and remove run directory
        if path_to_dir != path_to_run_dir:
            now = get_current_time()
            log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Moving log files of failed delineation to output directory...')
            move_files_errors = move_files_from_dir_to_dir(path_to_run_dir, path_to_dir)

            if len(move_files_errors) == 0:
                now = get_current_time()
                log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Log files moved successfully!')
                # Remove run directory
                log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Removing run directory...')
                shutil.rmtree(path_to_run_dir)
            else:
                now = get_current_time()
                log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Error! Failed to move log files from {path_to_run_dir} to {path_to_dir}')

            
        now = get_current_time()
        end = time.time()
        log_to_file(path_to_general_log, f'{now}: {nodename}: {thuc_id}: Did not finish (errors) {end-start} seconds')
        log_to_file(path_to_time_log, f'{thuc_id},{end-start}')

end = time.time()

print(f'Finished batch processing! Overall process took {(end-bigbang)/3600} hours')