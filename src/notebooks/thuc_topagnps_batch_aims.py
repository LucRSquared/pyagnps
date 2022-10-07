import sys, os
# sys.path.append('C:/Users/Luc/projects/pyagnps/src')
sys.path.append('src')
import geopandas as gpd
import pandas as pd

from pyagnps import topagnps
from pyagnps import utils
from pyagnps.utils import log_to_file, get_current_time, remove_all_files_from_dir_except_from_list

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

bigbang = time.process_time()

path_to_TOPAGNPS_bin = '/aims/luc/bins/TopAGNPS_v6.00.b.017_release_64-bit_Linux' # absolute or with respect to a sub directory in path_to_dir
path_to_thucs = '/aims/luc/data/tophuc_S_M_40000_closed_holes_with_container_thuc_merged_bbox_area_first_kept.gpkg'
root_dir = '/aims/luc/thuc_runs_40k_SM/'

path_to_log_dir = f'{root_dir}/LOGS/'
if not(os.path.exists(path_to_log_dir) and os.path.isdir(path_to_log_dir)):
    os.makedirs(path_to_log_dir)

path_to_qc_dir = f'{root_dir}/QualityControl/'
if not(os.path.exists(path_to_qc_dir) and os.path.isdir(path_to_qc_dir)):
    os.makedirs(path_to_qc_dir)

path_to_time_log = f'{path_to_log_dir}/batch_time_log.txt'
path_to_general_log = f'{path_to_log_dir}/batch_general_log.txt'

path_to_thuc_runlist = f'{path_to_log_dir}/lmrb.csv'
path_to_thuc_faillist = f'{path_to_log_dir}/fail_list.csv'

thucs = gpd.read_file(path_to_thucs) # GeoDataFrame containing the thucs and their geometry
thucs = thucs.sort_values(by=['bbox_area_sqkm'], ascending=True)

runlist = pd.read_csv(path_to_thuc_runlist, dtype=object)



runlist = runlist.iloc[:,0].to_list() # Get the list of thucs that need to be 
# log_to_file(path_to_time_log, 'thuc,time_s') # Initialize completion time log for thucs

for _, tuc in thucs.iterrows():

    thuc_id = tuc['tophucid']

    if thuc_id not in runlist:
        continue

    start = time.process_time()

    thucid_dir_name = f'thuc_{thuc_id}_40000_SM_res_10_buff_500'
    thuc_select = thucs[thucs['tophucid']==thuc_id]

    if os.path.exists(f'{root_dir}/{thucid_dir_name}') and os.path.isdir(f'{root_dir}/{thucid_dir_name}'):
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {thuc_id}: THUC previously computed: SKIPPING')
        continue

    path_to_dir = topagnps.create_topagnps_directory(root_dir, thucid_dir_name)

    try:
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {thuc_id}: Downloading DEM')
        dem, path_to_asc = topagnps.download_dem(thuc_select, path_to_dir, name=thucid_dir_name, resolution_m=10, buffer_m=500)
    except:
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {thuc_id}: Failed to download DEM')
        log_to_file(path_to_thuc_faillist, f'{thuc_id}')
        continue

    dem_filename = path_to_asc.rsplit('/',1)[-1] # Part of the string after the last / = "thuc_1173_rest_10_m.asc"

    topagnpsXML = {'DEMPROC': 2,
                   'FORMAT': 0,
                   'CSA': 10,
                   'MSCL': 250,
                   'KEEPFILES': 1,
                   'FILENAME': dem_filename}

    now = get_current_time()
    log_to_file(path_to_general_log, f'{now}: {thuc_id}: Creating TopAGNPS control file')

    topagnps.create_topagnps_xml_control_file(topagnpsXML, path_to_dir+'/TOPAGNPS.XML')

    try:
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {thuc_id}: TopAGNPS pre-processing step')
        topagnps.run_topagnps(path_to_dir, path_to_TOPAGNPS_bin)
    except:
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {thuc_id}: Error! Failed to run pre-processing')
        log_to_file(path_to_time_log, f'{thuc_id},0')
        log_to_file(path_to_thuc_faillist, f'{thuc_id}')
        continue

    try:
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {thuc_id}: Finding outlet...')
        xout, yout, rowout, colout, raster_crs = topagnps.find_outlet_uparea_shape_intersection(path_to_dir+'/UPAREA.ASC', thuc_select)

        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {thuc_id}: Outlet found! x = {xout}, y = {yout} ({raster_crs})')
    except:
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {thuc_id}: Error! Failed to find outlet')
        log_to_file(path_to_time_log, f'{thuc_id},0')
        log_to_file(path_to_thuc_faillist, f'{thuc_id}')
        continue

    try:
        topagnpsXML = {'DEMPROC': 0,
                       'FORMAT': 0,
                       'CSA': 10,
                       'MSCL': 250,
                       'KEEPFILES': 1,
                       'OUTROW': rowout,
                       'OUTCOL': colout,
                       'READOUT': 1,
                       'FILENAME': dem_filename}

        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {thuc_id}: Updating control file')

        topagnps.create_topagnps_xml_control_file(topagnpsXML, path_to_dir+'/TOPAGNPS.XML')

        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {thuc_id}: Finishing TopAGNPS processing with outlet information')

        topagnps.run_topagnps(path_to_dir, path_to_TOPAGNPS_bin)
    except:
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {thuc_id}: Error! Failed to finish TopAGNPS processing')
        log_to_file(path_to_time_log, f'{thuc_id},0')
        log_to_file(path_to_thuc_faillist, f'{thuc_id}')
        continue

    try:
        log_to_file(path_to_general_log, f'{now}: {thuc_id}: Computing quality control for THUC')
        path_to_cell_IDs_asc = f'{path_to_dir}/AnnAGNPS_Cell_IDs.asc'
        path_to_topagnps_wrn = f'{path_to_dir}/TopAGNPS_wrn.CSV'

        quality = topagnps.quality_control_areas_vs_boundary(path_to_cell_IDs_asc, thuc_select)
        quality['thuc'] = thuc_id

        touch_edges = topagnps.check_topagnps_wrn_log(path_to_topagnps_wrn)
        quality['touching_edges'] = touch_edges

        log_to_file(f'{path_to_qc_dir}/{thuc_id}.json', json.dumps(quality, indent=2))

    except:
        now = get_current_time()
        log_to_file(path_to_general_log, f'{now}: {thuc_id}: Error! Failed to compute quality control')
        log_to_file(path_to_time_log, f'{thuc_id},0')
        log_to_file(path_to_thuc_faillist, f'{thuc_id}')
        continue

    now = get_current_time()
    log_to_file(path_to_general_log, f'{now}: {thuc_id}: Deleting unnecessary files...')

    keep_files.append(f'{dem_filename}')
    keep_files.append(f'{dem_filename}'.replace('.asc','.prj'))
    file_del_errors = remove_all_files_from_dir_except_from_list(path_to_dir, keep_files)

    end = time.process_time()

    now = get_current_time()
    log_to_file(path_to_general_log, f'{now}: {thuc_id}: Finished normally in {end-start} seconds')
    log_to_file(path_to_time_log, f'{thuc_id},{end-start}')

print(f'Finished batch processing! Overall process took {(end-bigbang)/3600} hours')