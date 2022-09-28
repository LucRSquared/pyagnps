import sys
# sys.path.append('C:/Users/Luc/projects/pyagnps/src')
sys.path.append('src')
import geopandas as gpd

from pyagnps import topagnps

import time

path_to_TOPAGNPS_bin = 'C:/Users/Luc/projects/pyagnps/src/bins/TopAGNPS_v6.00.b.017_release_64-bit_Windows.exe' # absolute or with respect to a sub directory in path_to_dir
path_to_thucs = 'inputs/thucs/tophuc_S_40000.gpkg'
root_dir = 'D:/AIMS/WBD/TopAGNPS_Watershed_Delineation/thuc_tests_increasing_hucs/'

thucs = gpd.read_file(path_to_thucs)
thucs = thucs.drop_duplicates(subset='num_contained_hucs', keep='first').sort_values(by='num_contained_hucs')
# thucs = thucs.sort_values(by=['bbox_area_sqkm'], ascending=True)

counter = 0
for _, tuc in thucs.iterrows():
    counter +=1

    if counter > 30:
        break

    start = time.time()

    thuc_id = tuc['tophucid']

    thucid_dir_name = f'thuc_{thuc_id}_40000_res_10_buff_500'
    thuc_select = thucs[thucs['tophucid']==thuc_id]

    path_to_dir = topagnps.create_topagnps_directory(root_dir, thucid_dir_name)

    dem, path_to_asc = topagnps.download_dem(thuc_select, path_to_dir, name=thucid_dir_name, resolution_m=10, buffer_m=500)

    dem_filename = path_to_asc.rsplit('/',1)[-1] # Part of the string after the last / = "thuc_1173_rest_10_m.asc"

    topagnpsXML = {'DEMPROC': 2,
                'FORMAT': 0,
                'CSA': 10,
                'MSCL': 250,
                'KEEPFILES': 1,
                'FILENAME': dem_filename}

    topagnps.create_topagnps_xml_control_file(topagnpsXML, path_to_dir+'/TOPAGNPS.XML')

    try:

        topagnps.run_topagnps(path_to_dir, path_to_TOPAGNPS_bin)

        xout, yout, rowout, colout, raster_crs = topagnps.find_outlet_uparea_shape_intersection(path_to_dir+'/UPAREA.ASC', thuc_select)

        topagnpsXML = {'DEMPROC': 0,
                'FORMAT': 0,
                'CSA': 10,
                'MSCL': 250,
                'KEEPFILES': 1,
                'OUTROW': rowout,
                'OUTCOL': colout,
                'READOUT': 1,
                'FILENAME': dem_filename}

        topagnps.create_topagnps_xml_control_file(topagnpsXML, path_to_dir+'/TOPAGNPS.XML')
        topagnps.run_topagnps(path_to_dir, path_to_TOPAGNPS_bin)
    except:
        print(f'Error handling {thuc_id}')

        original_stdout = sys.stdout
        with open('D:/AIMS/WBD/TopAGNPS_Watershed_Delineation/thuc_tests/batch_log.txt', 'a') as batch_log:
            sys.stdout = batch_log
            print(f'{thuc_id},0')
            sys.stdout = original_stdout
    else:

        end = time.time()

        original_stdout = sys.stdout
        with open('D:/AIMS/WBD/TopAGNPS_Watershed_Delineation/thuc_tests/batch_log.txt', 'a') as batch_log:
            sys.stdout = batch_log
            print(f'{thuc_id},{end-start}')
            sys.stdout = original_stdout

print('finished!')