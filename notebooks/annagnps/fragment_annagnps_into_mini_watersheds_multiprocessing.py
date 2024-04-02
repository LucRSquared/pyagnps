from pathlib import Path
import copy
import pandas as pd

from multiprocessing.pool import Pool, ThreadPool
from functools import partial
from tqdm import tqdm

import pyagnps

def make_mini_watershed_reach_cell_data_section(reach_id, df_og_reaches, df_og_cells, df_soil, df_soil_layers, annagnps_master_template_dict, mini_watersheds_dir):
    print(f'Processing reach {reach_id}')
    mini_watershed = mini_watersheds_dir / f"reach_{reach_id:010.0f}"

    if mini_watershed.exists():
        print(f'Skipping reach {reach_id}')
        return
    else:
        mini_watershed.mkdir(exist_ok=True)
    
    df_og_reaches = df_og_reaches[df_og_reaches['length']!=0]
    
    dfr = df_og_reaches[df_og_reaches['reach_id'].eq(reach_id)]
    dfr_valid = pyagnps.annagnps.make_df_reaches_valid(dfr)

    df_contributing_cells = df_og_cells[df_og_cells['reach_id']==reach_id]

    path_to_cells = mini_watershed / 'cell_data_section.csv'
    path_to_reaches = mini_watershed / 'reach_data_section.csv'
    path_to_soil_data = mini_watershed / 'soil_data.csv'
    path_to_soil_layers_data = mini_watershed / 'soil_layers_data.csv'

    annagnps_master = copy.deepcopy(annagnps_master_template_dict)
    annagnps_master['Reach Data'] = './reach_data_section.csv'
    annagnps_master['Cell Data'] = './cell_data_section.csv'
    annagnps_master['Soil Data'] = './soil_data.csv'
    annagnps_master['Soil Layer Data'] = './soil_layers_data.csv'


    if (num_cells:=df_contributing_cells.shape[0]) == 3:
        # This is a source reach we don't need to do anything special
        pass

    elif num_cells <= 2:
        # Copy first cell and append at the end:
        first_cell = df_contributing_cells.iloc[0].copy(deep=True).to_frame().T
        first_cell['cell_id'] = 0
        first_cell['cell_area'] = 1e-4
        first_cell['mgmt_field_id'] = "Water"
        first_cell['reach_location_code'] = 0
        df_contributing_cells = pd.concat([df_contributing_cells, first_cell]).reset_index(drop=True)

    else:
        raise Exception(f"Invalid number of cells ({num_cells}) for reach {reach_id}")
    
    dfr_valid = dfr_valid.apply(pd.to_numeric, errors='ignore')
    df_contributing_cells = df_contributing_cells.apply(pd.to_numeric, errors='ignore')

    # Make specific soil data sections
    specific_soils = df_contributing_cells['soil_id'].to_list()
    df_soil_specific = df_soil[df_soil['Soil_ID'].isin(specific_soils)]
    df_soil_layers_specific = df_soil_layers[df_soil_layers['Soil_ID'].isin(specific_soils)]

    df_soil_specific.to_csv(path_to_soil_data, index=False, float_format='%1.5f')
    df_soil_layers_specific.to_csv(path_to_soil_layers_data, index=False, float_format='%1.5f')

    dfr_valid.to_csv(path_to_reaches, index=False, float_format='%1.5f')
    df_contributing_cells.to_csv(path_to_cells, index=False, float_format='%1.5f')

    df_master = pd.DataFrame({
        'Data Section ID': annagnps_master.keys(),
        'File Name': annagnps_master.values()})

    df_master.to_csv(mini_watershed / 'annagnps_master.csv', index=False)

    annagnps_fil = mini_watershed / 'AnnAGNPS.fil'
    annagnps_fil.write_text('annagnps_master.csv');


if __name__ == '__main__':
    print('Script called as main')

    print('Reading data...')

    # p_og_annagnps = Path("C:/Users/Luc/projects/AnnAGNPS/AnnAGNPS_simulations_by_AIMS/pelahatchie_creek_and_bay_NLDAS2_2006_2016")
    # p_og_annagnps = Path("C:/Users/Luc/projects/AnnAGNPS/AnnAGNPS_simulations_by_AIMS/TESTS_1055_NESTED/Iowa_1055_Reach_30366")
    # p_og_annagnps = Path("C:/Users/Luc/projects/AnnAGNPS/AnnAGNPS_simulations_by_AIMS/TESTS_1055_NESTED/Iowa_1055_Reach_29996")
    # p_og_annagnps = Path("C:/Users/Luc/projects/AnnAGNPS/AnnAGNPS_simulations_by_AIMS/TESTS_1055_NESTED/Iowa_1055_Reach_28857")
    p_og_annagnps = Path("C:/Users/Luc/projects/AnnAGNPS/AnnAGNPS_simulations_by_AIMS/TESTS_1055_NESTED/Iowa_1055_Reach_5465")
    p_og_master = p_og_annagnps / 'annagnps_master.csv'

    df_og_master = pd.read_csv(p_og_master)

    og_files = {key:p_og_annagnps / Path(val) for key, val in zip(df_og_master['Data Section ID'], df_og_master['File Name'])}

    # Things that need to be particularized
    df_og_reaches      = pd.read_csv(og_files['Reach Data'], low_memory=False)
    df_og_cells        = pd.read_csv(og_files['Cell Data'])

    # Things that could be shared if they are not too big (and they generally aren't)
    df_soil            = pd.read_csv(og_files['Soil Data'])
    df_soil_layers     = pd.read_csv(og_files['Soil Layer Data'])

    df_mgmt_field      = pd.read_csv(og_files['Management Field Data'])
    df_mgmt_oper       = pd.read_csv(og_files['Management Operation Data'])
    df_mgmt_schedule   = pd.read_csv(og_files['Management Schedule Data'])

    # Do some reformatting
    df_mgmt_oper = pyagnps.annagnps.format_mgmt_operation_for_output(df_mgmt_oper)
    df_mgmt_schedule = pyagnps.annagnps.format_mgmt_schedule_for_output(df_mgmt_schedule)

    df_crop            = pd.read_csv(og_files['Crop Data'])
    df_crop_growth     = pd.read_csv(og_files['Crop Growth Data'])
    df_non_crop        = pd.read_csv(og_files['Non-Crop Data'])

    df_roc             = pd.read_csv(og_files['Runoff Curve Number Data'])

    # Things that are shared
    df_globalfac       = pd.read_csv(og_files['Global IDs Factors and Flags Data'])
    df_out_opts_aa     = pd.read_csv(og_files['Output Options - AA'])
    df_out_opts_tbl    = pd.read_csv(og_files['Output Options - TBL'])
    df_out_opts_global = pd.read_csv(og_files['Output Options - Global'])
    df_watershed_data  = pd.read_csv(og_files['Watershed Data'])
    df_sim_period      = pd.read_csv(og_files['Simulation Period Data'])
    df_annaid          = pd.read_csv(og_files['AnnAGNPS ID'])

    p_og_climate_dir = p_og_annagnps / 'climate' # Path to original climate files

    print('Creating shared files')
    # master_output_dir = Path("C:/Users/Luc/projects/AnnAGNPS/AnnAGNPS_simulations_by_AIMS/pelahatchie_creek_and_bay_NLDAS2_2006_2016_fragmented")
    master_output_dir = Path("C:/Users/Luc/projects/AnnAGNPS/AnnAGNPS_simulations_by_AIMS/TESTS_1055_NESTED_fragmented") / p_og_annagnps.name
    master_output_dir.mkdir(exist_ok=True, parents=True)

    # Folder that will contain shared data
    master_output_dir_shared = master_output_dir / 'shared'
    master_output_dir_shared.mkdir(exist_ok=True)

    # Folder that will contain all the mini watersheds
    mini_watersheds_dir = master_output_dir / 'mini_watersheds'
    mini_watersheds_dir.mkdir(exist_ok=True)

    new_master_file_template = {}

    # Iterable containing DataFrame names and corresponding filenames
    dataframes_and_filenames = [
        ("AnnAGNPS ID",                       df_annaid,          "annaid.csv"),
        # ("Soil Data",                         df_soil,            "soil_data.csv"),
        # ("Soil Layer Data",                   df_soil_layers,     "soil_layers_data.csv"),
        ("Management Field Data",             df_mgmt_field,      "management_field.csv"),
        ("Management Operation Data",         df_mgmt_oper,       "management_oper.csv"),
        ("Management Schedule Data",          df_mgmt_schedule,   "management_schedule.csv"),
        ("Crop Data",                         df_crop,            "crop_data.csv"),
        ("Crop Growth Data",                  df_crop_growth,     "crop_growth.csv"),
        ("Non-Crop Data",                     df_non_crop,        "non_crop.csv"),
        ("Runoff curve Number Data",          df_roc,             "runoffcurve.csv"),
        ("Global IDs Factors and Flags Data", df_globalfac,       "globfac.csv"),
        ("Output Options - AA",               df_out_opts_aa,     "outopts_aa.csv"),
        ("Output Options - TBL",              df_out_opts_tbl,    "outopts_tbl.csv"),
        ("Output Options - Global",           df_out_opts_global, "outopts_global.csv"),
        ("Watershed Data",                    df_watershed_data,  "watershed_data.csv"),
        ("Simulation Period Data",            df_sim_period,      "sim_period.csv"),
    ]

    dummy_sim_dir = mini_watersheds_dir / 'dummy_reach'

    for master_key, df_name, filename in dataframes_and_filenames:
        df = df_name  # Access DataFrame using its name
        new_path = master_output_dir_shared / filename

        new_master_file_template[master_key] = "./" + str(pyagnps.utils.get_relative_path(dummy_sim_dir,new_path)).replace("\\","/")
        df.to_csv(new_path, index=False)

    print('Copying climate files...')
    # Copy climate files
    shared_climate_dir = master_output_dir_shared / 'climate'
    shared_climate_dir.mkdir(exist_ok=True)
    pyagnps.utils.copy_files_from_dir_to_dir(p_og_climate_dir, shared_climate_dir);
    new_master_file_template['CLIMATE DATA - DAILY'] = "./" + str(pyagnps.utils.get_relative_path(dummy_sim_dir,shared_climate_dir / 'climate_daily.csv')).replace("\\","/")
    new_master_file_template['CLIMATE DATA - STATION'] = "./" + str(pyagnps.utils.get_relative_path(dummy_sim_dir,shared_climate_dir / 'climate_station.csv')).replace("\\","/")

    # Filter the reaches
    df_og_reaches = df_og_reaches[df_og_reaches['length'] != 0]

    # Get the unique reach IDs
    reach_ids = df_og_reaches['reach_id'].unique()

    func = partial(make_mini_watershed_reach_cell_data_section, df_og_reaches=df_og_reaches, 
                   df_og_cells=df_og_cells, 
                   df_soil=df_soil, 
                   df_soil_layers=df_soil_layers, 
                   annagnps_master_template_dict=new_master_file_template, 
                   mini_watersheds_dir=mini_watersheds_dir)

    # Set the number of processes (adjust as needed)
    num_processes = 8 # mp.cpu_count()
    # Call the parallel processing function

    df_og_reaches = df_og_reaches[df_og_reaches['length']!=0]

    # with Pool() as pool:
    with ThreadPool() as pool:
        for result in pool.imap_unordered(func, reach_ids):
            pass