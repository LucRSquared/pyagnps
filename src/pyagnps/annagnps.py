from pathlib import Path
import copy
from pyagnps.utils import find_rows_containing_pattern, copy_files_from_dir_to_dir, get_relative_path, get_bounds_from_cells
from pyagnps import constants

from multiprocessing.pool import ThreadPool
from functools import partial

from tqdm.auto import tqdm

import pandas as pd

def make_df_reaches_valid(df_reaches):
    reaches = set(df_reaches['reach_id'])
    receiving_reaches = set(df_reaches['receiving_reach'])

    outlet_reach = list(receiving_reaches - reaches)[0]

    if outlet_reach == 'OUTLET':
        return df_reaches
    else:
        outlet_row = df_reaches[df_reaches['receiving_reach']==outlet_reach].copy()
        outlet_row['reach_id'] = outlet_reach
        outlet_row['receiving_reach'] = 'OUTLET'
        outlet_row['length'] = 0

        df_reaches_valid = pd.concat([outlet_row, df_reaches], ignore_index=True)
        return df_reaches_valid

def make_annagnps_inputs_dirs(output_folder=Path().cwd(), subdirs=['general', 'climate', 'simulation', 'watershed', 'GIS']):
    output_folder.mkdir(exist_ok=True, parents=True)
    subdirs_paths = []
    for subdir in subdirs:
        category_dir = output_folder / subdir
        category_dir.mkdir(exist_ok=True)
        subdirs_paths.append(category_dir)
    return subdirs_paths

def format_mgmt_operation_for_output(df):
    df['Effect_Code_01'] = df['Effect_Code_01'].astype('Int64')
    df['Effect_Code_02'] = df['Effect_Code_02'].astype('Int64')
    df['Effect_Code_03'] = df['Effect_Code_03'].astype('Int64')
    df['Effect_Code_04'] = df['Effect_Code_04'].astype('Int64')
    df['Effect_Code_05'] = df['Effect_Code_05'].astype('Int64')
    return df

def format_mgmt_schedule_for_output(df):
    df['Event_Year'] = df['Event_Year'].astype('Int64')
    df['Event_Month'] = df['Event_Month'].astype('Int64')
    df['Event_Day'] = df['Event_Day'].astype('Int64')

    return df

def check_cell_data(df):
    # Check if values in the sheet_flow_slope column are greater than 3 and if so set them to 3
    df.loc[df['sheet_flow_slope'] > 3, 'sheet_flow_slope'] = 3
    df.loc[df['sheet_flow_slope'] < 0.00001, 'sheet_flow_slope'] = 0.00001

    df.loc[df['shallow_conc_flow_slope'] > 3, 'shallow_conc_flow_slope'] = 3
    df.loc[df['shallow_conc_flow_slope'] < 0.00001, 'shallow_conc_flow_slope'] = 0.00001
    
    return df

def check_soil_layers(df):
    # Check if values in the Layer_Depth column are greater than 3000 and if so set them to 3000
    df.loc[df['Layer_Depth'] > 3000, 'Layer_Depth'] = 3000
    return df

# FUNCTIONS FOR POST PROCESSING ANNAGNPS OUTPUTS

def read_all_annagnps_output_files(output_folder, prepare_for_db=False, thuc_id=''):
    """
    Reads all .csv files in the output folder and returns dataframes

    Args:
        output_folder (str): Path to the output folder containing .out files.
        prepare_for_db (bool): Whether to prepare the data for the database. Defaults to False.

        
    # Read all .csv files in the CSV_Output_Files in the root folder
    """

    processed_outputs = {}

    for file in output_folder.glob('CSV_Output_Files/UA_RR_Output/*.csv'):
        name = file.with_suffix('').name

        if name.startswith('AnnAGNPS_AA_Sediment_erosion_UA_RR_Total_All') or \
           name.startswith('AnnAGNPS_AA_Sediment_erosion_UA_RR_Total_All') or \
           name.startswith('AnnAGNPS_AA_Sediment_erosion_UA_RR_Total_Clay') or \
           name.startswith('AnnAGNPS_AA_Sediment_erosion_UA_RR_Total_Gully') or \
           name.startswith('AnnAGNPS_AA_Sediment_erosion_UA_RR_Total_Lg_Agg') or \
           name.startswith('AnnAGNPS_AA_Sediment_erosion_UA_RR_Total_Sand') or \
           name.startswith('AnnAGNPS_AA_Sediment_erosion_UA_RR_Total_Silt') or \
           name.startswith('AnnAGNPS_AA_Sediment_erosion_UA_RR_Total_Sm_Agg') or \
           name.startswith('AnnAGNPS_AA_Sediment_erosion_UA_RR_Total_SnR_Gly_Pnd') or \
           name.startswith('AnnAGNPS_AA_Sediment_erosion_UA_RR_Total_SnR') or \
           name.startswith('AnnAGNPS_AA_Sediment_yield_UA_RR_Total_All') or \
           name.startswith('AnnAGNPS_AA_Sediment_yield_UA_RR_Total_Clay') or \
           name.startswith('AnnAGNPS_AA_Sediment_yield_UA_RR_Total_Gully') or \
           name.startswith('AnnAGNPS_AA_Sediment_yield_UA_RR_Total_Sand') or \
           name.startswith('AnnAGNPS_AA_Sediment_yield_UA_RR_Total_Silt') or \
           name.startswith('AnnAGNPS_AA_Sediment_yield_UA_RR_Total_SnR_Gly_Pnd') or \
           name.startswith('AnnAGNPS_AA_Sediment_yield_UA_RR_Total_SnR') or \
           name.startswith('AnnAGNPS_AA_Water_yield_UA_RR_Total'):
            
            processed_outputs[name] = pd.read_csv(file, index_col=False)

    # Read all .csv files in the root folder
    for file in output_folder.glob('*.csv'):
        name = file.with_suffix('').name

        if name == 'AnnAGNPS_AA':
            df_aa_n, df_aa_oc, df_aa_p = read_annagnps_aa_file(file)

            # Merge the dataframes

            df_n_oc_p = df_aa_n.merge(df_aa_oc).merge(df_aa_p)

            processed_outputs['AnnAGNPS_AA'] = df_n_oc_p
            
    # Concatenate all the dataframes that contain "AnnAGNPS_AA_Sediment_erosion_UA_RR_Total"
    df_tmp = pd.concat([processed_outputs.pop(x) for x in list(processed_outputs.keys()) if "Sediment_erosion_UA_RR_Total" in x])
    processed_outputs["AnnAGNPS_AA_Sediment_erosion_UA_RR_Total"] = df_tmp.copy(deep=True)

    # Concatenate all the dataframes that contain "AnnAGNPS_AA_Sediment_yield_UA_RR_Total"
    df_tmp = pd.concat([processed_outputs.pop(x) for x in list(processed_outputs.keys()) if "Sediment_yield_UA_RR_Total" in x])
    # Remove rows where "Source" is "Subtotal"
    # df_tmp = df_tmp[df_tmp["Source"] != "Subtotal"]
    
    processed_outputs["AnnAGNPS_AA_Sediment_yield_UA_RR_Total"] = df_tmp.copy(deep=True)

    # Concatenate all the dataframes that contain "AnnAGNPS_AA_Water_yield_UA_RR_Total"
    processed_outputs["AnnAGNPS_AA_Water_yield_UA_RR_Total"] = pd.concat([processed_outputs.pop(x) for x in list(processed_outputs.keys()) if "Water_yield_UA_RR_Total" in x])
    
    if prepare_for_db:
        for name, df in processed_outputs.items():
            if name == "AnnAGNPS_AA":
                rename_dict_aa = {
                    'Cell ID': 'cell_id',
                    'Receiving Reach ID': 'receiving_reach_id',
                    'Drainage Area [ha]': 'drainage_area_ha',
                    'Attached N [kg/ha/yr  ]': 'attached_n_kg_ha_yr',
                    'Dissolved N [kg/ha/yr  ]': 'dissolved_n_kg_ha_yr',
                    'Attached OC [kg/ha/yr  ]': 'attached_oc_kg_ha_yr',
                    'Dissolved OC [kg/ha/yr  ]': 'dissolved_oc_kg_ha_yr',
                    'Attached P [kg/ha/yr  ]': 'attached_p_kg_ha_yr',
                    'Dissolved P [kg/ha/yr  ]': 'dissolved_p_kg_ha_yr'
                }

                df = df.drop(columns=[x for x in df.columns if 'Subtotal' in x])
                df = df.rename(columns=rename_dict_aa)
                df = df.assign(thuc_id=thuc_id)

            elif name == "AnnAGNPS_AA_Sediment_yield_UA_RR_Total":
                rename_dict_sediment_yield = {
                    'Description':                                     'description',
                    'Total_Yield_by_Mass_Mg_per_yr':                   'total_yield_by_mass_mg_per_yr',
                    'Cell_ID':                                         'cell_id',
                    'Rank':                                            'rank',
                    'Cell_Drainage_Area_hectare':                      'cell_drainage_area_hectare',
                    'Accumulated_Contributing_Cell_Area_percent':      'accumulated_contributing_cell_area_percent',
                    'Cell_Yield_by_Unit_Area_Mg_per_hectare_per_year': 'cell_yield_by_unit_area_mg_per_hectare_per_year',
                    'Yield_Unit_Area_Ranking_Ratio':                   'yield_unit_area_ranking_ratio',
                    'Cell_Yield_by_Mass_Mg_per_yr':                    'cell_yield_by_mass_mg_per_yr',
                    'Cell_Yield_by_Mass_percent':                      'cell_yield_by_mass_percent',
                    'Accumulated_Contributing_Cell_Yield_percent':     'accumulated_contributing_cell_yield_percent'
                }

                df = df.drop(columns=['Reach_ID', 'Reach_Location'], errors='ignore')
                df = df.rename(columns=rename_dict_sediment_yield)
                df = df.assign(thuc_id=thuc_id)

            elif name == "AnnAGNPS_AA_Sediment_erosion_UA_RR_Total":
                rename_dict_sediment_erosion = {
                    'Description':                                       'description',
                    'Total_Erosion_by_Mass_Mg_per_yr':                   'total_erosion_by_mass_mg_per_yr',
                    'Cell_ID':                                           'cell_id',
                    'Rank':                                              'rank',
                    'Cell_Drainage_Area_hectare':                        'cell_drainage_area_hectare',
                    'Accumulated_Contributing_Cell_Area_percent':        'accumulated_contributing_cell_area_percent',
                    'Cell_Erosion_by_Unit_Area_Mg_per_hectare_per_year': 'cell_erosion_by_unit_area_mg_per_hectare_per_year',
                    'Erosion_Unit_Area_Ranking_Ratio':                   'erosion_unit_area_ranking_ratio',
                    'Cell_Erosion_by_Mass_Mg_per_yr':                    'cell_erosion_by_mass_mg_per_yr',
                    'Cell_Erosion_by_Mass_percent':                      'cell_erosion_by_mass_percent',
                    'Accumulated_Contributing_Cell_Erosion_percent':     'accumulated_contributing_cell_erosion_percent'
                }

                df = df.drop(columns=['Reach_ID', 'Reach_Location'], errors='ignore')
                df = df.rename(columns=rename_dict_sediment_erosion)
                df = df.assign(thuc_id=thuc_id)

            elif name == "AnnAGNPS_AA_Water_yield_UA_RR_Total":
                rename_dict_water_yield = {
                    'Description':                                 'description',
                    'Total_Yield_by_Mass_Mg_per_yr':               'total_yield_by_mass_mg_per_yr',
                    'Cell_ID':                                     'cell_id',
                    'Rank':                                        'rank',
                    'Cell_Drainage_Area_hectare':                  'cell_drainage_area_hectare',
                    'Accumulated_Contributing_Cell_Area_percent':  'accumulated_contributing_cell_area_percent',
                    'Cell_Yield_by_Unit_Area_mm_per_year':         'cell_yield_by_unit_area_mm_per_year',
                    'Yield_Unit_Area_Ranking_Ratio':               'yield_unit_area_ranking_ratio',
                    'Cell_Yield_by_Mass_Mg_per_yr':                'cell_yield_by_mass_mg_per_yr',
                    'Cell_Yield_by_Mass_percent':                  'cell_yield_by_mass_percent',
                    'Accumulated_Contributing_Cell_Yield_percent': 'accumulated_contributing_cell_yield_percent'
                }

                df = df.drop(columns=['Reach_ID', 'Reach_Location'], errors='ignore')
                df = df.rename(columns=rename_dict_water_yield)
                df = df.assign(thuc_id=thuc_id)

            df = df[df['cell_id'] != 0].copy()
            processed_outputs[name] = df.copy(deep=True)

    return processed_outputs
        

def read_annagnps_aa_file(aa_file):
    """
    Reads an AnnAGNPS AgroAssessment file and returns dataframes for nitrogen, organic carbon, and phosphorus.

    Args:
        aa_file (str): Path to the AnnAGNPS AgroAssessment file.

    Returns:
        tuple: A tuple of pandas DataFrames for nitrogen, organic carbon, and phosphorus.
    """
    # Find the rows with the header and total rows for each chemical
    header_rows = find_rows_containing_pattern(aa_file, 'Cell ID,Receiving Reach ID', skiprows=0)
    last_rows = find_rows_containing_pattern(aa_file, 'Watershed,Totals', skiprows=0)

    # Initialize the dataframes
    df_n = df_oc = df_p = None

    # Read the data from each chemical
    for chem, n_header, n_last in zip(['nitrogen', 'organic_carbon', 'phosphorus'], header_rows, last_rows):
        # Read the data for each chemical into separate dataframes
        match chem:
            case 'nitrogen':
                df_n = pd.read_csv(aa_file, skiprows=n_header, nrows=n_last-n_header-1, index_col=False)
                df_n  = df_n[df_n['Cell ID'] != 0].copy()
            case 'organic_carbon':
                df_oc = pd.read_csv(aa_file, skiprows=n_header, nrows=n_last-n_header-1, index_col=False)
                df_oc = df_oc[df_oc['Cell ID'] != 0].copy()
            case 'phosphorus':
                df_p = pd.read_csv(aa_file, skiprows=n_header, nrows=n_last-n_header-1, index_col=False)
                df_p  = df_p[df_p['Cell ID'] != 0].copy()        
    
    # Return the dataframes
    return df_n, df_oc, df_p

def compute_dominant_storm_type(scs_storm_types, bounds):
    
    bounds_scs = bounds.overlay(scs_storm_types)
    bounds_scs['area'] = bounds_scs.to_crs('epsg:3857').geometry.area

    main_storm_type = bounds_scs.loc[bounds_scs['area'].argmax(), 'SCS Zone Type']

    if main_storm_type == 'I':
        main_storm_type = 'Std. SCS Type I'
    elif main_storm_type == 'IA':
        main_storm_type = 'Std. SCS Type Ia'
    elif main_storm_type == 'II':
        main_storm_type = 'Std. SCS Type II'
    elif main_storm_type == 'III':
        main_storm_type = 'Std. SCS Type III'
    else:
        # default
        main_storm_type = 'Std. SCS Type II'

    return main_storm_type

def compute_weighted_precip_zones_parameters(precip_zones, bounds):

    bounds_precip = bounds.overlay(precip_zones)
    bounds_precip['area'] = bounds_precip.to_crs('epsg:3857').geometry.area

    weighted_R_fctr = (bounds_precip['area'] * bounds_precip['R_factor']).sum() / bounds_precip['area'].sum()
    weighted_10_year_EI = (bounds_precip['area'] * bounds_precip['10_year_EI']).sum() / bounds_precip['area'].sum()

    dominant_EI = bounds_precip.loc[bounds_precip['area'].argmax(), 'EI_Zone']

    if dominant_EI == 'default':
        dominant_EI = constants.DEFAULT_EI_NUMBER # 100 
    else:
        dominant_EI = int(dominant_EI.replace('US_',''))

    return weighted_R_fctr, weighted_10_year_EI, dominant_EI

def fragment_watershed(annagnps_dir, mini_watersheds_dir, **kwargs):
    """
    Reads an AnnAGNPS master file in a specified directory and reads the master file
    and creates mini watersheds.

    Args:
        annagnps_dir (str): Path to the directory containing the AnnAGNPS master file.
        mini_watersheds_dir (str): Path to the directory where the mini watersheds will be saved.
        shared_climate_dir : 
            - If left empty the original climate directory will be used.
            - (str): Path to the directory containing the shared climate files.
            - 'shared' : The shared climate directory will be used.

        num_processes (int): Number of processes to use. Defaults to 8. 
        cells_geometry (GeoDataFrame) : GeoDataFrame containing the cell geometries.
        share_global_watershed_climate_params (bool): Whether to recompute the climate parameters. 
            - Defaults to True. If True will recompute R_fctr, EI, 10_year_EI and storm type.
        precip_zones (GeoDataFrame) : GeoDataFrame containing the precip zone geometries.
        scs_storm_types (GeoDataFrame) : GeoDataFrame containing the SCS storm type geometries.
        
    """

    num_processes = kwargs.get('num_processes', 8)

    cells_geometry = kwargs.get('cells_geometry', None)
    precip_zones = kwargs.get('precip_zones', None)
    scs_storm_types = kwargs.get('scs_storm_types', None)
    share_global_watershed_climate_params = kwargs.get('share_global_watershed_climate_params', True)

    if (scs_storm_types is None) or (precip_zones is None):
        share_global_watershed_climate_params = True

    recompute_watershed_climate_parameters = not(share_global_watershed_climate_params)

    p_og_climate_dir = annagnps_dir / 'climate' # Path to original climate files
    shared_climate_dir = kwargs.get('shared_climate_dir', p_og_climate_dir)

    # Read the master file
    p_og_master = annagnps_dir / 'annagnps_master.csv'

    df_og_master = pd.read_csv(p_og_master)

    og_files = {key:annagnps_dir / Path(val) for key, val in zip(df_og_master['Data Section ID'], df_og_master['File Name'])}

    # Things that need to be particularized
    df_og_reaches = kwargs.get('df_og_reaches')
    if df_og_reaches is None:
        df_og_reaches = pd.read_csv(og_files['Reach Data'], low_memory=False)

    df_og_cells = kwargs.get('df_og_cells')
    if df_og_cells is None:
        df_og_cells = pd.read_csv(og_files['Cell Data'])

    df_soil = kwargs.get('df_soil')
    if df_soil is None:
        df_soil = pd.read_csv(og_files['Soil Data'])

    df_soil_layers = kwargs.get('df_soil_layers')
    if df_soil_layers is None:
        df_soil_layers = pd.read_csv(og_files['Soil Layer Data'])

    df_mgmt_field = kwargs.get('df_mgmt_field')
    if df_mgmt_field is None:
        df_mgmt_field = pd.read_csv(og_files['Management Field Data'])

    df_mgmt_oper = kwargs.get('df_mgmt_oper')
    if df_mgmt_oper is None:
        df_mgmt_oper = pd.read_csv(og_files['Management Operation Data'])

    df_mgmt_schedule = kwargs.get('df_mgmt_schedule')
    if df_mgmt_schedule is None:
        df_mgmt_schedule = pd.read_csv(og_files['Management Schedule Data'])

    df_crop = kwargs.get('df_crop')
    if df_crop is None:
        df_crop = pd.read_csv(og_files['Crop Data'])

    df_crop_growth = kwargs.get('df_crop_growth')
    if df_crop_growth is None:
        df_crop_growth = pd.read_csv(og_files['Crop Growth Data'])

    df_non_crop = kwargs.get('df_non_crop')
    if df_non_crop is None:
        df_non_crop = pd.read_csv(og_files['Non-Crop Data'])

    df_roc = kwargs.get('df_roc')
    if df_roc is None:
        df_roc = pd.read_csv(og_files['Runoff Curve Number Data'])

    df_globalfac = kwargs.get('df_globalfac')
    if df_globalfac is None:
        df_globalfac = pd.read_csv(og_files['Global IDs Factors and Flags Data'])

    df_out_opts_aa = kwargs.get('df_out_opts_aa')
    if df_out_opts_aa is None:
        df_out_opts_aa = pd.read_csv(og_files['Output Options - AA'])

    df_out_opts_tbl = kwargs.get('df_out_opts_tbl')
    if df_out_opts_tbl is None:
        df_out_opts_tbl = pd.read_csv(og_files['Output Options - TBL'])

    df_out_opts_global = kwargs.get('df_out_opts_global')
    if df_out_opts_global is None:
        df_out_opts_global = pd.read_csv(og_files['Output Options - Global'])

    df_watershed_data = kwargs.get('df_watershed_data')
    if df_watershed_data is None:
        df_watershed_data = pd.read_csv(og_files['Watershed Data'])

    df_sim_period = kwargs.get('df_sim_period')
    if df_sim_period is None:
        df_sim_period = pd.read_csv(og_files['Simulation Period Data'])

    df_annaid = kwargs.get('df_annaid')
    if df_annaid is None:
        # Enforce column Output_Units and Input_Units column are read as integers
        df_annaid = pd.read_csv(og_files['AnnAGNPS ID'], dtype={'Output_Units': int, 'Input_Units': int})

    # Do some reformatting
    df_mgmt_oper = format_mgmt_operation_for_output(df_mgmt_oper)
    df_mgmt_schedule = format_mgmt_schedule_for_output(df_mgmt_schedule)

    
    master_output_dir = annagnps_dir
    master_output_dir.mkdir(exist_ok=True, parents=True)

    # Folder that will contain shared data
    master_output_dir_shared = master_output_dir / 'shared'
    master_output_dir_shared.mkdir(exist_ok=True)

    # Folder that will contain all the mini watersheds
    mini_watersheds_dir.mkdir(exist_ok=True)

    # Handling the shared files
    new_master_file_template = {}

    # Iterable containing DataFrame names and corresponding filenames
    dataframes_and_filenames = [
        ("AnnAGNPS ID",                       df_annaid,          "annaid.csv",              False),
        ("Soil Data",                         df_soil,            "soil_data.csv",           True),
        ("Soil Layer Data",                   df_soil_layers,     "soil_layers_data.csv",    True),
        ("Management Field Data",             df_mgmt_field,      "management_field.csv",    False),
        ("Management Operation Data",         df_mgmt_oper,       "management_oper.csv",     False),
        ("Management Schedule Data",          df_mgmt_schedule,   "management_schedule.csv", False),
        ("Crop Data",                         df_crop,            "crop_data.csv",           False),
        ("Crop Growth Data",                  df_crop_growth,     "crop_growth.csv",         False),
        ("Non-Crop Data",                     df_non_crop,        "non_crop.csv",            False),
        ("Runoff curve Number Data",          df_roc,             "runoffcurve.csv",         False),
        ("Global IDs Factors and Flags Data", df_globalfac,       "globfac.csv",             recompute_watershed_climate_parameters),
        ("Output Options - AA",               df_out_opts_aa,     "outopts_aa.csv",          False),
        ("Output Options - TBL",              df_out_opts_tbl,    "outopts_tbl.csv",         False),
        ("Output Options - Global",           df_out_opts_global, "outopts_global.csv",      False),
        ("Watershed Data",                    df_watershed_data,  "watershed_data.csv",      False),
        ("Simulation Period Data",            df_sim_period,      "sim_period.csv",          recompute_watershed_climate_parameters),
    ]

    dummy_sim_dir = mini_watersheds_dir / 'dummy_reach'

    for master_key, df_name, filename, recompute in dataframes_and_filenames:
        df = df_name  # Access DataFrame using its name

        if df.empty: # if dataframe is empty no need to add it to the master file
            continue

        if recompute:
            new_path = filename
            new_master_file_template[master_key] = f"./{filename}"
        else:
            new_path = master_output_dir_shared / filename
            new_master_file_template[master_key] = "./" + str(get_relative_path(dummy_sim_dir,new_path)).replace("\\","/")

            df.to_csv(new_path, index=False)


    # CLIMATE

    
    # Copy climate files

    if shared_climate_dir == 'shared':
        shared_climate_dir = master_output_dir_shared / 'climate'
    elif shared_climate_dir is None:
        shared_climate_dir = p_og_climate_dir
    else:
        shared_climate_dir = Path(shared_climate_dir)


    shared_climate_dir.mkdir(exist_ok=True)

    if p_og_climate_dir != shared_climate_dir:
        print('Copying climate files...')
        copy_files_from_dir_to_dir(p_og_climate_dir, shared_climate_dir);
    
    new_master_file_template['CLIMATE DATA - DAILY'] = "./" + str(get_relative_path(dummy_sim_dir,shared_climate_dir / 'climate_daily.csv')).replace("\\","/")
    new_master_file_template['CLIMATE DATA - STATION'] = "./" + str(get_relative_path(dummy_sim_dir,shared_climate_dir / 'climate_station.csv')).replace("\\","/")

    # Filter the reaches
    df_og_reaches = df_og_reaches[df_og_reaches['length'] != 0]

    # Get the unique reach IDs
    reach_ids = df_og_reaches['reach_id'].unique()

    func = partial(make_mini_watershed_reach_cell_data_section, df_og_reaches=df_og_reaches, 
                   df_og_cells=df_og_cells, 
                   df_soil=df_soil, 
                   df_soil_layers=df_soil_layers,
                   df_sim_period=df_sim_period,
                   df_globalfac=df_globalfac,
                   annagnps_master_template_dict=new_master_file_template, 
                   mini_watersheds_dir=mini_watersheds_dir,
                   recompute_watershed_climate_parameters=recompute_watershed_climate_parameters,
                   gdf_cells=cells_geometry,
                   precip_zones=precip_zones,
                   scs_storm_types=scs_storm_types)

    # Call the parallel processing function

    # df_og_reaches = df_og_reaches[df_og_reaches['length']!=0]
    total_tasks = len(reach_ids)
    # with Pool() as pool:
    mini_watersheds = []
    with ThreadPool(processes=num_processes) as pool:
        for mini_watershed in tqdm(pool.imap_unordered(func, reach_ids), total=total_tasks, desc="Processing reaches"):
            mini_watersheds.append(mini_watershed)

    # Return list of Pathlib objects containing the mini_watersheds source files
    # Sort mini_watersheds by increasing reach_id knowing that the reach_ids are in the format f"reach_{reach_id:010.0f}"
    mini_watersheds = sorted(mini_watersheds, key=lambda x: int(x.name.split("_")[1]))

    return mini_watersheds

def make_mini_watershed_reach_cell_data_section(reach_id, df_og_reaches, df_og_cells, df_soil, df_soil_layers, 
                                                annagnps_master_template_dict, mini_watersheds_dir, 
                                                recompute_watershed_climate_parameters=False, 
                                                df_sim_period=None, df_globalfac=None,
                                                gdf_cells=None, precip_zones=None, scs_storm_types=None):
    """
    This function takes the dataframes of the AnnAGNPS simulation and generates mini watersheds
    for each reach under the mini_watersheds_dir directory
    """

    if (df_sim_period is None) or (df_globalfac is None) or (gdf_cells is None) or (precip_zones is None) or (scs_storm_types is None):
        recompute_watershed_climate_parameters = False

    
    # print(f'Processing reach {reach_id}')
    mini_watershed = mini_watersheds_dir / f"reach_{reach_id:010.0f}"

    if mini_watershed.exists():
        # print(f'Skipping reach {reach_id}')
        return mini_watershed
    else:
        mini_watershed.mkdir(exist_ok=True)
    
    df_og_reaches = df_og_reaches[df_og_reaches['length']!=0]
    
    dfr = df_og_reaches[df_og_reaches['reach_id'].eq(reach_id)]
    dfr_valid = make_df_reaches_valid(dfr)

    df_contributing_cells = df_og_cells[df_og_cells['reach_id']==reach_id]

    path_to_cells = mini_watershed / 'cell_data_section.csv'
    path_to_reaches = mini_watershed / 'reach_data_section.csv'
    path_to_soil_data = mini_watershed / 'soil_data.csv'
    path_to_soil_layers_data = mini_watershed / 'soil_layers_data.csv'

    # Modify master file and everything will be the 
    annagnps_master = copy.deepcopy(annagnps_master_template_dict)
    annagnps_master['Reach Data'] = './reach_data_section.csv'
    annagnps_master['Cell Data'] = './cell_data_section.csv'
    annagnps_master['Soil Data'] = './soil_data.csv'
    annagnps_master['Soil Layer Data'] = './soil_layers_data.csv'

    if recompute_watershed_climate_parameters:
        annagnps_master['Simulation Period Data'] = './sim_period.csv'
        annagnps_master['Global IDs Factors and Flags Data'] = './globfac.csv'

        gdf_cells = gdf_cells[gdf_cells['cell_id'].isin(df_contributing_cells['cell_id'].tolist())]

        bounds = get_bounds_from_cells(gdf_cells)

        main_storm_type = compute_dominant_storm_type(scs_storm_types, bounds)
        
        df_globalfac.loc[0,'Wshd_Storm_Type_ID'] = main_storm_type

        df_globalfac.to_csv(mini_watershed / 'globfac.csv', index=False)
        
        weighted_R_fctr, \
        weighted_10_year_EI, \
        dominant_EI = compute_weighted_precip_zones_parameters(precip_zones, bounds)

        df_sim_period.loc[0,'Rainfall_Fctr'] = weighted_R_fctr
        df_sim_period.loc[0,'10-Year_EI'] = weighted_10_year_EI
        df_sim_period.loc[0,'EI_Number'] = dominant_EI

        df_sim_period.to_csv(mini_watershed / 'sim_period.csv', index=False)

    # Modify the cell data section
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
    
    dfr_valid = dfr_valid.apply(safe_to_numeric)
    df_contributing_cells = df_contributing_cells.apply(safe_to_numeric)

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

    return mini_watershed

def safe_to_numeric(x):
    try:
        return pd.to_numeric(x)
    except (ValueError, TypeError):
        return x