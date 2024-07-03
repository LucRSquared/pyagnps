from pathlib import Path
from pyagnps.utils import find_rows_containing_pattern

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

# FUNCTIONS FOR POST PROCESSING ANNAGNPS OUTPUTS

def read_all_annagnps_output_files(output_folder):
    """
    Reads all .csv files in the output folder and returns dataframes

    Args:
        output_folder (str): Path to the output folder containing .out files.

        
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
            # processed_outputs[f"{name}_nitrogen"] = df_aa_n
            # processed_outputs[f"{name}_organic_carbon"] = df_aa_oc
            # processed_outputs[f"{name}_phosphorus"] = df_aa_p
        
        # THESE FILES ARE REDUNDANT
        
        # elif name.startswith('AnnAGNPS_AA_Water_yield_(unit-area)') or \
        #      name.startswith('AnnAGNPS_AA_Sediment_yield_(unit-area)') or \
        #      name.startswith('AnnAGNPS_AA_Sediment_erosion_(unit-area)'):
            
        #     header_row = find_rows_containing_pattern(file, 'Cell ID')[0]

        #     if name.startswith('AnnAGNPS_AA_Water_yield_(unit-area)'):
        #         last_row = find_rows_containing_pattern(file, 'Watershed Totals')[0]
        #     else:
        #         last_row = find_rows_containing_pattern(file, 'Watershed,Totals')[0]
                
        #     df_tmp = pd.read_csv(file, skiprows=header_row, 
        #                         nrows=last_row-header_row-1, 
        #                         index_col=False)
            
        #     if "Source" in df_tmp.columns:
        #         df_tmp = df_tmp[df_tmp["Source"] != "Subtotal"]

        #     processed_outputs[name] = df_tmp.copy(deep=True)
            
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
    header_rows = find_rows_containing_pattern(aa_file, 'Cell ID', skiprows=0)
    last_rows = find_rows_containing_pattern(aa_file, 'Watershed,Totals', skiprows=0)

    # Initialize the dataframes
    df_n = df_oc = df_p = None

    # Read the data from each chemical
    for chem, n_header, n_last in zip(['nitrogen', 'organic_carbon', 'phosphorus'], header_rows, last_rows):
        # Read the data for each chemical into separate dataframes
        match chem:
            case 'nitrogen':
                df_n = pd.read_csv(aa_file, skiprows=n_header, nrows=n_last-n_header-1, index_col=False)
            case 'organic_carbon':
                df_oc = pd.read_csv(aa_file, skiprows=n_header, nrows=n_last-n_header-1, index_col=False)
            case 'phosphorus':
                df_p = pd.read_csv(aa_file, skiprows=n_header, nrows=n_last-n_header-1, index_col=False)

    # Return the dataframes
    return df_n, df_oc, df_p  