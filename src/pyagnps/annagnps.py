from pathlib import Path

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
    output_folder.mkdir(exist_ok=True)
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