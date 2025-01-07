from pyagnps import annagnps, aims
from pyagnps.utils import log_to_file, upsert_dataframe
from pathlib import Path

import pandas as pd

import sys

from sqlalchemy import text as sql_text

import argparse

import traceback

def get_thuc_num_cells_in_db_output_table(table_name, thuc_id, engine):
    match table_name:
        case 'pre_runs_annagnps_aa':
            query = f"""
            SELECT COUNT(cell_id) FROM {table_name} WHERE thuc_id = '{thuc_id}'
            """
        case 'pre_runs_annagnps_aa_water_yield_ua_rr_total':
            query = f"""
            SELECT COUNT(cell_id) FROM {table_name} WHERE thuc_id = '{thuc_id}'
            """
        case 'pre_runs_annagnps_aa_sediment_yield_ua_rr_total':
            query = f"""
            SELECT COUNT(cell_id) FROM {table_name}
            WHERE thuc_id = '{thuc_id}'
            AND description = 'Sediment_Total_All_Sources'
            """
        case 'pre_runs_annagnps_aa_sediment_erosion_ua_rr_total':
            query = f"""
            SELECT COUNT(cell_id) FROM {table_name}
            WHERE thuc_id = '{thuc_id}'
            AND description = 'Erosion_Total_All_Sources'
            """
    
    df = pd.read_sql_query(sql=sql_text(query), con=engine.connect())

    num_cells = df.iloc[0].values[0]

    return num_cells

def get_thuc_num_cells_in_db_cell_data_section(thuc_id, engine):
    query = f"""
    SELECT COUNT(cell_id) FROM thuc_{thuc_id}_annagnps_cell_data_section
    """

    df = pd.read_sql_query(sql=sql_text(query), con=engine.connect())

    num_cells = df.iloc[0].values[0]

    return num_cells


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--credentials',   type=str, help='Path to the credentials JSON file')
    parser.add_argument('--thuc_id',       type=str, help='THUC ID to check outputs for')
    # parser.add_argument('--')
    parser.add_argument('--table',         type=str, help='AnnAGNPS AA table to use',                  default='pre_runs_annagnps_aa')
    parser.add_argument('--log_file',      type=str, help='Path to the log file',                      default="generate_watershed_files.log")
    # parser.add_argument('--success_thucs', type=str, help='Path to list of success thucs',                      default="generate_watershed_files.log")


    args = parser.parse_args()

    credentials               = Path(args.credentials)

    log_file_path = Path(args.log_file)
    # success_thucs = Path(args.success_thucs)

    thuc_id  = args.thuc_id

    db_table = args.annagnps_aa_table

    db_url = aims.create_db_url_object(credentials)

    engine = aims.create_engine(db_url)


    log_to_file(log_file_path, f"Checking pre run tables for {thuc_id}...", add_timestamp=True)

    try:
        # Getting number of cells in thuc cell data section
        num_cells_in_db_cell_data_section = get_thuc_num_cells_in_db_cell_data_section(thuc_id, engine)
        log_to_file(log_file_path, f"Number of cells in cell data section for {thuc_id}: {num_cells_in_db_cell_data_section}", add_timestamp=True)

    except Exception as e:
        log_to_file(log_file_path, f"Error getting number of cells in cell data section for {thuc_id}: {e}\n{traceback.format_exc()}", add_timestamp=True)
        sys.exit(1)


    try:
        # Getting number of cells in db table
        num_cells_in_db_table = get_thuc_num_cells_in_db_output_table(db_table, thuc_id, engine)
        log_to_file(log_file_path, f"Number of cells in {db_table} for {thuc_id}: {num_cells_in_db_table}", add_timestamp=True)

    except Exception as e:
        log_to_file(log_file_path, f"Error getting number of cells in {db_table} for {thuc_id}: {e}\n{traceback.format_exc()}", add_timestamp=True)
        sys.exit(1)

    if num_cells_in_db_table != num_cells_in_db_cell_data_section:
        log_to_file(log_file_path, f"Number of cells in {db_table} for {thuc_id} does not match number of cells in cell data section for {thuc_id}", add_timestamp=True)
        sys.exit(1)
    else:
        log_to_file(log_file_path, f"Number of cells in {db_table} for {thuc_id} matches number of cells in cell data section for {thuc_id}", add_timestamp=True)
        # # Append thuc_id to success list using pathlib library
        # with success_thucs.open('a') as f:
        #     f.write(f"{thuc_id}\n")

    sys.exit(0)


if __name__ == '__main__':
    main()


