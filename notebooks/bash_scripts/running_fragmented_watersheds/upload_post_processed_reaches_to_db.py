from pyagnps import annagnps, aims
from pyagnps.utils import log_to_file, upsert_dataframe
from pathlib import Path

import pandas as pd

import sys

import argparse

import traceback


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--credentials',                            type=str, help='Path to the credentials JSON file')
    parser.add_argument('--post_processing_dir',                    type=str, help='Path to the output folder to store the generated files', required=True)
    parser.add_argument('--delete_post_processed_files_on_success', type=str, help='Delete post processed files if successful upload', default="false")
    parser.add_argument('--log_file',                               type=str, help='Path to the log file',                         default="upload_post_processed_reaches_to_db.log")


    args = parser.parse_args()

    credentials            = Path(args.credentials)
    post_processing_dir    = Path(args.post_processing_dir)

    delete_post_processed_files = args.delete_post_processed_files

    log_file_path          = Path(args.log_file)

    if delete_post_processed_files.lower() in ["true", "yes", "y","oui", "1"]:
        delete_post_processed_files = True
    else:
        delete_post_processed_files = False

    db_table_unique_columns = {
            'pre_runs_annagnps_aa': ['thuc_id', 'cell_id'],
            'pre_runs_annagnps_aa_sediment_erosion_ua_rr_total': ['thuc_id', 'cell_id', 'description'],
            'pre_runs_annagnps_aa_sediment_yield_ua_rr_total': ['thuc_id', 'cell_id', 'description'],
            'pre_runs_annagnps_aa_water_yield_ua_rr_total': ['thuc_id', 'cell_id', 'description'],
        }

    try:
        log_to_file(log_file_path, f"Reading post-processing files from {post_processing_dir}...", add_timestamp=True)

        db_url = aims.create_db_url_object(credentials)
        engine = aims.create_engine(db_url)

        for table_folder in post_processing_dir.iterdir():
            if not table_folder.is_dir():
                continue

            table_name = table_folder.name

            parquet_files = list(table_folder.glob('*.parquet'))

            if not parquet_files:
                log_to_file(log_file_path, f"No parquet files found in {table_folder}", add_timestamp=True)
                continue

            df_list = []

            for file in parquet_files:
                df = pd.read_parquet(file)
                df_list.append(df)

            df = pd.concat(df_list)

            try:
                unique_columns = db_table_unique_columns.get(table_name, None)
                upsert_dataframe(engine, df, table_name, unique_columns=unique_columns)

                if delete_post_processed_files:
                    log_to_file(log_file_path, f"Deleting post-processed files from {table_folder}...", add_timestamp=True)
                    for file in df_list:
                        file.unlink()

            except Exception as e:
                log_to_file(log_file_path, f"Error uploading data to {table_name}: {e}", add_timestamp=True)
                sys.exit(1)

        sys.exit(0)

    except Exception as e:
        log_to_file(log_file_path, f"Error reading post-processing files: {e}\n{traceback.format_exc()}", add_timestamp=True)
        sys.exit(1)


if __name__ == '__main__':
    main()

