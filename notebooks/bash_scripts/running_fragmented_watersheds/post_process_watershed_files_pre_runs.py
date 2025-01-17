from pyagnps import annagnps, aims
from pyagnps.utils import log_to_file, upsert_dataframe
from pathlib import Path

import sys

import argparse

import traceback


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--credentials',                    type=str, help='Path to the credentials JSON file')
    parser.add_argument('--output_folder',                  type=str, help='Path to the output folder to store the generated files', default=Path.cwd())
    parser.add_argument('--thuc_id',                        type=str, help='THUC ID to generate files for')
    # parser.add_argument('--')
    parser.add_argument('--annagnps_aa_table',              type=str, help='AnnAGNPS AA table to use',                     default='pre_runs_annagnps_aa')
    parser.add_argument('--aa_water_yield_table',           type=str, help='AnnAGNPS AA water yield table to use',         default='pre_runs_annagnps_aa_water_yield_ua_rr_total')
    parser.add_argument('--aa_sediment_yield_table',        type=str, help='AnnAGNPS AA sediment yield table to use',      default='pre_runs_annagnps_aa_sediment_yield_ua_rr_total')
    parser.add_argument('--aa_sediment_erosion_table',      type=str, help='AnnAGNPS AA sediment erosion table to use',    default='pre_runs_annagnps_aa_sediment_erosion_ua_rr_total')
    parser.add_argument('--save_method',                    type=str, help='Where to place results (\'db\' or \'files\')', default='db')
    parser.add_argument('--log_file',                       type=str, help='Path to the log file',                         default="generate_watershed_files.log")


    args = parser.parse_args()

    credentials               = Path(args.credentials)
    output_folder             = Path(args.output_folder)

    log_file_path             = Path(args.log_file)

    thuc_id                   = args.thuc_id

    save_method               = args.save_method

    annagnps_aa_table         = args.annagnps_aa_table
    aa_water_yield_table      = args.aa_water_yield_table
    aa_sediment_yield_table   = args.aa_sediment_yield_table
    aa_sediment_erosion_table = args.aa_sediment_erosion_table

    data_output_labels = ['AnnAGNPS_AA','AnnAGNPS_AA_Sediment_erosion_UA_RR_Total','AnnAGNPS_AA_Sediment_yield_UA_RR_Total','AnnAGNPS_AA_Water_yield_UA_RR_Total']

    match save_method:
        case 'db':
            db_url = aims.create_db_url_object(credentials)
            engine = aims.create_engine(db_url)
        case 'files':
            save_results_folder = output_folder / '..' / 'post_processed_files'
            save_results_folder.mkdir(exist_ok=True, parents=True)
        case _:
            raise ValueError(f"Invalid save_method: {save_method}")


    try:
        # log_to_file(log_file_path, f"Reading AnnAGNPS files for {thuc_id}...", add_timestamp=True)

        try:
            data = annagnps.read_all_annagnps_output_files(output_folder, prepare_for_db=True, thuc_id=thuc_id)
        except Exception as e:
            log_to_file(log_file_path, f"Error reading AnnAGNPS files: {e}\n{traceback.format_exc()}", add_timestamp=True)
            sys.exit(1)

        db_table_unique_columns = {
            'pre_runs_annagnps_aa': ['thuc_id', 'cell_id'],
            'pre_runs_annagnps_aa_sediment_erosion_ua_rr_total': ['thuc_id', 'cell_id', 'description'],
            'pre_runs_annagnps_aa_sediment_yield_ua_rr_total': ['thuc_id', 'cell_id', 'description'],
            'pre_runs_annagnps_aa_water_yield_ua_rr_total': ['thuc_id', 'cell_id', 'description'],
        }

        for label, db_table in zip(data_output_labels, [annagnps_aa_table, aa_sediment_erosion_table, aa_sediment_yield_table, aa_water_yield_table]):
            df = data[label]

            match save_method:
                case 'db':
                    # Upload to DB using upsert method
                    try:
                        upsert_dataframe(engine, df, db_table, unique_columns=db_table_unique_columns[db_table])
                        # df.to_sql(db_table, engine.connect(), if_exists='append', index=False)
                    except Exception as e:
                        print(f"Error uploading {label} data to {db_table}: {e}")
                        log_to_file(log_file_path, f"Error uploading {label} data to {db_table}: {e}", add_timestamp=True)
                        sys.exit(1)
                case 'files':
                    reach_name = output_folder.absolute().name
                    table_output_folder = save_results_folder / db_table
                    table_output_folder.mkdir(exist_ok=True, parents=True)

                    df.to_parquet(table_output_folder / f"{thuc_id}_{db_table}_{reach_name}.parquet", engine='pyarrow', compression='snappy', index=False)
                case _:
                    raise ValueError(f"Invalid save_method: {save_method}")

        sys.exit(0)

    except Exception as e:
        log_to_file(log_file_path, f"Error reading AnnAGNPS files: {e}\n{traceback.format_exc()}", add_timestamp=True)
        sys.exit(1)


if __name__ == '__main__':
    main()


