from pyagnps import annagnps, aims
from pyagnps.utils import log_to_file
from pathlib import Path

import sys

from sqlalchemy import create_engine, text as sql_text

import argparse


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--credentials',                    type=str, help='Path to the credentials JSON file')
    parser.add_argument('--output_folder',                  type=str, help='Path to the output folder to store the generated files')
    parser.add_argument('--thuc_id',                        type=str, help='THUC ID to generate files for')
    # parser.add_argument('--')
    parser.add_argument('--annagnps_aa_table',              type=str, help='AnnAGNPS AA table to use',                  default='pre_runs_annagnps_aa')
    parser.add_argument('--aa_water_yield_table',           type=str, help='AnnAGNPS AA water yield table to use',      default='pre_runs_annagnps_aa_water_yield_ua_rr_total')
    parser.add_argument('--aa_sediment_yield_table',        type=str, help='AnnAGNPS AA sediment yield table to use',   default='pre_runs_annagnps_aa_sediment_yield_ua_rr_total')
    parser.add_argument('--aa_sediment_erosion_table',      type=str, help='AnnAGNPS AA sediment erosion table to use', default='pre_runs_annagnps_aa_sediment_erosion_ua_rr_total')

    args = parser.parse_args()

    credentials               = Path(args.credentials)
    output_folder             = Path(args.output_folder)

    thuc_id                   = args.thuc_id

    annagnps_aa_table         = args.annagnps_aa_table
    aa_water_yield_table      = args.aa_water_yield_table
    aa_sediment_yield_table   = args.aa_sediment_yield_table
    aa_sediment_erosion_table = args.aa_sediment_erosion_table


    data_output_labels = ['AnnAGNPS_AA','AnnAGNPS_AA_Sediment_erosion_UA_RR_Total','AnnAGNPS_AA_Sediment_yield_UA_RR_Total','AnnAGNPS_AA_Water_yield_UA_RR_Total']

    db_url = aims.create_db_url_object(credentials)

    engine = create_engine(db_url)


    try:
        data = annagnps.read_all_annagnps_output_files(output_folder, prepare_for_db=True, thuc_id=thuc_id)
        for label, db_table in zip(data_output_labels, [annagnps_aa_table, aa_sediment_erosion_table, aa_sediment_yield_table, aa_water_yield_table]):
            df = data[label]

            # Upload to DB using append method
            try:
                df.to_sql(db_table, engine.connect(), if_exists='append', index=False)
            except Exception as e:
                print(f"Error uploading {label} data to {db_table}: {e}")
                sys.exit(1)

        sys.exit(0)

    except Exception as e:
        print(f"Error reading AnnAGNPS files: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()


