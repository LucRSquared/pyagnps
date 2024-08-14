import argparse

import os
from pathlib import Path
from datetime import datetime

import itertools

import json, re

from tqdm.auto import tqdm


import numpy as np
import pandas as pd

from sqlalchemy import URL

from pyagnps import climate
# from pyagnps.utils import log_to_file, get_current_time


def main(START_DATE, END_DATE, coords, path_nldas_daily_files, path_to_creds, db_table_name, MAXITER_GLOBAL):
    """
    Parameters
    ----------
    START_DATE : str
        YYYY-MM-DD
    END_DATE : str
        YYYY-MM-DD
    coords : list
        List of tuple of coordinates [(lon, lat), ...]
    path_nldas_daily_files : str, Path
        Path to NLDAS-2 daily files
    path_to_creds : str
        Path to credentials json file
        * Must contain the entries:
        - 'user'
        - 'password'
        - 'host'
        - 'port'
        - 'database'
    """

    # DATABASE SETUP
    path_to_creds = Path(path_to_creds)

    # creds = {
    #     'aims': open_creds_dict(path_to_creds),
    #     # 'menderes': open_creds_dict(path_to_creds_menderes),
    #     'docker': {
    #             'user': 'postgres',
    #             'password': 'postgres_pass',
    #             'host': 'localhost',
    #             'port': '5432',
    #             'database': 'test_db'
    #         }
    # }

    creds = open_creds_dict(path_to_creds)


    db_url = URL.create(
                    "postgresql",
                    username=creds['user'],
                    password=creds['password'],
                    host=creds['host'],
                    port=creds['port'],
                    database=creds['database']
                    )
    # MAIN LOOP
    for iter_global in range(MAXITER_GLOBAL):
        
        try:

            nldas_clm_daily = climate.ClimateAnnAGNPSCoords(coords=coords,
                                                            start=START_DATE,
                                                            end=END_DATE,
                                                            date_mode="daily")

            nldas_clm_daily.read_nldas_daily_data(path_nldas_daily_files=path_nldas_daily_files)

            nldas_clm_daily.generate_annagnps_daily_climate_data_from_nldas_daily(
                                                            saveformat="database",
                                                            db_url=db_url,
                                                            db_table_name=db_table_name,
                                                            return_dataframes=False,
            )
            break

        except Exception as e:
            print(f"Failed to process {START_DATE} to {END_DATE}: Attempt {iter_global+1} of {MAXITER_GLOBAL}. Error: {e}")

    print(f"All done! for {START_DATE} to {END_DATE}")

def populate_parquet_to_db(**kwargs):
    """ This function explores all the files in a directory matching a certain pattern of climate station data
       and populates them to the database

    Parameters
    ----------
    - rods_dir : str, Path 
        Path to the folder where the chunks are stored (possibly in subdirectories)
    - path_to_creds : str, Path
        Path to the credentials json file
    - db_table_name : str
        Name of the DB table to update. Default is 'climate_nldas2'
    - pattern : str
        Pattern of the chunks to be explored with the glob function.
        Default is '**/climate_daily*chunk*.parquet'
    - delete_chunks_on_sucess : bool
        If True, the chunks will be deleted after being successfully processed.
        Default is False
    """

    rods_dir = Path(kwargs.get("rods_dir", None))
    path_to_creds = Path(kwargs.get("path_to_creds", None))
    db_table_name = kwargs.get("db_table_name", "climate_nldas2")
    pattern = kwargs.get("pattern", f"**/climate_daily_*_chunk*.parquet")
    delete_chunks_on_sucess = kwargs.get("delete_chunks_on_sucess", False)

    # DATABASE SETUP
    path_to_creds = Path(path_to_creds)

    creds = open_creds_dict(path_to_creds)

    db_url = URL.create(
                    "postgresql",
                    username=creds['user'],
                    password=creds['password'],
                    host=creds['host'],
                    port=creds['port'],
                    database=creds['database']
                    )

    print("Scanning for all available parquet files...") 
    all_chunks = list(rods_dir.glob(pattern))

    print("Gathering station information...")
    stations = {}
    for chunk in all_chunks:
        station_id = re.findall(r'climate_daily_(.*?)_', str(chunk))[0]
        if station_id in stations:
            stations[station_id].append(chunk)
        else:
            stations[station_id] = [chunk]

    engine_creator = lambda : climate.create_engine_with_pool(db_url, max_connections=100)

    engine = engine_creator()
    
    print("Writing files to database...")
    for station_id, chunk_files in tqdm(stations.items(), desc="Writing chunks to database", ascii=True):
    
        df_station = pd.concat([pd.read_parquet(file, engine='pyarrow') for file in chunk_files])

        lon, lat = df_station.iloc[0].lon.item(), df_station.iloc[0].lat.item()

        # Process for database insertion
        available_dates = climate.get_available_dates_for_station(station_id, engine, table=db_table_name)
        missing_dates = climate.get_missing_dates(available_dates, df_station.index.min(), df_station.index.max())
        continuous_periods = climate.find_continuous_periods(missing_dates)

        gdf_clm = None
        for period in continuous_periods:
            start, end = period[0], period[-1]
            df_period = df_station[(df_station.index >= start) & (df_station.index <= end)]
            
            if len(df_period) == 0:
                continue

            gdf_clm_period = climate.prepare_annagnps_climate_for_db(df_period, station_id, lon, lat)
            
            if gdf_clm is None:
                gdf_clm = gdf_clm_period
            else:
                gdf_clm = pd.concat([gdf_clm, gdf_clm_period])

        if gdf_clm is not None:
            try:
                climate.insert_climate_nldas2(gdf_clm, engine, table=db_table_name)
                if delete_chunks_on_sucess:
                    for chunk in chunk_files:
                        chunk.unlink()
            except Exception as e:
                print(f"Couldn't insert {station_id}: {e} \nContinuing")

    print("All done!")

def open_creds_dict(path_to_json_creds):
    with open(path_to_json_creds, "r") as f:
        credentials = json.load(f)
        return credentials

def str2date(date_str):
    """
    Convert a string representation of a date in the format "YYYY-MM-DD" to a `date` object.

    Args:
        date_str (str): A string representing a date in the format "YYYY-MM-DD".

    Returns:
        datetime.date: A `date` object representing the input date.

    Raises:
        ValueError: If the input string cannot be parsed as a date in the specified format.
    """
    return datetime.strptime(date_str, "%Y-%m-%d").date()

def delete_hyriver_cache():
    """
    Delete the HyRiver cache folder
    """
    try:
        path_to_cache = Path(os.environ["HYRIVER_CACHE_NAME"])
        print("Deleting HyRiver cache")
        path_to_cache.unlink(missing_ok=True)
    except KeyError:
        pass

def cli_call():

    parser = argparse.ArgumentParser(description="Query climate data from NLDAS2 and populate it to the AIMS database")

    # THUCS can be provided directly in the command line or if it's a path to a file with a .csv extension (without headers) it will read it
    parser.add_argument('--start_date',             help='Start date in YYYY-MM-DD format',                   type=str, default="2000-01-01")
    parser.add_argument('--end_date',               help='End date in YYYY-MM-DD format',                     type=str, default="2022-12-31")
    parser.add_argument('--path_nldas_daily_files', help='Path to the directory with daily aggregated files', type=str)
    parser.add_argument('--path_to_creds',          help='Path to the database credentials JSON file',        type=str)
    parser.add_argument('--maxiter_global',         help='Maximum number of attempts',                        type=int, default=10)
    parser.add_argument('--db_table_name',          help='Name of the DB table to update',                    type=str, default="climate_nldas2")
    parser.add_argument('--lats', '-lt', nargs='+',
                        help="List of latitudes to extract (to be paired with matching lons)")
    parser.add_argument('--lons', '-ln', nargs='+',
                        help="List of longitudes to extract (to be paired with matching lats)")

    args = parser.parse_args()

    lons = args.lons
    lats = args.lats

    if not(lons) or not(lats):
        coords = None
    else:
        if isinstance(lons[0], str) and lons[0].lower() == "all":
            lons = np.arange(-124.9375, -67.0625+0.125, 0.125)

        if isinstance(lats[0], str) and lats[0].lower() == "all":
            lats = np.arange(25.0625,    52.9375+0.125, 0.125)

        coords = [(float(lon), float(lat)) for lon, lat in itertools.product(lons, lats)]

    main(
        START_DATE=args.start_date,
        END_DATE=args.end_date,
        coords=coords,
        path_nldas_daily_files=args.path_nldas_daily_files,
        path_to_creds=args.path_to_creds,
        db_table_name=args.db_table_name,
        MAXITER_GLOBAL=args.maxiter_global,
    )

def cli_call_pop_parquet():
    parser = argparse.ArgumentParser(description="Looks for AnnAGNPS compatible parquet files and populates them to the database")

    # THUCS can be provided directly in the command line or if it's a path to a file with a .csv extension (without headers) it will read it
    parser.add_argument('--rods_dir',                help='Directory with daily aggregated files in parquet format',                 type=str)
    parser.add_argument('--path_to_creds',           help='Path to the database credentials JSON file',                              type=str)
    parser.add_argument('--db_table_name',           help='Name of the DB table to update',                                          type=str, default="climate_nldas2")
    parser.add_argument('--pattern',                 help='Pattern of the parquet files to be explored with the glob function.',     type=str, default="**/climate_daily*chunk*.parquet")
    parser.add_argument('--delete_chunks_on_sucess', help='If True, the chunks will be deleted after being successfully processed.', type=bool, default=False)
    args = parser.parse_args()

    populate_parquet_to_db(
        rods_dir=args.rods_dir,
        path_to_creds=args.path_to_creds,
        db_table_name=args.db_table_name,
        pattern=args.pattern,
        delete_chunks_on_sucess=args.delete_chunks_on_sucess
    )

if __name__ == "__main__":
    cli_call()