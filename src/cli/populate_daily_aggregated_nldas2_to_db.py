# import psycopg2
import argparse

import os, socket
from pathlib import Path
from datetime import datetime

import random
import itertools

import json
import time

import math

import numpy as np
import pandas as pd
import geopandas as gpd

from shapely.geometry import Point

from sqlalchemy import URL, create_engine, text as sql_text
# from sqlalchemy.orm import sessionmaker

# from geoalchemy2 import Geometry

from pyagnps import climate
from pyagnps.utils import log_to_file, get_current_time


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

if __name__ == "__main__":
    cli_call()