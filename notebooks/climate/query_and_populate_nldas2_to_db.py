import psycopg2


from pathlib import Path
from datetime import datetime

import json
import time

import math

import numpy as np
import pandas as pd
import geopandas as gpd

from shapely.geometry import Point, box
import folium


from sqlalchemy import URL, create_engine, inspect, text as sql_text
from sqlalchemy.orm import sessionmaker

from geoalchemy2 import Geometry

from pyagnps import climate


def main():

    # DATABASE SETUP
    path_to_creds_aims = Path("../../inputs/db_credentials.json")
    # path_to_creds_menderes = Path("../../inputs/db_credentials_menderes.json")

    creds = {
        'aims': open_creds_dict(path_to_creds_aims),
        # 'menderes': open_creds_dict(path_to_creds_menderes),
        'docker': {
                'user': 'postgres',
                'password': 'postgres_pass',
                'host': 'localhost',
                'port': '5432',
                'database': 'test_db'
            }
    }

    url_object = lambda db : URL.create(
                        "postgresql",
                        username=creds[db]['user'],
                        password=creds[db]['password'],
                        host=creds[db]['host'],
                        port=creds[db]['port'],
                        database=creds[db]['database'])

    db_url = url_object('aims')

    engine = create_engine(db_url)
    conn = engine.connect().execution_options(stream_results=True)


    # PARAMETERS
    thucs_to_process = set(['1148'])
    START_DATE = "2000-01-01"
    END_DATE = "2021-01-31"

    print(f"Processing stations in THUCS {thucs_to_process}")
    print(f"Processing stations between {START_DATE} and {END_DATE}")

    print(f"Reading THUC GPKG")
    # THUCS
    path_thucs = Path('../../inputs/thucs/tophuc_S_M_40000_closed_holes_with_container_thuc_merged_bbox_area_first_kept.gpkg')
    thucs = gpd.read_file(path_thucs)

    print(f"Reading NLDAS2 GPKG")
    # NLDAS2
    path_grid = Path('../../inputs/climate/NLDAS2_GRID_CENTROIDS_epsg4326.gpkg')
    nldas2_grid = gpd.read_file(path_grid)

    MAXITER_GLOBAL = 10
    MAXITER_SINGLE_STATION = 10

    for thuc_id in thucs_to_process:
        print(f"Performing Join Operation to get list of stations in THUCS {thuc_id}")
        # Perform the spatial join to get list of stations in the buffered THUC
        my_thuc = thucs.loc[thucs['tophucid']==thuc_id,:]
        buffered_geom = my_thuc.geometry.iloc[0].buffer(math.sqrt(2)/2*0.125) # We buffer by sqrt(2)/2 * (nldas_2 spacing) in degrees to get all the likely needed gird points

        buffered_thuc = gpd.GeoDataFrame({'geometry': buffered_geom}, index=[0], crs=my_thuc.crs)

        contained_stations = gpd.sjoin(nldas2_grid, buffered_thuc, how='inner', predicate='within')

        stations_to_process = contained_stations['nldas2_grid_ID'].map(str).to_list()

        # MAIN LOOP
        for iter_global in range(MAXITER_GLOBAL):

            incomplete_stations = {station_id: [] for station_id in stations_to_process}

            skipped_stations = 0
            for row in nldas2_grid.iterfeatures():
                station_id = str(row['properties']['nldas2_grid_ID'])
                if station_id not in stations_to_process:
                    skipped_stations += 1
                    continue

                x, y = row['geometry']['coordinates']

                available_dates = get_available_dates_for_station(station_id, engine, table="climate_nldas2")
                missing_dates = get_missing_dates(available_dates, START_DATE, END_DATE)

                if missing_dates:
                    incomplete_stations[station_id].append((x, y, missing_dates))
                else:
                    del incomplete_stations[station_id]

            if not incomplete_stations:
                break

            num_incomplete_stations = len(incomplete_stations)
            if num_incomplete_stations > 0:
                print(f"Number of NLDAS-2 stations with incomplete data for THUC {thuc_id}: {num_incomplete_stations}")
            else:
                print(f"All NLDAS-2 stations have complete data for THUC {thuc_id}")


            for station_id, data in incomplete_stations.items():
                for x, y, missing_dates in data:
                    for iter_station in range(MAXITER_SINGLE_STATION):
            
                        try:
                            # Query and concatenate climate data for each continuous period
                            gdf_clm = None
                            continuous_periods = find_continuous_periods(missing_dates)
                            
                            for period in continuous_periods:
                                start, end = period[0], period[-1]
                                clm_annagnps = climate.ClimateAnnAGNPSCoords(coords=(x, y), start=start, end=end, date_mode="local")
                                clm = clm_annagnps.query_nldas2_generate_annagnps_climate_daily(float_format='%.2f')
                                gdf_clm_period = prepare_annagnps_climate_for_db(clm, station_id, x, y)
                                
                                if gdf_clm is None:
                                    gdf_clm = gdf_clm_period
                                else:
                                    gdf_clm = pd.concat([gdf_clm, gdf_clm_period])

                            insert_climate_nldas2(gdf_clm, engine)

                            num_incomplete_stations -= 1

                            print(f"THUC {thuc_id}, [{iter_global+1}/{MAXITER_GLOBAL} global attempts], {station_id}, x = {x}, y = {y}: Populated after {iter_station+1} attempt(s), {num_incomplete_stations} station(s) remaining")

                            break
                        except Exception as e:
                            print(f"THUC {thuc_id}, [{iter_global+1}/{MAXITER_GLOBAL} global attempts], {station_id}, x = {x}, y = {y}: Failed with error: {e}: RETRYING ({iter_station+1}/{MAXITER_SINGLE_STATION})")
                            time.sleep(1)

    print(f"All done! for THUCS {thucs_to_process}")

def find_continuous_periods(missing_dates):
# Find continuous periods of dates
    continuous_periods = []
    current_period = []
    for date in missing_dates:
        if not current_period or (date - current_period[-1]).astype('timedelta64[D]').astype(int) == 1:
            current_period.append(date)
        else:
            continuous_periods.append(current_period)
            current_period = [date]
    if current_period:
        continuous_periods.append(current_period)

    return continuous_periods

def open_creds_dict(path_to_json_creds):
    with open(path_to_json_creds, "r") as f:
        credentials = json.load(f)
        return credentials

def prepare_annagnps_climate_for_db(clm, station_id, xgrid, ygrid):
    """
    Prepare climate data for insertion into the climate_nldas2 table
    * Inputs:
    - clm: pandas.DataFrame in AnnAGNPS format
    - station_id: str
    - xgrid: float longitude in EPSG:4326
    - ygrid: float latitude in EPSG:4326
    * Output:
    - gdf_clm: GeoDataFrame in EPSG:4326
    """
    clm.columns = clm.columns.str.lower()
    clm['station_id'] = station_id

    gdf_clm = gpd.GeoDataFrame(clm, geometry=[Point(xgrid, ygrid)] * len(clm), crs="EPSG:4326")
    gdf_clm.rename(columns={'geometry': 'geom'}, inplace=True)
    gdf_clm.index.name = "date"

    gdf_clm = gdf_clm.set_geometry('geom')

    return gdf_clm

def insert_climate_nldas2(gdf_clm, engine, table="climate_nldas2"):
    gdf_clm.to_postgis(table, engine, if_exists="append", index=True)

def climate_table_has_station(station_name, engine, table="climate_nldas2"):
    """This function checks if the table contains data for the given station. Also return the maximum and minimum date in the table
    Inputs:
    - station_name: str
    - engine: sqlalchemy.engine
    - table: str
    Outputs:
    - has_data: bool
    - min_date: datetime.date
    - max_date: datetime.date
    """
    query = f"""
        SELECT MIN(date) AS min_date, MAX(date) AS max_date
        FROM {table}
        WHERE station_id = '{station_name}'
        GROUP BY station_id
    """
    with engine.connect() as connection:
        result = connection.execute(sql_text(query))
        
        # Check if the query returned any rows
        if result.rowcount > 0:
            row = result.fetchone()
            min_date = row[0]
            max_date = row[1]
            return True, min_date, max_date
        else:
            return False, None, None
        
def get_available_dates_for_station(station_name, engine, table="climate_nldas2"):
    """This function returns a DataFrame with the available dates for the given station stored in the table
    Inputs:
    - station_name: str
    - engine: sqlalchemy.engine
    - table: str
    Outputs:
    - df: pandas.DataFrame
    """
    query = f"""
        SELECT DISTINCT date
        FROM {table}
        WHERE station_id = '{station_name}'
        ORDER BY date
    """
    with engine.connect() as connection:
        df = pd.read_sql(sql_text(query), connection)
        return df
    
def get_missing_dates(df, start_date, end_date):
    """ Creates a start_date to end_date data range and returns the dates that are missing from the dataframe"""
    
    # Ensure the 'date' column is in datetime format
    df['date'] = pd.to_datetime(df['date'])
    
    # Create a date range from start_date to end_date
    date_range = pd.date_range(start=start_date, end=end_date)
    
    # Get the set of dates in the dataframe
    df_dates = set(df['date'])
    
    # Find dates in the date_range that are not in df_dates
    missing_dates = [np.datetime64(date) for date in date_range if date not in df_dates]
    
    return missing_dates

def filter_climate_data(gdf_clm, missing_dates):
    """Returns a subset of gdf_clm that is outside the min_date and max_date interval

    gdf_clm: GeoDataFrame
    missing_dates: list of np.datetime64
    """    
    # Filter out rows where the date is within the min_date and max_date interval
    filtered_gdf = gdf_clm[gdf_clm.index.isin(missing_dates)].copy()

    return filtered_gdf

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

if __name__ == "__main__":
    main()