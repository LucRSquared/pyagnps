# COLLECTION OF AIMS SPECIFIC SCRIPTS AND METHODS
# This assumes that a connection to a database with properly configured
# tables for thucs is available

from pathlib import Path
import json

import pandas as pd
import geopandas as gpd

from sqlalchemy import create_engine, text as sql_text
from sqlalchemy import URL

class AIMSWatershed:
    def __init__(self, path_to_json_creds, **kwargs):
        credentials = open_creds_dict(path_to_json_creds)
        url_object = create_db_url_object(credentials)
        self.engine = create_engine(url_object)
        self.path_to_json_creds = Path(path_to_json_creds)

        self.thuc_id = kwargs.get("thuc_id", None)

        self.outlet_x = kwargs.get("outlet_x", None)
        self.outlet_y = kwargs.get("outlet_y", None)

        self.start_time = kwargs.get("start_time", None)
        self.end_time = kwargs.get("end_time", None)

        self.climate_method = kwargs.get("climate_method", "nldas2_database")

    def get_thuc_id_by_xy(self, x=None, y=None):
        """ Returns the thuc_id for which the coordinates belong to. """
        if x is None or y is None:
            x, y = self.outlet_x, self.outlet_y

        sql = sql_text(
            f"SELECT thuc_near_run_id_tr({x},{y})"
        )
        thuc = pd.read_sql(sql, self.engine)
        thuc_id = thuc.iloc[0].values[0]

        self.thuc_id = thuc_id
        return thuc_id
    
    def get_outlet_xy_from_thuc_id_and_reach_id(self, thuc_id, reach_id=2):
        """ If not specified the reach_id defaults to 2, the most downstream reach.
        This function returns the x and y coordinates of centroid of the reach_id in thuc_id."""
        sql = sql_text(
            f"""SELECT ST_X(ST_Centroid(ST_Collect(ST_Transform(geom, 4326)))) AS lon, 
                ST_Y(ST_Centroid(ST_Collect(ST_Transform(geom, 4326)))) AS lat 
            FROM thuc_{thuc_id}_annagnps_reach_ids
            WHERE dn = {reach_id}
            """
        )
        df = pd.read_sql(sql, self.engine)

        x, y = df["lon"].to_list()[0], df["lat"].to_list()[0]

        self.outlet_x = x
        self.outlet_y = y

        return x, y
    
    def set_reaches_for_output(self, output_reaches=['OUTLET']):
        """ Set the output_reaches for the simulation. """
        if output_reaches:
            self.output_reaches = [reach for reach in output_reaches]


def open_creds_dict(path_to_json_creds):
    with open(path_to_json_creds, "r") as f:
        credentials = json.load(f)
        return credentials
    
def create_db_url_object(credentials):
    url_object = URL.create(
        "postgresql",
        username=credentials["user"],
        password=credentials["password"],
        host=credentials["host"],
        port=credentials["port"],
        database=credentials["database"],
    )
    return url_object