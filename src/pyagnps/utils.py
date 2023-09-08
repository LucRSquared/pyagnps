import rasterio
from rasterio import features
from shapely.geometry import shape
from shapely.geometry.polygon import Polygon, LinearRing
import numpy as np
import geopandas as gpd
import shutil, os, glob, sys

import dateutil
from datetime import datetime, timezone


def polygonize_cell_reach_IDs_asc(
    path_to_cell_IDs_asc,
    path_to_cell_IDs_shp="cells.gpkg",
    outtype="GPKG",
    outepsg=4326,
    return_gdf=False,
    writefile=True,
):
    with rasterio.open(path_to_cell_IDs_asc) as src:
        band = src.read(1).astype("int32")
        nodata = src.nodata

    mask = band != nodata

    shapes = features.shapes(
        source=band, mask=mask, connectivity=8, transform=src.transform
    )

    data = {"DN": [], "geometry": []}
    for geom, value in shapes:
        data["DN"].append(int(value))
        data["geometry"].append(shape(geom))

    gdf = gpd.GeoDataFrame(data, crs=src.crs).to_crs(epsg=outepsg)

    if writefile:
        try:
            gdf.to_file(path_to_cell_IDs_shp, driver=outtype)
        except:
            gdf.to_file(path_to_cell_IDs_shp, driver="GPKG")

    if return_gdf:
        return gdf
    else:
        return


def get_geometry_boundary(geom):
    if isinstance(geom, LinearRing):
        return geom
    elif isinstance(geom, Polygon):
        return geom.exterior
    elif isinstance(geom, gpd.geoseries.GeoSeries):
        if geom.shape[0] == 1:
            return geom.values[0].exterior
        else:
            raise Exception("Geoseries has more than one object")
    else:
        raise Exception("Invalid geometry object")


def remove_all_files_from_dir_except_from_list(path_to_dir, keep_list=None):
    """
    path_to_dir : directory containing files to be deleted
    keep_list : list containing file names to keep.
                If keep_list is None it deletes everything
                If keep_list is 'all' then do nothing
    """

    delete_error_files = []

    if keep_list is None:
        shutil.rmtree(path_to_dir)
    elif keep_list == "all":
        return []
    else:
        all_files = os.listdir(path_to_dir)
        delete_list = [f for f in all_files if f not in keep_list]

        for elem in delete_list:
            try:
                os.remove(os.path.join(path_to_dir, elem))
            except:
                delete_error_files.append(elem)

    return delete_error_files


def move_files_from_dir_to_dir(path_to_dir, path_to_new_dir):
    """
    path_to_dir : directory containing files to be moved
    path_to_new_dir : directory where files will be moved
    """

    move_error_files = []

    all_files = glob.glob(f"{path_to_dir}/*")

    for elem in all_files:
        try:
            shutil.move(elem, path_to_new_dir)
        except:
            move_error_files.append(elem)

    return move_error_files


def copy_files_from_dir_to_dir(path_to_dir, path_to_new_dir):
    """
    path_to_dir : directory containing files to be copied
    path_to_new_dir : directory where files will be copied
    """

    copy_error_files = []

    all_files = glob.glob(f"{path_to_dir}/*")

    for elem in all_files:
        try:
            shutil.copy2(elem, path_to_new_dir)
        except:
            copy_error_files.append(elem)

    return copy_error_files


def log_to_file(filepath, message):
    original_stdout = sys.stdout
    with open(filepath, "a") as log:
        sys.stdout = log
        print(message)
        sys.stdout = original_stdout


def get_current_time(format="%Y-%m-%d-%H-%M-%S"):
    now = datetime.now()
    now_str = now.strftime(format)
    return now_str

def get_date_from_string(date_string, outputtype=np.datetime64):
    """Converts a string to a datetime object
       Converts to UTC time zone (naive datetime object)
       
       date_string: string in ISO 8601 (friendly-ish) format"""

    date = dateutil.parser.parse(date_string)

    if not(date.tzinfo is None or date.tzinfo == dateutil.tz.tzutc()):
        date = date.astimezone(timezone.utc)
        
    date = date.replace(tzinfo=None)

    if outputtype == datetime:
        return date
    elif outputtype == np.datetime64:
        return np.datetime64(date)
    else:
        raise TypeError("outputtype must be datetime or np.datetime64")
