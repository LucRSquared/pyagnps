from pathlib import Path
import rasterio
from rasterio import features
from shapely.geometry import shape
from shapely.geometry.polygon import Polygon, LinearRing
import numpy as np

from tqdm import tqdm

import pandas as pd
import geopandas as gpd
import shutil, os, glob, sys
import time

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
    
def get_bounds_from_cells(gdf, buffer=0.00001):
    bounds = gdf[[gdf.geometry.name]].copy(deep=True)
    bounds.geometry = bounds.geometry.buffer(buffer)
    bounds = bounds.dissolve()
    bounds = bounds.to_crs('epsg:4326')
    return bounds

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


# def move_files_from_dir_to_dir(path_to_dir, path_to_new_dir):
#     """
#     path_to_dir : directory containing files to be moved
#     path_to_new_dir : directory where files will be moved
#     """

#     move_error_files = []

#     all_files = glob.glob(f"{path_to_dir}/*")

#     for elem in all_files:
#         try:
#             shutil.move(elem, path_to_new_dir)
#         except:
#             move_error_files.append(elem)

#     return move_error_files

def move_files_from_dir_to_dir(source_dir: Path, destination_dir: Path) -> list[str]:
    """Moves files from one directory to another, handling both string and Path objects.

    Args:
        source_dir (Path or str): The source directory containing files to be moved.
        destination_dir (Path or str): The destination directory where files will be moved.

    Returns:
        list[str]: A list of files that failed to move.
    """

    source_dir = Path(source_dir)  # Convert source_dir to Path object
    destination_dir = Path(destination_dir)  # Convert destination_dir to Path object

    move_error_files = []

    for file in source_dir.iterdir():
        if file.is_file():  # Ensure it's a file
            try:
                shutil.move(file, destination_dir / file.name)  # Use pathlib for path handling
            except Exception as e:
                move_error_files.append(file.name)
                print(f"Error moving {file.name}: {e}")  # Print error message

    return move_error_files


# def copy_files_from_dir_to_dir(path_to_dir, path_to_new_dir):
#     """
#     path_to_dir : directory containing files to be copied
#     path_to_new_dir : directory where files will be copied
#     """

#     copy_error_files = []

#     all_files = glob.glob(f"{path_to_dir}/*")

#     for elem in all_files:
#         try:
#             shutil.copy2(elem, path_to_new_dir)
#         except:
#             copy_error_files.append(elem)

#     return copy_error_files

def copy_files_from_dir_to_dir(source_dir: Path, destination_dir: Path) -> list[str]:
    """Copies files from one directory to another, handling both string and Path objects.

    Args:
        source_dir (Path or str): The source directory containing files to be copied.
        destination_dir (Path or str): The destination directory where files will be copied.

    Returns:
        list[str]: A list of files that failed to copy.
    """

    source_dir = Path(source_dir)  # Convert source_dir to Path object
    destination_dir = Path(destination_dir)  # Convert destination_dir to Path object

    copy_error_files = []

    for file in source_dir.iterdir():
        if file.is_file():  # Ensure it's a file
            try:
                shutil.copy2(file, destination_dir / file.name)  # Use pathlib for copying
            except Exception as e:
                copy_error_files.append(file.name)
                print(f"Error copying {file.name}: {e}")  # Print error message

    return copy_error_files

def get_relative_path(source_path: Path, target_path: Path) -> Path:
  """
  Calculates the relative path from the source path to the target path.

  Args:
      source_path (Path): The path from which the relative path is calculated.
      target_path (Path): The target path for which the relative path is desired.

  Returns:
      Path: The relative path from source_path to target_path.
  """

  source_parts = source_path.parts
  target_parts = target_path.parts
  num_common_parts = 0

  # Find the number of common path components
  while num_common_parts < len(source_parts) and num_common_parts < len(target_parts) and source_parts[num_common_parts] == target_parts[num_common_parts]:
    num_common_parts += 1

  # Build the relative path
  relative_parts = [Path("..")] * (len(source_parts) - num_common_parts)
  relative_parts.extend(target_parts[num_common_parts:])
  return Path(*relative_parts)

def relative_input_file_path(output_folder, path_to_file):
    path_to_file = Path(path_to_file)
    output_folder = Path(output_folder)
    relative_path = str(path_to_file.relative_to(output_folder))
    relative_path = relative_path.replace('\\', '/')
    relative_path = f"./{relative_path}"
    return relative_path

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

    if not (date.tzinfo is None or date.tzinfo == dateutil.tz.tzutc()):
        date = date.astimezone(timezone.utc)

    date = date.replace(tzinfo=None)

    if outputtype == datetime:
        return date
    elif outputtype == np.datetime64:
        return np.datetime64(date)
    else:
        raise TypeError("outputtype must be datetime or np.datetime64")

def month_difference(from_date, to_date):
    """
    gets the number of months between from_date and to_date
    could be negative or zero
    example usage:
        >>> month_difference(datetime(2010,10,1), datetime(2010,11,1))
        1
    """
    return (to_date.year - from_date.year) * 12 + (to_date.month - from_date.month)


def write_csv_control_file_from_dict(data_dict, output_path='control.csv'):
    # Writes the contents of kwargs with the key as a column and the value as the value
    # at output_path
    output_path = Path(output_path)
    # Create a DataFrame concisely using a list of tuples
    df = pd.DataFrame([list(data_dict.values())], columns=data_dict.keys())

    # Write the DataFrame to a CSV file, ensuring header row
    df.to_csv(output_path, index=False, header=True)

def find_rows_containing_pattern(file, pattern, skiprows=0):
    # Returns the list of row numbers (0-indexed)for the lines that contains the pattern
    # does not return row numbers smaller than skiprows

    row_nums = []

    with open(file, 'r') as f:
        lines = f.readlines()

        for i, line in enumerate(lines):
            if i < skiprows:
                continue
            if pattern in line:
                row_nums.append(i)

    return row_nums

def download_files_from_url(session, urls, out_dir):
    """
    visit and download the files from the
    for the list of urls for the files we wish to retrieve
    """
    responses = []

    # display progress bar (thanks tqdm!) and download the files
    for url in tqdm(urls):
        # extract the filename from the url to be used when saving the file
        filename = Path(url).name

        maxattempts = 10
        attempt = 1
        while attempt <= maxattempts:
            try:
                # save the file
                out_file = (out_dir / filename).resolve()
                if out_file.exists():
                    # print(f'File {filename} already exists, skipping')
                    break
                else:
                    # submit the request using the session
                    response = session.get(url, stream=True)
                    responses.append(response)

                    # raise an exception in case of http errors
                    response.raise_for_status()

                    with out_file.open('wb') as fd:
                        for chunk in response.iter_content(chunk_size=1024 * 1024):
                            fd.write(chunk)
               

            # except requests.exceptions.HTTPError as e:
            #     print(f'HTTP Error {e.response.status_code} for url {url}, Retrying {attempt}/{maxattempts}')
            #     attempt += 1
            #     time.sleep(5)
            except:
                print(f'\nHTTP Error for url {url}, Retrying {attempt}/{maxattempts}')
                attempt += 1
                time.sleep(5)

            if attempt == maxattempts:
                print(f'Failed to download {url} after {maxattempts} attempts, appending to file {out_dir}/failed_urls.txt')
                with (out_dir / 'failed_urls.txt').open('a') as f:
                    f.write(f'{url}\n')

    return responses
