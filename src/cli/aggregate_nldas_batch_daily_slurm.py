import warnings
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=RuntimeWarning, module="dask.array.reductions")

import argparse
from pathlib import Path
import math
import time
from datetime import datetime
from dateutil.relativedelta import relativedelta
import xarray as xr
import pandas as pd
import numpy as np
import geopandas as gpd
from dask.distributed import Client, progress
from dask_jobqueue import SLURMCluster
from pyagnps.utils import get_date_from_string, download_simple_file
from pyagnps import climate

def calculate_optimal_chunk_size(total_dataset_size, available_memory, safety_factor=0.2):
    max_chunk_size = available_memory * (1 - safety_factor)
    return min(max_chunk_size, total_dataset_size * 0.01)  # Cap at 1% of total dataset size

def setup_slurm_dask_cluster(nodes=5, cores_per_node=16, memory_per_node="64GB"):
    cluster = SLURMCluster(
        queue='your_queue_name',
        project='your_project_name',
        cores=cores_per_node,
        memory=memory_per_node,
    )
    cluster.scale(nodes * cores_per_node)
    client = Client(cluster)
    return client

def get_gmt_offset(lon, lat, gmt_offset_dict):
    return gmt_offset_dict.get((lon, lat), None)

def main(**kwargs):
    begin_time = time.time()

    # Extract parameters
    from_date = kwargs.get('from_date')
    to_date = kwargs.get('to_date')
    files_dir = Path(kwargs.get('files_dir'))
    product = kwargs.get('product', 'NLDAS_FORA0125_H.2.0')
    nldas_grid_gmt_file = kwargs.get('nldas_grid_gmt_file')
    output_dir_agg_netcdf = Path(kwargs.get('output_dir_agg_netcdf') or files_dir / f'nldas-{from_date}-{to_date}_DAILY')
    chunk_frac_tot = kwargs.get('chunk_frac_tot', 1.0)

    netcdf_fileroot = product.replace('H.2.0', 'D.')

    if not files_dir:
        raise ValueError('Missing source files directory')

    output_dir_agg_netcdf.mkdir(parents=True, exist_ok=True)

    if nldas_grid_gmt_file is None:
        url = 'https://github.com/LucRSquared/pyagnps/raw/main/inputs/climate/NLDAS_GRID_LAT_LON_GMT_OFFSET.gpkg'
        dirname = Path().cwd() / 'ancillary_files'
        nldas_grid_gmt_file = download_simple_file(url, dirname)
    else:
        nldas_grid_gmt_file = Path(nldas_grid_gmt_file)

    print(f'Reading grid GMT offset file: {nldas_grid_gmt_file}')
    grid_gmt_offset = gpd.read_file(nldas_grid_gmt_file)

    gmt_offset_min = grid_gmt_offset['GMT_OFFSET'].min()
    gmt_offset_max = grid_gmt_offset['GMT_OFFSET'].max()

    gmt_offset_dict = {(row['lon'], row['lat']): row['GMT_OFFSET'] for _, row in grid_gmt_offset.iterrows()}

    start_date = get_date_from_string(from_date, outputtype=datetime)
    end_date = get_date_from_string(to_date, outputtype=datetime)

    total_duration = end_date - start_date
    frac_dur = chunk_frac_tot * total_duration
    time_processing_chunk = relativedelta(days=math.ceil(frac_dur.total_seconds()/(24*3600)))

    # Calculate optimal chunk size
    total_dataset_size = 1.5 * 1024 * 1024 * 1024 * 1024  # 1.5TB in bytes (estimate)
    available_memory = 64 * 1024 * 1024 * 1024  # 64GB in bytes
    optimal_chunk_size = calculate_optimal_chunk_size(total_dataset_size, available_memory)
    chunk_frac_tot = optimal_chunk_size / total_dataset_size

    # Set up SLURM-aware Dask cluster
    client = setup_slurm_dask_cluster()
    print("Dashboard address:", client.dashboard_link)

    curr_date = start_date
    end_date_proc = end_date

    while curr_date < end_date_proc:
        if (end_date_proc - curr_date).days < 365*time_processing_chunk.years + 30*time_processing_chunk.months + time_processing_chunk.days:
            time_processing_chunk = relativedelta(seconds=int((end_date_proc - curr_date).total_seconds()), microseconds=(end_date_proc - curr_date).microseconds)

        print(f"Processing data from {curr_date} to {curr_date + time_processing_chunk}")

        with progress():
            ds = climate.read_dataset(files_dir, curr_date+relativedelta(hours=-float(gmt_offset_max)),
                                      time_processing_chunk+relativedelta(hours=-float(gmt_offset_min)),
                                      product=product, augment=True)

            lon, lat = xr.broadcast(ds.lon, ds.lat)
            ds['GMT_OFFSET'] = xr.apply_ufunc(get_gmt_offset, lon, lat, kwargs={'gmt_offset_dict': gmt_offset_dict}, vectorize=True)
            ds['GMT_OFFSET'].attrs['units'] = 'hours'

            grouped = ds.groupby('GMT_OFFSET')

            grouped_agg_list = []
            for gmt_offset, group in grouped:
                group['time_local'] = group['time'].values + pd.to_timedelta(gmt_offset, unit='h')

                if 'FORA' in product:
                    agg = climate.aggregate_forcing_data_daily(group)
                # Add NOAH aggregation here when implemented

                grouped_agg_list.append(agg.unstack('stacked_lon_lat'))

            combined_agg = xr.concat(grouped_agg_list, dim='GMT_OFFSET').compute()

            for t in combined_agg.time.values:
                combined_agg_t = combined_agg.sel(time=t)
                daystamp = pd.to_datetime(combined_agg_t.time.values)
                filename = f'{netcdf_fileroot}A{daystamp.year}{daystamp.month:02d}{daystamp.day:02d}.020.nc'
                
                print(f"Writing {filename}")
                combined_agg_t.to_netcdf(output_dir_agg_netcdf / filename)

        curr_date += time_processing_chunk

    elapsed_time = time.time() - begin_time
    print(f'Elapsed time: {elapsed_time} seconds')

    client.close()

def cli_call():
    parser = argparse.ArgumentParser()
    parser.add_argument('--from_date', '-fd', help="From Date/Time Format (lower boundary included) YYYY-MM-DDTHHMM in local time", type=str)
    parser.add_argument('--to_date', '-td', help="To Date/Time (upper boundary excluded) Format YYYY-MM-DDTHHMM in local time", type=str)
    parser.add_argument('--product', '-prd', help="The NLDAS product to download, default is NLDAS_FORA0125_H.2.0")
    parser.add_argument('--files_dir', '-f', help="Source files directory", type=str)
    parser.add_argument('--nldas_grid_gmt_file', '-grd', help="NLDAS Grid GMT File (geopackage with lat, lon, GMT offset). If not provided, the file will be downloaded.", type=str)
    parser.add_argument('--output_dir_agg_netcdf', '-onetcdf', help="Output Directory for output aggregated netcdf files", type=str, required=False)
    parser.add_argument('--chunk_frac_tot', '-ch', help="Chunk fraction of total data", type=float, required=False)

    args = parser.parse_args()
    main(**vars(args))

if __name__ == '__main__':
    cli_call()