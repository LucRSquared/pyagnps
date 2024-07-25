import os
from pathlib import Path
import argparse
import xarray as xr
import pandas as pd
import numpy as np
from datetime import datetime
from dateutil.relativedelta import relativedelta
from dask.distributed import Client
from timezonefinder import TimezoneFinder
import pytz
import dask.dataframe as dd

from pyagnps.utils import get_date_from_string

# import pna  # Assuming this contains your aggregation function


def read_netcdf_files(files_dir, start_date, end_date, product):
    file_pattern = f"{files_dir}/{product}*.nc"
    return xr.open_mfdataset(file_pattern, combine='by_coords', parallel=True, chunks={'time': 24})

def create_timezone_mapping(lats, lons):
    tf = TimezoneFinder()
    timezone_mapping = {}
    for lat, lon in zip(lats, lons):
        timezone_str = tf.certain_timezone_at(lat=lat, lon=lon)
        if timezone_str is None:
            timezone_str = 'UTC'
        timezone_mapping[(lat, lon)] = pytz.timezone(timezone_str)
    return timezone_mapping

def group_by_local_date(ds, timezone_mapping):
    def to_local_date(time, lat, lon):
        utc_time = pd.Timestamp(time).to_pydatetime().replace(tzinfo=pytz.UTC)
        local_tz = timezone_mapping[(lat, lon)]
        local_time = utc_time.astimezone(local_tz)
        return local_time.date()

    local_dates = xr.apply_ufunc(
        to_local_date,
        ds.time,
        ds.lat,
        ds.lon,
        vectorize=True,
        dask='allowed'
    )
    
    ds['local_date'] = xr.DataArray(local_dates, dims=('time', 'lat', 'lon'))
    return ds

def process_and_save_data(ds, output_dir, latslons, output_formats):
    for lat, lon in latslons:
        df = ds.sel(lat=lat, lon=lon, method='nearest').to_dataframe()
        
        for output_format in output_formats:
            outfile = output_dir / f'nldas_{lat}_{lon}.{output_format}'
            if output_format == 'parquet':
                df.to_parquet(outfile)
            elif output_format == 'csv':
                df.to_csv(outfile, index=False)

def main():
    parser = argparse.ArgumentParser()
    # Add all your argument parsers here
    args = parser.parse_args()

    start_date = get_date_from_string(args.from_date, outputtype=datetime)
    end_date = get_date_from_string(args.to_date, outputtype=datetime)

    # Set up output directories
    output_dir = args.output_dir or (Path.cwd() / f'nldas-{start_date.strftime("%Y%m%d")}.to.{end_date.strftime("%Y%m%d")}')
    output_dir_agg_netcdf = args.output_dir_agg_netcdf or (output_dir / 'aggregated_netcdf')

    output_dir.mkdir(exist_ok=True)
    output_dir_agg_netcdf.mkdir(exist_ok=True)

    lats = [float(lat) for lat in args.lats]
    lons = [float(lon) for lon in args.lons]
    latslons = list(zip(lats, lons))

    with Client() as client:
        print(client)

        # Read NetCDF files
        ds = read_netcdf_files(args.files_dir, start_date, end_date, args.product)
        
        # Create timezone mapping
        timezone_mapping = create_timezone_mapping(lats, lons)
        
        # Group by local date
        ds = group_by_local_date(ds, timezone_mapping)
        
        # Aggregate data
        agg_ds = pna.aggregate_forcing_data_daily(ds)
        
        # Save aggregated data
        for t in agg_ds.time:
            agg_t = agg_ds.sel(time=t)
            filename = f'{args.product}_A{t.dt.strftime("%Y%m%d").values}.nc'
            agg_t.to_netcdf(output_dir_agg_netcdf / filename)

        # Process and save data for each lat/lon pair
        process_and_save_data(agg_ds, output_dir, latslons, args.output_formats)

    print(f'{datetime.now()} - Done')

if __name__ == '__main__':
    main()