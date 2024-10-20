import os, sys, shutil
import argparse
from pathlib import Path
import math

from datetime import datetime
from dateutil.relativedelta import relativedelta
import xarray as xr
import pandas as pd
import numpy as np

import geopandas as gpd

from dask.distributed import Client

from pyagnps.utils import get_date_from_string, download_simple_file
# import gesdisc_get.prepare_nldas_annagnps as pna

# This script is used to generate climate data from pre-aggregated netcdf files

def main():
    # Parse the command line arguments

    parser = argparse.ArgumentParser()

    parser.add_argument('--from_date', '-fd',
                        help="From Date/Time Format (lower boundary included) YYYY-MM-DDTHHMM in local time",
                        type=str)

    parser.add_argument('--to_date', '-td',
                        help="To Date/Time (upper boundary excluded) Format YYYY-MM-DDTHHMM in local time",
                        type=str)

    parser.add_argument('--lats', '-lt', nargs='+',
                        help="List of latitudes to extract (to be paired with matching lons)")

    parser.add_argument('--lons', '-ln', nargs='+',
                        help="List of longitudes to extract (to be paired with matching lats)")

    parser.add_argument('--product', '-prd',
                        help=("The NLDAS product to download, default is NLDAS_FORA0125_H.2.0"
                              "see also: 'https://hydro1.gesdisc.eosdis.nasa.gov/data/NLDAS/'"))

    parser.add_argument('--files_dir', '-f',
                        help=("Source files directory)"),
                        type=str)

    parser.add_argument('--nldas_grid_gmt_file', '-grd',
                        help=("NLDAS Grid GMT File (geopackage with lat, lon, GMT offset). If not provided, GMT offset will be assumed to be 0. If 'default', the default path for the file will be checked."),
                        type=str)

    parser.add_argument('--output_dir', '-o',
                        help=("Output Directory for climate files "
                            "(default ./nldas-<from_date>.to.<to_date>"),
                        type=str,
                        required=False)

    parser.add_argument('--output_dir_agg_netcdf', '-onetcdf',
                        help=("Output Directory for output aggregated netcdf files "
                            "(default ./nldas_netcdf-<from_date>.to.<to_date>"),
                        type=str,
                        required=False)

    parser.add_argument('--output_formats', '-of', nargs='+',
                        help=("List of Output Formats"
                            "(default parquet"),
                        required=False)

    parser.add_argument('--chunk_frac_tot', '-ch',
                        help=("Chunk fraction of total data"),
                        type=float,
                        required=False)


    args = parser.parse_args()

    if args.product is None:
        product = 'NLDAS_FORA0125_H.2.0'
    else:
        product = args.product

    netcdf_fileroot = product.replace('H.2.0', 'D.')
    product_daily = product.replace('H.2.0', 'D.2.0')

    if args.output_dir is None:
        default_output_dir = Path().cwd() / f'nldas-{start_date.year}{start_date.month:02d}{start_date.day:02d}-{end_date.year}{end_date.month:02d}{end_date.day:02d}'
        output_dir = default_output_dir
    else:
        output_dir = Path(args.output_dir)

    if args.output_dir_agg_netcdf is None and (output_dir == default_output_dir):
        output_dir_agg_netcdf = Path().cwd() / f'nldas-{start_date.year}{start_date.month:02d}{start_date.day:02d}-{end_date.year}{end_date.month:02d}{end_date.day:02d}_DAILY'
    elif args.output_dir_agg_netcdf is None:
        output_dir_agg_netcdf = output_dir / f'nldas-{start_date.year}{start_date.month:02d}{start_date.day:02d}-{end_date.year}{end_date.month:02d}{end_date.day:02d}_DAILY'
    else:
        output_dir_agg_netcdf = Path(args.output_dir_agg_netcdf)

    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)

    if not output_dir_agg_netcdf.exists():
        output_dir_agg_netcdf.mkdir(parents=True, exist_ok=True)

    if args.files_dir is None:
        raise ValueError('Missing source files directory')
    else:
        files_dir = Path(args.files_dir)

     
    if args.nldas_grid_gmt_file is None:
        # Download from : https://github.com/LucRSquared/pyagnps/raw/main/inputs/climate/NLDAS_GRID_LAT_LON_GMT_OFFSET.gpkg
        # and save in the same directory if the file doesn't exist there

        dirname = Path().cwd() / 'ancillary_files'

        nldas_grid_gmt_file = download_simple_file(args.nldas_grid_gmt_file, dirname)

        LOCAL_GMT_OFFSET = True
    else:
        nldas_grid_gmt_file = Path(args.nldas_grid_gmt_file)

    if args.output_formats is None:
        output_formats = ['parquet']
    else:
        output_formats = args.output_formats

    if args.chunk_frac_tot is None:
        chunk_frac_tot = 1.0
    else:
        chunk_frac_tot = args.chunk_frac_tot

    lats = args.lats
    lons = args.lons

    # Convert to lists of floats
    if lats is not None:
        lats = [float(lat) for lat in lats]
    if lons is not None:
        lons = [float(lon) for lon in lons]

    if (lats == 'all') and (lons == 'all'):
        lats = np.arange(-124.9375, -67.0625+0.125, 0.125)
        lons = np.arange(25.0625, 52.9375+0.125, 0.125)

    # Read grid GMT offset file if possible
    if LOCAL_GMT_OFFSET:
        try:
            print(f'Reading grid GMT offset file : {nldas_grid_gmt_file}')
            grid_gmt_offset = gpd.read_file(nldas_grid_gmt_file)
            LOCAL_GMT_OFFSET = True
        except:
            print('Error reading grid GMT offset file using GMT offset = 0')
            LOCAL_GMT_OFFSET = False
            GMT_OFFSET = 0
    else:
        GMT_OFFSET = 0

    start_date = get_date_from_string(args.from_date, outputtype=datetime)
    end_date = get_date_from_string(args.to_date, outputtype=datetime)

    if args.output_dir is None:
        output_dir = os.path.join(os.curdir,
                                f'nldas-{start_date.year}{start_date.month:02d}{start_date.day:02d}.to.{end_date.year}{end_date.month:02d}{end_date.day:02d}')
    else:
        output_dir = os.path.join(os.curdir, args.output_dir)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Make a temporary directory under output_dir to store temporary files
    tmp_dir = os.path.join(output_dir, 'tmp')
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # Processing chunk size (in days)
    total_duration = end_date - start_date
    frac_dur = chunk_frac_tot * total_duration
    time_processing_chunk = relativedelta(days=math.ceil(frac_dur.total_seconds()/(24*3600))) # rounded up to the nearest day


    ######################################################
    # Aggregate netCDF files based on the GMT offset and write to files if they don't exist

    if LOCAL_GMT_OFFSET:
        # Get the unique GMT offsets
        unique_gmt_offsets = grid_gmt_offset['GMT_OFFSET'].unique()
    else:
        unique_gmt_offsets = [GMT_OFFSET]

    latslons_input = list(zip(lats, lons))

    print(f'Inputs : {latslons_input}')

    # Loop through each GMT offset and aggregate the netCDF files
    for gmt_offset in unique_gmt_offsets:

        if LOCAL_GMT_OFFSET:
            # Get the grid cells with this GMT offset
            grid_gmt_offset_subset = grid_gmt_offset[grid_gmt_offset['GMT_OFFSET'] == gmt_offset]

            # Get the lat/lon values for the grid cells with this GMT offset
            lats_gmt = grid_gmt_offset_subset['lat'].values
            lons_gmt = grid_gmt_offset_subset['lon'].values

        else:
            lats_gmt = lats
            lons_gmt = lons

        latslons_gmt = list(zip(lats_gmt, lons_gmt))

        # Get lat/lon values for the grid cells with this GMT offset
        latslons = list(set(latslons_input).intersection(set(latslons_gmt)))

        if not(latslons):
            print(f'No grid cells with GMT offset = {gmt_offset}')
            continue
        
        with Client() as client:

            print(client)

            # get the files for each dataset type requested
            agg_list = []

            curr_date = start_date + relativedelta(hours=-float(gmt_offset)) # THE LOCAL GMT OFFSET IS ADDED HERE. I NEED TO DO MINUS TO EXPRESS IT IN UTC
            end_date_proc = end_date + relativedelta(hours=-float(gmt_offset)) # THE LOCAL GMT OFFSET IS ADDED HERE. I NEED TO DO MINUS TO EXPRESS IT IN UTC

            netcdf_file_list = []

            while curr_date < end_date_proc:

                if (end_date_proc - curr_date).days < 365*time_processing_chunk.years + 30*time_processing_chunk.months + time_processing_chunk.days:
                    time_processing_chunk = end_date_proc - curr_date
                    time_processing_chunk = relativedelta(seconds=int(time_processing_chunk.total_seconds()), 
                                                          microseconds=time_processing_chunk.microseconds)

                # next_date = curr_date + time_processing_chunk
                print(f'{datetime.now()} - Aggregating {curr_date.year}-{curr_date.month}-{curr_date.day} (UTC TIME, GMT OFFSET {gmt_offset})')
                
                # Read dataset
                ds = pna.read_dataset(files_dir, curr_date, time_processing_chunk, product=product, augment=True)

                if 'FORA' in product:

                    # Aggregate the forcing data
                    agg = pna.aggregate_forcing_data_daily(ds)
                
                elif 'NOAH' in product:

                    # Aggregate NOAH data
                    agg = pna.aggregate_noah_data_daily(ds)

                agg = agg.compute()

                # Write the aggregated forcing data
                for t in agg.time.values:

                    agg_t = agg.sel(time=t)

                    try:
                        daystamp = pd.to_datetime(agg_t.time.values[0])
                    except:
                        daystamp = pd.to_datetime(agg_t.time.values)

                    filename = f'{netcdf_fileroot}A{daystamp.year}{daystamp.month:02d}{daystamp.day:02d}.020.nc'
                    agg_t.to_netcdf(os.path.join(output_dir_agg_netcdf, filename), compute=True)

                    netcdf_file_list.append(os.path.join(output_dir_agg_netcdf, filename))

                agg_list.append(agg)
                
                curr_date += time_processing_chunk


            #####################################

            # Create dictionary of list of parquets files to later assemble for each coordinates
            parquet_files = {}
            oneday = relativedelta(days=1)

            curr_date = start_date
            chunk_counter = 1
            end_date_proc = end_date - oneday # -oneday to avoid including the next day

            # (Adjusted) Processing chunk size (in days)
            total_duration = end_date_proc - start_date
            frac_dur = chunk_frac_tot * total_duration
            time_processing_chunk = relativedelta(days=math.ceil(frac_dur.total_seconds()/(24*3600))) # rounded up to the nearest day

            while curr_date < end_date_proc:

                if (end_date_proc - curr_date).days < 365*time_processing_chunk.years + 30*time_processing_chunk.months + time_processing_chunk.days:
                    time_processing_chunk = end_date_proc - curr_date
                    time_processing_chunk = relativedelta(seconds=int(time_processing_chunk.total_seconds()), 
                                                            microseconds=time_processing_chunk.microseconds)

                # next_date = curr_date + time_processing_chunk
                print(f'{datetime.now()} - Processing {curr_date.year}-{curr_date.month}-{curr_date.day} - Chunk {chunk_counter}')
                
                file_list = pna.get_file_list(output_dir_agg_netcdf, curr_date, curr_date+time_processing_chunk, product=product_daily)

                print(f'Reading {len(file_list)} files...')

                agg_ds_chunk = xr.open_mfdataset(file_list, chunks={'time': 24}, 
                            parallel=False, combine='nested', concat_dim='time', engine='h5netcdf')
                
                agg_ds_chunk = agg_ds_chunk.persist()

                for lat, lon in latslons:
                    print(f'{datetime.now()} - Extracting data for lat: {lat}, lon: {lon}')
                    df = agg_ds_chunk.sel(lat=lat, lon=lon, method='nearest').to_dataframe()

                    # Save the data 
                    output_format = 'parquet'
                        
                    print(f'{datetime.now()} - Saving data for lat: {lat}, lon: {lon} in {output_format} format')
                    outfile = os.path.join(tmp_dir, f'chunk_{chunk_counter}_nldas_{product}_{lat}_{lon}_from_{start_date.year}-{start_date.month:02d}-{start_date.day:02d}_to_{end_date_proc.year}-{end_date_proc.month:02d}-{end_date.day:02d}.{output_format}')

                    parquet_files[(lat, lon)] = parquet_files.get((lat, lon), []) + [outfile]

                    df_out = pna.generate_climate_file_daily(df.copy(), output_filepath=outfile, saveformat=output_format)

                curr_date += time_processing_chunk
                chunk_counter += 1

            # print(parquet_files)

            # Assemble the parquet files
            for lat, lon in latslons:
                outfile = os.path.join(output_dir, f'nldas_{product}_{lat}_{lon}_from_{start_date.year}-{start_date.month:02d}-{start_date.day:02d}_to_{end_date_proc.year}-{end_date_proc.month:02d}-{end_date_proc.day:02d}.parquet')
                
                files = parquet_files[(lat, lon)]

                print(f'{datetime.now()} - Assembling parquet files for lat: {lat}, lon: {lon}')
                full_df = pd.concat([pd.read_parquet(f) for f in files])

                for output_format in output_formats:
                    if output_format == 'parquet':
                        print(f'{datetime.now()} - Saving data for lat: {lat}, lon: {lon} in {output_format} format')
                        full_df.to_parquet(outfile)
                    elif output_format == 'csv':
                        print(f'{datetime.now()} - Saving data for lat: {lat}, lon: {lon} in {output_format} format')
                        full_df.to_csv(outfile.replace('.parquet', '.'+output_format), index=False)
                    else:
                        raise ValueError(f'Unknown output format: {output_format}')

            client.close()

        # Remove temporary directory
    shutil.rmtree(tmp_dir)

    print(f'{datetime.now()} - Done')


if __name__ == '__main__':
    main()