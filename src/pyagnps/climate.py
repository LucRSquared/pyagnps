import warnings

import itertools

import psutil

from pathlib import Path
import pynldas2
import py3dep
import pydaymet

# import datetime
# import pytz
from timezonefinder import TimezoneFinder

import numpy as np
import pandas as pd
import geopandas as gpd

from sqlalchemy import URL, create_engine, text as sql_text

from shapely.geometry import Point

import xarray as xr

from tqdm import tqdm
import time
from datetime import datetime, timedelta, timezone
from dateutil.relativedelta import relativedelta

from .utils import get_date_from_string, month_difference

from .constants import _NLDAS_PRODUCTS_20, _NLDAS_PRODUCTS_002
# Do not create
import os, shutil

os.environ[
    "HYRIVER_CACHE_NAME"
] = "/tmp/climate_cache.sqlite"  # On Windows systems it will be under C:/tmp
# os.environ["HYRIVER_CACHE_DISABLE"] = "true"


class ClimateAnnAGNPSCoords:
    def __init__(
        self,
        coords: tuple = None,
        start: str = "2022-01-01",
        end: str = "2022-12-31",
        date_mode: str = "local",
    ):
        """
        ## Positional arguments:
        - coords: (longitude, latitude) in EPSG:4326 projection or
                  list of (longitude, latitude) [(lon1, lat1), (lon2, lat2), ...] in EPSG:4326 projection
        - start: Start date YYYY-MM-DD (assumes it starts at Midnight)
                 Can specify hour too: YYYY-MM-DDTHH
                 The data point returned will cover the time period YYYY-MM-DDTHH to YYYY-MM-DDT HH+1
        - end: End date, same specification as start

            - [Start date ------- Included Range ------- End date]

        ## Key-value arguments:
        - date_mode: "local" will assume that the start and end dates are provided
                             in the coords local time zone
                     "utc" will assume that the start and end dates are provided in
                           UTC timezone
        """

        self.coords = coords

        self.update_coords_start_end_dates(coords=coords,
                                           start=start,
                                           end=end,
                                           date_mode=date_mode)

        self.warnings = []
        self.clm = None  # DataFrame
        self.clm_resampled = None  # DataFrame resampled
        self.ds = None  # xarray.DataSet for CMIP6 outputs
        self.ds_select = (
            None  # xarray.DataSet for ds selection for coord and time bounds
        )

    def update_coords_start_end_dates(self,
                               coords = None,
                               start: str = None,
                               end: str = None,
                               date_mode: str = "local"):
        """
        ## Positional arguments:
        - coords: (longitude, latitude) in EPSG:4326 projection. Optional, if blank previous coords are kept
        - start: Start date YYYY-MM-DD (assumes it starts at Midnight)
                 Can specify hour too: YYYY-MM-DDTHH
                 The data point returned will cover the time period YYYY-MM-DDTHH to YYYY-MM-DDT HH+1
        - end: End date, same specification as start

            - [Start date ------- Included Range ------- End date]

        ## Key-value arguments:
        - date_mode: "local" will assume that the start and end dates are provided
                             in the coords local time zone
                     "utc" will assume that the start and end dates are provided in
                           UTC timezone
                     "daily" will assume that the times given are daily in local 24hours days
        """

        if not coords:
            coords = self.coords
            # If no coords were provided to begin with we just keep None for coords
            # and assume the dates are UTC
            if coords is None:
                LON_ALL = np.arange(-124.9375, -67.0625+0.125, 0.125)
                LAT_ALL = np.arange(25.0625,    52.9375+0.125, 0.125)
                coords = list(itertools.product(LON_ALL, LAT_ALL))
                date_mode = "utc"

        if not start:
            start = self.start
            date_mode = "utc"

        if not end:
            end = self.end
            date_mode = "utc"

        if date_mode.lower() == "utc":
            self.start = pd.Timestamp(start)
            self.end = pd.Timestamp(end)
            self.timezone = "UTC"

            self.start_input = self.start.tz_localize(None)
            self.end_input = self.end.tz_localize(None)
        
        elif date_mode.lower() == "daily":
            self.start = pd.Timestamp(start)
            self.end = pd.Timestamp(end)
            self.timezone = "local-24h-day"

            self.start_input = self.start.tz_localize(None)
            self.end_input = self.end.tz_localize(None)

        elif date_mode.lower() == "local":
            # Identify timezone:
            if isinstance(self.coords, list):
                if len(self.coords) > 1:
                    raise Exception(
                        "Date mode 'local' is only available for single point query"
                    )
                elif len(self.coords) == 1:
                    coords = self.coords[0]

            lon, lat = coords
            tf = TimezoneFinder()
            tmz_name = tf.timezone_at(lng=lon, lat=lat)
            self.timezone = tmz_name

            start_stamp_local, end_stamp_local = pd.Timestamp(
                start, tz=tmz_name
            ), pd.Timestamp(end, tz=tmz_name)

            self.start_input = start_stamp_local.tz_localize(None)
            self.end_input = end_stamp_local.tz_localize(None)

            start_stamp_utc = start_stamp_local.tz_convert("UTC").tz_localize(None)
            end_stamp_utc = end_stamp_local.tz_convert("UTC").tz_localize(None)

            self.start = start_stamp_utc
            self.end = end_stamp_utc
        else:
            raise Exception(f"Invalid ``date_mode`` {date_mode}")

        self.date_mode = date_mode

        self.coords = coords
        self.coords_actual = None  # Actual coordinates queried using nearest method
        

    def _query_nldas2_climate(
        self,
        source="netcdf",
        variables=[
            "prcp",
            "temp",
            "rsds",
            "pet",
            "humidity",
            "wind_u",
            "wind_v",
            "psurf",
        ],
    ):
        clm = pynldas2.get_bycoords(
            coords=self.coords,
            start_date=self.start,
            end_date=self.end,
            variables=variables,
            source=source,
            n_conn=4,
        )

        # Express dates in local mode
        if self.date_mode == "local":
            clm.index = clm.index.tz_convert(self.timezone).tz_localize(None)
        else:
            clm.index = clm.index.tz_localize(None)

        self.clm = clm

    def _resample(self, rule="1D"):
        clm = self.clm
        full_how_dict = {
            "prcp": "sum",
            "pet": "sum",
            "rsds": "mean",
            "wind_speed": "mean",
            "wind_direction": "mean",
            "tdew": "mean",
            "temp_min": "min",
            "temp_max": "max",
            "RH": "mean",
            "temp": "mean",
            "psurf": "mean",
            "vp (Pa)": "mean",
            "esat": "mean",
            "humidity": "mean",
            "wind_u": "mean",
            "wind_v": "mean",
        }

        # Prepare min/max temp columns
        clm["temp_min"] = clm["temp"]
        clm["temp_max"] = clm["temp"]

        # Create aggregation dict
        vars_how = {var: full_how_dict[var] for var in clm.columns}

        clm = clm.resample(rule=rule, origin="start").agg(vars_how)
        clm.index = clm.index.tz_localize(None)

        self.clm_resampled = clm

        return clm

    def _compute_additional_climate_variables_method_nldas2(self):
        self._compute_RH()
        self._compute_dew_point()
        self._compute_wind_direction()
        self._compute_wind_speed()

    def _resample_and_compute_additional_climate_variables_method_nldas2_and_daymet(
        self,
    ):
        # query vapor pressure from daymet
        daymet_vp_daily = pydaymet.get_bycoords(
            self.coords, dates=(self.start, self.end), variables="vp"
        )
        daymet_vp_daily = daymet_vp_daily.rename_axis("date")

        # self.clm['temp_min'] = self.clm['temp']
        # self.clm['temp_max'] = self.clm['temp']
        self._compute_wind_speed()
        self._compute_wind_direction()
        self._resample(rule="1D")

        clm_daily = self.clm_resampled
        clm_daily["esat"] = compute_esat(clm_daily["temp"], Tunit="K")
        # Merge with daymet data
        clm_daily = pd.merge(
            clm_daily, daymet_vp_daily, how="outer", right_on="date", left_index=True
        )

        # Compute RH
        clm_daily["RH"] = 100 * clm_daily["vp (Pa)"] / clm_daily["esat"]
        clm_daily["RH"] = clm_daily["RH"].where(
            (clm_daily["RH"] <= 100) | np.isnan(clm_daily["RH"]), 100
        )

        # Linearly interpolate possibly missing values
        clm_daily.interpolate(inplace=True)
        if clm_daily.index.name != "date":
            clm_daily.set_index("date", inplace=True)

        # Compute Dew Point
        clm_daily["tdew"] = compute_dew_point(clm_daily["RH"], clm_daily["temp"])

        self.clm_resampled = clm_daily

    def query_nldas2_generate_annagnps_climate_daily(self, **kwargs):
        """
        Generate climate_daily.csv AnnAGNPS file/DataFrame from NLDAS-2 (and DAYMET)

        ### Key-Value Arguments:
        - output_filepath : str, optional
               Path to write the output file, by default None
        - saveformat : str, optional
            Format to save the output file, by default 'csv', also accepts 'parquet'
        - float_format : str, optional, default= '%.3f' for printing csv file
        """

        if not self.coords:
            raise Exception("Coordinates are missing. Please provide coords!")

        self._query_nldas2_climate(source="netcdf")
        # Test if there are any NaN values
        if self.clm.isna().any().any():
            self.warnings.append('NaN values were found in the NLDAS-2 data, values supplemented by DAYMET data') 
            # Use alternative method using Day Met
            self._query_nldas2_climate(
                source="grib",
                variables=["prcp", "temp", "rsds", "pet", "wind_u", "wind_v"],
            )
            self._resample_and_compute_additional_climate_variables_method_nldas2_and_daymet()
        else:
            # do normal way
            self._compute_additional_climate_variables_method_nldas2()
            self._resample(rule="1D")

        self._keep_annagnps_columns_only()

        df_daily = self._generate_climate_file_daily(use_resampled=True, **kwargs)
        return df_daily

    def read_cmip6_generate_annagnps_climate_daily(self, cmip_files, **kwargs):
        """
        Generate climate_daily.csv AnnAGNPS file/DataFrame from CMIP6 output files

        ### Arguments:
        - cmip_files: Iterable (list), paths to CMIP6 files. The required variables are:
                      'pr'    : Precipitation (kg/(m².s))
                      'hurs'  : Relative Humidity (%)
                      'rsds'  : Downward Shortwave Radiation (W/m²)
                      'tas'   : Average temperature (K)
                      'tasmax': Max Daily temperature (K)
                      'tasmin': Min Daily temperature (K)
                      'uas'   : Eastward (Zonal) Wind component near surface (m/s)
                      'vas'   : Northward (Meridional) Wind component near surface (m/s)

        ### Key-Value Arguments:
        - output_filepath : str, optional
               Path to write the output file, by default None
        - saveformat : str, optional
            Format to save the output file, by default 'csv', also accepts 'parquet'
        - float_format : str, optional, default= '%.3f' for printing csv file
        """
        self.read_cmip6_data(cmip_files)

        df_daily = self.generate_annagnps_daily_climate_data_cmip6(**kwargs)

        return df_daily

    def read_cmip5_maca_generate_annagnps_climate_daily(self, cmip_files, **kwargs):
        """
        Generate climate_daily.csv AnnAGNPS file/DataFrame from CMIP5 MACA_v2_METDATA output files

        ### Arguments:
        - cmip_files: Iterable (list), paths to CMIP6 files. The required variables are:
                      'pr'     : Precipitation (mm = kg/m²)
                      'rsds'   : Downward Shortwave Radiation (W/m²)
                      'tasmax' : Max Daily temperature (K)
                      'tasmin' : Min Daily temperature (K)
                      'uas'    : Eastward (Zonal) Wind component near surface (m/s)
                      'vas'    : Northward (Meridional) Wind component near surface (m/s)
                      'vpd'    : Mean vapor pressure deficit (kPa)
                      --- VARIABLES BELOW ARE NOT NECESSARY ---
                      'huss'   : Specific humidity (kg/kg)
                      'rhsmin' : Daily Min Relative Humidity (%)
                      'rhsmax' : Daily Min Relative Humidity (%)

        ### Key-Value Arguments:
        - output_filepath : str, optional
               Path to write the output file, by default None
        - saveformat : str, optional
            Format to save the output file, by default 'csv', also accepts 'parquet'
        - float_format : str, optional, default= '%.3f' for printing csv file
        """
        self.read_cmip5_maca_data(cmip_files)
    
        df_daily = self.generate_annagnps_daily_climate_data_cmip5_maca(**kwargs)

        return df_daily

    def read_nldas_daily_generate_climate_daily(self, path_nldas_daily_files, augment=True, **kwargs):
        """
        Read and generate AnnAGNPS ready climate daily files from NLDAS-2 data

        # Arguments:
        - path_nldas_daily_files: str, path to the folder containing the daily files
        - augment: bool, whether to augment the data with the additional AnnAGNPS required variables
        """

        # ### Key-Value Arguments:
        # - output_dir : str, path optional
        #        Path to write the output file, by default current working directory
        # - saveformat : str, optional
        #     Format to save the output file, by default None, also accepts 'parquet' and 'csv'. 
        #     If None it will not write a file
        # - float_format : str, optional, default= '%.3f' for printing csv file
        # - return_dataframes : bool, optional default False

        # # Outputs
        # - results_df : returns a dictionary with dataframes if return_dataframes = True. The keys are (x,y) tuples
        #                default is False

        self.read_nldas_daily_data(path_nldas_daily_files, augment=augment)
        # self._select_nldas_file_coords_timeslice_data()
        results_df = self.generate_annagnps_daily_climate_data_from_nldas_daily(**kwargs)

        return results_df

    def read_cmip6_data(self, cmip_files):
        """
        Reads CMIP6 data files from the ScenarioMIP Activity and "day" table.

        ### Arguments:
        - cmip_files: Iterable (list), paths to CMIP6 files. The required variables are:
                      'pr'    : Precipitation (kg/(m².s))
                      'hurs'  : Relative Humidity (%)
                      'rsds'  : Downward Shortwave Radiation (W/m²)
                      'tas'   : Average temperature (K)
                      'tasmax': Max Daily temperature (K)
                      'tasmin': Min Daily temperature (K)
                      'uas'   : Eastward (Zonal) Wind component near surface (m/s)
                      'vas'   : Northward (Meridional) Wind component near surface (m/s)
        """
        # We read in chunks leveraging dask "lazy reading" as to not overcharge the memory
        self.ds = xr.open_mfdataset(
            cmip_files, compat="override", chunks={"time": "auto"}, parallel=True
        )

        loaded_variables = list(self.ds.keys())

        required_vars = ["pr", "hurs", "rsds", "tas", "tasmax", "tasmin", "uas", "vas"]

        for var in required_vars:
            if var not in loaded_variables:
                raise ValueError(
                    f"Variable {var} missing from input files, make sure all required input files are provided"
                )

    def read_cmip5_maca_data(self, cmip_files):
        """
        Reads CMIP5 MACA resampled data table (Tested with MACA_v2_METDATA).

        ### Arguments:
        - cmip_files: Iterable (list), paths to CMIP6 files. The required variables are:
                      'pr'     : Precipitation (mm = kg/m²)
                      'rsds'   : Downward Shortwave Radiation (W/m²)
                      'tasmax' : Max Daily temperature (K)
                      'tasmin' : Min Daily temperature (K)
                      'uas'    : Eastward (Zonal) Wind component near surface (m/s)
                      'vas'    : Northward (Meridional) Wind component near surface (m/s)
                      'vpd'    : Mean vapor pressure deficit (kPa)
                      --- VARIABLES BELOW ARE NOT NECESSARY
                      'huss'   : Specific humidity (kg/kg)
                      'rhsmin' : Daily Min Relative Humidity (%)
                      'rhsmax' : Daily Min Relative Humidity (%)
        """

        # Function to preprocess each file before merging the dataset
        def preprocess(ds):
            if "air_temperature" in ds:
                if "minimum" in ds["air_temperature"].cell_methods:
                    ds = ds.rename_vars({"air_temperature": "tasmin"})
                elif "maximum" in ds["air_temperature"].cell_methods:
                    ds = ds.rename_vars({"air_temperature": "tasmax"})
            elif "relative_humidity" in ds:  # Those are not necessary
                if "minimum" in ds["relative_humidity"].cell_methods:
                    ds = ds.rename_vars({"relative_humidity": "rhsmin"})
                elif "maximum" in ds["relative_humidity"].cell_methods:
                    ds = ds.rename_vars({"relative_humidity": "rhsmax"})
            elif "eastward_wind" in ds:
                ds = ds.rename_vars({"eastward_wind": "uas"})
            elif "northward_wind" in ds:
                ds = ds.rename_vars({"northward_wind": "vas"})
            elif "precipitation" in ds:
                ds = ds.rename_vars({"precipitation": "pr"})
            elif "surface_downwelling_shortwave_flux_in_air" in ds:
                ds = ds.rename_vars(
                    {"surface_downwelling_shortwave_flux_in_air": "rsds"}
                )
            elif "vpd" in ds:
                ds["vpd"] = ds["vpd"] * 1000  # Convert kPa --> Pa
                ds["vpd"].attrs["units"] = "Pa"
            return ds

        self.ds = xr.open_mfdataset(
            cmip_files, chunks={"time": "auto"}, preprocess=preprocess, parallel=True
        )

        loaded_variables = list(self.ds.keys())

        required_vars = ["pr", "rsds", "tasmax", "tasmin", "uas", "vas", "vpd"]

        for var in required_vars:
            if var not in loaded_variables:
                raise ValueError(
                    f"Variable {var} missing from input files, make sure all required input files are provided"
                )
            
    def read_nldas_daily_data(self, path_nldas_daily_files, augment=True):
        """ Read NLDAS-2 daily data from files downloaded from GESDISC. 

        These files have the same structure as the hourly files but were pre aggregated
        at a daily scale (i.e. one file per day).

        # Arguments:
        - path_nldas_daily_files: str, path to the folder containing the daily files
        - augment: bool, whether to augment the data with the additional AnnAGNPS required variables
        """

        def preprocess(ds):
            try:
                del ds['spatial_ref']
            except:
                pass
            return ds

        path_nldas_daily_files = Path(path_nldas_daily_files)

        print("Opening daily files...")

        # try:
        #     print('trying open_mfdataset with parallel and auto chunking')
        #     self.ds = xr.open_mfdataset(
        #         path_nldas_daily_files.glob("*.nc"),
        #         combine="nested",
        #         concat_dim="time",
        #         preprocess=preprocess,
        #         chunks={"time": "auto"}, parallel=True
        #     )
        # except:

        # self.ds = xr.open_mfdataset(
        #     path_nldas_daily_files.glob("*.nc"),
        #     combine="nested",
        #     concat_dim="time",
        #     preprocess=preprocess,
        #     chunks={"time": 1}, parallel=False)

        # self.ds = xr.open_mfdataset(
        #     path_nldas_daily_files.glob("*.nc"),
        #     combine="nested",
        #     concat_dim="time",
        #     preprocess=preprocess,
        #     chunks={"time": 30}, parallel=True)

        # Get list of files
        start_date = self.start
        end_date = self.end
        list_of_files_to_open = get_file_list(path_nldas_daily_files, start_date, end_date, product="NLDAS_FORA0125_D.2.0")

        datasets = []
        # for file in tqdm(list(path_nldas_daily_files.glob("*.nc"))):
        for file in tqdm(list_of_files_to_open):
            ds = xr.open_dataset(file, chunks={"time": 1})
            ds = preprocess(ds)
            datasets.append(ds)

        print('Concatenating')

        self.ds = xr.concat(datasets, dim="time")

        self.ds = self.ds.sortby('time')

        self.ds = self.ds.chunk({'time': 365})

        print('Finished Concatenating and rechunking')

    def generate_annagnps_daily_climate_data_cmip5_maca(self, **kwargs):
        ### Key-Value Arguments:
        """
        - output_filepath : str, optional
               Path to write the output file, by default None
        - saveformat : str, optional
            Format to save the output file, by default 'csv', also accepts 'parquet'
        - float_format : str, optional, default= '%.3f' for printing csv file
        """

        self._select_cmip_coords_timeslice_data()
        self._generate_cmip5_coords_timeslice_dataframe()

        df_daily = self._generate_climate_file_daily(use_resampled=False, **kwargs)
        return df_daily
    
    def generate_annagnps_daily_climate_data_cmip6(self, **kwargs):
        """
        Generate climate_daily.csv AnnAGNPS file/DataFrame from CMIP6 output files

        ### Key-Value Arguments:
        - output_filepath : str, optional
               Path to write the output file, by default None
        - saveformat : str, optional
            Format to save the output file, by default 'csv', also accepts 'parquet'
        - float_format : str, optional, default= '%.3f' for printing csv file
        """
        self._select_cmip_coords_timeslice_data()
        self._generate_cmip6_coords_timeslice_dataframe()

        df_daily = self._generate_climate_file_daily(use_resampled=False, **kwargs)

        return df_daily
    
    def _select_nldas_file_coords_timeslice_data(self):
        """
        Slices the NLDAS-2 xarray.DataSet for the coords
        """

        if not self.coords:
            raise Exception("Coordinates are missing. Please provide coords!")

        longitudes = np.unique(np.sort([lon for lon, _ in self.coords]))
        latitudes  = np.unique(np.sort([lat for _, lat in self.coords]))

        # variables_to_keep_and_rename = {"PotEvap": "Potential_ET",
        #                                 "Rainf": "Precip", 
        #                                 "SWdown": "Solar_Radiation",
        #                                 "Tdew": "Dew_Point", 
        #                                 "wind_speed": "Wind_Speed", 
        #                                 "wind_direction": "Wind_Direction", 
        #                                 "Tmin": "Min_Air_Temperature", 
        #                                 "Tmax": "Max_Air_Temperature"}
        
        # Throw an error if the self.start and self.end (included) are not in the dataset
        if self.ds.time.values.min() > pd.to_datetime(self.start):
            raise Exception(f"Start date {self.start} beyond dataset range")
        if self.ds.time.values.max() < pd.to_datetime(self.end):
            raise Exception(f"End date {self.end} beyond dataset range")

        time_mask = (self.ds.time >= self.start) & (self.ds.time <= self.end)
        self.ds_select = self.ds.isel(time=time_mask).sel(lat=latitudes, lon=longitudes, method="nearest")

        # self.ds_select = self.ds.sel(lat=latitudes, lon=longitudes, method="nearest").sel(
        #     time=slice(self.start, self.end)
        # )

        # Rename columns to match AnnAGNPS
        # self.ds_select = self.ds_select.rename_vars(variables_to_keep_and_rename)
        self.coords_actual = self.coords

    def generate_annagnps_daily_climate_data_from_nldas_daily(self, **kwargs):
        """
        Generate climate_daily.csv AnnAGNPS file/DataFrame from NLDAS output files

        ### Key-Value Arguments:
        - output_dir : str, path optional
            Path to write the output file, by default current working directory
        - saveformat : str, optional
            Format to save the output file, by default None, also accepts 'parquet' and 'csv'. 
            If None it will not write a file
        - engine : sqlalchemy.engine, optional
            SQLAlchemy engine to connect to the database, by default None
        - db_table_name : str, optional by default 'climate_nldas2'
        - MAXITER_SINGLE_STATION : int, optional by default 10. Number of attempts to upload data to database for a single station
        - float_format : str, optional, default= '%.3f' for printing csv file
        - return_dataframes : bool, optional default False

        # Outputs
        - results_df : returns a dictionary with dataframes if return_dataframes = True. 
                    The keys are the station_id of each NDLAS-2 cell
        """
        # Unpack kwargs
        output_dir = kwargs.get("output_dir", Path.cwd())
        saveformat = kwargs.get("saveformat", None)
        float_format = kwargs.get("float_format", "%.3f")
        return_dataframes = kwargs.get("return_dataframes", False)
        engine = kwargs.get("engine", None)
        db_table_name = kwargs.get("db_table_name", "climate_nldas2")
        MAXITER_SINGLE_STATION = kwargs.get("MAXITER_SINGLE_STATION", 10)

        variables_to_keep_and_rename = {
            "PotEvap": "Potential_ET", "Rainf": "Precip", 
            "SWdown": "Solar_Radiation", "Tdew": "Dew_Point", 
            "wind_speed": "Wind_Speed", "wind_direction": "Wind_Direction", 
            "Tmin": "Min_Air_Temperature", "Tmax": "Max_Air_Temperature"
        }

        if not output_dir.exists() and saveformat not in [None, "database"]:
            print(f'Creating output directory {output_dir}')
            output_dir.mkdir(parents=True, exist_ok=True)

        if saveformat in ['csv', 'parquet']:
                # Create temp directory
                output_dir_temp = output_dir / "temp"
                output_dir_temp.mkdir(parents=True, exist_ok=True)

                # Delete the contents of the temp directory
                for f in output_dir_temp.iterdir():
                    f.unlink()

        # Select and process all data at once
        self._select_nldas_file_coords_timeslice_data()
        ds_all = self.ds_select

        print("Computing additional climate variables...")
        ds_all = augment_forcing_data(ds_all)
        ds_all = ds_all[variables_to_keep_and_rename.keys()].rename(variables_to_keep_and_rename)

        # Convert to Dask DataFrame
        df_all = ds_all.to_dask_dataframe()

        # Apply transformations
        df_all = df_all.assign(
            Day=df_all.time.dt.day,
            Month=df_all.time.dt.month,
            Year=df_all.time.dt.year,
            Max_Air_Temperature=df_all.Max_Air_Temperature - 273.15,
            Min_Air_Temperature=df_all.Min_Air_Temperature - 273.15,
            Dew_Point=df_all.Dew_Point - 273.15,
            Sky_Cover=None,
            Storm_Type_ID=None,
            Actual_ET=None,
            Actual_EI=None,
            Input_Units_Code=1
        )

        # Reorder columns
        df_all = df_all[["lon","lat",
            "Month", "Day", "Year", "Max_Air_Temperature", "Min_Air_Temperature",
            "Precip", "Dew_Point", "Sky_Cover", "Wind_Speed", "Wind_Direction",
            "Solar_Radiation", "Storm_Type_ID", "Potential_ET", "Actual_ET",
            "Actual_EI", "Input_Units_Code"
        ]]

        # Optimize memory usage
        df_all = df_all.astype({
            "Month": "int8", "Day": "int8", "Year": "int16",
            "Max_Air_Temperature": "float32", "Min_Air_Temperature": "float32",
            "Precip": "float32", "Dew_Point": "float32",
            "Wind_Speed": "float32", "Wind_Direction": "float32",
            "Solar_Radiation": "float32", "Potential_ET": "float32",
            "Input_Units_Code": "int8"
        })

        results_df = {}
        all_stations = set()

        # Get total memory
        total_memory = psutil.virtual_memory().total

        # Process in chunks
        print("Repartitioning and processing in chunks...")
        df_all = df_all.repartition(partition_size=int(total_memory * 0.5)) # Using 50% of total memory
        chunks = df_all.partitions
        num_partitions = df_all.npartitions
        chunk_id = 0
        for chunk in tqdm(chunks, desc="Processing chunks", total=num_partitions):
            chunk_id += 1
            processed_chunk = chunk.compute()
            
            if saveformat == 'database':
                for (lon, lat), group in processed_chunk.groupby(['lon', 'lat']):
                    _, _, station_id = get_grid_position(lon, lat)

                    
                    # Get available dates already in the database for this station
                    available_dates = get_available_dates_for_station(station_id, engine, table=db_table_name)

                    # Get missing dates that need to be added
                    missing_dates = get_missing_dates(available_dates, self.start_input, self.end_input)

                    for iter_station in range(MAXITER_SINGLE_STATION):
                        try:
                            # Select and concatenate climate data for each continuous period
                            gdf_clm = None
                            continuous_periods = find_continuous_periods(missing_dates)

                            for period in continuous_periods:
                                start, end = period[0], period[-1]
                                df_period = group[(group.index >= start) & (group.index <= end)]
                                
                                if len(df_period) == 0:
                                    continue

                                gdf_clm_period = prepare_annagnps_climate_for_db(df_period, station_id, lon, lat)
                                
                                if gdf_clm is None:
                                    gdf_clm = gdf_clm_period
                                else:
                                    gdf_clm = pd.concat([gdf_clm, gdf_clm_period])

                            if gdf_clm is not None:
                                insert_climate_nldas2(gdf_clm, engine, table=db_table_name)
                            break

                        except Exception as e:
                            print(f"Station {station_id}, x = {lon}, y = {lat}: Failed with error: {e}: RETRYING ({iter_station+1}/{MAXITER_SINGLE_STATION})")
                            time.sleep(1)

            elif saveformat in ['csv', 'parquet']:
            # Write chunks to disk
                for (lon, lat), group in processed_chunk.groupby(['lon', 'lat']):
                    _, _, station_id = get_grid_position(lon, lat)
                    all_stations.add(station_id)

                    output_filepath = output_dir_temp / f"climate_daily_{station_id}_chunk_{chunk_id}.{saveformat}"
                    if saveformat == 'csv':
                        group.to_csv(output_filepath, index=False, float_format=float_format)
                    elif saveformat == 'parquet':
                        group.to_parquet(output_filepath, index=True)

            if return_dataframes:
                for (lon, lat), group in processed_chunk.groupby(['lon', 'lat']):
                    _, _, station_id = get_grid_position(lon, lat)

                    if 'climate_data' not in results_df[station_id]:
                        results_df[station_id] = {
                            'climate_data': group.copy(deep=True),
                            'coords': (lon, lat)
                        }
                    else:
                        results_df[station_id]['climate_data'] = pd.concat([results_df[station_id]['climate_data'], group.copy(deep=True)])#, ignore_index=True)

        if saveformat in ['csv', 'parquet']:
            # Loop through all stations and chunks in the temp directory and concatenate them and write them into a single file in the main directory
            for station_id in all_stations:
                output_filepath = output_dir / f"climate_daily_{station_id}.{saveformat}"
                chunks = Path(output_dir_temp).glob(f"climate_daily_{station_id}_chunk_*.{saveformat}")

                if saveformat == 'csv':
                    df = pd.concat([pd.read_csv(chunk) for chunk in chunks])
                    df.to_csv(output_filepath, index=False, float_format=float_format)
                elif saveformat == 'parquet':
                    df = pd.concat([pd.read_parquet(chunk) for chunk in chunks])
                    df.to_parquet(output_filepath, index=True)

            # Delete output_dir_temp and all its contents
            shutil.rmtree(output_dir_temp)


        return results_df if return_dataframes else None

    # def generate_annagnps_daily_climate_data_from_nldas_daily(self, **kwargs):
    #     """
    #     Generate climate_daily.csv AnnAGNPS file/DataFrame from NLDAS output files

    #     ### Key-Value Arguments:
    #     - output_dir : str, path optional
    #            Path to write the output file, by default current working directory
    #     - saveformat : str, optional
    #         Format to save the output file, by default None, also accepts 'parquet' and 'csv'. 
    #         If None it will not write a file
    #     - engine : sqlalchemy.engine, optional
    #         SQLAlchemy engine to connect to the database, by default None
    #     - db_table_name : str, optional by default 'climate_nldas2'
    #     - MAXITER_SINGLE_STATION : int, optional by default 10. Number of attempts to upload data to database for a single station
    #     - float_format : str, optional, default= '%.3f' for printing csv file
    #     - return_dataframes : bool, optional default False


    #     # Outputs
    #     - results_df : returns a dictionary with dataframes if return_dataframes = True. 
    #                    The keys are the station_id of each NDLAS-2 cell

    #     """
    #     # Unpack kwargs
    #     if "output_dir" in kwargs:
    #         output_dir = kwargs["output_dir"]
    #     else:
    #         output_dir = Path.cwd()

    #     if "saveformat" in kwargs:
    #         saveformat = kwargs["saveformat"]
    #     else:
    #         saveformat = None

    #     if "float_format" in kwargs:
    #         float_format = kwargs["float_format"]
    #     else:
    #         float_format = "%.3f"

    #     if "return_dataframes" in kwargs:
    #         return_dataframes = kwargs["return_dataframes"]
    #     else:
    #         return_dataframes = False

    #     if "engine" in kwargs:
    #         engine = kwargs["engine"]
    #     else:
    #         engine = None

    #     if "db_table_name" in kwargs:
    #         db_table_name = kwargs["db_table_name"]
    #     else:
    #         db_table_name = "climate_nldas2"

    #     if "MAXITER_SINGLE_STATION" in kwargs:
    #         MAXITER_SINGLE_STATION = kwargs["MAXITER_SINGLE_STATION"]
    #     else:
    #         MAXITER_SINGLE_STATION = 10

    #     self._select_nldas_file_coords_timeslice_data()

    #     variables_to_keep_and_rename = {"PotEvap": "Potential_ET",
    #                                     "Rainf": "Precip", 
    #                                     "SWdown": "Solar_Radiation",
    #                                     "Tdew": "Dew_Point", 
    #                                     "wind_speed": "Wind_Speed", 
    #                                     "wind_direction": "Wind_Direction", 
    #                                     "Tmin": "Min_Air_Temperature", 
    #                                     "Tmax": "Max_Air_Temperature"}


    #     if not output_dir.exists() and ((saveformat is not None) or (saveformat != "database")):
    #         print(f'Creating output directory {output_dir}')
    #         output_dir.mkdir(parents=True, exist_ok=True)

    #     results_df = {}

    #     for lon, lat in tqdm(list(itertools.product(self.ds_select.lon.values, self.ds_select.lat.values))):
            
    #         ds_curr = self.ds_select.sel(lon=lon, lat=lat, method="nearest")

    #         ds_curr = augment_forcing_data(ds_curr)
    #         ds_curr = ds_curr[variables_to_keep_and_rename.keys()]

    #         ds_curr = ds_curr.rename_vars(variables_to_keep_and_rename)

    #         ds_curr = ds_curr.sortby('time')

    #         df = ds_curr.to_dataframe()

    #         df = df.copy()

    #         # Readjust index so that it matches the input
    #         df.index = df.index - (df.index[0] - self.start_input)

    #         # Get the time series for the nearest grid cell
    #         df["Day"] = df.index.day
    #         df["Month"] = df.index.month
    #         df["Year"] = df.index.year

    #         # Convert temperatures to Celsius
    #         df["Max_Air_Temperature"] = df["Max_Air_Temperature"] - 273.15
    #         df["Min_Air_Temperature"] = df["Min_Air_Temperature"] - 273.15            
    #         df["Dew_Point"] = df["Dew_Point"] - 273.15

            

    #         # Add blank columns for the other variables
    #         df["Sky_Cover"] = None
    #         df["Storm_Type_ID"] = None
    #         df["Actual_ET"] = None
    #         df["Actual_EI"] = None

    #         # We use SI units because we respect ourserlves
    #         df["Input_Units_Code"] = 1


    #         # Reorder columns
    #         df = df[
    #             [
    #                 "Month",
    #                 "Day",
    #                 "Year",
    #                 "Max_Air_Temperature",
    #                 "Min_Air_Temperature",
    #                 "Precip",
    #                 "Dew_Point",
    #                 "Sky_Cover",
    #                 "Wind_Speed",
    #                 "Wind_Direction",
    #                 "Solar_Radiation",
    #                 "Storm_Type_ID",
    #                 "Potential_ET",
    #                 "Actual_ET",
    #                 "Actual_EI",
    #                 "Input_Units_Code",
    #             ]
    #         ]

    #         # Optimize memory usage
    #         df = df.astype(
    #             {
    #                 "Month": "int8",
    #                 "Day": "int8",
    #                 "Year": "int16",
    #                 "Max_Air_Temperature": "float32",
    #                 "Min_Air_Temperature": "float32",
    #                 "Precip": "float32",
    #                 "Dew_Point": "float32",
    #                 # 'Sky_Cover': 'float32',
    #                 "Wind_Speed": "float32",
    #                 "Wind_Direction": "float32",
    #                 "Solar_Radiation": "float32",
    #                 # 'Storm_Type_ID': 'int8',
    #                 "Potential_ET": "float32",
    #                 # 'Actual_ET': 'float32',
    #                 # 'Actual_EI': 'float32',
    #                 "Input_Units_Code": "int8",
    #             }
    #         )

    #         _, _, station_id = get_grid_position(lon, lat)

    #         if return_dataframes:
    #             results_df[station_id] = {
    #                                         'climate_data' : df.copy(deep=True),
    #                                         'coords' : (lon, lat)
    #                                     }


    #         output_filepath = output_dir / f"climate_daily_{station_id}.{saveformat}"
    #         if saveformat == "csv":
    #             df.to_csv(output_filepath, index=False, float_format=float_format)

    #         elif saveformat == "parquet":
    #             df.to_parquet(output_filepath, index=True)

    #         elif saveformat == "database":
    #             # Get available dates already in the database for this station
    #             available_dates = get_available_dates_for_station(station_id, engine, table=db_table_name)

    #             # Get missing dates that need to be added
    #             missing_dates = get_missing_dates(available_dates, self.start_input, self.end_input)

    #             for iter_station in range(MAXITER_SINGLE_STATION):
    #                 try:
    #                     # Select and concatenate climate data for each continuous period
    #                     gdf_clm = None
    #                     continuous_periods = find_continuous_periods(missing_dates)

    #                     for period in continuous_periods:
    #                         start, end = period[0], period[-1]

    #                         # Filter df to the current period
    #                         df_period = df[(df.index >= start) & (df.index <= end)]

    #                         if len(df_period) == 0:
    #                             continue

    #                         gdf_clm_period = prepare_annagnps_climate_for_db(df_period, station_id, lon, lat)
                            
    #                         if gdf_clm is None:
    #                             gdf_clm = gdf_clm_period
    #                         else:
    #                             gdf_clm = pd.concat([gdf_clm, gdf_clm_period])

    #                     if gdf_clm is None:
    #                         print(f"0 rows to insert for station {station_id}.")
    #                     else:
    #                         # print(f"Inserting {len(gdf_clm)} rows for station {station_id}.")
    #                         insert_climate_nldas2(gdf_clm, engine, table=db_table_name)
    #                     break

    #                 except Exception as e:
    #                     print(f"Station {station_id}, x = {lon}, y = {lat}: Failed with error: {e}: RETRYING ({iter_station+1}/{MAXITER_SINGLE_STATION})")
    #                     time.sleep(1)


    #                     break

    #         elif saveformat is None:
    #             pass
    #         else:
    #             raise ValueError(f"Invalid saveformat: {saveformat}")

    #     return results_df

    def _select_cmip_coords_timeslice_data(self):
        """
        Slices the CMIP6 or CMIP6-MACAv2-METDATA xarray.DataSet for the coords
        """

        if not self.coords:
            raise Exception("Coordinates are missing. Please provide coords!")

        longitude = (
            self.coords[0] + 360
        ) % 360  # For CMIP6 and CMIP5 the longitude needs to be in the [0, 360[ range
        latitude = self.coords[1]
        # self.ds_select = self.ds.sel(lat=latitude, lon=longitude, method="nearest").sel(
        #     time=slice(self.start.floor('D'), self.end.floor('D'))
        # )
        self.ds_select = self.ds.sel(lat=latitude, lon=longitude, method="nearest").sel(
            time=slice(self.start, self.end + pd.Timedelta('1D'))
        )

        self.coords_actual = (
            float(self.ds_select["lon"]) - 360,
            float(self.ds_select["lat"]),
        )

        self.ds_select = self.ds_select.drop(
            ["height", "time_bounds", "lat", "lon", "crs"], errors="ignore"
        )

    def _generate_cmip6_coords_timeslice_dataframe(self):
        """
        Generate a DataFrame from the selected xarray.DataSet
        """

        df = self.ds_select.to_dataframe()

        df["pr"] = (
            df["pr"] * 3600 * 24
        )  # compute total precipitation for one day based on the flux

        # We rename columns using the same variable names as the NLDAS-2 pynldas2 class so that we can
        # use the same methods
        df.rename(
            columns={
                "pr": "prcp",
                "hurs": "RH",
                "tas": "temp",
                "tasmin": "temp_min",
                "tasmax": "temp_max",
                "uas": "wind_u",
                "vas": "wind_v",
            },
            inplace=True,
            errors="ignore",
        )

        self.clm = df
        self.clm_resampled = df

        self._compute_dew_point()
        self._compute_wind_speed()
        self._compute_wind_direction()
        self._keep_annagnps_columns_only()

    def _generate_cmip5_coords_timeslice_dataframe(self):
        """
        Generate a DataFrame from the selected xarray.DataSet
        """

        df = self.ds_select.to_dataframe()

        # We rename columns using the same variable names as the NLDAS-2 pynldas2 class so that we can
        # use the same methods
        df.rename(
            columns={
                "pr": "prcp",
                "tasmin": "temp_min",
                "tasmax": "temp_max",
                "uas": "wind_u",
                "vas": "wind_v",
            },
            inplace=True,
            errors="ignore",
        )

        self.clm = df
        self.clm_resampled = df

        self._compute_dew_point_cmip5_maca()
        self._compute_wind_speed()
        self._compute_wind_direction()
        self._keep_annagnps_columns_only()

    def generate_cmip_lon_lat_secondary_climate_id(self, coords=None) -> int:
        """
        Generate a secondary climate ID based on provided coordinates and a climate dataset.

        ### Parameters:
        - coords : (lon, lat) tuple in (EPSG:4326) 
            A tuple containing latitude and longitude values for which the secondary climate ID is generated.

            Optional. If not provided, the self.coords values will be used

        ### Returns:
        - int
            A unique secondary climate ID calculated from the nearest indices of the provided latitude and longitude in the dataset.

        ### Notes:
        This function calculates the nearest indices for the given latitude and longitude inside the provided climate dataset 
        and generates a unique secondary climate ID using these indices.
        """

        ds = self.ds

        if not coords:
            coords = self.coords

        N_lon = ds.dims['lon']

        lon, lat = coords

        lon = (lon + 360) % 360  # For CMIP6 and CMIP5 the longitude needs to be in the [0, 360[ range

        # find nearest index inside ds for each lat lon
        idx_lat = abs(ds['lat']-lat).argmin().item()
        idx_lon = abs(ds['lon']-lon).argmin().item()

        # Generate unique id for this pair
        return idx_lon + idx_lat * N_lon
    
    def _compute_wind_speed(self):
        clm = self.clm
        clm["wind_speed"] = compute_wind_speed(clm["wind_u"], clm["wind_v"])

    def _compute_wind_direction(self, method="annagnps"):
        clm = self.clm
        clm["wind_direction"] = compute_wind_direction(
            clm["wind_u"], clm["wind_v"], method=method
        )

    def _compute_RH(self):
        clm = self.clm
        clm["RH"] = compute_RH(clm["psurf"], clm["temp"], clm["humidity"])

    def _compute_dew_point(self):
        clm = self.clm
        clm["tdew"] = compute_dew_point(clm["RH"], clm["temp"])

    def _compute_dew_point_cmip5_maca(self):
        clm = self.clm
        # Compute Tavg
        clm["tavg"] = (clm["temp_min"] + clm["temp_max"]) * 0.5
        # Compute esat
        clm["esat"] = compute_esat(clm["tavg"], Tunit="K")
        # Compute RH
        clm["RH"] = 100 * (1 - clm["vpd"] / clm["esat"])
        # Compute Tdew
        clm["tdew"] = compute_dew_point(clm["RH"], clm["tavg"], Tunit="K")

    def _keep_annagnps_columns_only(self):
        """
        Once extra variables have been computed (Tdew, Wind Speed, Wind Direction) some variables are no longer needed
        for AnnAGNPS climate file generation and can be dropped
        """
        columns_to_remove = [
            "RH",
            "psurf",
            "temp",
            "vp (Pa)",
            "esat",
            "humidity",
            "wind_u",
            "wind_v",
        ]

        if self.clm is not None:
            self.clm.drop(columns=columns_to_remove, errors="ignore", inplace=True)

        if self.clm_resampled is not None:
            self.clm_resampled.drop(
                columns=columns_to_remove, errors="ignore", inplace=True
            )

    def _generate_climate_file_daily(
        self,
        use_resampled=True,
        output_filepath=None,
        saveformat="csv",
        float_format="%.3f",
    ):
        """Generate a climate file for a given period. Returns a DataFrame and writes to csv if output_filepath is provided.

        Parameters
        ----------
        df : pandas.DataFrame with all necessary variables
        use_resampled: if True it will used the resampled dataframe if it exists, if not it will use the current one
        output_filepath : str, optional
            Path to write the output file, by default None
        saveformat : str, optional
            Format to save the output file, by default 'csv', also accepts 'parquet'
        float_format : str, optional default = '%.3f'

        Returns
        -------
        pandas.DataFrame
            DataFrame containing the climate data properly formatted
        """

        if use_resampled:
            df = self.clm_resampled.copy(deep=True)
        else:
            df = self.clm.copy(deep=True)

        # Columns to keep:
        columns_to_keep = set(
            [
                "Month",
                "Day",
                "Year",
                "Max_Air_Temperature",
                "Min_Air_Temperature",
                "Precip",
                "Dew_Point",
                "Sky_Cover",
                "Wind_Speed",
                "Wind_Direction",
                "Solar_Radiation",
                "Storm_Type_ID",
                "Potential_ET",
                "Actual_ET",
                "Actual_EI",
                "Input_Units_Code",
            ]
        )
        # Readjust index so that it matches the input
        df.index = df.index - (df.index[0] - self.start_input)

        # Get the time series for the nearest grid cell
        df["Day"] = df.index.day
        df["Month"] = df.index.month
        df["Year"] = df.index.year

        # Convert temperatures to Celsius
        df["temp_max"] = df["temp_max"] - 273.15
        df["temp_min"] = df["temp_min"] - 273.15
        # df['Tair'] = df['Tair'] - 273.15
        try:
            df["tdew"] = df["tdew"] - 273.15
        except KeyError:
            df["tdew"] = None
        # Total Solar Radiation (W/m2)
        df["Solar_Radiation"] = df["rsds"]

        # Add blank columns for the other variables
        df["Sky_Cover"] = None
        df["Storm_Type_ID"] = None
        df["Actual_ET"] = None
        df["Actual_EI"] = None

        # We use SI units because we respect ourserlves
        df["Input_Units_Code"] = 1

        # No need to convert precipitation (PotEvap) to mm/day, if we assume rhow = 1000 kg/m3 = 1 kg/L -> 1 mm = 1 L/m2 = 1 kg/m2
        if "pet" in df:
            df["pet"] = df["pet"].apply(lambda x: max(x, 0) if x is not None else None)
        else:
            df["pet"] = None

        # Rename columns
        df = df.rename(
            columns={
                "wind_speed": "Wind_Speed",
                "wind_direction": "Wind_Direction",
                "prcp": "Precip",
                "pet": "Potential_ET",
                "temp_min": "Min_Air_Temperature",
                "temp_max": "Max_Air_Temperature",
                "tdew": "Dew_Point",
            }
        )

        # Drop useless columns
        columns = set(df.columns)
        useless_columns = columns - columns.intersection(columns_to_keep)
        df = df.drop(columns=useless_columns)

        # Reorder columns
        df = df[
            [
                "Month",
                "Day",
                "Year",
                "Max_Air_Temperature",
                "Min_Air_Temperature",
                "Precip",
                "Dew_Point",
                "Sky_Cover",
                "Wind_Speed",
                "Wind_Direction",
                "Solar_Radiation",
                "Storm_Type_ID",
                "Potential_ET",
                "Actual_ET",
                "Actual_EI",
                "Input_Units_Code",
            ]
        ]

        # Optimize memory usage
        df = df.astype(
            {
                "Month": "int8",
                "Day": "int8",
                "Year": "int16",
                "Max_Air_Temperature": "float32",
                "Min_Air_Temperature": "float32",
                "Precip": "float32",
                "Dew_Point": "float32",
                # 'Sky_Cover': 'float32',
                "Wind_Speed": "float32",
                "Wind_Direction": "float32",
                "Solar_Radiation": "float32",
                # 'Storm_Type_ID': 'int8',
                "Potential_ET": "float32",
                # 'Actual_ET': 'float32',
                # 'Actual_EI': 'float32',
                "Input_Units_Code": "int8",
            }
        )

        if output_filepath is not None:
            if saveformat == "csv":
                df.to_csv(output_filepath, index=False, float_format=float_format)
            elif saveformat == "parquet":
                df.to_parquet(output_filepath, index=True)
            else:
                raise ValueError("Invalid saveformat")

        return df

    def generate_climate_station_file(
        self,
        output_filepath=Path("climate_station.csv"),
        climate_station_name="Station",
        beginning_climate_date=None,
        ending_climate_date=None,
        latitude=None,
        longitude=None,
        elevation="3dep",
        temperature_lapse_rate="",
        precipitation_n="",
        global_storm_type_id="",
        first_elevation_difference="",
        first_elevation_rain_factor="",
        second_elevation_difference="",
        second_elevation_rain_factor="",
        two_year_24_hour_precipitation="",
        calibration_or_areal_correction_coefficient="",
        calibration_or_areal_correction_exponent="",
        minimum_interception_evaporation="",
        maximum_interception_evaporation="",
        winter_storm_type_id="",
        spring_storm_type_id="",
        summer_storm_type_id="",
        autumn_storm_type_id="",
        version=6.00,
        input_units_code=1,
    ):
        """
        Generate the AnnAGNPS climate station file for this current query (unless overriden)
        Will generate a "climate_station.csv" file in the specified output_dir.

        ### Arguments:
        - output_dir (Path): The directory where to save the file. Defaults to the current working directory.
        - version (str or float): The version number. Defaults to 6.00.
        - input_units_code (int): The input units code. Defaults to 1 for SI units.
        - climate_station_name (str): The climate station name. Defaults to "Station".
        - beginning_climate_date (str): The beginning climate date in 'mm/dd/yyyy' format. Defaults to the queried start date
        - ending_climate_date (str): The ending climate date in 'mm/dd/yyyy' format. Defaults to queried end date
        - latitude (float): The latitude. Defaults to queried latitude.
        - longitude (float): The longitude. Defaults to queried longitude.
        - elevation (float/str):
            - (float): the elevation.
            - (str): '3dep' (default) queries the elevation at (longitude, latitude) querying the 3dep service
        - temperature_lapse_rate (float): The temperature lapse rate. Defaults to an empty string.
        - precipitation_n (float): Precipitation N. Defaults to an empty string.
        - global_storm_type_id (str): Global storm type ID. Defaults to an empty string.
        - first_elevation_difference (float): First elevation difference. Defaults to an empty string.
        - first_elevation_rain_factor (float): First elevation rain factor. Defaults to an empty string.
        - second_elevation_difference (float): Second elevation difference. Defaults to an empty string.
        - second_elevation_rain_factor (float): Second elevation rain factor. Defaults to an empty string.
        - two_year_24_hour_precipitation (float): 2 Yr 24 hr Precipitation. Defaults to an empty string.
        - calibration_or_areal_correction_coefficient (float): Calibration or Areal Correction Coefficient. Defaults to an empty string.
        - calibration_or_areal_correction_exponent (float): Calibration or Areal Correction Exponent. Defaults to an empty string.
        - minimum_interception_evaporation (float): Minimum Interception Evaporation. Defaults to an empty string.
        - maximum_interception_evaporation (float): Maximum Interception Evaporation. Defaults to an empty string.
        - winter_storm_type_id (str): Winter storm type ID. Defaults to an empty string.
        - spring_storm_type_id (str): Spring storm type ID. Defaults to an empty string.
        - summer_storm_type_id (str): Summer storm type ID. Defaults to an empty string.
        - autumn_storm_type_id (str): Autumn storm type ID. Defaults to an empty string.
        """

        df = generate_climate_station_file(
            output_filepath=output_filepath,
            climate_station_name=climate_station_name,
            beginning_climate_date=beginning_climate_date
            if beginning_climate_date
            else self.start.strftime("%m/%d/%Y"),
            ending_climate_date=ending_climate_date
            if ending_climate_date
            else self.end.strftime("%m/%d/%Y"),
            latitude=latitude if latitude is not None else self.coords[1],
            longitude=longitude if longitude is not None else self.coords[0],
            elevation=elevation,
            temperature_lapse_rate=temperature_lapse_rate,
            precipitation_n=precipitation_n,
            global_storm_type_id=global_storm_type_id,
            first_elevation_difference=first_elevation_difference,
            first_elevation_rain_factor=first_elevation_rain_factor,
            second_elevation_difference=second_elevation_difference,
            second_elevation_rain_factor=second_elevation_rain_factor,
            two_year_24_hour_precipitation=two_year_24_hour_precipitation,
            calibration_or_areal_correction_coefficient=calibration_or_areal_correction_coefficient,
            calibration_or_areal_correction_exponent=calibration_or_areal_correction_exponent,
            minimum_interception_evaporation=minimum_interception_evaporation,
            maximum_interception_evaporation=maximum_interception_evaporation,
            winter_storm_type_id=winter_storm_type_id,
            spring_storm_type_id=spring_storm_type_id,
            summer_storm_type_id=summer_storm_type_id,
            autumn_storm_type_id=autumn_storm_type_id,
            version=version,
            input_units_code=input_units_code,
        )

        return df


def generate_climate_station_file(
    output_filepath=Path("climate_station.csv"),
    climate_station_name="Station",
    beginning_climate_date=None,
    ending_climate_date=None,
    latitude=None,
    longitude=None,
    elevation="3dep",
    temperature_lapse_rate="",
    precipitation_n="",
    global_storm_type_id="",
    first_elevation_difference="",
    first_elevation_rain_factor="",
    second_elevation_difference="",
    second_elevation_rain_factor="",
    two_year_24_hour_precipitation="",
    calibration_or_areal_correction_coefficient="",
    calibration_or_areal_correction_exponent="",
    minimum_interception_evaporation="",
    maximum_interception_evaporation="",
    winter_storm_type_id="",
    spring_storm_type_id="",
    summer_storm_type_id="",
    autumn_storm_type_id="",
    version=6.00,
    input_units_code=1,
):
    """
    Generate the AnnAGNPS climate station file for this current query (unless overriden)
    Will generate a "climate_station.csv" file in the specified output_dir.

    ### Arguments:
    - output_dir (Path): The directory where to save the file. Defaults to the current working directory.
    - version (str or float): The version number. Defaults to 6.00.
    - input_units_code (int): The input units code. Defaults to 1 for SI units.
    - climate_station_name (str): The climate station name. Defaults to "Station".
    - beginning_climate_date (str): The beginning climate date in 'mm/dd/yyyy' format.
    - ending_climate_date (str): The ending climate date in 'mm/dd/yyyy' format.
    - latitude (float): The latitude.
    - longitude (float): The longitude.
    - elevation (float/str):
        - (float): the elevation.
        - (str): '3dep' (default) queries the elevation at (longitude, latitude) querying the 3dep service
    - temperature_lapse_rate (float): The temperature lapse rate. Defaults to an empty string.
    - precipitation_n (float): Precipitation N. Defaults to an empty string.
    - global_storm_type_id (str): Global storm type ID. Defaults to an empty string.
    - first_elevation_difference (float): First elevation difference. Defaults to an empty string.
    - first_elevation_rain_factor (float): First elevation rain factor. Defaults to an empty string.
    - second_elevation_difference (float): Second elevation difference. Defaults to an empty string.
    - second_elevation_rain_factor (float): Second elevation rain factor. Defaults to an empty string.
    - two_year_24_hour_precipitation (float): 2 Yr 24 hr Precipitation. Defaults to an empty string.
    - calibration_or_areal_correction_coefficient (float): Calibration or Areal Correction Coefficient. Defaults to an empty string.
    - calibration_or_areal_correction_exponent (float): Calibration or Areal Correction Exponent. Defaults to an empty string.
    - minimum_interception_evaporation (float): Minimum Interception Evaporation. Defaults to an empty string.
    - maximum_interception_evaporation (float): Maximum Interception Evaporation. Defaults to an empty string.
    - winter_storm_type_id (str): Winter storm type ID. Defaults to an empty string.
    - spring_storm_type_id (str): Spring storm type ID. Defaults to an empty string.
    - summer_storm_type_id (str): Summer storm type ID. Defaults to an empty string.
    - autumn_storm_type_id (str): Autumn storm type ID. Defaults to an empty string.
    """

    data = {
        "Version": str(version),
        "Input Units Code": input_units_code,
        "Climate Station Name": climate_station_name,
        "Beginning Climate Date 'mm/dd/yyyy'": beginning_climate_date,
        "Ending Climate Date 'mm/dd/yyyy'": ending_climate_date,
        "Latitude": latitude,
        "Longitude": longitude,
        "Elevation": elevation,
        "Temperature Lapse Rate": temperature_lapse_rate,
        "Precipitation N": precipitation_n,
        "Global Storm Type ID": global_storm_type_id,
        "1st Elevation Difference": first_elevation_difference,
        "1st Elevation Rain Factor": first_elevation_rain_factor,
        "2nd Elevation Difference": second_elevation_difference,
        "2nd Elevation Rain Factor": second_elevation_rain_factor,
        "2 Yr 24 hr Precipitation": two_year_24_hour_precipitation,
        "Calibration or Areal Correction Coefficient": calibration_or_areal_correction_coefficient,
        "Calibration or Areal Correction Exponent": calibration_or_areal_correction_exponent,
        "Minimum Interception Evaporation": minimum_interception_evaporation,
        "Maximum Interception Evaporation": maximum_interception_evaporation,
        "Winter Storm Type ID": winter_storm_type_id,
        "Spring Storm Type ID": spring_storm_type_id,
        "Summer Storm Type ID": summer_storm_type_id,
        "Autumn Storm Type ID": autumn_storm_type_id,
    }

    if elevation == "3dep":
        elevation = py3dep.elevation_bycoords(
            (data["Longitude"], data["Latitude"]), source="tep"
        )
        data["Elevation"] = elevation
    elif not (isinstance(elevation, float)):
        raise Exception(
            f"Invalid elevation: {elevation}, please provide an elevation in meters (float) or '3dep' to query the elevation automatically"
        )

    df = pd.DataFrame(data, index=[0])
    df.to_csv(output_filepath, float_format="%1.3f", index=False)

    # # Write header line
    # output_filepath.write_text(f"{','.join(data.keys())}\n")
    # # Write values:
    # with output_filepath.open(mode='a') as f:
    #     f.write(f"{','.join([str(x) for x in data.values()])}")

    return df


def compute_esat(Tair, Tunit="K"):
    """
    Compute saturation vapor pressure

    Parameters:
    -----------
    Tair: Air temperature

    Returns:
    -------
    esat : Saturation vapor pressure
    """
    # Constants
    A1 = 17.625
    B1 = 243.04  # degC
    C1 = 610.94  # Pa

    if Tunit == "K":
        d = 273.15
    elif Tunit == "degC":
        d = 0

    # Compute saturation vapor pressure
    esat = C1 * np.exp(
        A1 * (Tair - d) / (B1 + (Tair - d))
    )  # Magnus formula but not really it's actually Alduchov and Eskrige (1996)

    return esat


def compute_RH(Psurf, Tair, Qair, Tunit="K"):
    """
    Compute relative humidity from air temperature and specific humidity.
    This formula has < 0.4% error for -40 °C < Tair < 50 °C.

    Parameters
    ----------
    Psurf : xarray.DataArray
        Surface pressure [Pa]
    Tair : xarray.DataArray
        Air temperature [K]
    Qair : xarray.DataArray
        Specific humidity [kg/kg]

    Returns
    -------
    RH : xarray.DataArray
        Relative humidity [%]
    """

    es = compute_esat(Tair, Tunit)

    epsilon = 0.622

    # Compute vapor pressure
    e = (
        Qair * Psurf / epsilon
    )  # The real formula should be e = Qair * (Pair - es) / (epsilon * es)
    # Compute relative humidity
    RH = e / es * 100

    # replace values greater than 100 with 100
    RH = RH.where((RH <= 100) | np.isnan(RH), 100)

    return RH


def compute_dew_point(RH, Tair, Tunit="K"):
    """Compute dew point temperature from relative humidity and air temperature.
    This formula has < 0.4% error for -40 °C < Tair < 50 °C.

    (See Mark G. Lawrence "The Relationship between Relative Humidity and Dewpoint Temperature", 2005)

    Parameters
    ----------
    RH : xarray.DataArray
        Relative humidity [%]
    Tair : xarray.DataArray
        Air temperature [K]

    Returns
    -------
    Tdew : xarray.DataArray
        Dew point temperature [K]
    """
    # Constants
    A1 = 17.625
    B1 = 243.04  # degC
    # C1 = 610.94 # Pa

    if Tunit == "K":
        d = 273.15
    elif Tunit == "degC":
        d = 0

    magnus_exp_arg_ratio = A1 * (Tair - d) / (B1 + (Tair - d))

    # Set Relative Humidity to 1% if it ends up being an invalid value
    # If remaining NaN values subsist backward or forward fill where appropriate

    RH = RH.copy()  # Ensure we are working with a copy
    try: # If it's a DataFrame
        RH.loc[RH<0] = 1
        RH.loc[RH>100] = 100
    except:
        RH = RH.where(RH <= 0, 1)
        RH = RH.where(RH >= 100, 100)
    
    try: # If it's a DataFrame
        RH = RH.ffill().bfill()
    except:
        RH = RH.ffill(dim="time").bfill(dim="time")

    log_RH = np.log(RH / 100)

    Tdew = (
        B1 * (log_RH + magnus_exp_arg_ratio) / (A1 - log_RH - magnus_exp_arg_ratio) + d
    )

    return Tdew


def compute_wind_speed(Uwind, Vwind):
    """Compute wind speed from U and V wind components.

    Parameters
    ----------
    Uwind : xarray.DataArray
        U wind component [m/s]
    Vwind : xarray.DataArray
        V wind component [m/s]

    Returns
    -------
    wind_speed : xarray.DataArray
        Wind speed [m/s]
    """
    wind_speed = np.sqrt(Uwind**2 + Vwind**2)

    return wind_speed


def compute_wind_direction(Uwind, Vwind, method="annagnps"):
    """Compute wind direction from U and V wind components.
    Wind coming from North is 0°, East is 90°, South is 180°, West is 270°.

    Parameters
    ----------
    Uwind : xarray.DataArray
        U wind component [m/s]
    Vwind : xarray.DataArray
        V wind component [m/s]
    method : str, optional
        Method to compute wind direction, by default 'annagnps' (° Measured clock-wise from North)
        'typical' (Counterclockwise from South)

    Returns
    -------
    wind_direction : xarray.DataArray
        Wind direction [deg]
    """
    # wind_direction = (np.angle(Uwind + 1j * Vwind, deg=True) + 90) % 360 % For some reason xarray does not like this
    if method == "annagnps":
        wind_direction = np.arctan2(Uwind, Vwind) * 180 / np.pi % 360
    elif method == "typical":
        wind_direction = (np.arctan2(Vwind, Uwind) * 180 / np.pi + 90) % 360
    else:
        raise ValueError(f"Method {method} is not implemented")

    return wind_direction

# METHODS TO DOWNLOAD AND PROCESS NLDAS-2 DATA FROM FILES DOWNLOADED FROM GESDISC
def get_time_indices_from_time_bnds(ds, from_date, to_date):
    """Get time indices from time bounds."""
    time_bnds = ds['time_bnds'].values
    from_date = np.datetime64(from_date)
    to_date = np.datetime64(to_date)

    bounds_later_than_from = time_bnds[:, 0] >= from_date
    bounds_earlier_than_to = time_bnds[:, 1] <= to_date

    idx_from = np.argmax(bounds_later_than_from)
    idx_to = bounds_earlier_than_to.shape[0] - np.argmax(bounds_earlier_than_to[::-1]) - 1

    idx_from = int(idx_from)
    idx_to = int(idx_to)

    return (idx_from, idx_to)

def get_file_list(root_dir, from_date, to_date, product='NLDAS_FORA0125_H.2.0'):
    """Get a list of files for a given product and date range (in datetime format)"""

    target_dates = []
    file_list = []

    if '.002' in product:
        # GRIB files
        product_dict = _NLDAS_PRODUCTS_002
        file_extension = '.002.grb'
    elif '.2.0' in product:
        # NetCDF files
        product_dict = _NLDAS_PRODUCTS_20
        file_extension = '.020.nc'

    time_type = product_dict[product]['type']
    fileroot = product_dict[product]['fileroot']

    if time_type == 'hourly':
        onehour = timedelta(hours=1)
        # Hourly data

        curr_date = from_date

        while curr_date <= to_date:
            target_dates.append(curr_date)
            curr_date += onehour

        for dat in target_dates:
            file_list.append((f'{root_dir}/'
                              f'{fileroot}'
                              f'{dat.year}'
                              f'{dat.month:02d}'
                              f'{dat.day:02d}.'
                              f'{dat.hour:02d}00{file_extension}'
                                  ))

    elif time_type == 'daily':
        oneday = timedelta(days=1)
        curr_date = from_date

        while curr_date <= to_date:
            target_dates.append(curr_date)
            curr_date += oneday

        for dat in target_dates:
            file_list.append((f'{root_dir}/'
                              f'{fileroot}'
                              f'{dat.year}'
                              f'{dat.month:02d}'
                              f'{dat.day:02d}{file_extension}'))

    elif time_type == 'monthly':
        for i in range(month_difference(from_date, to_date) + 1):
            target_dates.append((from_date + relativedelta(months=i)))

        for i in range(len(target_dates)):
            file_list.append((f'{root_dir}/'
                              f'{fileroot}{target_dates[i].year}'
                              f'{target_dates[i].month:02d}{file_extension}'
                                  ))

    elif time_type == 'monthly_climatology':
        print(f'Data product is {product} which is a monthly climatology averaged dataset, from_date and to_date will be ignored and all months will be downloaded')
        for i in range(1, 12 + 1):
            file_list.append((f'{root_dir}/{fileroot}{i:02d}{file_extension}'))

    return file_list

def read_dataset(root_dir, start_date, end_date_or_timedelta, product='NLDAS_FORA0125_H.2.0', augment=True):
    """Read dataset (FORA or NOAH or...).
    
    Parameters
    ----------
    root_dir : str
        Root directory
    start_date : str
        Start date in friend ISO format (YYYY-MM-DD) (hours and timezone can be specified if necessary)
    end_date_or_timedelta : 
        str
        End date (YYYY-MM-DD)
        OR a timedelta object
    product : str, optional
        Product name, by default 'NLDAS_FORA0125_H.2.0'
    augment : bool, optional
        Augment forcing data with additional variables, by default True
    
    Returns
    -------
    forcing_data : xarray.Dataset
        Forcing data or NOAH data
    """

    # Get the list of files
    if isinstance(start_date, str):
        start_date = get_date_from_string(start_date, outputtype=datetime)

    if isinstance(end_date_or_timedelta, str):
        end_date = get_date_from_string(end_date_or_timedelta, outputtype=datetime)
    elif isinstance(end_date_or_timedelta, timedelta) or isinstance(end_date_or_timedelta, relativedelta):
        end_date = start_date + end_date_or_timedelta

    file_list = get_file_list(root_dir, start_date, end_date, product=product)
    # Read data
    data = xr.open_mfdataset(file_list, chunks={'time': 24}, 
                    parallel=False, combine='nested', concat_dim='time', engine='h5netcdf')

    # Because the way the files are timestamp, to get the data that is within the time range it's safer to make sure with
    # the time_bnds variable
    idx_from, idx_to = get_time_indices_from_time_bnds(data, start_date, end_date)
    data = data.isel(time=slice(idx_from, idx_to+1))

    data = data.rio.write_crs(4087)
    
    # Augment forcing data
    if augment and (product == 'NLDAS_FORA0125_H.2.0'):
        data = augment_forcing_data(data)
    
    return data

def augment_forcing_data(forcing_data):
    """Augment forcing data with additional variables.
    
    Parameters
    ----------
    forcing_data : xarray.Dataset
        Forcing data
    
    Returns
    -------
    forcing_data : xarray.Dataset
        Forcing data with additional variables
    """
    # Compute relative humidity
    forcing_data['RH'] = compute_RH(forcing_data['PSurf'], forcing_data['Tair'], forcing_data['Qair'])
    forcing_data['RH'].attrs['units'] = '%'
    forcing_data['RH'].attrs['long_name'] = 'Relative humidity at 2m [%]'
    forcing_data['RH'].attrs['standard_name'] = 'Relative Humidity'
    forcing_data['RH'].attrs['cell_methods'] = 'time: point'
    forcing_data['RH'].attrs['vmin'] = forcing_data['RH'].min().values
    forcing_data['RH'].attrs['vmax'] = forcing_data['RH'].max().values

    # Compute dew point temperature
    forcing_data['Tdew'] = compute_dew_point(forcing_data['RH'], forcing_data['Tair'])
    forcing_data['Tdew'].attrs['units'] = 'K'
    forcing_data['Tdew'].attrs['long_name'] = 'Dew point temperature at 2m [K]'
    forcing_data['Tdew'].attrs['standard_name'] = 'Dew Point Temperature'
    forcing_data['Tdew'].attrs['cell_methods'] = 'time: point'
    forcing_data['Tdew'].attrs['vmin'] = forcing_data['RH'].min().values
    forcing_data['Tdew'].attrs['vmax'] = forcing_data['RH'].max().values

    # Compute wind speed
    forcing_data['wind_speed'] = compute_wind_speed(forcing_data['Wind_E'], forcing_data['Wind_N'])
    forcing_data['wind_speed'].attrs['units'] = 'm/s'
    forcing_data['wind_speed'].attrs['long_name'] = '10-meter above ground Zonal wind speed [m/s]'
    forcing_data['wind_speed'].attrs['standard_name'] = 'Near surface magnitude wind speed'
    forcing_data['wind_speed'].attrs['cell_methods'] = 'time: point'
    forcing_data['wind_speed'].attrs['vmin'] = forcing_data['wind_speed'].min().values
    forcing_data['wind_speed'].attrs['vmax'] = forcing_data['wind_speed'].max().values

    # Compute wind direction
    # pipi = compute_wind_direction(forcing_data['Wind_E'], forcing_data['Wind_N'])
    forcing_data['wind_direction'] = compute_wind_direction(forcing_data['Wind_E'], forcing_data['Wind_N'])
    forcing_data['wind_direction'].attrs['units'] = '°'
    forcing_data['wind_direction'].attrs['long_name'] = 'Wind direction in degrees [°]'
    forcing_data['wind_direction'].attrs['standard_name'] = 'Wind direction'
    forcing_data['wind_direction'].attrs['cell_methods'] = 'time: point'
    forcing_data['wind_direction'].attrs['vmin'] = forcing_data['wind_direction'].min().values
    forcing_data['wind_direction'].attrs['vmax'] = forcing_data['wind_direction'].max().values
    
    return forcing_data

def variable_agg_func_forA(x):
    """Aggregate variables in a dataset"""
    # Loop over the variables
    if x.name in ['Tair', 'Qair', 'PSurf', 'Wind_E', 'Wind_N', 'LWdown', 'CAPE', 'SWdown', 'RH', 'Tdew', 'wind_speed', 'wind_direction']:
        return x.mean(dim='time')
    elif x.name in ['Rainf', 'CRainf_frac', 'PotEvap']:
        return x.sum(dim='time')
    else:
        warnings.warn(f"Unexpected variable name, doesn't know what to do with {x.name}, sum, mean, something else?")

def variable_agg_func_forNOAH(x):
    """Aggregate variables in a dataset"""
    # Loop over the variables
    if x.name in ['SWdown', 'LWdown', 'SWnet', 'LWnet', 'Qle', 'Qh', 'Qg', 'Qf', 'AvgSurfT', 'Albedo', 'SWE', 'SnowDepth', 'SnowFrac', \
                  'SoilT_0_10cm', 'SoilT_10_40cm', 'SoilT_40_100cm', 'SoilT_100_200cm', \
                  'SoilM_0_10cm', 'SoilM_10_40cm', 'SoilM_40_100cm', 'SoilM_100_200cm', 'SoilM_0_100cm', 'SoilM_0_200cm',
                  'RootMoist', 'SMLiq_0_10cm', 'SMLiq_10_40cm', 'SMLiq_40_100cm', 'SMLiq_100_200cm', 'SMAvail_0_100cm', 'SMAvail_0_200cm',
                  'PotEvap', 'ECanop', 'TVeg', 'ESoil', 'SubSnow', 'CanopInt', 'ACond', 'CCond', 'RCS', 'RCT', 'RCQ', 'RCSOL',
                  'RSmin', 'RSMacr', 'LAI', 'GVEG', 'Streamflow']:
        return x.mean(dim='time')
    elif x.name in ['Snowf', 'Rainf', 'Evap', 'Qs', 'Qsb', 'Qsm']:
        return x.sum(dim='time')
    else:
        warnings.warn(f"Unexpected variable name, doesn't know what to do with {x.name}, sum, mean, something else?")
    
def aggregate_forcing_data(ds, start_date=None, end_date_or_timedelta=None, agg_func=variable_agg_func_forA):
    """Aggregate forcing data for a given period.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset containing the forcing data
    start_date : str
        Start date of the period to aggregate
    end_date_or_timedelta : 
        str
        End date (YYYY-MM-DD)
        OR a timedelta object
    agg_func : function, optional
        Function to aggregate the variables, by default variable_agg_func_forA
    """

    if start_date is None:
        start_date = ds.time.min().values # Use the first date in the dataset

    if end_date_or_timedelta is None:
        end_date = ds.time.max().values # Use the last date in the dataset


    if isinstance(start_date, str):
        start_date = get_date_from_string(start_date, outputtype=datetime)

    if isinstance(end_date_or_timedelta, str):
        end_date = get_date_from_string(end_date_or_timedelta)
    elif isinstance(end_date_or_timedelta, timedelta):
        end_date = start_date + end_date_or_timedelta
    elif isinstance(end_date_or_timedelta, datetime):
        end_date = end_date_or_timedelta

    # Identify indices for the start and end date so that we can use time_bnds
    idx_from, idx_to = get_time_indices_from_time_bnds(ds, start_date, end_date)

    # Delete useless variables not needed for AnnAGNPS
    try:
        del ds['time_bnds'] # Remove time bounds
    except:
        pass

    mid_date = start_date + (end_date - start_date)/2

    # Convert all time bounds to np.datetime64
    start_date = np.datetime64(start_date, 'D')
    end_date = np.datetime64(end_date, 'D')
    mid_date = np.datetime64(mid_date, 'D')

    # mid_date = mid_date.astype('datetime64[D]') # Convert to datetime64[D] to avoid issues with time bounds
    dscopy = ds.copy()
    
    try:
        del dscopy['GMT_OFFSET']
    except:
        pass

    # Aggregate forcing data
    agg_data = dscopy.isel(time=slice(idx_from, idx_to + 1)).map(agg_func, keep_attrs=True)
    # agg_data = ds.sel(time=slice(start_date, end_date)).map(agg_func)

    # Add Tmin and Tmax
    agg_data['Tmin'] = ds['Tair'].sel(time=slice(start_date, end_date)).min(dim='time', skipna=True, keep_attrs=True)
    agg_data['Tmax'] = ds['Tair'].sel(time=slice(start_date, end_date)).max(dim='time', skipna=True, keep_attrs=True)

    agg_data['Tmin'].attrs['standard_name'] = 'Near surface 24h minimum temperature'
    agg_data['Tmin'].attrs['long_name'] = '2-meter above ground 24h minimum temperature [K]'

    agg_data['Tmax'].attrs['standard_name'] = 'Near surface 24h maximum temperature'
    agg_data['Tmax'].attrs['long_name'] = '2-meter above ground 24h maximum temperature [K]'

    try:
        del agg_data['Tmin'].attrs['vmin']
        del agg_data['Tmin'].attrs['vmax']
    except:
        pass

    try:
        del agg_data['Tmax'].attrs['vmin']
        del agg_data['Tmax'].attrs['vmax']
    except:
        pass

    agg_data = agg_data.expand_dims(dim={'time': [mid_date]}, axis=0) # Add a time dimension

    agg_data.attrs['time_definition'] = 'aggregated' 

    return agg_data

def aggregate_noah_data(ds, start_date=None, end_date_or_timedelta=None, agg_func=variable_agg_func_forNOAH):
    """Aggregate NOAH data for a given period.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset containing the forcing data
    start_date : str
        Start date of the period to aggregate
    end_date_or_timedelta : 
        str
        End date (YYYY-MM-DD)
        OR a timedelta object
    agg_func : function, optional
        Function to aggregate the variables, by default variable_agg_func_forNOAH
    """

    if start_date is None:
        start_date = ds.time.min().values # Use the first date in the dataset

    if end_date_or_timedelta is None:
        end_date = ds.time.max().values # Use the last date in the dataset


    if isinstance(start_date, str):
        start_date = get_date_from_string(start_date, outputtype=datetime)

    if isinstance(end_date_or_timedelta, str):
        end_date = get_date_from_string(end_date_or_timedelta)
    elif isinstance(end_date_or_timedelta, timedelta):
        end_date = start_date + end_date_or_timedelta
    elif isinstance(end_date_or_timedelta, datetime):
        end_date = end_date_or_timedelta

    # Identify indices for the start and end date so that we can use time_bnds
    idx_from, idx_to = get_time_indices_from_time_bnds(ds, start_date, end_date)

    # Delete useless variables not needed for AnnAGNPS
    try:
        del ds['time_bnds'] # Remove time bounds
    except:
        pass

    mid_date = start_date + (end_date - start_date)/2

    # Convert all time bounds to np.datetime64
    start_date = np.datetime64(start_date, 'D')
    end_date = np.datetime64(end_date, 'D')
    mid_date = np.datetime64(mid_date, 'D')

    # mid_date = mid_date.astype('datetime64[D]') # Convert to datetime64[D] to avoid issues with time bounds
    
    # Aggregate forcing data
    agg_data = ds.isel(time=slice(idx_from, idx_to + 1)).map(agg_func, keep_attrs=True)

    agg_data = agg_data.expand_dims(dim={'time': [mid_date]}, axis=0) # Add a time dimension

    agg_data.attrs['time_definition'] = 'aggregated' 

    return agg_data

def aggregate_forcing_data_daily(ds, start_date=None, end_date_or_timedelta=None, agg_func=variable_agg_func_forA):
    """Aggregate forcing data for a given period.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset containing the forcing data
    start_date : str
        Start date of the period to aggregate
    end_date_or_timedelta : 
        str
        End date (YYYY-MM-DD)
        OR a timedelta object
    agg_func : function, optional
        Function to aggregate the variables, by default variable_agg_func_forA
    """

    if start_date is None:
        start_date = ds.time.min().values # Use the first date in the dataset
        start_date = datetime.fromtimestamp(start_date.astype(datetime)*1e-9, tz=timezone.utc)

    if end_date_or_timedelta is None:
        end_date = ds.time.max().values # Use the last date in the dataset
        end_date = datetime.fromtimestamp(end_date.astype(datetime)*1e-9, tz=timezone.utc)

    if isinstance(start_date, str):
        start_date = get_date_from_string(start_date, outputtype=datetime)

    if isinstance(end_date_or_timedelta, str):
        end_date = get_date_from_string(end_date_or_timedelta)
    elif isinstance(end_date_or_timedelta, timedelta) or isinstance(end_date_or_timedelta, relativedelta):
        end_date = start_date + end_date_or_timedelta
    elif isinstance(end_date_or_timedelta, datetime):
        end_date = end_date_or_timedelta
    

    # forcing_ds['time_bnds'].values == get_date_from_string('1992-11-01T01:00:00') # To work out later
    
    oneday = timedelta(days=1)

    if end_date - start_date < oneday:
        oneday = end_date - start_date # If the period is less than one day, then use the period length

    # List of datasets
    ds_list = []

    curr_date = start_date

    while curr_date < end_date:
        # Aggregate forcing data
        agg_data_tmp = aggregate_forcing_data(ds.copy(), curr_date, curr_date+oneday)

        try:
            del agg_data_tmp['time_bnds'] # Remove time bounds
        except:
            pass

        ds_list.append(agg_data_tmp)

        curr_date += oneday

    agg_data = xr.concat(ds_list, dim='time')

    agg_data.attrs['time_definition'] = 'aggregated_daily' 

    return agg_data

def get_grid_position(lon, lat, 
                      lon_min=-124.9375, lon_max=-67.0625, 
                      lat_min=25.0625, lat_max=52.9375, 
                      total_columns=464, total_rows=224):
    """
    Calculate the grid position of a given longitude and latitude within a specified range.
    Defaults to NLDAS-2 grid

    Args:
        lon (float): The longitude value.
        lat (float): The latitude value.
        lon_min (float, optional): The minimum longitude value of the range. Defaults to -124.9375.
        lon_max (float, optional): The maximum longitude value of the range. Defaults to -67.0625.
        lat_min (float, optional): The minimum latitude value of the range. Defaults to 25.0625.
        lat_max (float, optional): The maximum latitude value of the range. Defaults to 52.9375.
        total_columns (int, optional): The total number of columns in the grid. Defaults to 464.
        total_rows (int, optional): The total number of rows in the grid. Defaults to 224.

    Returns:
        tuple: A tuple containing the column, row, and unique ID of the grid position.
            - column (int): The column index of the grid position.
            - row (int): The row index of the grid position.
            - station_id (int): The unique ID of the grid position.
    """
    
    # Calculate the step size for longitude and latitude
    lon_step = (lon_max - lon_min) / (total_columns - 1)
    lat_step = (lat_max - lat_min) / (total_rows - 1)
    
    # Calculate the column and row
    column = int((lon - lon_min) / lon_step) + 1
    row = int((lat - lat_min) / lat_step) + 1

    # Calculate the unique ID
    station_id = (row - 1) * total_columns + column
    
    return column, row, station_id

def lon_lat_from_station_id(station_id, 
                            lon_min=-124.9375, lon_max=-67.0625, 
                            lat_min=25.0625, lat_max=52.9375, 
                            total_columns=464, total_rows=224):

    # Calculate the step size for longitude and latitude
    lon_step = (lon_max - lon_min) / (total_columns - 1)
    lat_step = (lat_max - lat_min) / (total_rows - 1)

    row = max((station_id // total_columns), 1)
    column = total_columns if (station_id % total_columns) == 0 else (station_id % total_columns)

    # Calculate the longitude and latitude
    lon = lon_min + (column - 1) * lon_step
    lat = lat_min + (row - 1) * lat_step

    return lon, lat


# METHODS TO INTERACT WITH THE DATABASE
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

    clm = clm[["month",
               "Day",
               "Year",
               "max_air_temperature",
               "min_air_temperature",
               "precip",
               "dew_point",
               "sky_cover",
               "wind_speed",
               "wind_direction",
               "solar_radiation",
               "storm_type_id",
               "potential_et",
               "actual_et",
               "actual_ei",
               "input_units_code"]]

    clm = clm.assign(station_id=station_id)

    gdf_clm = gpd.GeoDataFrame(clm, geometry=[Point(xgrid, ygrid)] * len(clm), crs="EPSG:4326")
    gdf_clm.rename(columns={'geometry': 'geom'}, inplace=True)
    gdf_clm.index.name = "date"

    gdf_clm = gdf_clm.set_geometry('geom')

    return gdf_clm

def insert_climate_nldas2(gdf_clm, engine, table="climate_nldas2", if_exists="append"):
    gdf_clm.to_postgis(table, engine, if_exists=if_exists, index=True)

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