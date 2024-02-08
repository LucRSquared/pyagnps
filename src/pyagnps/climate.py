from pathlib import Path
import pynldas2
import py3dep
import pydaymet

# import datetime
# import pytz
from timezonefinder import TimezoneFinder

import numpy as np
import pandas as pd

import xarray as xr

# Do not create
import os

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
        - coords: (longitude, latitude) in EPSG:4326 projection
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
                               coords: tuple = None,
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
        """

        if not coords:
            coords = self.coords
            # If no coords were provided to begin with we just keep None for coords
            # and assume the dates are UTC
            if coords is None:
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

        elif date_mode.lower() == "local":
            # Identify timezone:
            lon, lat = coords
            tf = TimezoneFinder()
            tmz_name = tf.timezone_at(lng=lon, lat=lat)
            self.timezone = tmz_name

            start_stamp_local, end_stamp_local = pd.Timestamp(
                start, tz=tmz_name
            ), pd.Timestamp(end, tz=tmz_name)
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
        self.ds_select = self.ds.sel(lat=latitude, lon=longitude, method="nearest").sel(
            time=slice(self.start, self.end)
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

    RH.loc[RH<0] = 1
    RH.loc[RH>100] = 100
    RH = RH.ffill().bfill()

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
