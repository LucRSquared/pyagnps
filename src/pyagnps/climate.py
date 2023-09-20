from pathlib import Path
import pynldas2
# import pydaymet

import datetime
# import pytz
from timezonefinder import TimezoneFinder

import numpy as np
import pandas as pd

# Do not create 
import os
os.environ["HYRIVER_CACHE_NAME"] = "/tmp/climate_cache.sqlite" # On Windows systems it will be under C:/tmp
# os.environ["HYRIVER_CACHE_DISABLE"] = "true"


class clm_annagnps_coords():
    def __init__(self, 
                 coords: tuple, 
                 start: str = "2022-01-01", 
                 end: str ="2022-12-31", 
                 date_mode: str ="local",
                 variables:list[str] = ['prcp','temp', 'rsds', 'pet', 
                                        'humidity', 'wind_u', 'wind_v', 
                                        'psurf']):
        '''
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
        - variables: str or list of str, optional
                Variables to download. If None, all variables are downloaded.
                Valid variables are: ``prcp`` (precipitation in kg/m² ~ mm), ``pet`` (potential evaporation in kg/m² ~ mm),
                                     ``temp`` (temperature in K), 
                                     ``wind_u``, ``wind_v`` (10-m above ground zonal and meridional winds respectively, in m/s)
                                     ``rlds`` (surface downward longwave radiation flux in W/m²), 
                                     ``rsds`` (surface downward shortwave radiation flux in W/m²),
                                     ``humidity`` (2-m above ground specific humidity kg/kg)
                                     and ``psurf`` (surface pressure in Pa if provided from netcdf NLDAS-2 data rods source``)
        '''

        if date_mode == "utc":
            self.start = pd.Timestamp(start)
            self.end = pd.Timestamp(end)
            self.timezone = "UTC"
        elif date_mode == "local":
            # Identify timezone:
            lon, lat = coords
            tf = TimezoneFinder()
            tmz_name = tf.timezone_at(lng=lon, lat=lat)
            self.timezone = tmz_name

            start_stamp_local, end_stamp_local = pd.Timestamp(start, tz=tmz_name), pd.Timestamp(end, tz=tmz_name)
            start_stamp_utc = start_stamp_local.tz_convert("UTC").tz_localize(None)
            end_stamp_utc = end_stamp_local.tz_convert("UTC").tz_localize(None)
            self.start = start_stamp_utc
            self.end = end_stamp_utc
        else:
            raise Exception(f"Invalid ``date_mode`` {date_mode}")
        
        self.date_mode = date_mode
        
        self.coords = coords
        self.variables = variables

        self.clm = None # DataFrame
        self.clm_resampled = None

    def query_nldas2_climate(self):

        clm = pynldas2.get_bycoords(
        coords=self.coords,
        start_date=self.start,
        end_date=self.end,
        variables = self.variables,
        source='netcdf',
        n_conn=4)

        # Express dates in local mode
        if self.date_mode == "local":
            clm.index = clm.index.tz_convert(self.timezone).tz_localize(None)
        else:
            clm.index = clm.index.tz_localize(None)

        
        self.clm = clm

    def compute_wind_speed(self):
        clm = self.clm
        clm["wind_speed"] = compute_wind_speed(clm["wind_u"], clm["wind_v"])

    def compute_wind_direction(self, method='annagnps'):
        clm = self.clm
        clm["wind_direction"] = compute_wind_direction(clm["wind_u"], clm["wind_v"], method=method)

    def compute_RH(self):
        clm = self.clm
        clm["RH"] = compute_RH(clm["psurf"], clm["temp"], clm["humidity"])

    def compute_dew_point(self):
        clm = self.clm
        clm["tdew"] = compute_dew_point(clm["RH"], clm["temp"])

    def keep_annagnps_columns_only(self):
        """
        Once extra variables have been computed (Tdew, Wind Speed, Wind Direction) some variables are no longer needed
        for AnnAGNPS climate file generation and can be dropped 
        """
        columns_to_remove = ["RH", "psurf", "vp (Pa)", "esat", "humidity", "wind_u", "wind_v"]

        if self.clm is not None:
            self.clm.drop(columns=columns_to_remove, errors="ignore", inplace=True)

        if self.clm_resampled is not None:
            self.clm_resampled.drop(columns=columns_to_remove, errors="ignore", inplace=True)

    def resample(self, rule="1D"):
        clm = self.clm
        full_how_dict = {
                'prcp': 'sum',
                'pet': 'sum', 
                'rsds': 'mean',
                'wind_speed': 'mean',
                'wind_direction': 'mean',
                'tdew': 'mean', 
                'temp_min': 'min',
                'temp_max': 'max',
                'RH': 'mean',
                'temp': 'mean',
                'psurf': 'mean',
                'vp (Pa)': 'mean',
                'esat': 'mean',
                'humidity': 'mean',
                'wind_u': 'mean',
                'wind_v': 'mean'}
        
        # Prepare min/max temp columns
        clm["temp_min"] = clm["temp"]
        clm["temp_max"] = clm["temp"]
        
        # Create aggregation dict 
        vars_how = {var:full_how_dict[var] for var in clm.columns}

        clm = clm.resample(rule=rule).agg(vars_how)
        clm.index = clm.index.tz_localize(None)

        self.clm_resampled = clm

        return clm
    
    def compute_additional_climate_variables(self):
        self.compute_RH()
        self.compute_dew_point()
        self.compute_wind_direction()
        self.compute_wind_speed()

    def query_nldas2_generate_annagnps_climate_daily(self, **kwargs):
        """
        Generate climate_daily.csv AnnAGNPS file/DataFrame

        ### Key-Value Arguments:
        - output_filepath : str, optional
               Path to write the output file, by default None
        - saveformat : str, optional
            Format to save the output file, by default 'csv', also accepts 'parquet'
        - float_format : str, optional, default= '%.3f' for printing csv file
        """
        self.query_nldas2_climate()
        self.compute_additional_climate_variables()
        self.keep_annagnps_columns_only()
        self.resample(rule="1D")
        
        df_daily = self.generate_climate_file_daily(use_resampled=True, **kwargs)
        return df_daily

    def generate_climate_file_daily(self, use_resampled=True, output_filepath=None, saveformat='csv', float_format='%.3f'):
        """ Generate a climate file for a given period. Returns a DataFrame and writes to csv if output_filepath is provided.
        
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
        columns_to_keep = set(['Month', 'Day', 'Year', 'Max_Air_Temperature', 'Min_Air_Temperature', 'Precip', 'Dew_Point',
                'Sky_Cover', 'Wind_Speed', 'Wind_Direction', 'Solar_Radiation', 'Storm_Type_ID',
                'Potential_ET', 'Actual_ET', 'Actual_EI', 'Input_Units_Code'])

        # Get the time series for the nearest grid cell
        df['Day'] = df.index.day
        df['Month'] = df.index.month
        df['Year'] = df.index.year

        # Convert temperatures to Celsius
        df['temp_max'] = df['temp_max'] - 273.15
        df['temp_min'] = df['temp_min'] - 273.15
        # df['Tair'] = df['Tair'] - 273.15
        df['tdew'] = df['tdew'] - 273.15

        # Total Solar Radiation (W/m2)
        df['Solar_Radiation'] = df['rsds']

        # Add blank columns for the other variables
        df['Sky_Cover'] = None
        df['Storm_Type_ID'] = None
        df['Actual_ET'] = None
        df['Actual_EI'] = None

        # We use SI units because we respect ourserlves
        df['Input_Units_Code'] = 1

        # No need to convert precipitation (PotEvap) to mm/day, if we assume rhow = 1000 kg/m3 = 1 kg/L -> 1 mm = 1 L/m2 = 1 kg/m2
        df["pet"] = df["pet"].apply(lambda x: max(x,0))

        # Rename columns
        df = df.rename(columns={'wind_speed': 'Wind_Speed',
                                'wind_direction': 'Wind_Direction',
                                'prcp': 'Precip',
                                'pet': 'Potential_ET',
                                'temp_min': 'Min_Air_Temperature',
                                'temp_max': 'Max_Air_Temperature',
                                'tdew': 'Dew_Point'})
        
        # Drop useless columns
        columns = set(df.columns)
        useless_columns = columns - columns.intersection(columns_to_keep)
        df = df.drop(columns=useless_columns)


        # Reorder columns
        df = df[['Month', 'Day', 'Year', 'Max_Air_Temperature', 'Min_Air_Temperature', 'Precip', 'Dew_Point',
                'Sky_Cover', 'Wind_Speed', 'Wind_Direction', 'Solar_Radiation', 'Storm_Type_ID',
                'Potential_ET', 'Actual_ET', 'Actual_EI', 'Input_Units_Code']]

        # Optimize memory usage
        df = df.astype({'Month': 'int8',
                        'Day': 'int8',
                        'Year': 'int16',
                        'Max_Air_Temperature': 'float32',
                        'Min_Air_Temperature': 'float32',
                        'Precip': 'float32',
                        'Dew_Point': 'float32',
                        # 'Sky_Cover': 'float32',
                        'Wind_Speed': 'float32',
                        'Wind_Direction': 'float32',
                        'Solar_Radiation': 'float32',
                        # 'Storm_Type_ID': 'int8',
                        'Potential_ET': 'float32',
                        # 'Actual_ET': 'float32',
                        # 'Actual_EI': 'float32',
                        'Input_Units_Code': 'int8'})

        if output_filepath is not None:
            if saveformat == 'csv':
                df.to_csv(output_filepath, index=False, float_format=float_format)
            elif saveformat == 'parquet':
                df.to_parquet(output_filepath, index=True)
            else:
                raise ValueError('Invalid saveformat')

        return df
    
    def generate_climate_station_file(self, output_dir=Path().cwd(), **kwargs):
        """
        Generate the AnnAGNPS climate station file for this current query (unless overriden)
        Will generate a "climate_station.csv" file in the specified output_dir.
        """
        # columns = [
        #     "Version",
        #     "Input_Units_Code",
        #     "Climate_Station_Name",
        #     "Beginning_Climate_Date",
        #     "Ending_Climate_Date",
        #     "Latitude",
        #     "Longitude",
        #     "Elevation",
        #     "Temperature_Lapse_Rate",
        #     "Precipitation_N",
        #     "Global_Storm_Type_ID",
        #     "1st_Elevation_Difference",
        #     "1st_Elevation_Rain_Factor",
        #     "2nd_Elevation_Difference",
        #     "2nd_Elevation_Rain_Factor",
        #     "2_Yr_24_hr_Precipitation",
        #     "Calibration_or_Areal_Correction_Coefficient",
        #     "Calibration_or_Areal_Correction_Exponent",
        #     "Minimum_Interception_Evaporation",
        #     "Maximum_Interception_Evaporation",
        #     "Winter_Storm_Type_ID",
        #     "Spring_Storm_Type_ID",
        #     "Summer_Storm_Type_ID",
        #     "Autumn_Storm_Type_ID"
        # ]

        columns_defaults = {
            "Version": 6.00,
            "Input Units Code": 1,
            "Climate Station Name": "Station",
            "Beginning Climate Date 'mm/dd/yyyy'": self.start.strftime("%m/%d/%Y"),
            "Ending Climate Date 'mm/dd/yyyy'": self.end.strftime("%m/%d/%Y"),
            "Latitude": self.coords[1],
            "Longitude": self.coords[0],
            "Elevation": 0,
            "Temperature Lapse Rate": "",
            "Precipitation N": "",
            "Global Storm Type ID": "",
            "1st Elevation Difference": "",
            "1st Elevation Rain Factor": "",
            "2nd Elevation Difference": "",
            "2nd Elevation Rain Factor": "",
            "2 Yr 24 hr Precipitation": "",
            "Calibration or Areal Correction Coefficient": "",
            "Calibration or Areal Correction Exponent": "",
            "Minimum Interception Evaporation": "",
            "Maximum Interception Evaporation": "",
            "Winter Storm Type ID": "",
            "Spring Storm Type ID": "",
            "Summer Storm Type ID": "",
            "Autumn Storm Type ID": ""
        }


    


def compute_RH(Psurf, Tair, Qair, Tunit='K'):
    """Compute relative humidity from air temperature and specific humidity.
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
    # Constants
    A1 = 17.625
    B1 = 243.04 # degC
    C1 = 610.94 # Pa

    if Tunit == 'K':
        d = 273.15
    elif Tunit == 'degC':
        d = 0

    epsilon = 0.622

    # Compute saturation vapor pressure
    es = C1 * np.exp(A1 * (Tair - d) / (B1 + (Tair - d))) # Magnus formula but not really it's actually Alduchov and Eskrige (1996)

    # Compute vapor pressure
    e = Qair * Psurf / epsilon # The real formula should be e = Qair * (Pair - es) / (epsilon * es)
    # Compute relative humidity
    RH = e / es * 100

    # replace values greater than 100 with 100
    RH = RH.where((RH <= 100) | np.isnan(RH), 100)

    return RH

def compute_dew_point(RH, Tair, Tunit='K'):
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
    B1 = 243.04 # degC
    # C1 = 610.94 # Pa

    if Tunit == 'K':
        d = 273.15
    elif Tunit == 'degC':
        d = 0

    magnus_exp_arg_ratio = A1 * (Tair - d) / (B1 + (Tair - d))
    log_RH = np.log(RH / 100)

    Tdew = B1 * (log_RH + magnus_exp_arg_ratio) / (A1 - log_RH - magnus_exp_arg_ratio) + d
    
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

def compute_wind_direction(Uwind, Vwind, method='annagnps'):
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
    if method=='annagnps':
        wind_direction = np.arctan2(Uwind, Vwind) * 180 / np.pi % 360
    elif method == 'typical':
        wind_direction = (np.arctan2(Vwind, Uwind)*180/np.pi + 90) % 360
    else:
        raise ValueError(f"Method {method} is not implemented")
    
    return wind_direction