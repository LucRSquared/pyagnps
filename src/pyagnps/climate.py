import pynldas2
import pydaymet

import datetime
# import pytz
from timezonefinder import TimezoneFinder

import numpy as np
import pandas as pd


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
            self.start = start
            self.end = end
        elif date_mode == "local":
            # Identify timezone:
            lon, lat = coords
            tf = TimezoneFinder()
            tmz_name = tf.timezone_at(lng=lon, lat=lat)

            # Other method
            start_stamp_local, end_stamp_local = pd.Timestamp(start, tz=tmz_name), pd.Timestamp(end, tz=tmz_name)
            start_stamp_utc = start_stamp_local.tz_convert("UTC").tz_localize(None)
            end_stamp_utc = end_stamp_local.tz_convert("UTC").tz_localize(None)
            self.start = start_stamp_utc
            self.end = end_stamp_utc
        else:
            raise Exception(f"Invalid ``date_mode`` {date_mode}")
        
        self.coords = coords
        self.variables = variables

        self.clm = None # DataFrame

    def query_nldas2_climate(self):
        clm = pynldas2.get_bycoords(
        coords=self.coords,
        start_date=self.start,
        end_date=self.end,
        variables = self.variables,
        source='netcdf',
        n_conn=4)
        
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

    def resample(self, rule="1D"):
        clm = self.clm
        how_dict = {
                'prcp': 'sum',
                'pet': 'sum', 
                'rsds': 'mean',
                'wind_speed': 'mean',
                'wind_direction': 'mean',
                'tdew': 'mean', 
                'temp_min': 'min',
                'temp_max': 'max'}
        pass

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