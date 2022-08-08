import sys, os
sys.path.append('..')

from pyagnps import soil_data_market as sdm

path_to_NITA_exe = '../../src/bins/NITA.exe' # With respect to the nita_files
# county_codes = ['OH107','OH011']

# combine_list = [1, 1]
# units_out = [1, 1]

county_codes = ['MS107']
combine_list = [1]
units_out = [1]

filefolder = './outputs/soil_MS107'

if not(os.path.isdir(filefolder)):
    os.mkdir(filefolder)

sdm.run_batch_write_files(county_codes, outpath=filefolder)

sdm.run_nita(filefolder, path_to_NITA_exe, combine_list, units_out)