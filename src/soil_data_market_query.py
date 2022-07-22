import requests
import pandas as pd
import glob, os, subprocess
from pathlib import Path

def run_one_query(county_code):

    url = "https://SDMDataAccess.sc.egov.usda.gov/Tabular/SDMTabularService/post.rest"

    sQuery = f''' select
    sa.saverest,
    l.areasymbol,
    l.areaname,
    mu.musym,
    hydgrp,
    kwfact,
    albedodry_r,
    (SELECT CASE when min(resdept_r) is null then '>200' else cast(min(resdept_r) as varchar) END

    from component left outer join corestrictions on component.cokey = corestrictions.cokey where component.cokey = c.cokey and reskind is not null) as restrictiondepthr,
    partdensity,
    c.compname,
    texdesc,
    hzdepb_r,
    dbovendry_r,
    claytotal_r,
    silttotal_r,
    sandtotal_r,
    (select sum(cf.fragvol_r) as fragvol FROM chfrags cf WHERE cf.chkey = ch.chkey ) as fragvol,
    sandvf_r,
    caco3_r,
        ksat_r,
        wthirdbar_r,
        wfifteenbar_r,
        om_r,
        ph1to1h2o_r,
        c.comppct_r
        FROM
        legend l INNER JOIN mapunit mu ON mu.lkey = l.lkey
        LEFT OUTER JOIN sacatalog sa ON sa.areasymbol = l.areasymbol
        LEFT OUTER JOIN component c ON c.mukey = mu.mukey
        LEFT OUTER JOIN chorizon ch ON ch.cokey = c.cokey
        LEFT OUTER JOIN chtexturegrp ct ON ch.chkey=ct.chkey
        WHERE l.areasymbol = '{county_code}' and ct.rvindicator = 'yes'
        Order by l.areasymbol, musym, hzdepb_r'''

    dRequest = dict()
    dRequest["FORMAT"] = "JSON+COLUMNNAME"
    dRequest["QUERY"] = sQuery
    
    resp = requests.post(url, data=dRequest).json()
    
    df = pd.DataFrame(resp['Table'])
    df.columns = resp['Table'][0]
    df = df.drop(labels=[0], axis=0)

    return df

def run_batch_write_files(county_codes, outpath=''):

    for county_code in county_codes:
        df = run_one_query(county_code)
        outfile = Path(outpath,f'{county_code}_nasis.csv')
        df.to_csv(outfile,index=False,sep=',')

def run_nita(filefolder, path_to_NITA_exe, combine_list=None, units_out=None):

    list_of_files = glob.glob(str(Path(filefolder,"*_nasis.csv")))

    # The files are not combined by default
    if combine_list is None:
        combine_list = [0 for _ in list_of_files]
    
    # SI Units by default
    if units_out is None:
        units_out = [1 for _ in list_of_files]

    with open(Path(filefolder,'NITA_CONTROL.csv'), 'w') as control:
        control.write('FILENAME,UNITS_OUT,COMBINE\n')
        for i, file in enumerate(list_of_files):
            filename = Path(file).name
            control.write(f'{filename},{units_out[i]},{combine_list[i]}\n')
    
    command = path_to_NITA_exe+' /f:'+ 'NITA_CONTROL.csv'

    os.chdir(filefolder)
    subprocess.call(command)

def main(county_codes, filefolder, combine_list=None, units_out=None, path_to_NITA_exe='NITA.exe'):

    if not(os.path.isdir(filefolder)):
        os.mkdir(filefolder)
    
    run_batch_write_files(county_codes, outpath=filefolder)

    run_nita(filefolder, path_to_NITA_exe, combine_list, units_out)

        

if __name__ == "__main__":
    # This code assumes that in your current working directory
    # - a filefolder subdirectory will be created 
    # - path_to_NITA_exe is the relative path to the executable with respect to filefolder

    path_to_NITA_exe = '../../src/bins/NITA.exe' # With respect to the nita_files
    county_codes = ['OH107','OH011']

    combine_list = [1, 1]
    units_out = [1, 1]

    filefolder = './outputs/soil'
    main(county_codes,filefolder,combine_list,units_out,path_to_NITA_exe)

