import requests
import pandas as pd
import geopandas as gpd
from shapely import wkt
import glob, os, subprocess
from pathlib import Path

def download_soil_geodataframe(bbox=None):
    # bbox = (minlon,minlat, maxlon, maxlat) In EPSG:4326 CRS
    # e.g. bbox = (-89.94724,34.22708,-89.76632,34.31553)

    if bbox is None:
        raise Exception('Please provide a bounding box in EPSG:4326 coordinate format (minlon,minlat, maxlon, maxlat)')

    # Prepare the JSON query
    body = {
    "format": "JSON",
    "query": f"select Ma.*, M.mupolygonkey, M.areasymbol, M.nationalmusym, M.mupolygongeo from mupolygon M, muaggatt Ma where M.mupolygonkey in \
    (select * from SDA_Get_Mupolygonkey_from_intersection_with_WktWgs84('polygon(\
        ({bbox[0]} {bbox[1]}, {bbox[2]} {bbox[1]}, {bbox[2]} {bbox[3]}, {bbox[0]} {bbox[3]}, {bbox[0]} {bbox[1]})\
    )')) and M.mukey=Ma.mukey",
    }

    url = "https://sdmdataaccess.sc.egov.usda.gov/TABULAR/post.rest"

    # Send the query and collect the response
    soil_response = requests.post(url, json=body).json()

    # Reshaping Data
    data = {"musym":[],
        "muname":[],
        "mustatus":[],
        "slopegraddcp":[], 
        "slopegradwta":[], 
        "brockdepmin":[],
        "wtdepannmin":[],
        "wtdepaprjunmin":[],
        "flodfreqdcd":[],
        "flodfreqmax":[],
        "pondfreqprs":[],
        "aws025wta":[],
        "aws050wta":[],
        "aws0100wta":[],
        "aws0150wta":[],
        "drclassdcd":[],
        "drclasswettest":[],
        "hydgrpdcd":[],
        "iccdcd":[],
        "iccdcdpct":[],
        "niccdcd":[],
        "niccdcdpct":[],
        "engdwobdcd":[],
        "engdwbdcd":[],
        "engdwbll":[],
        "engdwbml":[],
        "engstafdcd":[],
        "engstafll":[],
        "engstafml":[],
        "engsldcd":[],
        "engsldcp":[],
        "englrsdcd":[],
        "engcmssdcd":[],
        "engcmssmp":[],
        "urbrecptdcd":[],
        "urbrecptwta":[],
        "forpehrtdcp":[],
        "hydclprs":[],
        "awmmfpwwta":[],
        "mukey":[],
        "mupolygonkey":[],
        "areasymbol":[],
        "nationalmusym":[],
        "geometry":[]
}

    for d in soil_response['Table']:
        for i, kv in enumerate(data.items()):
            data[kv[0]].append(d[i])

    df = pd.DataFrame(data)

    df['geometry'] = df['geometry'].apply(wkt.loads)
    gdf = gpd.GeoDataFrame(df, crs='epsg:4326')

    gdf = gpd.clip(gdf, bbox)

    return gdf



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