import requests
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.mask import mask
from rasterio.features import geometry_window
import xarray as xr
import rioxarray
from pyproj import CRS
from rasterstats import zonal_stats
from shapely import wkt
from shapely.geometry import box
import glob, os, subprocess
from pathlib import Path
from tqdm import tqdm


def download_soil_geodataframe(bbox=None):
    # bbox = (minlon,minlat, maxlon, maxlat) In EPSG:4326 CRS
    # e.g. bbox = (-89.94724,34.22708,-89.76632,34.31553)

    if bbox is None:
        raise Exception(
            "Please provide a bounding box in EPSG:4326 coordinate format (minlon,minlat, maxlon, maxlat)"
        )

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

    if not soil_response:
        return None

    # Reshaping Data
    data = {
        "musym": [],
        "muname": [],
        "mustatus": [],
        "slopegraddcp": [],
        "slopegradwta": [],
        "brockdepmin": [],
        "wtdepannmin": [],
        "wtdepaprjunmin": [],
        "flodfreqdcd": [],
        "flodfreqmax": [],
        "pondfreqprs": [],
        "aws025wta": [],
        "aws050wta": [],
        "aws0100wta": [],
        "aws0150wta": [],
        "drclassdcd": [],
        "drclasswettest": [],
        "hydgrpdcd": [],
        "iccdcd": [],
        "iccdcdpct": [],
        "niccdcd": [],
        "niccdcdpct": [],
        "engdwobdcd": [],
        "engdwbdcd": [],
        "engdwbll": [],
        "engdwbml": [],
        "engstafdcd": [],
        "engstafll": [],
        "engstafml": [],
        "engsldcd": [],
        "engsldcp": [],
        "englrsdcd": [],
        "engcmssdcd": [],
        "engcmssmp": [],
        "urbrecptdcd": [],
        "urbrecptwta": [],
        "forpehrtdcp": [],
        "hydclprs": [],
        "awmmfpwwta": [],
        "mukey": [],
        "mupolygonkey": [],
        "areasymbol": [],
        "nationalmusym": [],
        "geometry": [],
    }

    for d in soil_response["Table"]:
        for i, kv in enumerate(data.items()):
            data[kv[0]].append(d[i])

    df = pd.DataFrame(data)

    df["geometry"] = df["geometry"].apply(wkt.loads)
    gdf = gpd.GeoDataFrame(df, crs="epsg:4326")

    gdf = gpd.clip(gdf, bbox)

    return gdf


def download_soil_geodataframe_tiles(bbox=None, tile_size=0.1, explode_geometries=True):
    # bbox = (minlon,minlat, maxlon, maxlat) In EPSG:4326 CRS
    # e.g. bbox = (-89.94724,34.22708,-89.76632,34.31553)
    # tile_size in degrees of latitude or longitude, defaults to 0.1
    # explode_geometries : Boolean (True default), to explode multi-part geometries into single geometries

    if bbox is None:
        raise Exception(
            "Please provide a bounding box in EPSG:4326 coordinate format (minlon,minlat, maxlon, maxlat)"
        )

    min_lon, min_lat, max_lon, max_lat = bbox
    if min_lon > max_lon or min_lat > max_lat:
        raise Exception("Invalid bounding box")

    # calculate the number of chunks to divide the area into
    n_lon_chunks = int((max_lon - min_lon) / tile_size) + 1
    n_lat_chunks = int((max_lat - min_lat) / tile_size) + 1

    # divide the bounding box into smaller sub-boxes
    sub_boxes = []
    for i in range(n_lon_chunks):
        for j in range(n_lat_chunks):
            sub_box = (
                min_lon + i * tile_size,
                min_lat + j * tile_size,
                min_lon + (i + 1) * tile_size,
                min_lat + (j + 1) * tile_size,
            )
            sub_boxes.append(sub_box)

    # download the data for each sub-box and concatenate into a single geodataframe
    gdf_list = []
    for sub_box in tqdm(sub_boxes):
        sub_gdf = download_soil_geodataframe(sub_box)
        gdf_list.append(sub_gdf)

    gdf_list = [gdf_tmp for gdf_tmp in gdf_list if gdf_tmp is not None]

    gdf = gpd.GeoDataFrame(pd.concat(gdf_list, ignore_index=True), crs="epsg:4326")
    gdf = gpd.clip(gdf, bbox)

    columns = gdf.columns

    # Get rid of tile junctions by dissolving the geometries for identical mukeys
    gdf = gdf.dissolve(by=["mukey"], as_index=False)
    gdf = gdf[columns]

    if explode_geometries:
        # Explode in order to seperate MULTIPOLYGONs into simple POLYGONs
        gdf = gdf.explode(ignore_index=True, index_parts=False)

    return gdf

def calculate_zonal_stats(geometry, raster_dataset, agg_method='majority'):
    # Extract the first band of the raster data
    raster_values = raster_dataset.values[0]
    stats = zonal_stats(geometry, raster_values, affine=raster_dataset.rio.transform(), stats=agg_method, all_touched=True)
    return stats[0][agg_method] if len(stats) > 0 else None

def assign_attr_zonal_stats_raster_layer(
    geo_bins, raster_layer, agg_method='majority', attr="mukey"):

    if not isinstance(geo_bins, gpd.GeoDataFrame):
        geo_bins = gpd.read_file(geo_bins)

    UTM_CRS = geo_bins.estimate_utm_crs()

    geo_bins = geo_bins.to_crs(UTM_CRS)

    bbox = geo_bins.total_bounds
    bbox_geometry = box(*bbox)

    raster_data = rioxarray.open_rasterio(raster_layer, masked=True)

    # Reproject the bounding box to the raster CRS
    bbox_gdf = gpd.GeoDataFrame({'geometry': [bbox_geometry]}, crs=UTM_CRS)
    bbox_gdf = bbox_gdf.to_crs(raster_data.rio.crs)
    bbox_transformed = bbox_gdf.total_bounds

    cropped_data = raster_data.rio.clip_box(*bbox_transformed)
    cropped_data = cropped_data.rio.reproject(UTM_CRS.to_string())

    zonal_func = lambda x: calculate_zonal_stats(x, cropped_data, agg_method=agg_method)

    geo_bins[attr] = geo_bins.geometry.apply(zonal_func)

    return geo_bins


# def assign_attr_zonal_stats_raster_layer(
#     geo_bins, raster_layer, agg_method='majority', attr="mukey"):
    
#     # Assign attributes to a GeoDataFrame based on the plurality of a given attribute in a GeoDataFrame
#     # - geo_bins: GeoDataFrame of the cells to perform zonal statistics on
#     # - raster_layer: path to raster file
#     # - attr: attribute of the geo_attributes_layer that should be attributed to geo_bins

#     # Returns geo_bins with a new column called attr

#     if not isinstance(geo_bins, gpd.GeoDataFrame):
#         geo_bins = gpd.read_file(geo_bins)

#     UTM_CRS = geo_bins.estimate_utm_crs()

#     geo_bins = geo_bins.to_crs(UTM_CRS)

#     bbox = geo_bins.total_bounds
#     bbox_geometry = box(*bbox)

#     with rasterio.open(raster_layer) as src:
#         raster_crs = CRS.from_string(src.crs.to_string())
#         bbox_geometry_reprojected = gpd.GeoSeries(bbox_geometry, crs=UTM_CRS).to_crs(raster_crs).iloc[0]

#         window = rasterio.windows.from_bounds(*bbox_geometry_reprojected.bounds, transform=src.transform)
#         cropped_data = src.read(window=window)
#         cropped_transform = src.window_transform(window)

#         # write_cropped_raster_to_file(cropped_data, cropped_transform, src, 'cropped_raster.tif')

#         reprojected_data = np.empty(cropped_data.shape, dtype=cropped_data.dtype)

#         dst_transform = rasterio.transform.from_bounds(*geo_bins.total_bounds, 
#                                                        reprojected_data.shape[1], 
#                                                        reprojected_data.shape[2])
        
#         nodata_value = src.nodata

#         rasterio.warp.reproject(
#             source=cropped_data,
#             src_transform=cropped_transform,
#             src_crs=src.crs,
#             destination=reprojected_data,
#             dst_transform=dst_transform,
#             dst_crs=UTM_CRS,
#             resampling=rasterio.warp.Resampling.bilinear,
#             src_nodata=nodata_value,
#             dst_nodata=nodata_value
#         )

#         # zonal_func = lambda x: calculate_zonal_stats(x, reprojected_data[0], affine=dst_transform, 
#         #                                             agg_method=agg_method)
#         zonal_func = lambda x: calculate_zonal_stats(x, src, affine=cropped_transform, 
#                                                     agg_method=agg_method)

#         geo_bins[attr] = geo_bins.geometry.apply(zonal_func)

#     return geo_bins


def assign_attr_plurality_vector_layer(
    geo_bins, geo_attributes_layer, attr="mukey", bin_id="DN"
):
    # Assign attributes to a GeoDataFrame based on the plurality of a given attribute in a GeoDataFrame
    # - geo_bins: GeoDataFrame of the bins
    # - geo_attributes_layer: GeoDataFrame of the attributes
    # - attr: attribute of the geo_attributes_layer that should be attributed to geo_bins
    # - bin_id: id of the bin column in geo_bins

    # Returns geo_bins with a new column called attr

    if not isinstance(geo_bins, gpd.GeoDataFrame):
        geo_bins = gpd.read_file(geo_bins)

    if not isinstance(geo_attributes_layer, gpd.GeoDataFrame):
        geo_attributes_layer = gpd.read_file(geo_attributes_layer)

    UTM_CRS = geo_bins.estimate_utm_crs()
    geo_bins = geo_bins.to_crs(UTM_CRS)
    geo_attributes_layer = geo_attributes_layer.to_crs(UTM_CRS)

    # Compute intersection of cells with soil types
    gdf_overlay = gpd.overlay(
        geo_bins, geo_attributes_layer, how="intersection", keep_geom_type=False
    )

    gdf_overlay["area"] = gdf_overlay.geometry.area

    # Add the area of each soil group within the cells
    grouped_data = gdf_overlay.groupby([bin_id, attr])["area"].sum().reset_index()

    # Compute for each cell what is the soil type that has the maximum area
    majority_data = grouped_data.loc[
        grouped_data.groupby(bin_id)["area"].idxmax(), [bin_id, attr]
    ]

    # Merge data with majority soil type for each cell
    geo_bins = pd.merge(
        geo_bins, majority_data, on=bin_id, how="left", suffixes=("_del", None)
    )
    # cell_data_section = cell_data_section.drop(columns=['Soil_ID_del'])

    return geo_bins


def assign_soil_to_annagnps_cells(
    cell_data_section,
    cells_geometry,
    soil_data,
    attr="mukey",
    outpath_cell_data_section=None,
    write_csv=False,
):
    # - cell_data_section: GeoDataFrame or path to csv AnnAGNPS_Cell_Data_Section.csv
    # - cell_geometry: GeoDataFrame or path to shapefile AnnAGNPS_Cell_IDs.shp
    # - soil_data: GeoDataFrame of path to SSURGO shapefile
    # - attr : attribute of the soil_data that should be attributed to
    # - outpath_cell_data_section: path to modified AnnAGNPS_Cell_Data_Section.csv (if different from cell_data_section)

    if not isinstance(cell_data_section, pd.DataFrame):
        cell_data_section = pd.read_csv(cell_data_section)
        outpath_cell_data_section = cell_data_section
    elif outpath_cell_data_section is None:
        outpath_cell_data_section = "AnnAGNPS_Cell_Data_Section.csv"  # Default name and will write the file in the current directory

    if not isinstance(cells_geometry, gpd.GeoDataFrame):
        cells_geometry = gpd.read_file(cells_geometry)

    if not isinstance(soil_data, gpd.GeoDataFrame):
        soil_data = gpd.read_file(soil_data)  # SSURGO data

    UTM_CRS = cells_geometry.estimate_utm_crs()
    cells_geometry = cells_geometry.to_crs(UTM_CRS)
    soil_data = soil_data.to_crs(UTM_CRS)

    # Compute intersection of cells with soil types
    gdf_overlay = gpd.overlay(
        cells_geometry, soil_data, how="intersection", keep_geom_type=False
    )

    gdf_overlay["area"] = gdf_overlay.geometry.area

    # Add the area of each soil group within the cells
    grouped_data = gdf_overlay.groupby(["DN", "musym"])["area"].sum().reset_index()

    # Compute for each cell what is the soil type that has the maximum area
    majority_data = grouped_data.loc[
        grouped_data.groupby("DN")["area"].idxmax(), ["DN", "musym"]
    ]

    # Rename columns
    majority_data = majority_data.rename(columns={"DN": "Cell_ID", "musym": "Soil_ID"})

    # Store original columns order
    columns = cell_data_section.columns

    # Merge data with majority soil type for each cell
    cell_data_section = pd.merge(
        cell_data_section,
        majority_data,
        on="Cell_ID",
        how="left",
        suffixes=("_del", None),
    )
    cell_data_section = cell_data_section.drop(columns=["Soil_ID_del"])

    # Reorder columns according to original column order
    cell_data_section = cell_data_section[columns]

    if write_csv:
        cell_data_section.to_csv(outpath_cell_data_section, sep=",", index=False)

    return cell_data_section


def assign_soil_to_annagnps_cells_OLD(
    cell_data_section,
    cells_geometry,
    soil_data,
    outpath_cell_data_section=None,
    write_csv=False,
):
    # - cell_data_section: GeoDataFrame or path to csv AnnAGNPS_Cell_Data_Section.csv
    # - cell_geometry: GeoDataFrame or path to shapefile AnnAGNPS_Cell_IDs.shp
    # - soil_data: GeoDataFrame of path to SSURGO shapefile
    # - outpath_cell_data_section: path to modified AnnAGNPS_Cell_Data_Section.csv (if different from cell_data_section)

    if ~isinstance(cell_data_section, pd.DataFrame):
        cell_data_section = pd.read_csv(cell_data_section)
        outpath_cell_data_section = cell_data_section
    elif outpath_cell_data_section is None:
        outpath_cell_data_section = "AnnAGNPS_Cell_Data_Section.csv"  # Default name and will write the file in the current directory

    if ~isinstance(cells_geometry, gpd.GeoDataFrame):
        cells_geometry = gpd.read_file(cells_geometry)

    if ~isinstance(soil_data, gpd.GeoDataFrame):
        soil_data = gpd.read_file(soil_data)  # SSURGO data

    utm_crs = cells_geometry.estimate_utm_crs()

    soil_data = soil_data.to_crs(utm_crs)

    for _, cell in cells_geometry.iterrows():
        gdf_cell = gpd.GeoDataFrame(
            pd.DataFrame(data={"DN": [cell[0]], "geometry": [cell[1]]}),
            crs=cells_geometry.crs,
        )
        gdf_cell = gdf_cell.to_crs(utm_crs)

        intersection = soil_data.overlay(
            gdf_cell, how="intersection", keep_geom_type=False
        )
        intersection["area"] = intersection.geometry.area
        intersection = pd.pivot_table(
            intersection, index=["musym"], values=["area"], aggfunc="sum"
        )  # !!! TEST THIS (to make sure if there are multiple intersections with the same soil group it's counted as one)
        musymmajority = intersection["area"].idxmax()
        # del intersection

        # musymmajority = intersection.loc[intersection['area'].idxmax(),'musym']
        cell_data_section.loc[
            cell_data_section["Cell_ID"] == cell[0], "Soil_ID"
        ] = musymmajority

    if write_csv:
        cell_data_section.to_csv(outpath_cell_data_section, sep=",", index=False)

    return cell_data_section


def run_one_query(county_code):
    url = "https://SDMDataAccess.sc.egov.usda.gov/Tabular/SDMTabularService/post.rest"

    # Updated query after Ron's Email 2022-11-21
    sQuery = f"""select 
            sa.saverest, 
            l.areasymbol, 
            l.areaname,
            mu.musym, 
            hydgrp,
            kwfact,
            albedodry_r,
            (SELECT CASE when min(resdept_r) is null then '>200' else cast(min(resdept_r) as 
            varchar) END
                        from component left outer join corestrictions on component.cokey = 
            corestrictions.cokey where component.cokey = c.cokey and reskind is not null) as 
            restrictiondepthr,
            partdensity,
            c.compname, 
            texdesc, 
            hzdepb_r, 
            dbovendry_r, 
            claytotal_r,
            silttotal_r,
            sandtotal_r,
            (select sum(cf.fragvol_r) as fragvol  FROM chfrags cf WHERE cf.chkey = ch.chkey 
            ) as fragvol,
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
            LEFT OUTER JOIN component c ON c.mukey = mu.mukey and c.cokey = (SELECT TOP 1 component.cokey FROM component WHERE 
            component.mukey=mu.mukey ORDER BY component.comppct_r DESC)
            LEFT OUTER JOIN chorizon ch ON ch.cokey = c.cokey 
            LEFT OUTER JOIN chtexturegrp ct ON ch.chkey=ct.chkey 

            WHERE l.areasymbol = {county_code} and ct.rvindicator = 'yes'

            Order by l.areasymbol, musym, compname, hzdepb_r"""

    # sQuery = f''' select
    # sa.saverest,
    # l.areasymbol,
    # l.areaname,
    # mu.musym,
    # hydgrp,
    # kwfact,
    # albedodry_r,
    # (SELECT CASE when min(resdept_r) is null then '>200' else cast(min(resdept_r) as varchar) END

    # from component left outer join corestrictions on component.cokey = corestrictions.cokey where component.cokey = c.cokey and reskind is not null) as restrictiondepthr,
    # partdensity,
    # c.compname,
    # texdesc,
    # hzdepb_r,
    # dbovendry_r,
    # claytotal_r,
    # silttotal_r,
    # sandtotal_r,
    # (select sum(cf.fragvol_r) as fragvol FROM chfrags cf WHERE cf.chkey = ch.chkey ) as fragvol,
    # sandvf_r,
    # caco3_r,
    #     ksat_r,
    #     wthirdbar_r,
    #     wfifteenbar_r,
    #     om_r,
    #     ph1to1h2o_r,
    #     c.comppct_r
    #     FROM
    #     legend l INNER JOIN mapunit mu ON mu.lkey = l.lkey
    #     LEFT OUTER JOIN sacatalog sa ON sa.areasymbol = l.areasymbol
    #     LEFT OUTER JOIN component c ON c.mukey = mu.mukey
    #     LEFT OUTER JOIN chorizon ch ON ch.cokey = c.cokey
    #     LEFT OUTER JOIN chtexturegrp ct ON ch.chkey=ct.chkey
    #     WHERE l.areasymbol = '{county_code}' and ct.rvindicator = 'yes'
    #     Order by l.areasymbol, musym, hzdepb_r'''

    dRequest = dict()
    dRequest["FORMAT"] = "JSON+COLUMNNAME"
    dRequest["QUERY"] = sQuery

    resp = requests.post(url, data=dRequest).json()

    df = pd.DataFrame(resp["Table"])
    df.columns = resp["Table"][0]
    df = df.drop(labels=[0], axis=0)

    return df


def run_batch_write_files(county_codes, outpath=""):
    for county_code in county_codes:
        df = run_one_query(county_code)
        outfile = Path(outpath, f"{county_code}_nasis.csv")
        df.to_csv(outfile, index=False, sep=",")


def run_nita(filefolder, path_to_NITA_exe, combine_list=None, units_out=None):
    list_of_files = glob.glob(str(Path(filefolder, "*_nasis.csv")))

    # The files are not combined by default
    if combine_list is None:
        combine_list = [0 for _ in list_of_files]

    # SI Units by default
    if units_out is None:
        units_out = [1 for _ in list_of_files]

    with open(Path(filefolder, "NITA_CONTROL.csv"), "w") as control:
        control.write("FILENAME,UNITS_OUT,COMBINE\n")
        for i, file in enumerate(list_of_files):
            filename = Path(file).name
            control.write(f"{filename},{units_out[i]},{combine_list[i]}\n")

    command = path_to_NITA_exe + " /f:" + "NITA_CONTROL.csv"

    os.chdir(filefolder)
    subprocess.call(command)
