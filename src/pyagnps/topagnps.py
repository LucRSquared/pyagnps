"""
TopAGNPS related functions
"""
# from osgeo import gdal
# from importlib.resources import path
import os
import subprocess
from lxml import etree as ET
import geopandas as gpd
import numpy as np
import rasterio
import rioxarray
import py3dep
from osgeo import gdal
import pyagnps.utils as utils

# import logging


def create_topagnps_directory(root_dir, topagnps_dir):
    """
    Create the directory for the topagnps simulation
    """
    path_to_dir = root_dir + "/" + topagnps_dir

    if not os.path.exists(path_to_dir):
        os.mkdir(path_to_dir)

    return path_to_dir


def create_topagnps_xml_control_file(dico, savepath):
    """
    Create the xml control file for TopAGNPS
    """
    root = ET.Element("TOPAGNPS")

    for key, value in dico.items():
        line = ET.SubElement(root, "KEYWORD")
        line.text = str(value)
        line.set("ID", str(key))

    tree = ET.ElementTree(root)
    tree.write(savepath, pretty_print=True)


def read_topagnps_xml_control_file(path_to_xml):
    """
    Read the xml control file for TopAGNPS
    """
    tree = ET.parse(path_to_xml)
    root = tree.getroot()

    dico = {}
    for child in root:
        dico[child.attrib["ID"]] = child.text

    return dico


# def run_topagnps(path_to_xml_dir, path_to_bin):
#     """
#     Call TopAGNPS binary and run by passing it the XML control file
#     """
#     previousdir = os.getcwd()
#     os.chdir(path_to_xml_dir)

#     try:
#         with open('command_line_output.txt', 'w') as clo:
#             print("Running TopAGNPS")
#             print(path_to_xml_dir)
#             print(f"{path_to_bin} /f:TOPAGNPS.XML")
#             subprocess.call([path_to_bin, " /f:TOPAGNPS.XML"], stdout=clo)
#     except:
#         os.chdir(previousdir)
#     os.chdir(previousdir)


def run_topagnps(
    path_to_xml_dir,
    path_to_bin,
    memtrack=False,
    output_memtrack_filename="valgrind_output.txt",
):
    """
    Call TopAGNPS binary and run by passing it the XML control file
    """

    previousdir = os.getcwd()
    os.chdir(path_to_xml_dir)

    try:
        with open("command_line_output.txt", "w") as clo:
            if memtrack:
                subprocess.call(
                    [
                        "valgrind",
                        "--tool=massif",
                        "--stacks=yes",
                        f"--massif-out-file={output_memtrack_filename}",
                        path_to_bin,
                        " /f:TOPAGNPS.XML",
                    ],
                    stdout=clo,
                )
            else:
                subprocess.call([path_to_bin, " /f:TOPAGNPS.XML"], stdout=clo)

    except:
        os.chdir(previousdir)
    os.chdir(previousdir)


def download_dem(
    gdf, path_to_dir, name="dem_thuc", resolution_m=10, buffer_m=500, keeptif=False
):
    """This function downloads the DEM from USGS PY3DEP service directly"""
    # gdf : GeoDataFrame of the area (containing one polygon)
    # path_to_dir : path to the directory containing the topagnps simulations

    if not (path_to_dir.endswith("/")):
        path_to_dir = path_to_dir + "/"

    if isinstance(gdf, str):
        gdf = gpd.read_file(gdf)

    if gdf.shape[0] > 1:
        raise Exception("Please provide a GeoDataFrame containing only one line")

    path_to_asc = f"{path_to_dir}{name}_res_{resolution_m}_m.asc"
    path_to_tif = f"{path_to_dir}{name}_res_{resolution_m}_m.tif"

    output_crs = gdf.estimate_utm_crs()
    gdf_buffered = gdf.to_crs(output_crs).buffer(buffer_m).to_crs("epsg:4326")

    bbox = gdf_buffered.bounds

    dem = py3dep.get_map(
        "DEM",
        (
            bbox["minx"].values[0],
            bbox["miny"].values[0],
            bbox["maxx"].values[0],
            bbox["maxy"].values[0],
        ),
        resolution=resolution_m,
        geo_crs="epsg:4326",
        crs="epsg:4326",
    )

    dem.name = "dem"
    dem.attrs["units"] = "meters"
    dem = dem.rio.reproject(output_crs)

    dem.rio.to_raster(
        path_to_tif
    )  # Write first to tif file because rasterio cannot directly write to

    opts = gdal.WarpOptions(
        format="AAIGrid",
        srcNodata="nan",
        dstNodata=-999,
        xRes=resolution_m,
        yRes=resolution_m,
        resampleAlg=gdal.GRA_NearestNeighbour,
    )

    with open(path_to_asc, "w"):
        gdal.Warp(path_to_asc, path_to_tif, options=opts)

    if os.path.exists(path_to_tif) and not (keeptif):  # delete temporary
        os.remove(path_to_tif)
    else:
        print(f"The file {path_to_tif} does not exist")

    return dem, path_to_asc


def get_dem_from_raster(
    gdf,
    path_to_input_raster,
    path_to_dir,
    name="dem_thuc",
    resolution_m=10,
    buffer_m=500,
    keeptif=False,
):
    """This function retrieves the DEM from a provided raster that contains the shape specified by a GeoDataFrame"""
    # gdf : GeoDataFrame of the area (containing one polygon)
    # path_to_dir : path to the directory containing the topagnps simulations

    if not (path_to_dir.endswith("/")):
        path_to_dir = path_to_dir + "/"

    if isinstance(gdf, str):
        gdf = gpd.read_file(gdf)

    if gdf.shape[0] > 1:
        raise Exception("Please provide a GeoDataFrame containing only one line")

    path_to_asc = f"{path_to_dir}{name}_res_{resolution_m}_m.asc"
    path_to_tif = f"{path_to_dir}{name}_res_{resolution_m}_m.tif"

    output_crs = gdf.estimate_utm_crs()

    with rioxarray.open_rasterio(path_to_input_raster, masked=True) as raster_dem:
        raster_crs = raster_dem.rio.crs
        encoded_nodata = raster_dem.rio.encoded_nodata

        gdf_buffered = gdf.to_crs(output_crs).buffer(buffer_m).to_crs(raster_crs)

        bbox = gdf_buffered.bounds

        dem = raster_dem.rio.clip_box(
            minx=bbox["minx"].values[0],
            miny=bbox["miny"].values[0],
            maxx=bbox["maxx"].values[0],
            maxy=bbox["maxy"].values[0],
        )

        # nodata = dem._FillValue
        # dem = dem.rio.write_nodata(nodata, encoded=True)
        # dem = dem.where(dem != dem.rio.nodata)
        # dem = dem.rio.write_nodata(np.nan)

    if "hydrodem" in path_to_input_raster:
        # Data was stored in centimeters
        dem = dem / 100.0

    dem = dem.rio.write_nodata(encoded_nodata)

    # Check if nodata is gone or not

    dem.name = "dem"
    dem.attrs["units"] = "meters"

    # dem = dem.rio.reproject(output_crs)

    dem.rio.to_raster(
        path_to_tif
    )  # Write first to tif file because rasterio cannot directly write to

    with rioxarray.open_rasterio(path_to_tif, masked=True) as tmp_tif:
        nodataval = tmp_tif.rio.nodata

    opts = gdal.WarpOptions(
        format="AAIGrid",
        srcNodata=nodataval,
        dstNodata=-999,
        xRes=resolution_m,
        yRes=resolution_m,
        resampleAlg=gdal.GRA_NearestNeighbour,
    )

    with open(path_to_asc, "w"):
        gdal.Warp(path_to_asc, path_to_tif, options=opts)

    if os.path.exists(path_to_tif) and not (keeptif):  # delete temporary
        os.remove(path_to_tif)

    return dem, path_to_asc


def process_hydrodem_for_topagnps(
    path_to_raster_dir,
    path_to_write_dir=None,
    wall_threshold_m=4800,
    unwall=True,
    unfill=False,
    applyshift=True,
):
    """
    INPUTS:
    path_to_raster_dir : this folder needs to contain the files hydrodem.tif, elev_cm.tif, and filldepth.tif
    path_to_write_dir  : default is path_to_raster_dir, if specified: directory to write processed dem from
    wall_threshold_m   : if elevation is greater than wall_threshold_m then replace with elevation from elev_cm.tif
    unfill             : True = remove the filling applied
    unwall             : True = removes walls
    applyshift         : computes and apply a shift equal to abs(min(hydrodem)) + 1

    OUTPUTS:
    hydrodem[_unfilled][_shifted][_unwalled]_topagnps.tif in METERS

    """

    path_to_hydrodem = f"{path_to_raster_dir}/hydrodem.tif"
    path_to_elev_cm = f"{path_to_raster_dir}/elev_cm.tif"
    path_to_filldepth = f"{path_to_raster_dir}/filldepth.tif"

    if path_to_write_dir is None:
        path_to_write_dir = path_to_raster_dir

    output_raster = f"{path_to_write_dir}/hydrodem"

    if unfill:
        output_raster = f"{output_raster}_unfilled"

    if applyshift:
        output_raster = f"{output_raster}_shifted"

    if unwall:
        output_raster = f"{output_raster}_unwalled_{int(wall_threshold_m)}_m"

    output_raster = f"{output_raster}_topagnps.tif"

    with rioxarray.open_rasterio(path_to_hydrodem, masked=True) as hydrodem:
        encoded_nodata = hydrodem.rio.encoded_nodata

        if unfill:
            print("removing depression filling")
            with rioxarray.open_rasterio(path_to_filldepth, masked=True) as filldepth:
                hydrodem = hydrodem - filldepth

        if unwall:
            print("knocking down walls!")
            with rioxarray.open_rasterio(path_to_elev_cm, masked=True) as elevcm:
                mask = hydrodem - elevcm < int(wall_threshold_m * 100)

                hydrodem = hydrodem.where(mask, elevcm)

        if applyshift:
            print("applying shift")
            try:
                shift = abs(hydrodem.STATISTICS_MINIMUM) + 100
                shift = int(shift.values)
            except:
                shift = int(abs(hydrodem.min()) + 100)

            hydrodem = hydrodem + shift

            output_raster = output_raster.replace("_shifted", f"_shifted_{shift}_cm")

        else:
            shift = int(0)

        # hydrodem = hydrodem/100.0
        # if hydrodem.rio.nodata is None:
        #     print('Replacing lost nodata')
        hydrodem = hydrodem.rio.write_nodata(encoded_nodata)

        print("writing to raster")
        hydrodem.rio.to_raster(
            output_raster,
            num_threads="all_cpus",
            dtype="int32",
            driver="GTiff",
            windowed=True,
        )

    return hydrodem, shift


def find_outlet_uparea_shape_intersection(path_to_uparea_asc, boundary):
    # if boundary is string: read from file
    # if boundary is GeoDataFrame or GeoSeries use it as such
    # return (x, y, row, col)

    if isinstance(boundary, str):
        shape = gpd.read_file(boundary)
        shape = shape.geometry
    elif isinstance(boundary, gpd.GeoDataFrame):
        shape = boundary.geometry
    elif isinstance(boundary, gpd.GeoSeries):
        shape = boundary

    with rasterio.open(path_to_uparea_asc) as src:
        shape = shape.to_crs(src.crs)

        # Transform shape into raster of the same size as the uparea with True=inside
        rasterized_boundary = rasterio.features.geometry_mask(
            geometries=shape,
            out_shape=(src.height, src.width),
            transform=src.transform,
            all_touched=True,
            invert=True,
        )
        # Extract the values inside the mask
        upareas = np.extract(rasterized_boundary, src.read())

        # Find maximum value inside the mask
        maxuparea = max(upareas)

        # Find the indices inside the mask where the raster reaches its maximum value
        idxmaxuparea = np.argwhere(
            maxuparea == src.read()
        )  # Indices where uparea inside the boundary is equal to the maximum value

        numoutlets = idxmaxuparea.shape[0]
        if numoutlets > 1:
            raise Exception(
                f"Ambiguous outlet location, found {numoutlets} eligible candidates! Review hyrological boundary and/or UPAREA file"
            )
        elif numoutlets == 0:
            raise Exception(
                "This is weird, no outlets found ‽, this shouldn"
                "t happen, if you"
                "re seeing this error message you"
                "re in trouble"
            )
        elif numoutlets == 1:  # The normal case that should happen normally
            # Get row col indices
            rowout, colout = (
                idxmaxuparea[0, 1] + 1,
                idxmaxuparea[0, 2] + 1,
            )  # +1 because TopAGNPS indices start at 1
            # Get x, y
            xout, yout = src.xy(rowout - 1, colout - 1)

            return (xout, yout, rowout, colout, src.crs)


def check_topagnps_wrn_log(path_to_topagnps_wrn_log):
    """
    Returns a list of DEM directions where the cells touch the edge of the map
    """

    if "TopAGNPS_wrn.CSV" not in path_to_topagnps_wrn_log:
        raise Exception(f"Wrong log file: {path_to_topagnps_wrn_log}")

    with open(path_to_topagnps_wrn_log, "r") as topwarn:
        line_touch_edges = [
            line for line in topwarn if "The watershed boundary touches" in line
        ]

    touch_edges = []
    directions = ["North", "West", "South", "East"]

    for line in line_touch_edges:
        for dir in directions:
            if dir in line:
                touch_edges.append(dir)

    return touch_edges


def quality_control_areas_vs_boundary(cells, boundary_shape):
    if isinstance(cells, str):
        if "AnnAGNPS_Cell_IDs.asc" in cells:
            cells = utils.polygonize_cell_reach_IDs_asc(
                cells, outepsg=4326, return_gdf=True, writefile=False
            )
        else:
            # try to read as a GeoDataFrame object
            cells = gpd.read_file(cells)
    else:
        raise Exception(
            "The TopAGNPS cells should be provided as a GeoDataFrame object or a path to a shape object or AnnAGNPS_Cell_IDs.asc raster file"
        )

    if isinstance(boundary_shape, str):
        boundary_shape = gpd.read_file(boundary_shape)
    elif not (isinstance(boundary_shape, gpd.GeoDataFrame)):
        raise Exception(
            "The boundary shape should be either the path to a shape file or a GeoDataFrame object directly"
        )

    utm = cells.estimate_utm_crs()

    cells_bounds = (
        cells.to_crs(utm).buffer(0.01).unary_union
    )  # Buffer by 1cm to prevent floating point errors

    boundary_shape = boundary_shape.to_crs(
        utm
    )  # Make sure the boundary and the cells are on the same utm

    boundary_shape_bounds = boundary_shape.geometry.iloc[0]

    boundary_minus_cells = boundary_shape_bounds.difference(
        cells_bounds
    )  # Area inside boundary not covered by the cells
    cells_minus_boundary = cells_bounds.difference(
        boundary_shape_bounds
    )  # Cells area that are not in the original boundary
    intersection = boundary_shape_bounds.intersection(
        cells_bounds
    )  # Area where they overlap

    fraction_of_boundary_covered = intersection.area / boundary_shape_bounds.area
    fraction_of_boundary_missed = 1 - fraction_of_boundary_covered

    results = {}

    results["total_cells_area_sqm"] = cells_bounds.area
    results["total_boundary_area_sqm"] = boundary_shape_bounds.area

    results["cells_missed_area_sqm"] = boundary_minus_cells.area
    results["beyond_boundary_cells_area_sqm"] = cells_minus_boundary.area

    results["intersection_area_sqm"] = intersection.area
    results["fraction_boundary_covered"] = fraction_of_boundary_covered
    results["fraction_boundary_missed"] = fraction_of_boundary_missed

    return results
