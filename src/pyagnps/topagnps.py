# from osgeo import gdal
from lxml import etree as ET
import geopandas as gpd
import numpy as np
import rasterio

def create_topagnps_dictionary(**kwargs):
    pass

def create_topagnps_xml_control_file(dico, savepath):
    pass

def run_topagnps(topbinpath, path_to_control_file):
    pass

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
        rasterized_boundary = rasterio.features.geometry_mask(geometries=shape, 
                                                              out_shape=(src.height, src.width),
                                                              transform=src.transform,
                                                              all_touched=True,
                                                              invert=True)
        # Extract the values inside the mask
        upareas = np.extract(rasterized_boundary, src.read())

        # Find maximum value inside the mask
        maxuparea = max(upareas)

        # Find the indices inside the mask where the raster reaches its maximum value
        idxmaxuparea = np.argwhere(maxuparea==src.read()) # Indices where uparea inside the boundary is equal to the maximum value

        numoutlets = idxmaxuparea.shape[0]
        if numoutlets > 1:
            raise Exception(f'Ambiguous outlet location, found {numoutlets} eligible candidates! Review hyrological boundary and/or UPAREA file')
        elif numoutlets == 0:
            raise Exception('This is weird, no outlets found ‽, this shouldn''t happen, if you''re seeing this error message you''re in trouble')
        elif numoutlets == 1: # The normal case that should happen normally
            # Get row col indices
            rowout, colout = idxmaxuparea[0,1]+1, idxmaxuparea[0,2]+1 # +1 because TopAGNPS indices start at 1
            # Get x, y 
            xout, yout = src.xy(rowout-1, colout-1)

            return (xout, yout, rowout, colout, src.crs)


