import rasterio
from rasterio import features
from shapely.geometry import shape
from shapely.geometry.polygon import Polygon, LinearRing
import geopandas as gpd

def polygonize_cell_reach_IDs_asc(path_to_cell_IDs_asc, path_to_cell_IDs_shp, outepsg=4326, return_gdf=False):

    with rasterio.open(path_to_cell_IDs_asc) as src:
        band = src.read(1).astype('int32')
        nodata = src.nodata

    mask = band != nodata

    shapes = features.shapes(source=band, mask=mask, connectivity=8, transform=src.transform)

    data = {'DN':[], 'geometry':[]}
    for geom, value in shapes:
        data['DN'].append(int(value))
        data['geometry'].append(shape(geom))

    gdf = gpd.GeoDataFrame(data, crs=src.crs).to_crs(epsg=outepsg)

    gdf.to_file(path_to_cell_IDs_shp, driver='ESRI Shapefile')

    if return_gdf:
        return gdf
    else:
        return

def get_geometry_boundary(geom):
    if isinstance(geom, LinearRing):
        return geom
    elif isinstance(geom, Polygon):
        return geom.exterior
    elif isinstance(geom, gpd.geoseries.GeoSeries):
        if geom.shape[0] == 1:
            return geom.values[0].exterior
        else:
            raise Exception('Geoseries has more than one object')
    else:
        raise Exception('Invalid geometry object')