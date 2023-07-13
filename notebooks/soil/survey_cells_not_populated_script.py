import json
from pathlib import Path

from tqdm import tqdm

import geopandas as gpd
import pandas as pd

from sqlalchemy import create_engine, text as sql_text

# DATABASE SETUP
credentials = Path("../../inputs/db_credentials.json")
with open(credentials, "r") as f:
    credentials = json.load(f)

user = credentials["user"]
password = credentials["password"]
host = credentials["host"]
port = credentials["port"]
database = credentials["database"]

engine = create_engine(f"postgresql://{user}:{password}@{host}:{port}/{database}")

# path_to_thucs = Path(
#     "D:/AIMS/Datasets/THUCS_TopAGNPS_Delineations/40k_SM/tophuc_S_M_40000_closed_holes_with_container_thuc_merged_bbox_area_first_kept.gpkg"
# )

path_to_thucs = Path(
    "/aims-nas/luc/data/tophuc_S_M_40000_closed_holes_with_container_thuc_merged_bbox_area_first_kept.gpkg"
)

# path_to_thucs = Path('../../inputs/thucs/tophuc_S_M_40000_closed_holes_with_container_thuc_merged_bbox_area_first_kept.gpkg')

outpath_cells = Path('../../outputs/soil_data_market/db_population/cells_without_soil_ssurgo.gpkg')


thucs = gpd.read_file(
    path_to_thucs
)  
# GeoDataFrame containing the thucs and their geometry
thucs = thucs.sort_values(by=["bbox_area_sqkm"], ascending=True)

runlist = thucs['tophucid'].to_list()
# runlist = ['1957', '1958']

cells_without_soil = []

for thuc_id in tqdm(runlist):
    
    # Query for cells where soil is not populated
    query_cds = f"SELECT cell_id FROM thuc_{thuc_id}_annagnps_cell_data_section WHERE soil_id is NULL"

    try:
        with engine.connect() as conn:
            cell_data_section = pd.read_sql(sql=sql_text(query_cds), con=conn)

            unpopulated_cells = cell_data_section['cell_id'].to_list()

            if unpopulated_cells:

                query_geom = "SELECT * FROM thuc_{}_annagnps_cell_ids WHERE dn in ({})".format(thuc_id, ','.join(str(x) for x in unpopulated_cells))


                cells = gpd.read_postgis(sql=sql_text(query_geom), con=conn, geom_col="geom")
                cells = cells.to_crs('epsg:4326')
                cells['thuc'] = thuc_id

                del cells['fid']

                cells_without_soil.append(cells)

    except Exception as e:
        print(e)
        continue

empty_cells = gpd.GeoDataFrame(pd.concat(cells_without_soil, ignore_index=True), crs='epsg:4326', geometry="geom")

empty_cells.to_file(outpath_cells, driver='GPKG')
