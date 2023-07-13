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

outpath_cells = Path('../../outputs/soil_data_market/db_population/cells_unknown_soil_gNATSGO.gpkg')


thucs = gpd.read_file(
    path_to_thucs
)  
# GeoDataFrame containing the thucs and their geometry
thucs = thucs.sort_values(by=["bbox_area_sqkm"], ascending=True)

runlist = thucs['tophucid'].to_list()
# runlist = ['1957', '1958']

# Fetch the list of valid Soil_ID values from usa_valid_soil_data
soil_ids_query = """SELECT DISTINCT "Soil_ID" FROM usa_valid_soil_data"""
soil_ids_df = pd.read_sql(soil_ids_query, engine)
soil_ids = soil_ids_df['Soil_ID'].tolist()


rows_invalid_soil_id = []

for thuc_id in tqdm(runlist):
    
    # Query cells that have soil_id populated
    query_cds = f"SELECT cell_id, soil_id FROM thuc_{thuc_id}_annagnps_cell_data_section WHERE soil_id IS NOT NULL"

    try:
        with engine.connect() as conn:
            # query all populated cells
            cell_data_section = pd.read_sql(sql=sql_text(query_cds), con=conn)
            # filter out the cells that have a soil_id in in usa_valid_soil_data
            cell_data_section = cell_data_section[~cell_data_section['soil_id'].isin(soil_ids)]

            populated_invalid_cells = cell_data_section['cell_id'].to_list()
            populated_invalid_soils = cell_data_section['soil_id'].to_list()

            if populated_invalid_cells:

                query_geom = "SELECT dn, geom FROM thuc_{}_annagnps_cell_ids WHERE dn in ({})"\
                    .format(thuc_id, ','.join(str(x) for x in populated_invalid_cells))


                cells = gpd.read_postgis(sql=sql_text(query_geom), con=conn, geom_col="geom")
                cells = cells.to_crs('epsg:4326')
                cells['thuc'] = thuc_id

                del cells['fid']

                invalid_rows = cells.merge(cell_data_section, left_on='dn', right_on='cell_id')
                invalid_rows = invalid_rows[['geom', 'thuc', 'cell_id', 'soil_id']]

                rows_invalid_soil_id.append(invalid_rows)

    except Exception as e:
        print(e)
        continue

invalid_cells = gpd.GeoDataFrame(pd.concat(rows_invalid_soil_id, ignore_index=True), crs='epsg:4326', geometry="geom")

invalid_cells.to_file(outpath_cells, driver='GPKG')
