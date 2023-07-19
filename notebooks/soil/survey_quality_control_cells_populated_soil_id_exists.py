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

# outpath_cells = Path('../../outputs/soil_data_market/db_population/cells_unknown_soil_gNATSGO.gpkg')
outpath_cells = Path('/aims-nas/luc/thuc_soil_survey/cells_NITA_unvalid_soils_reason_gNATSGO.gpkg')


print('Reading Original THUC file')

thucs = gpd.read_file(
    path_to_thucs
)  
# GeoDataFrame containing the thucs and their geometry
thucs = thucs.sort_values(by=["bbox_area_sqkm"], ascending=False)

runlist = thucs['tophucid'].to_list()

# soil_root_path =  Path('D:/AIMS/Datasets/Soil/DATABASE_POPULATION_TASKS/SDM_QUERY_AND_NITA_PROCESSING/ALL_US_v2_SSURGO_STATSGO2_RSS/')
soil_root_path = Path('/aims-nas/data/datasets/SDM_QUERY_AND_NITA_PROCESSING/ALL_US_v2_SSURGO_STATSGO2_RSS/')

# Fetch the list of valid Soil_ID values from usa_valid_soil_data
# soil_ids_query = """SELECT DISTINCT "Soil_ID" FROM usa_valid_soil_data"""
# soil_ids_df = pd.read_sql(soil_ids_query, engine)
# soil_ids = soil_ids_df['Soil_ID'].tolist()

# Get soil_ids directly from NITA processed files:
path_to_NITA_processed_soil_data = soil_root_path / 'all_valid_soil_data.parquet'
df_soil_data = pd.read_parquet(path_to_NITA_processed_soil_data)
valid_soil_ids = df_soil_data['Soil_ID'].tolist()

# Get list of discarded records by NITA
path_to_NITA_discarded =  soil_root_path / 'nita_excluded_soil_ids_mukey.parquet'
df_nita_discarded = pd.read_parquet(path_to_NITA_discarded)
df_nita_discarded = df_nita_discarded.drop_duplicates(subset='soil_id', keep='first')

# Get all other records that were not queried in the first place
path_to_3soil_sources_rvindicator_yes = soil_root_path / 'raw_query_data' / 'soil_data_ALL_SOURCES_COMBINED_with_mukey.csv'
path_to_3soil_sources_rvindicator_no_and_NULL = soil_root_path / 'raw_query_data' / 'ALL_SOURCES_rvindicator_no_and_NULL.csv'

df_raw_all = pd.concat([pd.read_csv(path_to_3soil_sources_rvindicator_yes),
                        pd.read_csv(path_to_3soil_sources_rvindicator_no_and_NULL)])

# Generate df with soil_id and compname
df_soil_id_compname = df_raw_all[['mukey', 'compname']].copy(deep=True)
df_soil_id_compname.rename(columns={'mukey': 'soil_id'}, inplace=True)
df_soil_id_compname['soil_id'] = df_soil_id_compname['soil_id'].astype(str)

rows_invalid_soil_id = []

for thuc_id in tqdm(runlist):
    
    # Query cells that have soil_id populated
    query_cds = f"SELECT cell_id, soil_id FROM thuc_{thuc_id}_annagnps_cell_data_section WHERE soil_id IS NOT NULL"

    try:
        with engine.connect() as conn:
            # query all populated cells
            cell_data_section = pd.read_sql(sql=sql_text(query_cds), con=conn)
            # filter out the cells that have a soil_id in in usa_valid_soil_data
            cell_data_section = cell_data_section[~cell_data_section['soil_id'].isin(valid_soil_ids)] # the row that does the job

            populated_invalid_cells = cell_data_section['cell_id'].to_list()

            if populated_invalid_cells:

                query_geom = "SELECT dn, geom FROM thuc_{}_annagnps_cell_ids WHERE dn in ({})"\
                    .format(thuc_id, ','.join(str(x) for x in populated_invalid_cells))

                cells = gpd.read_postgis(sql=sql_text(query_geom), con=conn, geom_col="geom")
                cells = cells.to_crs('epsg:4326')
                cells['thuc'] = thuc_id

                invalid_rows = cells.merge(cell_data_section, left_on='dn', right_on='cell_id')
                invalid_rows = invalid_rows[['geom', 'thuc', 'cell_id', 'soil_id']]
                invalid_rows['soil_id'] = invalid_rows['soil_id'].astype(str)

                # Merge with discarded NITA files
                invalid_rows = invalid_rows.merge(df_nita_discarded, how='left', on='soil_id')

                # Add companame (soil description)
                invalid_rows = invalid_rows.merge(df_soil_id_compname, how='left', on='soil_id')

                # If problems value is NaN then it means it wasn't queried in the first place
                category = 'not in original query'
                if category in invalid_rows['problems'].cat.categories:
                    invalid_rows.loc[invalid_rows['problems'].isna(), 'problems'] = category
                else:
                    # Add the new category to the existing categories
                    invalid_rows['problems'] = invalid_rows['problems'].cat.add_categories(category)
                    invalid_rows.loc[invalid_rows['problems'].isna(), 'problems'] = category

                rows_invalid_soil_id.append(invalid_rows)
    except Exception as e:
        print(e)
        continue

print('Generating MEGA concatainated file of all problematic cells')

invalid_cells = gpd.GeoDataFrame(pd.concat(rows_invalid_soil_id, ignore_index=True), crs='epsg:4326', geometry="geom")
invalid_cells['problems'] = invalid_cells['problems'].astype(str)

print('Writing file to disk')
invalid_cells.to_file(outpath_cells, driver='GPKG')

print('Done!')