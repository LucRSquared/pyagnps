from pathlib import Path
import json
import pandas as pd
import geopandas as gpd
from tqdm import tqdm
from sqlalchemy import URL, create_engine, text as sql_text

credentials = Path('../../inputs/db_credentials.json')
with open(credentials, 'r') as f:
    credentials = json.load(f)

user     = credentials['user']
password = credentials['password']
host     = credentials['host']
port     = credentials['port']
database = credentials['database']

url_object = URL.create(
    "postgresql",
    username=user,
    password=password,
    host=host,
    port=port,
    database=database
)

# create a SQLAlchemy engine object
engine = create_engine(url_object)

path_to_thucs = Path(
    "/aims-nas/luc/data/tophuc_S_M_40000_closed_holes_with_container_thuc_merged_bbox_area_first_kept.gpkg"
)

print('Reading THUCs...')
thucs = gpd.read_file(
    path_to_thucs
)  # GeoDataFrame containing the thucs and their geometry
thucs = thucs.sort_values(by=["bbox_area_sqkm"], ascending=False)

runlist = thucs["tophucid"].to_list()
# runlist = ['1311']


# Get list of valid Soil_ID once 
valid_soil_ids = []
with engine.connect() as connection:
    result = connection.execute(sql_text("""SELECT "Soil_ID" FROM usa_valid_soil_data"""))
    valid_soil_ids = [row[0] for row in result]

for _, tuc in tqdm(thucs.iterrows(), total=thucs.shape[0]):
    thuc_id = tuc["tophucid"]

    if thuc_id not in runlist:
        continue
   
    with engine.connect() as connection:
        try:
            # Check if soil_id is in prefetched list 
            query = f"""
            ALTER TABLE thuc_{thuc_id}_annagnps_cell_data_section
            ADD COLUMN soil_id_annagnps_valid INT;

            UPDATE thuc_{thuc_id}_annagnps_cell_data_section 
            SET soil_id_annagnps_valid = CASE WHEN soil_id IN {tuple(valid_soil_ids)} THEN 1 ELSE 0 END
            """
            
            connection.execute(sql_text(query))
            connection.commit()
        except Exception as e:
            print(f'Error for THUC {thuc_id}') 
            # print(e)
            connection.rollback()
