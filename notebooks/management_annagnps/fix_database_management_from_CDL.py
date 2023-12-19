import socket
import time
import json
from pathlib import Path

from tqdm import tqdm

import geopandas as gpd
import pandas as pd

from pyagnps import soil_data_market as sdm
from pyagnps.utils import log_to_file, get_current_time

from sqlalchemy import URL, create_engine, text as sql_text
from sqlalchemy.orm import sessionmaker

# DATABASE SETUP
credentials = Path("../../inputs/db_credentials.json")
with open(credentials, "r") as f:
    credentials = json.load(f)

user = credentials["user"]
password = credentials["password"]
host = credentials["host"]
port = credentials["port"]
database = credentials["database"]

url_object = URL.create(
    "postgresql",
    username=user,
    password=password,
    host=host,
    port=port,
    database=database,
)

# create a SQLAlchemy engine object
engine = create_engine(url_object)

# PATHS SETUP
bigbang = time.time()
nodename = socket.gethostname()

path_to_thucs = Path(
    "/aims-nas/luc/data/tophuc_S_M_40000_closed_holes_with_container_thuc_merged_bbox_area_first_kept.gpkg"
)

raster_CDL_path = Path("/aims-nas/data/datasets/Management/CDL_Annual/CDL_2022.tif")

track_files_dir = Path("/aims-nas/luc/thuc_field_ID_CDL_db_fix/")
path_to_management_class_names = Path(
    "/aims-nas/data/datasets/Management/CDL_Annual/definition_of_CDL_Class_Names_for_use_in_AIMS/CDL_All_Codes.csv"
)

log_dir = track_files_dir / "LOGS"

general_log = log_dir / f"{nodename}_batch_field_ID_CDL_population_general_log.txt"
fail_list = log_dir / f"{nodename}_fail_list.txt"

print("Reading and initializing files...")
thucs = gpd.read_file(
    path_to_thucs
)  # GeoDataFrame containing the thucs and their geometry
thucs = thucs.sort_values(by=["bbox_area_sqkm"], ascending=False)

df_cdl = pd.read_csv(path_to_management_class_names)
dico = df_cdl[['Value','Class_Name']].set_index('Value').to_dict(orient='dict')['Class_Name']
dico = {key: '' if isinstance(value,float) else value for key, value in dico.items()}

classes = set(dico.values())
# Remove classes that need 
already_valid_classes = classes - set(["", "Clouds_No_Data", "Developed"])
reversed_already_valid_dico = {dico[val]: val for val in dico if dico[val] in already_valid_classes}
invalid_field_ids_to_fix = ("'NaN'", "'Clouds_No_Data'", "'Developed'")

# runlist = thucs["tophucid"].to_list()
runlist = ["0596", "1004"]

track_files_dir.mkdir(parents=True, exist_ok=True)
log_dir.mkdir(parents=True, exist_ok=True)

for _, tuc in tqdm(thucs.iterrows(), total=thucs.shape[0]):
    goodsofar = True

    thuc_id = tuc["tophucid"]

    if thuc_id not in runlist:
        continue

    start = time.time()

    # Make sure the path exists
    thucid_dir_name = f"thuc_{thuc_id}_CDL"
    thuc_dir = track_files_dir / thucid_dir_name

    if thuc_dir.exists():
        now = get_current_time()
        log_to_file(general_log, f"{now}: {nodename}: {thuc_id}: SKIPPING")
        continue
    else:
        thuc_dir.mkdir(parents=True)

    # Update Table and columns for valid existing categories
    try:
        now = get_current_time()
        log_to_file(
            general_log,
            f"{now}: {nodename}: {thuc_id}: Adding columns and updating columns for valid mgmt_field_id...",
        )

        with engine.connect() as connection:
            # Make sure that mgmt_field_id column is indeed text 
            # create a column for cdl value 2022
            # if the mgmt field values are NULL, NaN, Clouds_No_Data, or Developed the cdl_value is set to -1
            query = f"""
            ALTER TABLE thuc_{thuc_id}_annagnps_cell_data_section
            ALTER COLUMN mgmt_field_id SET DATA TYPE text USING CAST(mgmt_field_id AS text);

            ALTER TABLE thuc_{thuc_id}_annagnps_cell_data_section
            ADD COLUMN IF NOT EXISTS cdl_value_2022 INT;

            UPDATE thuc_{thuc_id}_annagnps_cell_data_section 
            SET cdl_value_2022 = CASE 
                WHEN mgmt_field_id IS NULL OR mgmt_field_id IN ({', '.join(invalid_field_ids_to_fix)}) THEN -1 
                ELSE NULL 
            END
            """

            connection.execute(sql_text(query))
            connection.commit()

    except Exception as e:
        connection.rollback()
        goodsofar = False
        now = get_current_time()
        log_to_file(general_log, f"{now}: {nodename}: {thuc_id}: {e}")

    # Find mgmt_field_id that have a cdl_value_2022 different than -1 and we apply the mapping
    if goodsofar:
        
        now = get_current_time()
        log_to_file(
            general_log,
            f"{now}: {nodename}: {thuc_id}: Applying mapping to already valid cells",
        )

        query = f"SELECT mgmt_field_id FROM thuc_{thuc_id}_annagnps_cell_data_section WHERE cdl_value_2022 <> -1 OR cdl_value_2022 IS NULL"
        df = pd.read_sql_query(sql=sql_text(query), con=engine.connect())
        df = df.drop_duplicates().reset_index(drop=True) # Find the unique classes present in that THUC

        with engine.connect() as connection:
            try:

                query = f"""
                UPDATE thuc_{thuc_id}_annagnps_cell_data_section
                SET cdl_value_2022 = CASE WHEN mgmt_field_id = :Mgmt_Field_ID THEN :CDL_Value_2022 ELSE -1 END
                WHERE mgmt_field_id = :Mgmt_Field_ID
                """
            
                for _, class_name in df.iterrows():
                    name = class_name['mgmt_field_id']
                    connection.execute(sql_text(query), {"CDL_Value_2022": reversed_already_valid_dico[name], "Mgmt_Field_ID": name})

                connection.commit()

            except Exception as e:
                print(f"Error for THUC {thuc_id}")
                print(e)
                connection.rollback()
                goodsofar = False

    else:
        pass


    # Query remaining cells and apply zonal startistic only to those
    if goodsofar:
        
        now = get_current_time()
        log_to_file(
            general_log,
            f"{now}: {nodename}: {thuc_id}: Performing plurality analysis on remaining invalid cells",
        )

        try:

            # Get cell ids to reprocess
            query = f"SELECT cell_id FROM thuc_{thuc_id}_annagnps_cell_data_section WHERE cdl_value_2022 = -1 OR cdl_value_2022 IS NULL"
            df = pd.read_sql_query(sql=sql_text(query), con=engine.connect())
            cells_ids_to_reprocess = tuple(df['cell_id'].to_list())

            # Collect cells geometries from database
            query = f"SELECT * FROM thuc_{thuc_id}_annagnps_cell_ids WHERE dn in {cells_ids_to_reprocess}"

            with engine.connect() as conn:
                cells = gpd.read_postgis(sql=sql_text(query), con=conn, geom_col="geom")
                utm = cells.estimate_utm_crs()
                cells = cells.to_crs(utm)

            # Perform the plurality analysis
            cells = sdm.assign_attr_zonal_stats_raster_layer(cells, raster_CDL_path, agg_method='majority', attr='CDL_Value')
            try:
                cells['CDL_Value'] = cells['CDL_Value'].astype('Int32')
            except Exception as e:
                log_to_file(
                general_log,
                f"{now}: {nodename}: {thuc_id}: {e}",
            )
                
            cells['Mgmt_Field_ID'] = cells['CDL_Value'].map(dico)
            cells = cells.rename(columns={"dn": "cell_id"})

            data_to_update = cells[["cell_id", "Mgmt_Field_ID", "CDL_Value"]].to_dict(orient="records")
        except:
            goodsofar = False

    else:
        pass

    # Populate data
    if goodsofar:

        now = get_current_time()
        log_to_file(
            general_log,
            f"{now}: {nodename}: {thuc_id}: Populating cells that were reprocessed",
        )

        try:
            # create a session factory
            Session = sessionmaker(bind=engine)
            # create a new session
            session = Session()
            # create a transaction
            transaction = session.begin()

            # execute your update query here
            query = f"""UPDATE thuc_{thuc_id}_annagnps_cell_data_section 
            SET mgmt_field_id = :Mgmt_Field_ID,
                cdl_value_2022 = :CDL_Value
            WHERE cell_id = :cell_id"""
            
            session.execute(sql_text(query), data_to_update)
            # commit the transaction
            transaction.commit()

        except Exception as e:
            # rollback the transaction on error
            transaction.rollback()
            goodsofar = False

        finally:
            # close the session
            session.close()

    else:
        pass

    if goodsofar:
        now = get_current_time()
        log_to_file(general_log, f"{now}: {nodename}: {thuc_id}: Done!")

    else:
        end = time.time()
        log_to_file(
            general_log,
            f"{now}: {nodename}: {thuc_id}: Did not finish (errors) {end-start} seconds (deleting dir)",
        )
        thuc_dir.rmdir()
        log_to_file(fail_list, f"{thuc_id}")

end = time.time()

print(f"Finished batch processing! Overall process took {(end-bigbang)/3600} hours")
