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
from sqlalchemy.exc import SQLAlchemyError

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
    database=database
)

# create a SQLAlchemy engine object
engine = create_engine(url_object)

# PATHS SETUP
bigbang = time.time()
nodename = socket.gethostname()

path_to_thucs = Path(
    "/aims-nas/luc/data/tophuc_S_M_40000_closed_holes_with_container_thuc_merged_bbox_area_first_kept.gpkg"
)

path_to_raster = Path("/aims-nas/data/datasets/Management/CDL_Annual/CDL_2022.tif")

track_files_dir = Path("/aims-nas/luc/thuc_field_ID_CDL/")
path_to_management_class_names = Path('/aims-nas/data/datasets/Management/CDL_Annual/definition_of_CDL_Class_Names_for_use_in_AIMS/CDL_Field_ID_dictionary.csv')

log_dir = track_files_dir / "LOGS"

general_log = log_dir / f"{nodename}_batch_field_ID_CDL_population_general_log.txt"
fail_list = log_dir / f"{nodename}_fail_list.txt"

print('Reading and initializing files...')
thucs = gpd.read_file(
    path_to_thucs
)  # GeoDataFrame containing the thucs and their geometry
thucs = thucs.sort_values(by=["bbox_area_sqkm"], ascending=True)

df_cdl = pd.read_csv(path_to_management_class_names)
dico = df_cdl[['CDL_Value','Modified_CDL_Category']].set_index('CDL_Value').to_dict(orient='dict')['Modified_CDL_Category']

runlist = thucs["tophucid"].to_list()
# runlist = ["1002", "1004"]
# runlist = ["1002"]

track_files_dir.mkdir(parents=True, exist_ok=True)
log_dir.mkdir(parents=True, exist_ok=True)

# runlist = pd.read_csv(path_to_thuc_runlist, dtype=object)
# runlist = runlist.iloc[:,0].to_list() # Get the list of thucs that need to be

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


    # Collect thuc cells geometry from database
    try:
        now = get_current_time()
        log_to_file(
            general_log,
            f"{now}: {nodename}: {thuc_id}: Querying cells from database...",
        )
        query = f"SELECT * FROM thuc_{thuc_id}_annagnps_cell_ids"

        with engine.connect() as conn:
            cells = gpd.read_postgis(sql=sql_text(query), con=conn, geom_col="geom")

        utm = cells.estimate_utm_crs()
        cells = cells.to_crs(utm)

    except Exception as e:
        goodsofar = False
        now = get_current_time()
        log_to_file(general_log, f"{now}: {nodename}: {thuc_id}: {e}")



    # Apply plurality analysis
    if goodsofar:
        try:
            now = get_current_time()
            log_to_file(
                general_log,
                f"{now}: {nodename}: {thuc_id}: Performing plurality analysis",
            )
            cells = sdm.assign_attr_zonal_stats_raster_layer(cells, path_to_raster, agg_method='majority', attr='CDL_Value')

            # Reformat the CDL Value column and map with the correct value
            cells['CDL_Value'] = cells['CDL_Value'].astype('Int32')
            cells.loc[cells['CDL_Value']==0, 'CDL_Value'] = 81 # Set 0 value to 81 = Cloud_No_Data
            cells['Mgmt_Field_ID'] = cells['CDL_Value'].map(dico)
            # this function uses rasterstats.zonal_stats and the "majority" function actually does the plurality operation by selecting the most common value

            cells = cells.rename(columns={"dn": "cell_id"})

            data_to_update = cells[["cell_id", "Mgmt_Field_ID"]].to_dict(orient="records")

        except Exception as e:
            goodsofar = False
            now = get_current_time()
            log_to_file(general_log, f"{now}: {nodename}: {thuc_id}: {e}")
    else:
        pass

    # Update table schema and make sure that mgmt_field_id is of type text
    if goodsofar:
        with engine.connect() as connection:
            try:
                query = f"ALTER TABLE thuc_{thuc_id}_annagnps_cell_data_section ALTER COLUMN mgmt_field_id TYPE TEXT"
                connection.execute(sql_text(query))
                # Commit the transaction explicitly
                connection.commit()

            except Exception as e:
                goodsofar = False
                
                now = get_current_time()
                log_to_file(general_log, f"{now}: {nodename}: {thuc_id}: {e} (rolling back)")

                # Rollback the transaction in case of an error
                connection.rollback()
    else:
        pass

    # Update database
    if goodsofar:
        now = get_current_time()
        log_to_file(general_log, f"{now}: {nodename}: {thuc_id}: Populating Field_ID with CDL data...")

        try:
            # create a session factory
            Session = sessionmaker(bind=engine)
            # create a new session
            session = Session()
            # create a transaction
            transaction = session.begin()

            # execute your update query here
            query = f"UPDATE thuc_{thuc_id}_annagnps_cell_data_section SET Mgmt_Field_ID = :Mgmt_Field_ID WHERE cell_id = :cell_id"
            session.execute(sql_text(query), data_to_update)
            # commit the transaction
            transaction.commit()

        except Exception as e:
            goodsofar = False
            # rollback the transaction on error
            transaction.rollback()
            now = get_current_time()
            log_to_file(
                general_log,
                f"{now}: {nodename}: {thuc_id}: Failed to update DB, rolling back... : {e}",
            )

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
