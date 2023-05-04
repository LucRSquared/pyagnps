import socket
import time
import json
from pathlib import Path

from tqdm import tqdm

import geopandas as gpd
import pandas as pd

from pyagnps import soil_data_market as sdm
from pyagnps.utils import log_to_file, get_current_time

from sqlalchemy import create_engine, text as sql_text
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

BUFFER = 100

# create a SQLAlchemy engine object
engine = create_engine(f"postgresql://{user}:{password}@{host}:{port}/{database}")

# PATHS SETUP
bigbang = time.time()
nodename = socket.gethostname()

path_to_thucs = Path(
    "/aims-nas/luc/data/tophuc_S_M_40000_closed_holes_with_container_thuc_merged_bbox_area_first_kept.gpkg"
)

path_to_gssurgo_gdb = Path("/aims-nas/data/datasets/gSSURGO/gSSURGO_CONUS_202210.gdb")

ssurgo_files_dir = Path("/aims-nas/luc/thuc_soils_shapefiles/")

log_dir = ssurgo_files_dir / "LOGS"

general_log = log_dir / f"{nodename}_batch_ssurgo_population_general_log.txt"
fail_list = log_dir / f"{nodename}_fail_list.txt"

thucs = gpd.read_file(
    path_to_thucs
)  # GeoDataFrame containing the thucs and their geometry
thucs = thucs.sort_values(by=["bbox_area_sqkm"], ascending=False)

runlist = thucs['tophucid'].to_list()
# runlist = ["1002", "1004"]
# runlist = ["1002"]

ssurgo_files_dir.mkdir(parents=True, exist_ok=True)
log_dir.mkdir(parents=True, exist_ok=True)

# runlist = pd.read_csv(path_to_thuc_runlist, dtype=object)
# runlist = runlist.iloc[:,0].to_list() # Get the list of thucs that need to be

for _, tuc in tqdm(thucs.iterrows(), total=thucs.shape[0]) :
    goodsofar = True

    thuc_id = tuc["tophucid"]

    if thuc_id not in runlist:
        continue

    start = time.time()

    # Make sure the path exists
    thucid_dir_name = f"thuc_{thuc_id}_ssurgo"
    thuc_dir = ssurgo_files_dir / thucid_dir_name

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


    # Select thuc for soil population
    cells_buffer = cells.copy(deep=True)
    cells_buffer['geom'] = cells_buffer.geometry.buffer(BUFFER)

    # Get SSURGO polygon geometry
    if goodsofar:
        try:
            now = get_current_time()
            log_to_file(
                general_log,
                f"{now}: {nodename}: {thuc_id}: Getting gSSURGO data",
            )

            geo_soil = gpd.read_file(path_to_gssurgo_gdb, driver='OpenFileGDB', layer='MUPOLYGON', bbox=cells_buffer)
            geo_soil = geo_soil.to_crs(utm)

        except Exception as e:
            goodsofar = False
            now = get_current_time()
            log_to_file(general_log, f"{now}: {nodename}: {thuc_id}: {e}")

    else:
        continue

    now = get_current_time()
    log_to_file(
        general_log,
        f"{now}: {nodename}: {thuc_id}: Checking if cells are covered by gSSURGO polygon",
    )

    # Check that the soil geometry covers all the cells
    if not (cells.within(geo_soil.unary_union.envelope).all()):
        now = get_current_time()
        log_to_file(
            general_log,
            f"{now}: {nodename}: {thuc_id}: WARNING /!\ Not all cells are covered by soil geometry!",
        )

    else:
        pass
        #print(f"{now}: {nodename}: {thuc_id}: All cells are covered with soil data")

    # Apply plurality analysis
    if goodsofar:
        try:
            now = get_current_time()
            log_to_file(
                general_log,
                f"{now}: {nodename}: {thuc_id}: Performing plurality analysis",
            )
            cells = sdm.assign_attr_plurality_vector_layer(
                cells, geo_soil, attr="MUKEY", bin_id="dn"
            )
            cells = cells.rename(columns={"dn": "cell_id", "MUKEY": "soil_id"})

        except Exception as e:
            goodsofar = False
            now = get_current_time()
            log_to_file(general_log, f"{now}: {nodename}: {thuc_id}: {e}")

    else:
        continue

    # Update database
    if goodsofar:
        now = get_current_time()
        log_to_file(general_log, f"{now}: {nodename}: {thuc_id}: Populating Soil_ID...")

        data_to_update = cells[["cell_id", "soil_id"]].to_dict(orient="records")

        # create a session factory
        Session = sessionmaker(bind=engine)
        # create a new session
        session = Session()
        # create a transaction
        transaction = session.begin()

        try:
            # execute your update query here
            query = f"UPDATE thuc_{thuc_id}_annagnps_cell_data_section SET soil_id = :soil_id WHERE cell_id = :cell_id"
            session.execute(sql_text(query), data_to_update)
            # commit the transaction
            transaction.commit()

        except SQLAlchemyError as e:
            goodsofar = False
            # rollback the transaction on error
            transaction.rollback()
            now = get_current_time()
            log_to_file(
                general_log,
                f"{now}: {nodename}: {thuc_id}: Failed to update DB, rolling back...",
            )

        finally:
            # close the session
            session.close()

    else:
        continue

    now = get_current_time()
    if goodsofar:
        log_to_file(general_log, f"{now}: {nodename}: {thuc_id}: Done!")

    else:
        end = time.time()
        log_to_file(
            general_log,
            f"{now}: {nodename}: {thuc_id}: Did not finish (errors) {end-start} seconds",
        )
        log_to_file(fail_list, f"{thuc_id}")

end = time.time()

print(f"Finished batch processing! Overall process took {(end-bigbang)/3600} hours")
