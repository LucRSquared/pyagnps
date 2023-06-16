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

# create a SQLAlchemy engine object
engine = create_engine(f"postgresql://{user}:{password}@{host}:{port}/{database}")

# PATHS SETUP
bigbang = time.time()
nodename = socket.gethostname()

path_to_thucs = Path(
    "/aims-nas/luc/data/tophuc_S_M_40000_closed_holes_with_container_thuc_merged_bbox_area_first_kept.gpkg"
)

path_to_gnatsgo_raster = Path("/aims-nas/data/datasets/gNATSGO/gNATSGO-mukey.tif")

path_to_gnatsgo_sapolygon = Path("/aims-nas/data/datasets/gNATSGO/gNATSGO_SAPOLYGON.gpkg")

gNATSGO_files_dir = Path("/aims-nas/luc/thuc_gNATSGO/")

log_dir = gNATSGO_files_dir / "LOGS"

general_log = log_dir / f"{nodename}_batch_ssurgo_population_general_log.txt"
fail_list = log_dir / f"{nodename}_fail_list.txt"

thucs = gpd.read_file(
    path_to_thucs
)  # GeoDataFrame containing the thucs and their geometry
thucs = thucs.sort_values(by=["bbox_area_sqkm"], ascending=False)

runlist = thucs["tophucid"].to_list()
# runlist = ["1002", "1004"]
# runlist = ["1002"]

gNATSGO_files_dir.mkdir(parents=True, exist_ok=True)
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
    thucid_dir_name = f"thuc_{thuc_id}_gNATSGO"
    thuc_dir = gNATSGO_files_dir / thucid_dir_name

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

    # Get SSURGO polygon geometry
    if goodsofar:
        try:
            now = get_current_time()
            log_to_file(
                general_log,
                f"{now}: {nodename}: {thuc_id}: Getting cells boundary",
            )

            boundary = cells.copy(deep=True)
            boundary['geom'] = boundary['geom'].buffer(0)
            boundary = boundary.unary_union
            boundary = gpd.GeoDataFrame(geometry=[boundary], crs=utm)

            now = get_current_time()
            log_to_file(
                general_log,
                f"{now}: {nodename}: {thuc_id}: Selecting gNATSGO SAPOLYGON in the boundary")

            sapolygon = gpd.read_file(path_to_gnatsgo_sapolygon, rows=0)
            crs_sapolygon = sapolygon.crs

            sapolygon = gpd.read_file(path_to_gnatsgo_sapolygon, bbox=boundary.to_crs(crs_sapolygon))

        except Exception as e:
            goodsofar = False
            now = get_current_time()
            log_to_file(general_log, f"{now}: {nodename}: {thuc_id}: {e}")

    else:
        pass

    # Only keep the cells covered by SAPOLYGON polygon
    if goodsofar:
        try:
            now = get_current_time()
            log_to_file(
                general_log,
                f"{now}: {nodename}: {thuc_id}: Restricting cells to local boundary not covered by SSURGO",
            )

            cells_to_update= cells.overlay(sapolygon.to_crs(utm), how='intersection')

            cells_to_update = cells_to_update.loc[cells_to_update['SOURCE'] != 'SSURGO', ['dn','geometry']]
            cells_to_update = cells_to_update.drop_duplicates()
            cells_to_update = cells_to_update.rename(columns={"geometry": "geom"})
            cells_to_update = cells_to_update.set_geometry('geom')

            if cells_to_update.empty:
                raise Exception("No overlap of cells and gSSURGO data")

        except Exception as e:
            goodsofar = False
            now = get_current_time()
            log_to_file(general_log, f"{now}: {nodename}: {thuc_id}: {e}")
    else:
        pass



    # Apply plurality analysis
    if goodsofar:
        try:
            now = get_current_time()
            log_to_file(
                general_log,
                f"{now}: {nodename}: {thuc_id}: Performing plurality analysis",
            )
            cells = sdm.assign_attr_zonal_stats_raster_layer(cells_to_update, path_to_gnatsgo_raster, agg_method='majority')
            # this function uses rasterstats.zonal_stats and the "majority" function actually does the plurality operation by selecting the most common value

            cells = cells.rename(columns={"dn": "cell_id", "MUKEY": "soil_id"})

        except Exception as e:
            goodsofar = False
            now = get_current_time()
            log_to_file(general_log, f"{now}: {nodename}: {thuc_id}: {e}")

    else:
        pass

    # Update database
    if goodsofar:
        now = get_current_time()
        log_to_file(general_log, f"{now}: {nodename}: {thuc_id}: Populating missing Soil_ID with gNATSGO data...")

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
        pass

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
