#!/bin/bash

# STRUCTURE
# ROOT_DIR
# ├── thuc_list_to_run.csv
# ├── thuc_0001
# │   ├── (bunch of regular input files)
# │   └── mini_watersheds
# │       ├── dir_list.csv
# │       ├── reach_000000002
# │       ├── reach_000000003
# │       ├── reach_000000004
# │       └── etc.
# │
# ├──  thuc_0002
#     ├── (bunch of regular input files)
#     └── mini_watersheds
#         ├── dir_list.csv
#         ├── reach_000000002
#         ├── reach_000000003
#         ├── reach_000000004
#         └── etc.

# Set the root directory (can be changed before script execution). This directory should contain thuc_list_to_run.csv and a directory for each watershed that needs to be processed
ROOT_DIR="/aims-nas/luc/annagnps_pre_runs_2000-01-01_2022-12-31_unforced_potet/"  # Needs to use absolute path

PY_BASH_DIR="/aims-nas/luc/code/pyagnps/notebooks/bash_scripts/running_fragmented_watersheds/" # the location of the python scripts are defined with respect to this

LOG_FILE="${ROOT_DIR}/annagnps_pre_runs_2000-01-01_2022-12-31.log"
FAILED_THUCS="${ROOT_DIR}/failed_thucs.csv"
chmod 666 "${FAILED_THUCS}"


# Set boolean parameters to generate and fragment, simulate, post-process results (and populate them to the database)
generate_main_files="true"


fragment_watershed="true"
share_global_watershed_parameters_with_mini_watersheds="true"


simulate_thuc="true"
force_simulate="false" # if true will overwrite existing simulation outputs

post_process="true"

climate_method="nldas2_database"
climate_table="climate_nldas2"

start_date="2000-01-01"
end_date="2022-12-31"
# end_date="2002-12-31"

pyagnps_dir="/aims-nas/luc/code/pyagnps" # the location of the python scripts are defined with respect to this

path_to_db_credentials="/aims-nas/luc/code/pyagnps/inputs/db_credentials.json"

path_to_nldas2_centroids="/aims-nas/data/datasets/CLIMATE/NLDAS2/NLDAS2_GRID_CENTROIDS_epsg4326.gpkg"
path_to_scs_storm_types="/aims-nas/data/datasets/TR-55/scs_storm_types.gpkg"
path_to_precip_zones="/aims-nas/data/datasets/RUSLE2/Climate/precip_zones_RUSLE2_cleaned_manually_extrapolated_pchip_linear_US_units.gpkg"

# partition="aims-highperf-oversubscribe,aims-default-oversubscribe"
partition="aims-highperf-oversubscribe"

# Nodes to exclude
exclude="aims-node6,aims-node7,aims-node11"


# Batch size for job simulations submissions
batch_size=500
maxiter=1000
num_processes=32

# Print parameters to the log file
{
  echo "ROOT_DIR: $ROOT_DIR"
  echo "PY_BASH_DIR: $PY_BASH_DIR"
  echo "LOG_FILE: $LOG_FILE"
  echo "generate_main_files: $generate_main_files"
  echo "fragment_watershed: $fragment_watershed"
  echo "share_global_watershed_parameters_with_mini_watersheds: $share_global_watershed_parameters_with_mini_watersheds"
  echo "simulate_thuc: $simulate_thuc"
  echo "force_simulate: $force_simulate"
  echo "post_process: $post_process"
  echo "climate_method: $climate_method"
  echo "climate_table: $climate_table"
  echo "start_date: $start_date"
  echo "end_date: $end_date"
  echo "pyagnps_dir: $pyagnps_dir"
  echo "path_to_db_credentials: $path_to_db_credentials"
  echo "path_to_nldas2_centroids: $path_to_nldas2_centroids"
  echo "path_to_scs_storm_types: $path_to_scs_storm_types"
  echo "path_to_precip_zones: $path_to_precip_zones"
  echo "partition: $partition"
  echo "exclude: $exclude"
  echo "batch_size: $batch_size"
  echo "maxiter: $maxiter"
  echo "num_processes: $num_processes"
} | tee -a "$LOG_FILE"

# Loop through thuc list
# For each thuc:
#   - check if AnnAGNPS.fil exists (if yes it means the regular input files were already generated)
#   - check if dir_list.csv exists within the mini_watersheds directory
#   - if not generate the files and fragment them using the python script "fragment_annagnps_into_mini_watersheds.py" with parameters startDate, endDate, thuc_ID
#   - if yes, use the existing dir_list.csv
#   - check that AnnAGNPS.log does not exist
#   - if not, run the simulation using the distributed method
#   - read results and upload to db
#   - cleanup files that are not needed
#   - log results

# Calculate the total number of jobs based on directory count
num_jobs=0

if [ -f "${ROOT_DIR}/thuc_list_to_run.csv" ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Reusing existing thuc_list_to_run.csv file" | tee -a "$LOG_FILE"
    num_jobs=$(wc -l < "${ROOT_DIR}/thuc_list_to_run.csv")
fi

echo "$(date '+%Y-%m-%d %H:%M:%S') - Found $num_jobs watersheds" | tee -a "$LOG_FILE"

# Loop to process thucs one by one
for ((thuc_index = 1; thuc_index <= num_jobs; thuc_index += 1)); do
  START_ALL=$(date +%s)

  thuc_id=$(sed -n "${thuc_index}p" "${ROOT_DIR}/thuc_list_to_run.csv" | tr -d '\r\n') # The tr command removes trailing new line

  THUC_LOG_FILE="${LOG_FILE%.*}_${thuc_id}.log"

  echo "$(date '+%Y-%m-%d %H:%M:%S') - Processing thuc $thuc_id" | tee -a "$THUC_LOG_FILE"

  # Check if regular input files exist (AnnAGNPS.fil) exists
  if [ -f "${ROOT_DIR}/thuc_${thuc}/AnnAGNPS.fil" ] || [ "$generate_main_files" = "false" ]; then
      echo "$(date '+%Y-%m-%d %H:%M:%S') - Skipping thuc $thuc_id, input files already exist" | tee -a "$THUC_LOG_FILE"
      generate_main_files="false"
  else
      generate_main_files="true"
  fi

  # Check if dir_list.csv exists
  if [ -f "${ROOT_DIR}/thuc_${thuc}/mini_watersheds/dir_list.csv" ] || [ "$fragment_watershed" = "false" ]; then
      fragment_watershed="false"
  else
      fragment_watershed="true"
  fi

  #---------------------------------------
  START_GENERATE=$(date +%s)
  if [[ "$generate_main_files" == "true" ]]; then
  
      
      echo "$(date '+%Y-%m-%d %H:%M:%S') - Generating input files for thuc $thuc_id and fragmenting it into mini watersheds" | tee -a "$THUC_LOG_FILE"

      "${PY_BASH_DIR}/generate_annagnps_files_full_and_fragment_watershed.sh" \
        --pyagnps_dir "$pyagnps_dir" \
        --py_bash_dir "$PY_BASH_DIR" \
        --simulation_dir "${ROOT_DIR}/thuc_${thuc_id}" \
        --credentials "$path_to_db_credentials" \
        --nldas2_centroids "$path_to_nldas2_centroids" \
        --scs_storm_types "$path_to_scs_storm_types" \
        --precip_zones "$path_to_precip_zones" \
        --start_date "$start_date" \
        --end_date "$end_date" \
        --climate_method "$climate_method" \
        --climate_table "$climate_table" \
        --thuc_id "$thuc_id" \
        --reach_id "$reach_id" \
        --generate_main_files "$generate_main_files" \
        --fragment_watershed "$fragment_watershed" \
        --share_global_watershed_parameters_with_mini_watersheds "$share_global_watershed_parameters_with_mini_watersheds" \
        --num_processes "$num_processes" \
        --log_file "$THUC_LOG_FILE" || { # what to do if generation of files fails
          cd "${ROOT_DIR}" ; 
          echo "$(date '+%m-%d %H:%M:%S') - Error: Failed to generate input files for thuc $thuc_id" | tee -a "$THUC_LOG_FILE"
          echo "$(date '+%Y-%Y-%m-%d %H:%M:%S'),$thuc_id,failed_generation" | tee -a "$FAILED_THUCS"
          continue
        }
  fi

  END_GENERATE=$(date +%s)
  echo "$(date '+%Y-%m-%d %H:%M:%S') - [$thuc_index/$num_jobs] - Finished generating thuc $thuc_id, elapsed time: $(($END_GENERATE - $START_GENERATE)) seconds" | tee -a "$LOG_FILE"
  
  START_SIMULATE=$(date +%s)
  if [[ "$simulate_thuc" == "true" ]]; then
      echo "$(date '+%Y-%m-%d %H:%M:%S') - Simulating thuc $thuc_id using distributed method" | tee -a "$THUC_LOG_FILE"
      cd "${ROOT_DIR}/thuc_${thuc_id}"
      "${PY_BASH_DIR}/process_one_fragmented_watershed_run_sim.sh" \
        --mini_watersheds_dir "./mini_watersheds" \
        --py_bash_dir "$PY_BASH_DIR" \
        --thuc_id "$thuc_id" \
        --batch_size "$batch_size" \
        --force_simulate "$force_simulate" \
        --maxiter "$maxiter" \
        --partition "$partition" \
        --exclude "$exclude" \
        --log_file "$THUC_LOG_FILE" \
        --failed_log_file "$FAILED_THUCS" || { # what to do if simulation fails
          cd "${ROOT_DIR}" ;
          echo "$(date '+%Y-%m-%d %H:%M:%S'),$thuc_id,failed_simulation" | tee -a "$FAILED_THUCS" 
          echo "$(date '+%Y-%m-%d %H:%M:%S') - Error: Simulation failed for thuc $thuc_id, continuing" | tee -a "$THUC_LOG_FILE"
          continue
        }
      cd "${ROOT_DIR}"
  fi

  END_SIMULATE=$(date +%s)
  echo "$(date '+%Y-%m-%d %H:%M:%S') - [$thuc_index/$num_jobs] - Finished simulating thuc $thuc_id, elapsed time: $(($END_SIMULATE - $START_SIMULATE)) seconds" | tee -a "$LOG_FILE"
  
  START_POSTPROCESS=$(date +%s)
  if [[ "$post_process" == "true" ]]; then
      echo "$(date '+%Y-%m-%d %H:%M:%S') - Post processing thuc $thuc_id and uploading results to db" | tee -a "$THUC_LOG_FILE"
      cd "${ROOT_DIR}/thuc_${thuc_id}"
      "${PY_BASH_DIR}/post_process_one_fragmented_watershed.sh" \
        --mini_watersheds_dir "./mini_watersheds" \
        --thuc_id "$thuc_id" \
        --pyagnps_dir "$pyagnps_dir" \
        --py_bash_dir "$PY_BASH_DIR" \
        --credentials "$path_to_db_credentials" \
        --partition "$partition" \
        --exclude "$exclude" \
        --batch_size 20 \
        --log_file "$THUC_LOG_FILE" \
        --failed_log_file "$FAILED_THUCS" || { # what to do if post processing fails
          cd "${ROOT_DIR}" ; 
          echo "$(date '+%Y-%m-%d %H:%M:%S'),$thuc_id,failed_postprocessing" | tee -a "$FAILED_THUCS"
          echo "$(date '+%Y-%m-%d %H:%M:%S') - Error: Post processing failed for thuc $thuc_id, continuing" | tee -a "$THUC_LOG_FILE"
          continue
        }
      cd "${ROOT_DIR}"
  fi

  END_POSTPROCESS=$(date +%s)
  echo "$(date '+%Y-%m-%d %H:%M:%S') - [$thuc_index/$num_jobs] - Finished post processing thuc $thuc_id, elapsed time: $(($END_POSTPROCESS - $START_POSTPROCESS)) seconds" | tee -a "$LOG_FILE"

  END_ALL=$(date +%s)
  echo "$(date '+%Y-%m-%d %H:%M:%S') - [$thuc_index/$num_jobs] - Finished processing thuc $thuc_id, elapsed time: $(($END_ALL - $START_ALL)) seconds" | tee -a "$LOG_FILE"
done


