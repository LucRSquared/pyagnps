#!/bin/bash

# Define function to handle arguments
parse_arguments() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --simulation_dir)
        SIMULATION_DIR="$2"
        shift 2
        ;;
      --pyagnps_dir)
        PYAGNPS_DIR="$2"
        shift 2
        ;;
      --py_bash_dir)
        PY_BASH_DIR="$2"
        shift 2
        ;;
      --credentials)
        path_to_db_credentials="$2"
        shift 2
        ;;
      --nldas2_centroids)
        path_to_nldas2_centroids="$2"
        shift 2
        ;;
      --scs_storm_types)
        path_to_scs_storm_types="$2"
        shift 2
        ;;
      --precip_zones)
        path_to_precip_zones="$2"
        shift 2
        ;;
      --start_date)
        start_date="$2"
        shift 2
        ;;
      --end_date)
        end_date="$2"
        shift 2
        ;;
      --climate_method)
        climate_method="$2"
        shift 2
        ;;
      --climate_table)
        climate_table="$2"
        shift 2
        ;;
      --thuc_id)
        thuc_id="$2"
        shift 2
        ;;
      --reach_id)
        reach_id="$2"
        shift 2
        ;;
      --generate_main_files)
        generate_main_files="$2"
        shift 2
        ;;
      --fragment_watershed)
        fragment_watershed="$2"
        shift 2
        ;;
      --share_global_watershed_parameters_with_mini_watersheds)
        share_global_watershed_parameters_with_mini_watersheds="$2"
        shift 2
        ;;
      --num_processes)
        num_processes="$2"
        shift 2
        ;;
      --log_file)
        LOG_FILE="$2"
        shift 2
        ;;
      *)
        echo "Invalid argument: $1"
        exit 1
        ;;
    esac
  done
}

# Call the function to parse arguments
parse_arguments "$@"

# Check if required arguments are provided
if [ -z "$SIMULATION_DIR" -o -z "$path_to_db_credentials" ]; then
  echo "$(date '+%Y-%m-%d %H:%M:%S') - Error: Missing required arguments: simulation_dir and credentials" | tee -a "$LOG_FILE"
  exit 1
fi

# Set defaults for optional arguments
if [ -z "$reach_id" ]; then
  reach_id=2
fi

if [ -z "$climate_method" ]; then
  climate_method="nldas2_database"
fi

if [ -z "$climate_table" ]; then
  climate_table="climate_nldas2"
fi

if [ -z "$fragment_watershed" ]; then
  fragment_watershed="true"
fi

if [ -z "$share_global_watershed_parameters_with_mini_watersheds" ]; then
  share_global_watershed_parameters_with_mini_watersheds="true"
fi

if [ -z "$generate_main_files" ]; then
  generate_main_files="true"
fi

if [ -z "$num_processes" ]; then
  num_processes=16
fi

if [ -z "$PYAGNPS_DIR" ]; then
  PYAGNPS_DIR="/aims-nas/luc/code/pyagnps/"
fi

# Make a default value of the LOG_FILE in case it is not specified so that it doesn't log to a file
if [ -z "$LOG_FILE" ]; then
  LOG_FILE="/dev/null"
fi

# Activate the virtual environment
source "$PYAGNPS_DIR/venv/bin/activate" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - Failed to activate virtual environment" | tee -a "$LOG_FILE"; exit 1; }

# export PYTHONUNBUFFERED=TRUE

# Run the script
python -u "$PY_BASH_DIR/generate_annagnps_files.py"
  --credentials "$path_to_db_credentials" \
  --nldas2_centroids "$path_to_nldas2_centroids" \
  --scs_storm_types "$path_to_scs_storm_types" \
  --precip_zones "$path_to_precip_zones" \
  --output_folder "$SIMULATION_DIR" \
  --start_date "$start_date" \
  --end_date "$end_date" \
  --climate_method "$climate_method" \
  --climate_table "$climate_table" \
  --thuc_id "$thuc_id" \
  --reach_id "$reach_id" \
  --generate_main_files "$generate_main_files" \
  --fragment_watershed "$fragment_watershed" \
  --num_processes "$num_processes" \
  --log_file "$LOG_FILE"