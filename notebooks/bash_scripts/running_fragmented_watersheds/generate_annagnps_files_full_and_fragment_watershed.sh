#!/bin/bash

# Define function to handle arguments
parse_arguments() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --simulation_dir)
        SIMULATION_DIR="$2"
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
      --num_processes)
        num_processes="$2"
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
  echo "Error: Missing required arguments: simulation_dir and credentials"
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

if [ -z "$generate_main_files" ]; then
  generate_main_files="true"
fi

if [ -z "$num_processes" ]; then
  num_processes=16
fi

# Rest of the script remains the same...

# Activate the virtual environment
source "$PYAGNPS_DIR/venv/bin/activate" || { echo "Failed to activate virtual environment"; exit 1; }

export PYTHONUNBUFFERED=TRUE

# Run the script
python generate_watershed_files_pre_runs.py \
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
  --num_processes "$num_processes"