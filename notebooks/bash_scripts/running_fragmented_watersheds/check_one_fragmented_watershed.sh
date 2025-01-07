#!/bin/bash

# Idea check against those queries and the number of cells in the original watershed
# located in watershed/cell_data_section.csv in the root directory of the thuc

# SELECT COUNT(cell_id) FROM pre_runs_annagnps_aa WHERE thuc_id = '0594'

# SELECT COUNT(cell_id) FROM pre_runs_annagnps_aa_sediment_erosion_ua_rr_total
# WHERE thuc_id = '0594'
# AND description = 'Erosion_Total_All_Sources'

# SELECT COUNT(cell_id) FROM pre_runs_annagnps_aa_sediment_yield_ua_rr_total
# WHERE thuc_id = '0594'
# AND description = 'Sediment_Total_All_Sources'

# SELECT COUNT(cell_id) FROM pre_runs_annagnps_aa_water_yield_ua_rr_total
# WHERE thuc_id = '0594'

# Define function to handle arguments
parse_arguments() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --py_bash_dir)
        PY_BASH_DIR="$2"
        shift 2
        ;;
      --pyagnps_dir)
        PYAGNPS_DIR="$2"
        shift 2
        ;;
      --thuc_id)
        thuc_id="$2"
        shift 2
        ;;
      --credentials)
        credentials="$2"
        shift 2
        ;;
      --check_tables)
        check_tables="$2"  # List of tables to check that are space separated
        shift 2
        ;;
      --log_file)
        LOG_FILE="$2"
        shift 2
        ;;
      --failed_log_file)
        FAILED_THUCS="$2"
        shift 2
        ;;
      --success_thucs)
        SUCCESS_THUCS="$2"
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

# Make a default value of the LOG_FILE in case it is not specified so that it doesn't log to a file
if [ -z "$LOG_FILE" ]; then
  LOG_FILE="/dev/null"
fi

# Make a default value of the FAILED_THUCS in case it is not specified so that it doesn't log to a file
if [ -z "$FAILED_THUCS" ]; then
  FAILED_THUCS="/dev/null"
fi

# Make a default value of the SUCCESS_THUCS in case it is not specified so that it doesn't log to a file
if [ -z "$SUCCESS_THUCS" ]; then
  SUCCESS_THUCS="/dev/null"
fi


if [ -z "$PY_BASH_DIR" ]; then
  echo "$(date '+%Y-%m-%d %H:%M:%S') - Error: Missing required argument: --py_bash_dir" | tee -a "$LOG_FILE"
  exit 1
fi

if [ -z "$PYAGNPS_DIR" ]; then
  echo "$(date '+%Y-%m-%d %H:%M:%S') - Error: Missing required argument: --pyagnps_dir" | tee -a "$LOG_FILE"
  exit 1
fi

if [ -z "$check_tables" ]; then
  check_tables="pre_runs_annagnps_aa pre_runs_annagnps_aa_sediment_erosion_ua_rr_total pre_runs_annagnps_aa_sediment_yield_ua_rr_total pre_runs_annagnps_aa_water_yield_ua_rr_total"
fi

if [ -z "$thuc_id" ]; then
  thuc_id=""
fi

if [ -z "$credentials" ]; then
  path_to_db_credentials=""
else
  path_to_db_credentials="$credentials"
fi


# Activate the virtual environment
source "$PYAGNPS_DIR/venv/bin/activate" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - Failed to activate virtual environment" | tee -a "$LOG_FILE"; exit 1; }

# export PYTHONUNBUFFERED=TRUE

# Parse check_tables (space separated) into an array
IFS=' ' read -r -a check_tables_array <<< "$check_tables"

# Check that all of them return successfully with 0 and if that's the case append the thuc_id at the end of a file
success=true
for table in "${check_tables_array[@]}"; do

  echo "$(date '+%Y-%m-%d %H:%M:%S') - Checking table: $table" | tee -a "$LOG_FILE"
  python -u "${PY_BASH_DIR}/check_one_fragmented_watershed.py" \
    --credentials "$path_to_db_credentials" \
    --thuc_id "$thuc_id" \
    --table "$table" \
    --log_file "$LOG_FILE"

  if [ $? -ne 0 ]; then
    success=false
  fi
done

if [ $success = true ]; then
  echo "$thuc_id" >> "$SUCCESS_THUCS"
fi

echo "$(date '+%Y-%m-%d %H:%M:%S') - Check finished with success status: $success" | tee -a "$LOG_FILE"

# If all the checks were successful, exit with 0
if [ $success = true ]; then
  exit 0
fi


