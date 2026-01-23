#!/bin/bash

# Define function to handle arguments
parse_arguments() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --mini_watersheds_dir)
        MINI_WATERSHEDS_DIR="$2"
        shift 2
        ;;
      --thuc_id)
        thuc_id="$2"
        shift 2
        ;;
      --note)
        note="$2"
        shift 2
        ;;
      --credentials)
        path_to_db_credentials="$2"
        shift 2
        ;;
      --csv_file)
        csv_file="$2"
        shift 2
        ;;
      --pyagnps_dir)
        pyagnps_dir="$2"
        shift 2
        ;;
      --py_bash_dir)
        PY_BASH_DIR="$2"
        shift 2
        ;;
      --log_file)
        LOG_FILE="$2"
        shift 2
        ;;
      --use_local_log)
        USE_LOCAL_LOG="$2"
        shift 2
        ;;
      --failed_log_file)
        FAILED_THUCS="$2"
        shift 2
        ;;
      --save_method)
        save_method="$2"
        shift 2
        ;;
      --annagnps_aa_table)
        annagnps_aa_table="$2"
        shift 2
        ;;
      --aa_water_yield_table)
        aa_water_yield_table="$2"
        shift 2
        ;;
      --aa_sediment_yield_table)
        aa_sediment_yield_table="$2"
        shift 2
        ;;
      --aa_sediment_erosion_table)
        aa_sediment_erosion_table="$2"
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

# Make a default value of the LOG_FILE in case it is not specified so that it doesn't log to a file
if [ -z "$USE_LOCAL_LOG" ]; then
  USE_LOCAL_LOG="false"
fi


# Make a default value of the FAILED_THUCS in case it is not specified so that it doesn't log to a file
if [ -z "$FAILED_THUCS" ]; then
  FAILED_THUCS="/dev/null"
fi

# Set default value for MINI_WATERSHEDS_DIR if not provided
if [ -z "$MINI_WATERSHEDS_DIR" ]; then
    MINI_WATERSHEDS_DIR=$(pwd)
fi

# Set default value for csv_file if not provided
if [ -z "$csv_file" ]; then
    csv_file="${MINI_WATERSHEDS_DIR}/dir_list.csv"
fi

if [ -z "$PYAGNPS_DIR" ]; then
  PYAGNPS_DIR="/aims-nas/luc/code/pyagnps/"
fi

if [ -z "$PY_BASH_DIR" ]; then
  echo "$(date '+%Y-%m-%d %H:%M:%S') - Error: Missing required argument: --py_bash_dir" | tee -a "$LOG_FILE"
  exit 1
fi

if [ -z "$save_method" ]; then
  save_method="db" # 'db' or 'files'
fi

if [ -z "$thuc_id" ]; then
  echo "$(date '+%Y-%m-%d %H:%M:%S') - Error: Missing required argument: --thuc_id" | tee -a "$LOG_FILE"
  exit 1
fi

if [ -z "$note" ]; then
  note="unforced_potet"
fi

if [ -z "$annagnps_aa_table" ]; then
  annagnps_aa_table="pre_runs_annagnps_aa"
fi

if [ -z "$aa_water_yield_table" ]; then
  aa_water_yield_table="pre_runs_annagnps_aa_water_yield_ua_rr_total"
fi

if [ -z "$aa_sediment_yield_table" ]; then
  aa_sediment_yield_table="pre_runs_annagnps_aa_sediment_yield_ua_rr_total"
fi

if [ -z "$aa_sediment_erosion_table" ]; then
  aa_sediment_erosion_table="pre_runs_annagnps_aa_sediment_erosion_ua_rr_total"
fi

# Store the global log path for failure reporting, but don't write to it routinely
GLOBAL_LOG_FILE="$LOG_FILE"
GLOBAL_FAILED_LOG="$FAILED_THUCS"

# Get the index of the directory to process from the job array
dir_index=$((SLURM_ARRAY_TASK_ID))

echo "$(date '+%Y-%m-%d %H:%M:%S') - Reading list of directories: $csv_file" # | tee -a "$LOG_FILE"

# Read the contents of dir_list.csv into an array
readarray -t dir_list < "$csv_file"

cd "$MINI_WATERSHEDS_DIR"

# Activate the virtual environment
source "$PYAGNPS_DIR/venv/bin/activate" || { echo "$(date '+%Y-%m-%d %H:%M:%S') - Failed to activate virtual environment" | tee -a "$LOG_FILE"; exit 1; }

# Check if the directory index is valid
if [ $dir_index -ge 0 ] && [ $dir_index -lt "${#dir_list[@]}" ]; then
    
    job_name=$(basename "${dir_list[$dir_index]}")

    # If the special flag is passed, create a local log inside the reach folder
    if [ "$USE_LOCAL_LOG" == "true" ]; then
        LOCAL_LOG_FILE="${dir_list[$dir_index]}/postproc_${thuc_id}.log"
    else
        # Fallback if run manually
        LOCAL_LOG_FILE="$LOG_FILE"
    fi
    
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Post Processing directory: $job_name" | tee -a "$LOCAL_LOG_FILE"
    
    cd "${dir_list[$dir_index]}" || exit 1

    # DO post processing here. --output_folder is the current working directory because we cdd there
    python -u "$PY_BASH_DIR/post_process_watershed_files_pre_runs.py" \
        --thuc_id "$thuc_id" \
        --note "$note" \
        --credentials "$path_to_db_credentials" \
        --output_folder "." \
        --save_method "$save_method" \
        --annagnps_aa_table "$annagnps_aa_table" \
        --aa_water_yield_table "$aa_water_yield_table" \
        --aa_sediment_yield_table "$aa_sediment_yield_table" \
        --aa_sediment_erosion_table "$aa_sediment_erosion_table" \
        --log_file "$LOCAL_LOG_FILE"
    exit_status=$?
    cd ..

    if [ $exit_status -eq 0 ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') - SUCCESS: ${job_name}" >> "$LOCAL_LOG_FILE"
        # Optional: You could rm the directory here if save_method=db
    else
        # ------------------------------------------------------------------------------
        # FAILURE HANDLING - Only touch global files if something breaks
        # ------------------------------------------------------------------------------
        if [ "$GLOBAL_LOG_FILE" != "/dev/null" ] && [ "$GLOBAL_LOG_FILE" != "USE_LOCAL_LOG" ]; then
             echo "${dir_list[$dir_index]}" >> "${GLOBAL_LOG_FILE%.*}_failed_post_process.log"
        fi
        
        if [ "$GLOBAL_FAILED_LOG" != "/dev/null" ]; then
             echo "$(date '+%Y-%m-%d %H:%M:%S'),$thuc_id,failed_postprocessing_${job_name}" >> "$GLOBAL_FAILED_LOG"
        fi
        
        echo "$(date '+%Y-%m-%d %H:%M:%S') - FAILED: ${job_name}" >> "$LOCAL_LOG_FILE"
    fi


else
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Invalid directory index: $dir_index, could not do post processing" | tee -a "$LOG_FILE"
    exit 1
fi

# Return to root directory
# cd "$MINI_WATERSHEDS_DIR"
# cd ..
# THUC_LOG_FILE="${LOG_FILE%.*}_${thuc_id}.log"
