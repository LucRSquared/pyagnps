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

# Get the index of the directory to process from the job array
dir_index=$((SLURM_ARRAY_TASK_ID))

# Read the contents of dir_list.csv into an array
readarray -t dir_list < "$csv_file"

cd "$MINI_WATERSHEDS_DIR"

# Check if the directory index is valid
if [ $dir_index -ge 0 ] && [ $dir_index -lt "${#dir_list[@]}" ]; then
    job_name=$(basename "${dir_list[$dir_index]}")
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Post Processing directory: $job_name" #| tee -a "$LOG_FILE"
    cd "${dir_list[$dir_index]}" || exit 1
    
    # DO post processing here. --output_folder is the current working directory because we cdd there
    python -u "$PY_BASH_DIR/post_process_watershed_files_pre_runs.py" \
        --thuc_id "$thuc_id" \
        --credentials "$path_to_db_credentials" \
        --output_folder "." \
        --annagnps_aa_table "$annagnps_aa_table" \
        --aa_water_yield_table "$aa_water_yield_table" \
        --aa_sediment_yield_table "$aa_sediment_yield_table" \
        --aa_sediment_erosion_table "$aa_sediment_erosion_table" & # \
        # --log_file "$LOG_FILE" &

    # cd "$MINI_WATERSHEDS_DIR" # This didn't work because MINI_WATERSHEDS_DIR = ./mini_watersheds which won't work once we are in the reach folder
    cd ..

    # If the python code finished successfully then delete the current directory
    if [ $? -eq 0 ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Post processing finished successfully, deleting directory ${dir_list[$dir_index]}" #| tee -a "$LOG_FILE"
        rm -rf "${dir_list[$dir_index]}"
    else
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Post processing failed, not deleting directory ${dir_list[$dir_index]}" | tee -a "$LOG_FILE"
    fi


else
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Invalid directory index: $dir_index, could not do post processing" | tee -a "$LOG_FILE"
    exit 1
fi

# Return to root directory
# cd "$MINI_WATERSHEDS_DIR"
# cd ..

