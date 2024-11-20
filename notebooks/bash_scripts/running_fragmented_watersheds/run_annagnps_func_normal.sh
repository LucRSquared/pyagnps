#!/bin/bash

# Set the root directory (can be changed before script execution)

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
      --force_simulate)
        force_simulate="$2"
        shift 2
        ;;
      --csv_file)
        csv_file="$2"
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

# Set default value for MINI_WATERSHEDS_DIR if not provided
if [ -z "$MINI_WATERSHEDS_DIR" ]; then
    MINI_WATERSHEDS_DIR=$(pwd)
fi

# Set default value for csv_file if not provided
if [ -z "$csv_file" ]; then
    csv_file="${MINI_WATERSHEDS_DIR}/dir_list.csv"
fi

if [ -z "$force_simulate" ]; then
    force_simulate="false"
fi

if [-z "$thuc_id" ]; then
    thuc_id=""
fi

# if [ -z "$PYAGNPS_DIR" ]; then
#   PYAGNPS_DIR="/aims-nas/luc/code/pyagnps/"
# fi

# Get the index of the directory to process from the job array
dir_index=$((SLURM_ARRAY_TASK_ID))

# Read the contents of dir_list.csv into an array
readarray -t dir_list < "$csv_file"

cd "$MINI_WATERSHEDS_DIR"

# Check if the directory index is valid
if [ $dir_index -ge 0 ] && [ $dir_index -lt "${#dir_list[@]}" ]; then
    job_name=$(basename "${dir_list[$dir_index]}")
    echo "Processing directory: $job_name"
    cd "${dir_list[$dir_index]}" || exit 1

    # If force_simulate == true then run annagnps otherwse check if AnnAGNPS.log exists
    if [ "$force_simulate" == "true" ] || [ ! -e ./AnnAGNPS.log ]; then
        # Run annagnps
        if ! annagnps; then
            ERROR_LOG_FILE="${LOG_FILE%.*}_failed_process.log"
            echo "${dir_list[$dir_index]}" | tee -a "$ERROR_LOG_FILE"

            FAILED_THUCS="${LOG_FILE%/*}/failed_thucs.csv"
            echo "$(date '+%Y-%m-%d %H:%M:%S'),$thuc_id,failed_processing" | tee -a "$FAILED_THUCS"
        fi
    else
        # echo "Skipping directory: $job_name" | tee -a "$LOG_FILE"
        cd ..
        exit 0
    fi
    
    # Run annagnps
    # annagnps

else
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Invalid directory index: $dir_index" | tee -a "$LOG_FILE"
    exit 1
fi

cd ..

