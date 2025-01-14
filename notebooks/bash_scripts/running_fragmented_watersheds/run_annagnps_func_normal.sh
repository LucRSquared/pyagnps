#!/bin/bash

echo "Hello with SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID" >&2

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
      --failed_log_file)
        FAILED_THUCS="$2"
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

# Make a default value of the LOG_FILE in case it is not specified so that it doesn't log to a file
if [ -z "$LOG_FILE" ]; then
  LOG_FILE="/dev/null"
fi

# Make a default value of the FAILED_THUCS in case it is not specified so that it doesn't log to a file
if [ -z "$FAILED_THUCS" ]; then
  FAILED_THUCS="/dev/null"
fi

# Set default value for csv_file if not provided
if [ -z "$csv_file" ]; then
    csv_file="${MINI_WATERSHEDS_DIR}/dir_list.csv"
fi

if [ -z "$force_simulate" ]; then
    force_simulate="false"
fi

if [ -z "$thuc_id" ]; then
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

echo "Hello from dir_index=$dir_index and reach ${dir_list[$dir_index]}, the array has length ${#dir_list[@]}" | tee -a "$LOG_FILE"


echo "[DEBUG] Starting job with:" >&2
echo "  dir_index: $dir_index" >&2
echo "  array length: ${#dir_list[@]}" >&2
echo "  current dir: $(pwd)" >&2
echo "  target dir: ${dir_list[$dir_index]:-NONE}" >&2

# Check if the directory index is valid
if [ $dir_index -ge 0 ] && [ $dir_index -lt "${#dir_list[@]}" ]; then
    job_name=$(basename "${dir_list[$dir_index]}")
    echo "Processing directory: $job_name" | tee -a "$LOG_FILE"
    
    if [ ! -d "${dir_list[$dir_index]}" ]; then
        echo "ERROR: Directory ${dir_list[$dir_index]} does not exist" >&2
        exit 1
    fi
    
    if ! cd "${dir_list[$dir_index]}"; then
        echo "ERROR: Could not cd to directory: $job_name" >&2
        exit 1
    fi

    # If force_simulate == true then run annagnps otherwise check if AnnAGNPS.log exists
    if [ "$force_simulate" = "true" ] || [ ! -e ./AnnAGNPS.log ]; then
        # Verify annagnps command exists
        if ! command -v annagnps &> /dev/null; then
            echo "ERROR: annagnps command not found" >&2
            exit 1
        fi
        
        # Run annagnps
        if ! annagnps; then
            ERROR_LOG_FILE="${LOG_FILE%.*}_failed_process.log"
            echo "${dir_list[$dir_index]}" | tee -a "$ERROR_LOG_FILE"
            echo "$(date '+%Y-%m-%d %H:%M:%S'),$thuc_id,failed_processing" | tee -a "$FAILED_THUCS"
            exit 1
        fi
    else
        echo "Skipping directory: $job_name (AnnAGNPS.log exists)" | tee -a "$LOG_FILE"
    fi
else
    echo "ERROR: Invalid directory index: $dir_index" >&2
    exit 1
fi

cd ".." || exit 1


# # Check if the directory index is valid
# if [ $dir_index -ge 0 ] && [ $dir_index -lt "${#dir_list[@]}" ]; then
#     job_name=$(basename "${dir_list[$dir_index]}")
#     echo "Processing directory: $job_name" | tee -a "$LOG_FILE"
#     cd "${dir_list[$dir_index]}" || { echo "Could not cd to directory: $job_name" | tee -a "$LOG_FILE" && exit 1 }

#     # If force_simulate == true then run annagnps otherwse check if AnnAGNPS.log exists
#     if [ "$force_simulate" == "true" ] || [ ! -e ./AnnAGNPS.log ]; then
#         # Run annagnps
#         if ! annagnps; then
#             ERROR_LOG_FILE="${LOG_FILE%.*}_failed_process.log"
#             echo "${dir_list[$dir_index]}" | tee -a "$ERROR_LOG_FILE"

#             echo "$(date '+%Y-%m-%d %H:%M:%S'),$thuc_id,failed_processing" | tee -a "$FAILED_THUCS"
#         fi
#     else
#         echo "Skipping directory: $job_name" | tee -a "$LOG_FILE"
#         cd ..
#         exit 0
#     fi
    
#     # Run annagnps
#     # annagnps

# else
#     echo "$(date '+%Y-%m-%d %H:%M:%S') - Invalid directory index: $dir_index" | tee -a "$LOG_FILE"
#     exit 1
# fi

# cd ..

