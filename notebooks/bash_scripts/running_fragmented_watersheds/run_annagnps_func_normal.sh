#!/bin/bash

# Set the root directory (can be changed before script execution)
MINI_WATERSHEDS_DIR="$1"
csv_file = "$2"

# Set default value for MINI_WATERSHEDS_DIR if not provided
if [ -z "$MINI_WATERSHEDS_DIR" ]; then
    MINI_WATERSHEDS_DIR=$(pwd)
fi

# Set default value for csv_file if not provided
if [ -z "$csv_file" ]; then
    csv_file="${MINI_WATERSHEDS_DIR}/dir_list.csv"
fi

# Get the index of the directory to process from the job array
dir_index=$((SLURM_ARRAY_TASK_ID))

# Read the contents of dir_list.csv into an array
readarray -t dir_list < csv_file

cd "$MINI_WATERSHEDS_DIR"

# Check if the directory index is valid
if [ $dir_index -ge 0 ] && [ $dir_index -lt "${#dir_list[@]}" ]; then
    job_name=$(basename "${dir_list[$dir_index]}")
    echo "Processing directory: $job_name"
    cd "${dir_list[$dir_index]}" || exit
    
    # Run annagnps
    annagnps

else
    echo "Invalid directory index: $dir_index"
    exit 1
fi

# # Using [[ ]]
# if [[ -e ./AnnAGNPS.log ]]; then
#     echo "Skipping directory: $job_name"
#     exit 1
# else
#     echo "Processing directory: $job_name"
# fi

# Return to root directory
cd "$MINI_WATERSHEDS_DIR"

