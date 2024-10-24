#!/bin/bash

# Define function to handle arguments
parse_arguments() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --mini_watsheds_dir)
        MINI_WATERSHEDS_DIR="$2"
        shift 2
        ;;
      --batch_size)
        batch_size="$2"
        shift 2
        ;;
      --maxiter)
        maxiter="$2"
        shift 2
        ;;
      --partition)
        partition="$2"
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
if [ -z "$MINI_WATERSHEDS_DIR" ]; then
  echo "Error: Missing required argument: --mini_watsheds_dir"
  exit 1
fi

# Set defaults for optional arguments
if [ -z "$batch_size" ]; then
  batch_size=1000
fi

if [ -z "$maxiter" ]; then
  maxiter=1000
fi

if [ -z "$partition" ]; then
  partition="aims-highperf-oversubscribe,aims-default-oversubscribe"
fi

# Calculate the total number of jobs based on directory count
num_jobs=0

if [ -f "${MINI_WATERSHEDS_DIR}/dir_list.csv" ]; then
    echo "Using existing dir_list.csv file"
    num_jobs=$(wc -l < "${ROOT_DIR}/dir_list.csv")
else
    echo "Looking for jobs to submit..."
    rm -f "${MINI_WATERSHEDS_DIR}/dir_list.csv"  # Remove the CSV file if it exists
    for subdir in "${ROOT_DIR}"/*/ ; do
        [[ -d $subdir ]] && echo "${subdir%/}" >> "${ROOT_DIR}/dir_list.csv" && ((num_jobs++))
    done
fi

echo "Found $num_jobs jobs, submitting them in batches of size $batch_size"

# Loop to submit jobs in batches
for ((start_index = 0; start_index < num_jobs; start_index += batch_size)); do
  end_index=$((start_index + batch_size - 1))
  
  # Ensure end_index doesn't exceed num_jobs
  [[ $end_index -ge $num_jobs ]] && end_index=$((num_jobs-1))

  echo "Submitting jobs ${start_index} to ${end_index}..."
  
  # Submit the job with the adjusted array range
  sbatch --oversubscribe \
         --array="${start_index}-${end_index}" \
         --partition="$partition" \
         --job-name="anna_${start_index}-${end_index}" \
         --output="annagnps_${start_index}-${end_index}_%A_%a_%N.out" \
         ./run_annagnps_func_normal.sh "$MINI_WATERSHEDS_DIR" &
  sleep 5

  num_running_jobs=$(squeue --noheader | wc -l)

#   # Optional delay between batch submissions
  # Check the number of currently running jobs
  iteration_count=0
  while [[ $num_running_jobs -gt 30 ]]; do
    ((iteration_count++))
    echo "Too many jobs already running, sleeping and retrying later... ($iteration_count/$maxiter)"
    sleep 5
    num_running_jobs=$(squeue --noheader | wc -l)
    # Exit with an error message if the maximum iterations are reached
    if [[ $iteration_count -eq $maxiter ]]; then
        echo "Error: Maximum iterations of $maxiter reached in the while loop."
        exit 1
    fi
  done

done
