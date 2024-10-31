#!/bin/bash

# Define function to handle arguments
parse_arguments() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --mini_watersheds_dir)
        MINI_WATERSHEDS_DIR="$2"
        shift 2
        ;;
      --py_bash_dir)
        PY_BASH_DIR="$2"
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

# Make a default value of the LOG_FILE in case it is not specified so that it doesn't log to a file
if [ -z "$LOG_FILE" ]; then
  LOG_FILE="/dev/null"
fi

# Check if required arguments are provided
if [ -z "$MINI_WATERSHEDS_DIR" ]; then
  echo "$(date '+%Y-%m-%d %H:%M:%S') - Error: Missing required argument: --mini_watersheds_dir" | tee -a "$LOG_FILE"
  exit 1
fi

if [ -z "$PY_BASH_DIR" ]; then
  echo "$(date '+%Y-%m-%d %H:%M:%S') - Error: Missing required argument: --py_bash_dir" | tee -a "$LOG_FILE"
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
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Using existing dir_list.csv file" | tee -a "$LOG_FILE"
    num_jobs=$(wc -l < "${MINI_WATERSHEDS_DIR}/dir_list.csv")
else
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Looking for jobs to submit..." | tee -a "$LOG_FILE"
    rm -f "${MINI_WATERSHEDS_DIR}/dir_list.csv"  # Remove the CSV file if it exists
    for subdir in "${MINI_WATERSHEDS_DIR}"/*/ ; do
        [[ -d $subdir ]] && echo "${subdir%/}" >> "${MINI_WATERSHEDS_DIR}/dir_list.csv" && ((num_jobs++))
    done
fi

echo "$(date '+%Y-%m-%d %H:%M:%S') - Found $num_jobs mini watersheds, submitting them in batches of size $batch_size" | tee -a "$LOG_FILE"
echo "$(date '+%Y-%m-%d %H:%M:%S') - Current directory: $PWD" | tee -a "$LOG_FILE"
echo "$(date '+%Y-%m-%d %H:%M:%S') - Mini watersheds directory: ${MINI_WATERSHEDS_DIR}" | tee -a "$LOG_FILE"
# Loop to submit jobs in batches
for ((start_index = 0; start_index < num_jobs; start_index += batch_size)); do
  end_index=$((start_index + batch_size - 1))
  
  # Ensure end_index doesn't exceed num_jobs
  [[ $end_index -ge $num_jobs ]] && end_index=$((num_jobs-1))

  echo "$(date '+%Y-%m-%d %H:%M:%S') - Submitting mini watersheds ${start_index} to ${end_index}..." | tee -a "$LOG_FILE"
  
  # Submit the job with the adjusted array range
  sbatch --oversubscribe \
         --requeue \
         --array="${start_index}-${end_index}" \
         --partition="$partition" \
         --job-name="anna_${start_index}-${end_index}" \
         --output="annagnps_${start_index}-${end_index}_%A_%a_%N.out" \
         "${PY_BASH_DIR}/run_annagnps_func_normal.sh" \
         --mini_watersheds_dir "$MINI_WATERSHEDS_DIR" \
         --log_file "$LOG_FILE" & #\
        #  --pyagnps_dir "$PYAGNPS_DIR" & # Not necessary but here it is anyway in case some python script is needed later
  sleep 5

  num_running_jobs=$(squeue --noheader | wc -l)

#   # Optional delay between batch submissions
  # Check the number of currently running jobs
  iteration_count=0
  while [[ $num_running_jobs -gt 30 ]]; do
    ((iteration_count++))
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Too many jobs already running, sleeping and retrying later... ($iteration_count/$maxiter)" | tee -a "$LOG_FILE"
    
    if (( iteration_count % 10 == 0 )); then
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Running release_requeue.sh script after $iteration_count iterations." | tee -a "$LOG_FILE"
        bash "${PY_BASH_DIR}/release_requeue.sh"
    fi

    sleep 5
    num_running_jobs=$(squeue --noheader | wc -l)

    if [[ $iteration_count -eq $maxiter ]]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Error: Maximum iterations of $maxiter reached in the while loop." | tee -a "$LOG_FILE"
        exit 1
    fi
  done

done
