#!/bin/bash
# monitor_held_jobs.sh

# Parse input arguments
parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --log_file)
                LOG_FILE="$2"
                shift
                ;;
            *)
                echo "Unknown option: $1"
                exit 1
                ;;
        esac
        shift
    done
}

# Define log file if argument is not provided
if [ -z "$LOG_FILE" ]; then
    LOG_FILE="/dev/null"
fi

# Function to log messages with timestamp
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOG_FILE"
}

# Function to handle a single held job
handle_held_job() {
    local jobid=$1
    local reason=$2
    
    log_message "Found held job $jobid with reason: $reason"
    
    # First release
    if scontrol release $jobid; then
        log_message "Successfully released job $jobid"
        sleep 2  # Give it a moment
        
        # Then requeue
        if scontrol requeue $jobid; then
            log_message "Successfully requeued job $jobid"
        else
            log_message "Failed to requeue job $jobid"
        fi
    else
        log_message "Failed to release job $jobid"
    fi
}

# Function to check for held jobs
check_held_jobs() {
    squeue -h -t PD -o "%i %r" | \
    grep -E "launch failed requeued held|launch failed held|dependency never satisfied|QOSMaxGRPJobsLimit|AssocGrpCPUMinutesLimit" || true
}

# Initialize counter for no-jobs-found iterations
no_jobs_counter=0
check_interval=5  # Check every 5 seconds
max_empty_checks=$((30 / check_interval))  # Number of checks to reach 30 seconds

log_message "Starting monitoring of held jobs..."

while true; do
    # Store the held jobs in a variable
    held_jobs=$(check_held_jobs)
    
    if [ -z "$held_jobs" ]; then
        # No held jobs found
        ((no_jobs_counter++))
        log_message "No held jobs found. Check $no_jobs_counter of $max_empty_checks"
        
        if [ $no_jobs_counter -ge $max_empty_checks ]; then
            log_message "No held jobs found for 30 seconds. Terminating monitor."
            exit 0
        fi
    else
        # Reset counter if jobs are found
        no_jobs_counter=0
        
        # Process each held job
        echo "$held_jobs" | while read -r jobid reason; do
            handle_held_job "$jobid" "$reason"
        done
    fi
    
    sleep $check_interval
done