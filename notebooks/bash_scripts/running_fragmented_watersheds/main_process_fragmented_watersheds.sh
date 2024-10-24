#!/bin/sbash

# STRUCTURE
# ROOT_DIR
# ├── thuc_list_to_run.csv
# ├── thuc_0001
# │   ├── (bunch of regular input files)
# │   └── mini_watersheds
# │       ├── dir_list.csv
# │       ├── reach_000000002
# │       ├── reach_000000003
# │       ├── reach_000000004
# │       └── etc.
# │
# ├──  thuc_0002
#     ├── (bunch of regular input files)
#     └── mini_watersheds
#         ├── dir_list.csv
#         ├── reach_000000002
#         ├── reach_000000003
#         ├── reach_000000004
#         └── etc.


# Set the root directory (can be changed before script execution). This directory should contain thuc_list_to_run.csv and a directory for each watershed that needs to be processed
ROOT_DIR="."  # Default value (current directory, assumes this script is located in ROOT_DIR)

# Set boolean parameters to generate and fragment, simulate, post-process results (and populate them to the database)
generate_main_files="true"
fragment_watershed="true"
simulate_thuc="true"
post_process="true"

path_to_db_credentials="/aims-nas/luc/code/pyagnps/inputs/db_credentials.json"

path_to_nldas2_centroids="/aims-nas/data/datasets/CLIMATE/NLDAS2/NLDAS2_GRID_CENTROIDS_epsg4326.gpkg"
path_to_scs_storm_types="/aims-nas/data/datasets/TR-55/scs_storm_types.gpkg"
path_to_precip_zones="/aims-nas/data/datasets/RUSLE2/Climate/precip_zones_RUSLE2_cleaned_manually_extrapolated_pchip_linear_US_units.gpkg"

# Batch size for job submissions
batch_size=1000
maxiter=1000
num_processes=32

# Loop through thuc list
# For each thuc:
#   - check if AnnAGNPS.fil exists (if yes it means the regular input files were already generated)
#   - check if dir_list.csv exists within the mini_watersheds directory
#   - if not generate the files and fragment them using the python script "fragment_annagnps_into_mini_watersheds.py" with parameters startDate, endDate, thuc_ID
#   - if yes, use the existing dir_list.csv
#   - check that AnnAGNPS.log does not exist
#   - if not, run the simulation using the distributed method
#   - read results and upload to db
#   - cleanup files that are not needed
#   - log results



# Calculate the total number of jobs based on directory count
num_jobs=0

if [ -f "${ROOT_DIR}/thuc_list_to_run.csv" ]; then
    echo "Reusing existing thuc_list_to_run.csv file"
    num_jobs=$(wc -l < "${ROOT_DIR}/thuc_list_to_run.csv")

echo "Found $num_jobs jobs, submitting them in batches of size $batch_size"

# Loop to process thucs one by one
for ((thuc_index = 0; thuc_index < num_jobs; thuc_index += 1)); do
  thuc_id=$(sed -n "${thuc_index}p" "${ROOT_DIR}/thuc_list_to_run.csv")
  echo "Processing thuc $thuc_id"

  # Check if regular input files exist (AnnAGNPS.fil) exists
  if [ -f "${ROOT_DIR}/thuc_${thuc}/AnnAGNPS.fil" ] || [ "$generate_main_files" = "false" ]; then
      echo "Skipping thuc $thuc_id, input files already exist"
      generate_main_files="false"
  else
      generate_main_files="true"
  fi

  # Check if dir_list.csv exists
  if [ -f "${ROOT_DIR}/thuc_${thuc}/mini_watersheds/dir_list.csv" ] || [ "$fragment_watershed" = "false" ]; then
      fragment_watershed="false"
  else
      fragment_watershed="true"
  fi

  #---------------------------------------

  if [ "$generate_main_files" == "true"]; then
      echo "Generating input files for thuc $thuc_id and fragmenting it into mini watersheds"
      ./generate_annagnps_files_full_and_fragment_watershed.sh \
        --simulation_dir "${ROOT_DIR}/thuc_${thuc}" \
        --credentials "$path_to_db_credentials" \
        --nldas2_centroids "$path_to_nldas2_centroids" \
        --scs_storm_types "$path_to_scs_storm_types" \
        --precip_zones "$path_to_precip_zones" \
        --start_date "$start_date" \
        --end_date "$end_date" \
        --climate_method "$climate_method" \
        --climate_table "$climate_table" \
        --thuc_id "$thuc_id" \
        --reach_id "$reach_id" \
        --generate_main_files "$generate_main_files" \
        --fragment_watershed "$fragment_watershed" \
        --num_processes "$num_processes"
  fi
     
  echo "Done generating input files for thuc $thuc_id"
  
  if [ "$simulate_thuc" == "true"]; then
      echo "Simulating thuc $thuc_id using distributed method"
      cd "${ROOT_DIR}/thuc_${thuc}"
      ./process_one_fragmented_watershed_run_sim.sh \
        --mini_watersheds_dir "./mini_watersheds" \
        --batch_size "$batch_size" \
        --maxiter "$maxiter" \
        --partition aims-highperf-oversubscribe,aims-default-oversubscribe
      cd "${ROOT_DIR}"
  fi

  echo "Done simulating thuc $thuc_id"

  if [ "$post_process" == "true"]; then
      echo "Post processing thuc $thuc_id"
      ./post_process_one_fragmented_watershed.sh \
        --simulation_dir "${ROOT_DIR}/thuc_${thuc}" \
        --credentials "$path_to_db_credentials" \
        --nldas2_centroids "$path_to_nldas2_centroids" \
        --scs_storm_types "$path_to_scs_storm_types" \
        --precip_zones "$path_to_precip_zones" \
        --start_date "$start_date" \
        --end_date "$end_date" \
        --climate_method "$climate_method" \
        --climate_table "$climate_table" \
        --thuc_id "$thuc_id" \
        --generate_main_files "$generate_main_files"
  fi
done



for ((start_index = 0; start_index < num_jobs; start_index += batch_size)); do
  end_index=$((start_index + batch_size - 1))
  
  # Ensure end_index doesn't exceed num_jobs
  [[ $end_index -ge $num_jobs ]] && end_index=$((num_jobs-1))

  echo "Submitting jobs ${start_index} to ${end_index}..."
  
  # Submit the job with the adjusted array range
  sbatch --oversubscribe --array="${start_index}-${end_index}" --partition=aims-highperf,aims-default ./run_annagnps_func_normal_csv.sh "$ROOT_DIR" &
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
