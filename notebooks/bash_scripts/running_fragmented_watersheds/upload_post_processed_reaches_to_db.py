from pyagnps import annagnps, aims
from pyagnps.utils import log_to_file, upsert_dataframe
from pathlib import Path
import pandas as pd
import sys
import argparse
import traceback
import gc  # Garbage Collector interface

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--credentials',                            type=str, help='Path to the credentials JSON file')
    parser.add_argument('--post_processing_dir',                    type=str, help='Path to the output folder to store the generated files', required=True)
    parser.add_argument('--file_batch_size',                        type=int, help='Number of files to merge in memory at once', default=500)
    parser.add_argument('--delete_post_processed_files_on_success', type=str, help='Delete post processed files if successful upload', default="false")
    parser.add_argument('--log_file',                               type=str, help='Path to the log file',                         default="upload_post_processed_reaches_to_db.log")

    args = parser.parse_args()

    credentials            = Path(args.credentials)
    post_processing_dir    = Path(args.post_processing_dir)
    log_file_path          = Path(args.log_file)
    file_batch_size        = args.file_batch_size

    # Robust boolean parsing
    delete_post_processed_files = args.delete_post_processed_files_on_success.lower() in ["true", "yes", "y", "oui", "1"]

    # Map tables to their unique constraints (Required for Upsert)
    db_table_unique_columns = {
            'pre_runs_annagnps_aa': ['thuc_id', 'cell_id', 'note'],
            'pre_runs_annagnps_aa_sediment_erosion_ua_rr_total': ['thuc_id', 'cell_id', 'description', 'note'],
            'pre_runs_annagnps_aa_sediment_yield_ua_rr_total': ['thuc_id', 'cell_id', 'description', 'note'],
            'pre_runs_annagnps_aa_water_yield_ua_rr_total': ['thuc_id', 'cell_id', 'description', 'note'],
        }

    try:
        log_to_file(log_file_path, f"Starting upload from {post_processing_dir}", add_timestamp=True)

        db_url = aims.create_db_url_object(credentials)

        # Connection args to prevent timeouts on long running uploads
        connect_args = {
            'keepalives': 1,
            'keepalives_idle': 60,
            'keepalives_interval': 10,
            'keepalives_count': 5
        }
        
        # Engine with pre-ping to handle connection drops gracefully
        engine = aims.create_engine(db_url, connect_args=connect_args,
                                            pool_pre_ping=True,
                                            pool_recycle=300)

        # Iterate over each table type folder
        for table_folder in post_processing_dir.iterdir():
            if not table_folder.is_dir():
                continue

            table_name = table_folder.name
            
            # Sort files to ensure deterministic order (important for restartability)
            parquet_files = sorted(list(table_folder.glob('*.parquet')))

            if not parquet_files:
                continue
            
            log_to_file(log_file_path, f"Found {len(parquet_files)} files for table {table_name}", add_timestamp=True)

            # --- DOUBLE BATCHING STRATEGY ---
            # Outer Loop: Process files in manageable chunks (e.g., 500 files at a time)
            # This ensures we never load all 50,000 files into RAM at once.
            for i in range(0, len(parquet_files), file_batch_size):
                
                batch_files = parquet_files[i : i + file_batch_size]
                df_list = []
                
                # 1. Read files into memory
                for file in batch_files:
                    try:
                        df = pd.read_parquet(file)
                        df_list.append(df)
                    except Exception as e:
                        log_to_file(log_file_path, f"Corrupt file skipped {file.name}: {e}")

                if not df_list:
                    continue

                # 2. Concatenate into one DataFrame for this batch
                df_batch = pd.concat(df_list)
                
                # 3. Upload (Upsert)
                # The upsert_dataframe function handles row-chunking internally
                try:
                    unique_columns = db_table_unique_columns.get(table_name, None)
                    if unique_columns is None:
                        log_to_file(log_file_path, f"Warning: No unique columns for {table_name}, skipping.")
                        break

                    upsert_dataframe(engine, df_batch, table_name, unique_columns=unique_columns)

                    # 4. Incremental Cleanup (Delete files ONLY after successful upload)
                    if delete_post_processed_files:
                        for file in batch_files:
                            try:
                                file.unlink()
                            except FileNotFoundError:
                                pass # Already gone

                    log_to_file(log_file_path, f"Uploaded batch {i//file_batch_size + 1} for {table_name} ({len(df_batch)} rows)", add_timestamp=True)

                except Exception as e:
                    # Log error and Exit.
                    # Because we use incremental cleanup, we can just restart the job later 
                    # and it will pick up exactly where it left off.
                    log_to_file(log_file_path, f"Critical Error uploading batch for {table_name}: {e}", add_timestamp=True)
                    sys.exit(1)
                
                # 5. Force Memory Cleanup
                # Python doesn't always release memory immediately. 
                # In high-throughput jobs, this prevents "memory creep".
                del df_batch
                del df_list
                gc.collect()

        log_to_file(log_file_path, "Upload process completed successfully.", add_timestamp=True)
        sys.exit(0)

    except Exception as e:
        log_to_file(log_file_path, f"Global Error: {e}\n{traceback.format_exc()}", add_timestamp=True)
        sys.exit(1)


if __name__ == '__main__':
    main()