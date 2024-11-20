# CLI commands to download and process NLDAS-2 data

## Download NLDAS-2 hourly netCDF files from GESDISC

Once installed the CLI command `download-nldas2` can be used to download NLDAS-2 netCDF files from GESDISC.

```bash
download-nldas2 --from_date 2024-01-01 --to_date 2024-03-01 --products NLDAS_FORA0125_H.2.0 --output_dir Path/to/NLDAS2/ --username myusername --password mypassword
```

## Aggregate NLDAS-2 netCDF files to daily netCDF files

The CLI command `aggregate-nldas2` can be used to aggregate NLDAS-2 netCDF files to daily netCDF files.

```bash
aggregate-nldas2 --product NLDAS_FORA0125_H.2.0 --files_dir path/to/netcdf_hourly/ --output_dir_agg_netcdf path/to/netcdf_daily/ --from_date 2020-01-02 --to_date 2020-01-04 --chunk_frac_tot 0.2
```

The chunk_frac_tot parameter is downloads the data in batches of that specified fraction of the total data. Adjust according to your specific configuration (lower value if downloads fail).
