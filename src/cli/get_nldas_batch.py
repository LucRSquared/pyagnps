import argparse
from pathlib import Path

import requests

from dateutil.relativedelta import relativedelta
import datetime


from pyagnps.utils import month_difference, download_files_from_url

from pyagnps.constants import _BASE_URL_NLDAS, _NLDAS_PRODUCTS_20, _NLDAS_PRODUCTS_002

# Example use from command line
# download-nldas2 --from_date 2024-01-01 --to_date 2024-03-1 --products NLDAS_FORA0125_H.002 --output_dir Path/to/NLDAS2/ --username myusername --password mypassword

class SessionWithHeaderRedirection(requests.Session):
    """
    adapted from
    https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
    """

    AUTH_HOST = 'urs.earthdata.nasa.gov'

    def __init__(self, username, password):
        super().__init__()
        self.auth = (username, password)

    # Overrides from the library to keep headers when redirected to or from
    # the NASA auth host.

    def rebuild_auth(self, prepared_request, response):
        headers = prepared_request.headers
        url = prepared_request.url

        if 'Authorization' in headers:
            original_parsed = requests.utils.urlparse(response.request.url)
            redirect_parsed = requests.utils.urlparse(url)

            if (original_parsed.hostname != redirect_parsed.hostname) and \
                    redirect_parsed.hostname != self.AUTH_HOST and \
                    original_parsed.hostname != self.AUTH_HOST:
                del headers['Authorization']
        return
    
def get_nldas_url_list(from_date, to_date, product='NLDAS_FORA0125_H.2.0', base_url=_BASE_URL_NLDAS, product_dict=_NLDAS_PRODUCTS_20):
    """
    create a list of URLs to download based on the date range
    and the product
    Appends the xml file url for each GRIB/netCDF file
    """
    date_diff = to_date - from_date
    target_dates = []
    download_urls = []

    if '.002' in product:
        # GRIB files
        product_dict = _NLDAS_PRODUCTS_002
        file_extension = '.002.grb'

    elif '.2.0' in product:
        # NetCDF files
        product_dict = _NLDAS_PRODUCTS_20
        file_extension = '.020.nc'

    time_type = product_dict[product]['type']
    fileroot = product_dict[product]['fileroot']

    if time_type == 'hourly':
        for i in range(date_diff.days + 1):
            target_dates.append((from_date + datetime.timedelta(days=i)))

        for i in range(len(target_dates)):
            for j in range(0, 24):
                download_urls.append((f'{_BASE_URL_NLDAS}/{product}/{target_dates[i].year}/'
                                      f'{target_dates[i].timetuple().tm_yday:03d}/'
                                      f'{fileroot}{target_dates[i].year}'
                                      f'{target_dates[i].month:02d}'
                                      f'{target_dates[i].day:02d}.{j:02d}00{file_extension}'
                                      ))
                download_urls.append((f'{_BASE_URL_NLDAS}/{product}/{target_dates[i].year}/'
                                      f'{target_dates[i].timetuple().tm_yday:03d}/'
                                      f'{fileroot}{target_dates[i].year}'
                                      f'{target_dates[i].month:02d}'
                                      f'{target_dates[i].day:02d}.{j:02d}00{file_extension}.xml'
                                      ))

    elif time_type == 'monthly':
        for i in range(month_difference(from_date, to_date) + 1):
            target_dates.append((from_date + relativedelta(months=i)))

        for i in range(len(target_dates)):
            download_urls.append((f'{_BASE_URL_NLDAS}/{product}/{target_dates[i].year}/'
                                  f'{fileroot}{target_dates[i].year}'
                                  f'{target_dates[i].month:02d}{file_extension}'
                                  ))
            download_urls.append((f'{_BASE_URL_NLDAS}/{product}/{target_dates[i].year}/'
                                  f'{fileroot}{target_dates[i].year}'
                                  f'{target_dates[i].month:02d}{file_extension}.xml'
                                  ))

    elif time_type == 'monthly_climatology':
        print(f'Data product is {product} which is a monthly climatology averaged dataset, from_date and to_date will be ignored and all months will be downloaded')
        for i in range(1, 12 + 1):
            download_urls.append((f'{_BASE_URL_NLDAS}/{product}/'
                                  f'{fileroot}{i:02d}{file_extension}'
                                  ))
            download_urls.append((f'{_BASE_URL_NLDAS}/{product}/'
                                  f'{fileroot}{i:02d}{file_extension}.xml'
                                  ))

    return download_urls

def main():
    # Parse the command line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('--from_date', '-fd',
                        help="From Date/Time Format (inclusive) YYYY-MM-DD",
                        type=datetime.date.fromisoformat)

    parser.add_argument('--to_date', '-td',
                        help="To Date/Time (inclusive) Format YYYY-MM-DD",
                        type=datetime.date.fromisoformat)

    parser.add_argument('--products', '-prd',
                        help=("The NLDAS product to download, default is NLDAS_NOAH0125_H.002"
                              "see also: 'https://hydro1.gesdisc.eosdis.nasa.gov/data/NLDAS/'"))

    parser.add_argument('--username', '-u',
                        help=("Username to access data from "
                            "https://disc.gsfc.nasa.gov/data-access"))

    parser.add_argument('--password', '-p',
                        help=("Password to access data from "
                            "https://disc.gsfc.nasa.gov/data-access"))

    parser.add_argument('--output_dir', '-o',
                        help=("Output Directory "
                            "(default ./nldas-<from_date>.to.<to_date>"),
                        required=False)

    parser.add_argument('--hourly', '-H', default=False,
                        action='store_const', const=True, help='Download Hourly Forcing Data')

    parser.add_argument('--monthly', '-M', default=False,
                        action='store_const', const=True, help='Download Monthly Forcing Data')

    args = parser.parse_args()
    data_sets = []

    # # TODO: verify user meant to use both hourly and daily
    # if args.hourly:
    #     data_sets.append(FORCINGS_DATA_SETS['hourly'])
    # if args.monthly:
    #     data_sets.append(FORCINGS_DATA_SETS['monthly'])

    # validate the arguments
    if args.to_date < args.from_date:
        raise ValueError("from_date must be less than or equal to to_date")

    if args.output_dir is None:
        output_dir = Path.cwd() / f'nldas-{args.from_date}.to.{args.to_date}'
    else:
        output_dir = Path(args.output_dir)

    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)

    # create session with the user credentials that will be used to authenticate
    # access to the data
    session = SessionWithHeaderRedirection(args.username, args.password)

    # get the files for each dataset type requested
    products = args.products.split(',')

    for data_set in products:
        urls = get_nldas_url_list(args.from_date, args.to_date, data_set)
        file_count = len(urls)
        # download_size = GRIB_FILE_SIZE * file_count

        print(f'Proceeding to download {file_count} files...')
        # print(f'Downloading approx. {download_size} MB (assume {GRIB_FILE_SIZE} MB per file)')

        download_files_from_url(session, urls, output_dir)

        print(f'Downloaded files to {output_dir}')


if __name__ == '__main__':
    main()