import os
import platform
import subprocess
import sys
from urllib.request import urlretrieve

def get_gdal_wheel_url(python_version, gdal_version):
    arch = 'win_amd64' if platform.architecture()[0] == '64bit' else 'win32'
    base_url = "https://download.lfd.uci.edu/pythonlibs/archived/"

    wheel_name = f"GDAL-{gdal_version}-cp{python_version.replace('.', '')}-cp{python_version.replace('.', '')}-{arch}.whl"
    return base_url + wheel_name

def get_gdal_version_linux():
    try:
        gdal_config = os.environ.get('GDAL_CONFIG', 'gdal-config')
        return subprocess.check_output([gdal_config, '--version']).decode('utf-8').strip()
    except (OSError, subprocess.CalledProcessError):
        raise RuntimeError('GDAL must be installed and gdal-config executable must be in your PATH')


def install_gdal():
    if sys.platform.startswith('linux'):
        try:
            gdal_version_output = subprocess.check_output(['gdal-bin', '--version']).decode().strip()
            gdal_version = gdal_version_output.split(',')[0].split()[1]
            print(f"Detected GDAL version: {gdal_version}")
            subprocess.check_call([sys.executable, '-m', 'pip', 'install', f'GDAL=={gdal_version}'])
        except FileNotFoundError:
            print("gdal-bin is not installed. Please install gdal-bin and libgdal-dev on your system.")
            sys.exit(1)
        except subprocess.CalledProcessError as e:
            print(f"Error checking GDAL version: {e}")
            sys.exit(1)
    elif sys.platform.startswith('win'):
        python_version = f"{sys.version_info.major}.{sys.version_info.minor}"
        gdal_version = "3.4.3"  # Specify a default or fallback GDAL version if needed
        gdal_url = get_gdal_wheel_url(python_version, gdal_version)
        wheel_path, _ = urlretrieve(gdal_url)
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', wheel_path])
    else:
        raise OSError("Unsupported operating system")

if __name__ == "__main__":
    print('Installing GDAL...')
    install_gdal()
    print("Successfully installed GDAL, try for yourself: python -c 'from osgeo import gdal'")
