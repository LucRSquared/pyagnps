# PyAGNPS

PyAGNPS is a series of scripts and notebooks that helps in the creation of input files for TopAGNPS and AnnAGNPS

This README outlines the steps to install the `pyagnps` package along with its GDAL dependencies.

# Installation

If you don't intend to use the `topagnps` and `subannagnps` submodules, then GDAL is not necessary and you can ignore the rest of this installation tutorial. Just run:

```Bash
pip install git+https://github.com/LucRSquared/pyagnps.git
```

And you're good to go! Otherwise, see below

# Installation of GDAL

## TL;DR

if you can execute without errors the following line of code within your environment then go to step 4
```bash
python -c "from osgeo import gdal"
```

## 1. Installing GDAL for Linux/MacOS

If you're using Linux/macOS, you need to install the system GDAL libraries before installing the GDAL Python bindings

The installation method for GDAL varies depending on your operating system and package manager. Here are some examples:

**Ubuntu/Debian:**

```bash
    sudo apt update  # Update package lists
    sudo apt install -y libgdal-dev python3-gdal  # Install GDAL libraries
```

 **Fedora/CentOS:**

```bash
    sudo dnf update  # Update package lists
    sudo dnf install -y gdal python3-gdal
```

 **Homebrew (macOS):**

```bash
    brew update
    brew install gdal
```

## 2. Installing GDAL for Windows

The `install_gdal.py` will attempt to find a suitable wheel to install in your python environment. You can skip executing this script if you already have a prefered way to install gdal in your python environment.


## 3. Running `install_gdal.py`

- For Linux/macOS this will install the python bindings for GDAL
- For Windows this will download and install a GDAL wheel from [Christoph Gohlke's website](https://www.cgohlke.com/). If it doesn't work then God be with you, you're on your own, use your favorite search engine to figure out how to install GDAL refer to the TL;DR, that's what you want to achieve.

- Run the following command to download and execute the `install_gdal.py` script:

```bash
   curl -fsSL https://raw.githubusercontent.com/LucRSquared/pyagnps/master/install_gdal.py| python3 -
```
    
This command downloads the script from the `pyagnps` GitHub repository and pipes it to the Python interpreter for execution.


## 4. Once GDAL is installed in your `python` environment:

Now you're ready to install `pyagnps` ! Just run:

```Bash
pip install git+https://github.com/LucRSquared/pyagnps.git
```

## 5. Check installation

If the following command does not return any error then the installation is a success!

```bash
python -c "import pyagnps"
```

## Additional Notes:

Ensure you have activated your virtual environment (if you're using one) before running these commands.
