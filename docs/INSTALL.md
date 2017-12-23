# HYDROID: Installation instructions

## Quick installation options
### Using Anaconda Python distribution (Mac, Linux, Windows(only HYDROIDexp module))
Install Miniconda with Python2.7 for your platform from [https://conda.io/miniconda.html](https://conda.io/miniconda.html).
Open terminal (Anaconda terminal on Windows) and run:
```
conda install -c hydroid hydroid
```

Test HYDROID:

```
HYDROID_test_exp #Tests exeprimental data analysis module
HYDROID_test_pred #Tests molecular structure analysis module
```

### Using PyPI and pip with native Python
Have Python2.7 and [pip](https://pypi.python.org/pypi/pip) installed.
Open terminal and run:
```
pip install hydroid
```
To enable functionality of HYDROIDpred module (currently for Mac and Linux) compile and install [freesasa library](http://freesasa.github.io) v2.0.2.

Test HYDROID:

```
HYDROID_test_exp #Tests exeprimental data analysis module
HYDROID_test_pred #Tests molecular structure analysis module
```

## Generic instructions for advanced users
### Requirements
- Python 2.7
- Python modules specified in [requirements.txt](../requirements.txt). May be installed with `pip install -r requirements.txt`
- [freesasa library](http://freesasa.github.io) v2.0.2 should be compiled with python bindings and accessible in python (only needed for HYDROIDpred module).
### Manual installation
After installing the required components, download and unarchive the [latest release of HYDROID](https://github.com/ncbi/HYDROID/releases).
Run:
```
python setup.py install
```

## From scratch copy-paste instructions on various platforms (using virtual environments)
Below we provide several copy-paste examples for different operating systems using native Python and separating the install into a virtual environment. 

## On Ubuntu Linux (v16.04 x64) with native Python
Open bash terminal an execute following commands.
~~~~
#prepare package manager and core packages
sudo apt-get install software-properties-common
sudo apt-add-repository universe
sudo apt-get update
sudo apt-get -y install python-pip
sudo apt-get -y install python-tk

#Enable virtualenv
pip install virtualenv
virtualenv venv
source venv/bin/activate
pip install --upgrade pip

#Install HYDROID
pip install https://github.com/ncbi/HYDROID/archive/master.tar.gz

#Install FREESASA (optional, only for HYDROIDpred)
pip install Cython
wget https://github.com/mittinatten/freesasa/releases/download/2.0.2/freesasa-2.0.2.tar.gz
mkdir freesasa
tar -zxf freesasa-2.0.2.tar.gz -C freesasa --strip-components=1
cd freesasa
./configure --enable-python-bindings --disable-json --disable-xml CFLAGS="-fPIC -O2" --prefix=`pwd`
make; make install

#Test that it's working
HYDROID_test_exp #Tests exeprimental data analysis module
HYDROID_test_pred #Tests molecular structure analysis module
~~~~

## On MacOS with native Python:
Open Terminal an execute following commands.
~~~~
#Enable virtualenv
pip install virtualenv
virtualenv venv
source venv/bin/activate
pip install --upgrade pip
pip install wget

#Install HYDROID
pip install https://github.com/ncbi/HYDROID/archive/master.tar.gz

#Install FREESASA (optional, only for HYDROIDpred)
pip install Cython
python -m wget https://github.com/mittinatten/freesasa/releases/download/2.0.2/freesasa-2.0.2.tar.gz
mkdir freesasa
tar -zxf freesasa-2.0.2.tar.gz -C freesasa --strip-components=1
cd freesasa
./configure --enable-python-bindings --disable-json --disable-xml CFLAGS="-fPIC -O2" --prefix=`pwd`
make; make install
cd ..

#MacOS specific stuff
deactivate
export PYTHONHOME=`pwd`/venv

#Test that it's working
python
from hydroid.command_line import test_exp, test_pred
test_exp()
test_pred()
~~~~

## On MacOS with Continuum Anaconda Python

First, install Miniconda with Python2.7 from [https://conda.io/miniconda.html](https://conda.io/miniconda.html).

Then open terminal an execute following commands.
~~~~
conda install wget

#download HYDROID
wget --no-check-certificate https://github.com/ncbi/HYDROID/archive/master.tar.gz
tar -zxf master.tar.gz
cd HYDROID-master

#Create environments and install packages
conda env create -f conda_env.yml
source activate hydroid

#Install HYDROID
python setup.py install

#Install FREESASA (optional, only for HYDROIDpred)
conda install Cython
wget --no-check-certificate https://github.com/mittinatten/freesasa/releases/download/2.0.2/freesasa-2.0.2.tar.gz
mkdir freesasa
tar -zxf freesasa-2.0.2.tar.gz -C freesasa --strip-components=1
cd freesasa
./configure --enable-python-bindings --disable-json --disable-xml CFLAGS="-fPIC -O2" --prefix=`pwd`
make; make install
cd ..

#Test that it's working
HYDROID_test_exp #Tests exeprimental data analysis module
HYDROID_test_pred #Tests molecular structure analysis module
~~~~
