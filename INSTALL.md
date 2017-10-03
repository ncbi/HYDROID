# Installation examples for HYDROID
## General comments
There are many different ways of installing Python together with its modules.
Installing Python and its modules is usually system specific. Using `virtualenv` together with `pip` or other python package managers (Continuum Anaconda, Enthought Canopy, WinPython) to match the exact versions of the modules is recommended. We provide environment files for `pip` ([requirements.txt](requirements.txt)) and Anaconda ([conda_env.yml](conda_env.yml)).

The most system specific requirement of HYDROID is the [matplotlib](http://matplotlib.org/users/installing.html) library and its graphical backends.

Below we provide several examples for different operating systems that can get you started. 

## On Ubuntu Linux (v16.04) with native Python
Open terminal an execute following commands.
~~~~
#prepare package manager and core packages
sudo apt-get install software-properties-common
sudo apt-add-repository universe
sudo apt-get update
sudo apt-get -y install python-pip
sudo apt-get -y install python-tk

#download HYDROID
wget https://github.com/ncbi/HYDROID/archive/v0.0.2.tar.gz
tar -zxf v0.0.2.tar.gz
mv HYDROID-0.0.2 HYDROID
cd HYDROID

#Enable virtualenv
pip install virtualenv
virtualenv venv
source venv/bin/activate

#Install packages
pip install --upgrade pip
pip install -r requirements.txt

#Install FREESASA (optional, only for HYDROIDpred)
pip install Cython
wget https://github.com/mittinatten/freesasa/releases/download/2.0.1/freesasa-2.0.1.tar.gz
mkdir freesasa
tar -zxf freesasa-2.0.1.tar.gz -C freesasa --strip-components=1
cd freesasa
./configure --enable-python-bindings --disable-json --disable-xml CFLAGS="-fPIC -O2" --prefix=`pwd`
make; make install
cd ..

#Test that it's working
cd example1
python exp_s2_assign_peaks.py
~~~~

## On Ubuntu Linux (v16.04) with Continuum Anaconda Python

First, install Anaconda with Python2.7 from [https://www.continuum.io](https://www.continuum.io)
Do not forget to add Anaconda to your PATH. Then open terminal an execute following commands.
~~~~
#download HYDROID
wget https://github.com/ncbi/HYDROID/archive/v0.0.2.tar.gz
tar -zxf v0.0.2.tar.gz
mv HYDROID-0.0.2 HYDROID
cd HYDROID

#Create environments and install packages
conda env create -f conda_env.yml
source activate hydroid

#Install FREESASA (optional, only for HYDROIDpred)
pip install Cython
wget https://github.com/mittinatten/freesasa/releases/download/2.0.1/freesasa-2.0.1.tar.gz
mkdir freesasa
tar -zxf freesasa-2.0.1.tar.gz -C freesasa --strip-components=1
cd freesasa
./configure --enable-python-bindings --disable-json --disable-xml CFLAGS="-fPIC -O2" --prefix=`pwd`
make; make install
cd ..

#Test that it's working
cd example1
python exp_s2_assign_peaks.py
~~~~


## On MacOS with native Python:
Open terminal an execute following commands.
~~~~
#download HYDROID
wget --no-check-certificate https://github.com/ncbi/HYDROID/archive/v0.0.2.tar.gz
tar -zxf v0.0.2.tar.gz
mv HYDROID-0.0.2 HYDROID
cd HYDROID

#Enable virtualenv
pip install virtualenv
virtualenv venv
source venv/bin/activate

#Install packages
pip install --upgrade pip
pip install -r requirements.txt

#Install FREESASA (optional, only for HYDROIDpred)
pip install Cython
wget --no-check-certificate https://github.com/mittinatten/freesasa/releases/download/2.0.1/freesasa-2.0.1.tar.gz
mkdir freesasa
tar -zxf freesasa-2.0.1.tar.gz -C freesasa --strip-components=1
cd freesasa
./configure --enable-python-bindings --disable-json --disable-xml CFLAGS="-fPIC -O2" --prefix=`pwd`
make; make install
cd ..

#MacOS specific stuff
deactivate
export PYTHONHOME=`pwd`/venv

#Test that it's working
cd example1
python exp_s2_assign_peaks.py
~~~~

## On MacOS with Continuum Anaconda Python

First, install Anaconda with Python2.7 from [https://www.continuum.io](https://www.continuum.io).
Then open terminal an execute following commands.
~~~~
#download HYDROID
wget --no-check-certificate https://github.com/ncbi/HYDROID/archive/v0.0.2.tar.gz
tar -zxf v0.0.2.tar.gz
mv HYDROID-0.0.2 HYDROID
cd HYDROID

#Create environments and install packages
conda env create -f conda_env.yml
source activate hydroid

#Install FREESASA (optional, only for HYDROIDpred)
pip install Cython
wget --no-check-certificate https://github.com/mittinatten/freesasa/releases/download/2.0.1/freesasa-2.0.1.tar.gz
mkdir freesasa
tar -zxf freesasa-2.0.1.tar.gz -C freesasa --strip-components=1
cd freesasa
./configure --enable-python-bindings --disable-json --disable-xml CFLAGS="-fPIC -O2" --prefix=`pwd`
make; make install
cd ..

#Test that it's working
cd example1
python exp_s2_assign_peaks.py
~~~~


## On Windows with Continuum Anaconda Python (currently only for HYDROIDexp part):
First, install Anaconda with Python2.7 from [https://www.continuum.io](https://www.continuum.io)
Do not forget to add Anaconda to your PATH.
Download HYDROID manually from https://github.com/ncbi/HYDROID/archive/v0.0.2.zip

Unzip the contents and change the directory name to `HYDROID`.
Open command line, chande directory to `HYDROID` folder and execute following commands.
~~~~
conda env create -f conda_env.yml
activate hydroid

#Test that it's working
cd example1
python exp_s2_assign_peaks.py
~~~~