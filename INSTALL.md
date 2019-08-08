### Install with anaconda

To install the code in your computer, you need first to install the anaconda
python (https://www.anaconda.com/download).

**You will have to use the Python 3.7 distribution.**

If you have a previous Python distribution installed, it is advisable to start from scracth
with a new Anaconda download of Python 3.7.

It is a good practice to update all your installed packages to have the latest versions:

conda update --all -y

The packages *lmfit*, *reproject*, and *fitsio* have to be installed from other channels:

conda install -c cprescher lmfit -y

conda install -c conda-forge reproject -y

conda install -c conda-forge fitsio -y

Now you are ready to install sospex:

conda install -c darioflute sospex -y

### Updating with anaconda

If there is a new release, you will have only to update the package:

conda update -c darioflute sospex -y

Alternatively, add the channel to your anaconda:

conda config --add channels darioflute

and simply do:

conda update sospex

### Install with pipy  (Currently the version on pipy is pretty old)

pip install sospex

to update:

pip install --upgrade sospex
