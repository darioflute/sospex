### Install with anaconda

To install the code in your computer, you need first to install the anaconda
python (https://www.anaconda.com/download).

**You will have to use the Python 3.7 distribution.**

If you have a previous Python distribution installed, it is advisable to start from scratch
with a new Anaconda download of Python 3.7.

It is a good practice to update all your installed packages to have the latest versions:

conda update --all -y

The packages *lmfit*, *reproject*, and *fitsio* have to be installed from other channels:

conda install -c conda-forge lmfit -y

conda install -c conda-forge reproject -y

conda install -c conda-forge fitsio -y

Now you are ready to install sospex:

conda install -c darioflute sospex -y

### Using an environment

It is possible to install it inside an environment.
To create it use the command:

conda create -n sospex python=3.7 anaconda

Then activate it:

conda activate sospex

and install sospex by following the installation commands in the previous section.
To exit the environment, type:

conda deactivate

### Note for windows users

There is currently no distribution available for fitsio for windows.
It is still possible to install sospex without the fitsio library, since the default to open fits is the astropy library.

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
