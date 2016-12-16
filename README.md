# sospex
 (SOFIA SPectrum EXplorer)

This code displays FIFI-LS spectral cubes and allows interactions with them.
In particular, the cube is shown as a 2D image (spatial image obtained as
average along the wavelength dimension) and a spectrum (sum of spatial pixels
of the original cube).
The user can define the ellipse to compute the spectrum on the 2D image.
It is possible also to crop the cube over the spatial image or cut it along
the wavelength dimension.
Lines and continuum can be fitted on the spectrum.
External spectra and images can be uploaded and overplotted on the FIFI-LS
cube.

The code runs with Python 2.x since its uses wxpython (not
yet officially ported to Python 3).
The code has been developed on Ubuntu Linux. It runs also on Mac OS-X.

INSTALLATION NOTES

To install the code in your computer, you need first to install the anaconda
python (https://www.continuum.io/downloads).
You will have to use the Python 2.x distribution since the code
uses wxpython which is still not ported to Python 3.

You have also to install the conda package lmfit:
conda install -c astropy lmfit

Then clone the repository:
git clone https://github.com/darioflute/sospex.git

Build the conda package:
conda build sospex

And install it locally:

conda install --use-local sospex

At this point you can start the code everywhere by
typing:

sospex

since the executable is in the ~/anaconda/bin directory.
