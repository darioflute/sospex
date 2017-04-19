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

### INSTALLATION NOTES

To install the code in your computer, you need first to install the anaconda
python (https://www.continuum.io/downloads).
You will have to use the Python 2.x distribution since the code
uses wxpython which is still not ported to Python 3.

You have also to install the conda package lmfit:
conda install -c astropy lmfit

At this point you can install sospex:

conda install -c darioflute sospex

### Start the program

To start simply type:

sospex

This will open a window with a few instructions and a series of buttons.
Each button has a tooltip help, i.e. hovering over the button will give
a succinct explanation of what the button is supposed to do.


