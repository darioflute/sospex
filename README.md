# sospex

## SOFIA SPectrum EXplorer

This code displays FIFI-LS and GREAT spectral cubes and allows interactions with them.
In particular, the cube is shown as a 2D image (spatial image obtained as
average along the wavelength dimension) and a spectrum (sum of spatial pixels
of the original cube).
The user can define apertures to compute the spectrum on the 2D image.
It is possible also to crop the cube over the spatial image or cut it along
the wavelength dimension.
External images can be uploaded and overplotted on the spectral cube.

The code runs with Python 3.
The code has been developed on Mac OS-X 10.11.6 (El Capitan) and Ubuntu Linux (14.04 LTS, Trusty Tahr).

### Current version

The current version released on anaconda is: 0.26.beta 

Release date:  April 30, 2018

### How to install the program

To install the code in your computer, you need first to install the anaconda
python (https://www.continuum.io/downloads).
You will have to use the Python 3.x distribution.

The packages *lmfit* and *reproject* have to be installed from the astropy channel:

conda install -c astropy lmfit

conda install -c astropy reproject

Now you are ready to install sospex:

conda install -c darioflute sospex

If there is a new release, you will have only to update the package:

conda update -c darioflute sospex


### How to start the program

To start simply type:

sospex

This will open a window with a few instructions and a series of buttons.
Each button has a tooltip help, i.e. hovering over the button will make appear
a succinct explanation of what the button is supposed to do.

To open a FITS file containing a spectral cube from FIFI-LS, you have to
press the folder icon. This will pop-up a open window to navigate the directory
tree and select a file.

To quit the program, click on the icon shaped as an exiting person.

The question mark icon open a series of tutorials on your web browser.
