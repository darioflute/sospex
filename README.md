# <img alt="SoSpEx" src="sospex/icons/sospexlogo.png" height="100">

**sospex** displays FIFI-LS, PACS, and GREAT spectral cubes and allows interactions with them.
In particular, the cube is shown as a 2D image (spatial image obtained as
average along the wavelength dimension) and a spectrum (sum of spatial pixels
of the original cube).
The user can define apertures to compute the spectrum on the 2D image.
It is possible also to crop the cube over the spatial image or cut it along
the wavelength dimension.
External images and cubes can be uploaded and overplotted on the spectral cube.

Current anaconda version: 0.37.beta  (May 6, 2019)

Developed and tested on Mac OS-X 10.11.6 (El Capitan) and 
Ubuntu Linux (14.04 LTS, Trusty Tahr). An untested version is provided for Windows.

- **Source:** https://github.com/darioflute/sospex
- **Bug reports:** https://github.com/darioflute/sospex/issues
- **How to install:** https://github.com/darioflute/sospex/blob/master/INSTALL.md
- **Tutorials:** https://nbviewer.jupyter.org/github/darioflute/sospex/blob/master/sospex/help/tutorials.html

Once installed, start it by typing in a window:

sospex

This will open a window with a few instructions and a series of buttons.
Each button has a tooltip help, i.e. hovering over the button will make appear
a succinct explanation of what the button is supposed to do.

To open a FITS file containing a spectral cube from FIFI-LS, press the folder icon. 
A pop-up window will allow you to navigate the directory tree and select a file.

To quit the program, click on the icon shaped as an exiting person.