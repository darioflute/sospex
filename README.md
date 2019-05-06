# <img alt="SoSpEx" src="sospex/icons/sospexlogo.png" height="100">

SoSpEx displays FIFI-LS, PACS, and GREAT spectral cubes and allows interactions with them.
In particular, the cube is shown as a 2D image (spatial image obtained as
average along the wavelength dimension) and a spectrum (sum of spatial pixels
of the original cube).
The user can define apertures to compute the spectrum on the 2D image.
It is possible also to crop the cube over the spatial image or cut it along
the wavelength dimension.
External images and cubes can be uploaded and overplotted on the spectral cube.

- **Source:** https://github.com/darioflute/sospex
- **Bug reports:** https://github.com/darioflute/sospex/issues
- **Install instructions:** https://github.com/darioflute/sospex/blob/master/INSTALL.md

### Current version

The current version released on anaconda is: 0.37.beta  

Release date:  May 6, 2019

It has been developed and tested on Mac OS-X 10.11.6 (El Capitan) and 
Ubuntu Linux (14.04 LTS, Trusty Tahr).
A version is provided for Windows with limited testing.


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

### Tutorials

[Load a spectral cube](https://nbviewer.jupyter.org/github/darioflute/sospex/blob/master/sospex/help/start.ipynb)

[Zoom and Pan](https://nbviewer.jupyter.org/github/darioflute/sospex/blob/master/sospex/help/zoom.ipynb)

[Change image intensity levels](https://nbviewer.jupyter.org/github/darioflute/sospex/blob/master/sospex/help/intensity.ipynb)

[Using the spectral panel](https://nbviewer.jupyter.org/github/darioflute/sospex/blob/master/sospex/help/specpanel.ipynb)

[Sliding along the spectral dimension](https://nbviewer.jupyter.org/github/darioflute/sospex/blob/master/sospex/help/slider.ipynb)

[Select a slice of the cube](https://nbviewer.jupyter.org/github/darioflute/sospex/blob/master/sospex/help/slice.ipynb)

[Crop and cut the cubes](https://nbviewer.jupyter.org/github/darioflute/sospex/blob/master/sospex/help/cutcrop.ipynb)

[Select and modify apertures](https://nbviewer.jupyter.org/github/darioflute/sospex/blob/master/sospex/help/apertures.ipynb)

[Load external images](https://nbviewer.jupyter.org/github/darioflute/sospex/blob/master/sospex/help/extimages.ipynb)

[Blink between images](https://nbviewer.jupyter.org/github/darioflute/sospex/blob/master/sospex/help/blink.ipynb)

[Overlap contours](https://nbviewer.jupyter.org/github/darioflute/sospex/blob/master/sospex/help/contours.ipynb)

[Save images and spectra](https://nbviewer.jupyter.org/github/darioflute/sospex/blob/master/sospex/help/save.ipynb)

[Fit continuum and compute moments](https://nbviewer.jupyter.org/github/darioflute/sospex/blob/master/sospex/help/moments.ipynb)

[Erase part of the cube](https://nbviewer.jupyter.org/github/darioflute/sospex/blob/master/sospex/help/erase.ipynb)
