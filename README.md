# <img alt="SoSpEx" src="sospex/icons/sospexlogo.png" height="100">

[![Anaconda-Server Badge](https://anaconda.org/darioflute/sospex/badges/version.svg)](https://anaconda.org/darioflute/sospex)
[![Anaconda-Server Badge](https://anaconda.org/darioflute/sospex/badges/latest_release_date.svg)](https://anaconda.org/darioflute/sospex)
[![PyPI version](https://badge.fury.io/py/sospex.svg)](https://badge.fury.io/py/sospex)
[![GitHub version](https://badge.fury.io/gh/darioflute%2Fsospex.svg)](https://badge.fury.io/gh/darioflute%2Fsospex)
[![Anaconda-Server Badge](https://anaconda.org/darioflute/sospex/badges/license.svg)](https://anaconda.org/darioflute/sospex)
[![Anaconda-Server Badge](https://anaconda.org/darioflute/sospex/badges/platforms.svg)](https://anaconda.org/darioflute/sospex)

**sospex** displays FIFI-LS, PACS, and GREAT spectral cubes and allows interactions with them.
The cube is shown as a 2D image (spatial image at a given wavelength) and a spectrum (plot of single
spatial pixel). The user can define apertures to compute the spectrum on the 2D image or a wavelength
range to show the spatial average emission. The cube can be manipolated (cropped and trimmed) and
moments of the emission can be computed and displayed. 
External images and cubes can be uploaded and overplotted on the spectral cube.

- **Source:** https://github.com/darioflute/sospex
- **Bug reports:** https://github.com/darioflute/sospex/issues
- **Anaconda:** https://anaconda.org/darioflute/sospex
- **How to install:** https://github.com/darioflute/sospex/blob/master/INSTALL.md
- **Tutorials:** https://nbviewer.jupyter.org/github/darioflute/sospex/blob/master/sospex/help/tutorials.ipynb

Once installed, start it by typing in a terminal window:
```bash
    sospex
```
This will open a window with a few instructions and a series of buttons.
Each button has a tooltip help, i.e. hovering over the button will make appear
a succinct explanation of what the button is supposed to do.

- To open a FITS file containing a spectral cube from FIFI-LS, press the folder icon. 
A pop-up window will allow you to navigate the directory tree and select a file.
- To quit the program, click on the icon shaped as an exiting person.