# <img alt="SoSpEx" src="sospex/icons/sospexlogo.png" height="100">

[![Anaconda-Server Badge](https://anaconda.org/darioflute/sospex/badges/version.svg?branch=master&kill_cache=1&service=github)](https://anaconda.org/darioflute/sospex)
[![Anaconda-Server Badge](https://anaconda.org/darioflute/sospex/badges/latest_release_date.svg?branch=master&kill_cache=1&service=github)](https://anaconda.org/darioflute/sospex)
[![PyPI version](https://badge.fury.io/py/sospex.svg?branch=master&kill_cache=1)](https://badge.fury.io/py/sospex)
[![Anaconda-Server Badge](https://anaconda.org/darioflute/sospex/badges/license.svg)](https://anaconda.org/darioflute/sospex)
[![Anaconda-Server Badge](https://anaconda.org/darioflute/sospex/badges/platforms.svg)](https://anaconda.org/darioflute/sospex)

**sospex** is a GUI tool to display and analyse [SOFIA](https://www.sofia.usra.edu) ([FIFI-LS](https://www.sofia.usra.edu/science/instruments/fifi-ls) and [GREAT](https://www.sofia.usra.edu/science/instruments/great)) and [Herschel](http://sci.esa.int/herschel/) ([PACS](https://www.cosmos.esa.int/web/herschel/pacs-overview)) spectral cubes.

- **Source:** https://github.com/darioflute/sospex
- **Bug reports:** https://github.com/darioflute/sospex/issues
- **Anaconda:** https://anaconda.org/darioflute/sospex
- **How to install:** https://github.com/darioflute/sospex/blob/master/INSTALL.md
- **Tutorials:** https://nbviewer.jupyter.org/github/darioflute/sospex/blob/master/sospex/help/tutorials.ipynb

Once installed, type in a terminal window:
```bash
    sospex
```
This opens a GUI with a few instructions and a series of buttons.
Buttons have tooltip helps, i.e. hovering over a button makes appear
a succinct explanation of what the button is supposed to do.

Cubes are shown as images (spatial images at given wavelengths) and spectra (plots in a sky direction). Cubes can be explored by moving the cursor on the image and a bar on the spectrum. Spectra can be extracted by defining apertures on the image.
By selecting a range on the spectrum, the average spatial emission in the wavelength range is displayed. The cube can be manipolated (cropped and trimmed) and moments of the emission can be computed and displayed. 
External images and cubes can be uploaded and their contours can be overplotted on the image tabs.
