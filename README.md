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

### How to install the program

To install the code in your computer, you need first to install the anaconda
python (https://www.continuum.io/downloads).
You will have to use the Python 2.x distribution since the code
uses wxpython which is still not ported to Python 3.

You have also to install the conda package lmfit:
conda install -c astropy lmfit

At this point you can install sospex:

conda install -c darioflute sospex

### How to start the program

To start simply type:

sospex

This will open a window with a few instructions and a series of buttons.
Each button has a tooltip help, i.e. hovering over the button will make appear
a succinct explanation of what the button is supposed to do.

To open a FITS file containing a spectral cube from FIFI-LS, you have to
press the double arrow icon. This will pop-up a open window to navigate the directory
tree and select a file.

To quit the program, click on the icon shaped as a red cross.

The question mark icon open this README file.

The recycling icon allows one to reload the same file in the case the original file
has been cut or cropped.

### How to interact with the cube

The main window is divided in two panels.

**Left panel**
On the left, an image is displayed. Clicking on the right button of the mouse it is
possible to select which image to display (flux, uncorrected flux, or exposure).
The name of the selected cube is then printed as title.

An ellipse will also appear centered on the point with the highest flux in the image.
It is possible to move and modify the ellipse with the mouse.

Two sliders on the bottom of the image allows one to change the limits of dynamical
intensity range displayed.

Using the mouse wheel, it is possible to zoom on the image. The central button allows one
to center the image on the position of the mouse.
Finally, to go back to the original limits, click on the yellow arrow icon.

Clicking on the contour icon, it is possible to display contours of one of the cubes (flux,
unc. flux, or exposure) over the displayed image.

The camera icon allows one to take a snapshot of the image and save it as a file in several formats.

The header icon simply displays the header of the selected cube.

**Right panel**
Also here, clicking on the right button of the mouse, you can choose to display several
curves (flux, uncorr. flux, exposure, etc.)
By default, the left panel shows the average of the entire spectrum for each spatial module
or spaxel.
You can choose to average only part of the spectrum by clicking and dragging the mouse
over the wavelength range of interest.

To pan or zoom over the plot, use the four arrows and the lens icons.
The yellow arrow send back to the original limits.
Do not forget to click again on the pan and zoom icons to exit from the panning/zooming mode.

Also here, a snapshot of the plot can be taken using the camera icon.

### Modifying the cube

Sometimes, one can be interested in exploring or saving only part of a cube.
For this reason, we included two functions:

**cropping**  on the left panel, the icon with the cutter allows one to crop the spatial part
of the cube. After clicking the icon, on the image click and drag to select a rectangle.
Then a dialog will appear to confirm the cropping selection.

**cutting**  on the right panel, one can select a wavelength region by clicking a dragging the
mouse over the spectrum. Then, clicking on the scissor icon, a dialog will appear to confirm the
selection of wavelengths made.



