import numpy as np
import matplotlib
import os
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT 
from matplotlib.ticker import ScalarFormatter
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Polygon,FancyArrowPatch
#import matplotlib.ticker as ticker

# Matplotlib parameters
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family']='STIXGeneral'
rcParams['font.size']=13
rcParams['mathtext.fontset']='stix'
rcParams['legend.numpoints']=1

from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.widgets import SpanSelector
import matplotlib.transforms as transforms

from astropy.wcs.utils import proj_plane_pixel_scales as pixscales
from astropy.coordinates import SkyCoord
from astropy import units as u

from PyQt5.QtWidgets import (QVBoxLayout, QHBoxLayout, QLineEdit, QSizePolicy, QInputDialog, 
                             QPushButton,QLabel,QMessageBox,QScrollArea,QWidget)
from PyQt5.QtGui import QFont
from PyQt5.QtCore import Qt, QSize, pyqtSignal
from PyQt5.QtTest import QTest

# https://het.as.utexas.edu/HET/Software/Astropy-1.0/_modules/astropy/visualization/stretch.html
from astropy.visualization import (LinearStretch, SqrtStretch, SquaredStretch, SinhStretch,AsinhStretch,
                                   LogStretch, ImageNormalize, PowerStretch)

class ScrollMessageBox(QMessageBox):
    def __init__(self, l, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.list = l
        font = QFont()
        font.setFamily("Monaco")
        font.setPointSize(10)
        self.setWindowTitle("Header")
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        self.content = QWidget()
        scroll.setWidget(self.content)
        vLayout = QVBoxLayout(self.content)
        hLayout = QHBoxLayout()
        self.filter = QPushButton("filter", self)
        self.filter.clicked.connect(self.filterClicked)
        hLayout.addWidget(self.filter)
        self.lineEdit = QLineEdit(self)
        hLayout.addWidget(self.lineEdit)
        vLayout.addLayout(hLayout)
        for item in l:
            label = QLabel(item, self)
            label.setFont(font)
            vLayout.addWidget(label)
        self.layout().addWidget(scroll, 0, 0, 1, self.layout().columnCount())
        self.setStyleSheet("QScrollArea{min-width:600 px; min-height: 400px}")
    
    def filterClicked(self):
        import re
        filter_text = str(self.lineEdit.text()).strip()
        for l in self.list:
            if bool(re.search(filter_text, l)):
                print('String ', l)


def ds9cmap():
    """Adding DS9 colormap.
    Adapted from http://nbviewer.jupyter.org/gist/adonath/c9a97d2f2d964ae7b9eb"""
    from matplotlib.cm import register_cmap, cmap_d
    ds9b = {'red': lambda v : 4 * v - 1, 
            'green': lambda v : 4 * v - 2,
            'blue': lambda v : np.select([v < 0.25, v < 0.5, v < 0.75, v <= 1],
                                         [4 * v, -4 * v + 2, 0, 4 * v - 3])}
    ds9b_r = {'red': lambda v : (4 * v - 1)[::-1], 
            'green': lambda v : (4 * v - 2)[::-1],
            'blue': lambda v : np.select([v < 0.25, v < 0.5, v < 0.75, v <= 1],
                                         [4 * v, -4 * v + 2, 0, 4 * v - 3][::-1])}    
    # Note that this definition slightly differs from ds9cool, but make more sense to me...
    ds9cool = {'red': lambda v : 2 * v - 1, 
               'green': lambda v : 2 * v - 0.5,
               'blue': lambda v : 2 * v}
    ds9cool_r = {'red': lambda v : (2 * v - 1)[::-1], 
               'green': lambda v : (2 * v - 0.5)[::-1],
                 'blue': lambda v : (2 * v)[::-1]}    
    ds9a = {'red': lambda v : np.interp(v, [0, 0.25, 0.5, 1],
                                        [0, 0, 1, 1]),
            'green': lambda v : np.interp(v, [0, 0.25, 0.5, 0.77, 1],
                                          [0, 1, 0, 0, 1]),
            'blue': lambda v : np.interp(v, [0, 0.125, 0.5, 0.64, 0.77, 1],
                                         [0, 0, 1, 0.5, 0, 0])}
    ds9a_r = {'red': lambda v : np.interp(v, [0, 0.25, 0.5, 1],
                                        [0, 0, 1, 1])[::-1],
            'green': lambda v : np.interp(v, [0, 0.25, 0.5, 0.77, 1],
                                          [0, 1, 0, 0, 1])[::-1],
            'blue': lambda v : np.interp(v, [0, 0.125, 0.5, 0.64, 0.77, 1],
                                         [0, 0, 1, 0.5, 0, 0])[::-1]}    
    ds9i8 = {'red': lambda v : np.where(v < 0.5, 0, 1), 
             'green': lambda v : np.select([v < 1/8., v < 0.25, v < 3/8., v < 0.5,
                                            v < 5/8., v < 0.75, v < 7/8., v <= 1],
                                           [0, 1, 0, 1, 0, 1, 0, 1]),
             'blue': lambda v : np.select([v < 1/8., v < 0.25, v < 3/8., v < 0.5,
                                           v < 5/8., v < 0.75, v < 7/8., v <= 1],
                                          [0, 0, 1, 1, 0, 0, 1, 1])}    
    ds9i8_r = {'red': lambda v : np.where(v > 0.5, 0, 1), 
             'green': lambda v : np.select([v < 1/8., v < 0.25, v < 3/8., v < 0.5,
                                            v < 5/8., v < 0.75, v < 7/8., v <= 1],
                                           [0, 1, 0, 1, 0, 1, 0, 1][::-1]),
             'blue': lambda v : np.select([v < 1/8., v < 0.25, v < 3/8., v < 0.5,
                                           v < 5/8., v < 0.75, v < 7/8., v <= 1],
                                          [0, 0, 1, 1, 0, 0, 1, 1][::-1])}    
    ds9aips0 = {'red': lambda v : np.select([v < 1/9., v < 2/9., v < 3/9., v < 4/9., v < 5/9.,
                                             v < 6/9., v < 7/9., v < 8/9., v <= 1],
                                            [0.196, 0.475, 0, 0.373, 0, 0, 1, 1, 1]), 
                'green': lambda v : np.select([v < 1/9., v < 2/9., v < 3/9., v < 4/9., v < 5/9.,
                                               v < 6/9., v < 7/9., v < 8/9., v <= 1],
                                              [0.196, 0, 0, 0.655, 0.596, 0.965, 1, 0.694, 0]),
                'blue': lambda v : np.select([v < 1/9., v < 2/9., v < 3/9., v < 4/9., v < 5/9.,
                                              v < 6/9., v < 7/9., v < 8/9., v <= 1],
                                             [0.196, 0.608, 0.785, 0.925, 0, 0, 0, 0, 0])}    
    ds9aips0_r = {'red': lambda v : np.select([v < 1/9., v < 2/9., v < 3/9., v < 4/9., v < 5/9.,
                                             v < 6/9., v < 7/9., v < 8/9., v <= 1],
                                            [0.196, 0.475, 0, 0.373, 0, 0, 1, 1, 1][::-1]), 
                'green': lambda v : np.select([v < 1/9., v < 2/9., v < 3/9., v < 4/9., v < 5/9.,
                                               v < 6/9., v < 7/9., v < 8/9., v <= 1],
                                              [0.196, 0, 0, 0.655, 0.596, 0.965, 1, 0.694, 0][::-1]),
                'blue': lambda v : np.select([v < 1/9., v < 2/9., v < 3/9., v < 4/9., v < 5/9.,
                                              v < 6/9., v < 7/9., v < 8/9., v <= 1],
                                             [0.196, 0.608, 0.785, 0.925, 0, 0, 0, 0, 0][::-1])}
    ds9rainbow = {'red': lambda v : np.interp(v, [0, 0.2, 0.6, 0.8, 1], [1, 0, 0, 1, 1]),
                  'green': lambda v : np.interp(v, [0, 0.2, 0.4, 0.8, 1], [0, 0, 1, 1, 0]),
                  'blue': lambda v : np.interp(v, [0, 0.4, 0.6, 1], [1, 1, 0, 0])}    
    ds9rainbow_r = {'red': lambda v : np.interp(v, [0, 0.2, 0.6, 0.8, 1], [1, 0, 0, 1, 1])[::-1],
                  'green': lambda v : np.interp(v, [0, 0.2, 0.4, 0.8, 1], [0, 0, 1, 1, 0])[::-1],
                  'blue': lambda v : np.interp(v, [0, 0.4, 0.6, 1], [1, 1, 0, 0])[::-1]}
    # This definition seems a bit strange...
    ds9he = {'red': lambda v : np.interp(v, [0, 0.015, 0.25, 0.5, 1],
                                         [0, 0.5, 0.5, 0.75, 1]),
             'green': lambda v : np.interp(v, [0, 0.065, 0.125, 0.25, 0.5, 1],
                                           [0, 0, 0.5, 0.75, 0.81, 1]),
             'blue': lambda v : np.interp(v, [0, 0.015, 0.03, 0.065, 0.25, 1],
                                          [0, 0.125, 0.375, 0.625, 0.25, 1])}    
    ds9heat = {'red': lambda v : np.interp(v, [0, 0.34, 1], [0, 1, 1]),
               'green': lambda v : np.interp(v, [0, 1], [0, 1]),
               'blue': lambda v : np.interp(v, [0, 0.65, 0.98, 1], [0, 0, 1, 1])}    
    ds9he_r = {'red': lambda v : np.interp(v, [0, 0.015, 0.25, 0.5, 1],
                                         [0, 0.5, 0.5, 0.75, 1])[::-1],
             'green': lambda v : np.interp(v, [0, 0.065, 0.125, 0.25, 0.5, 1],
                                           [0, 0, 0.5, 0.75, 0.81, 1])[::-1],
             'blue': lambda v : np.interp(v, [0, 0.015, 0.03, 0.065, 0.25, 1],
                                          [0, 0.125, 0.375, 0.625, 0.25, 1])[::-1]}    
    ds9heat_r = {'red': lambda v : np.interp(v, [0, 0.34, 1], [0, 1, 1])[::-1],
               'green': lambda v : np.interp(v, [0, 1], [0, 1])[::-1],
               'blue': lambda v : np.interp(v, [0, 0.65, 0.98, 1], [0, 0, 1, 1])[::-1]}
    # Additional GAIA map
    gaiareal = {'red': lambda v : np.interp(v, [0,  0.5, 1],[0, 1, 1]),
            'green': lambda v : np.interp(v, [0, 0.5, 1],[0, 0.5, 1]),
            'blue': lambda v : np.interp(v, [0, 0.5,  1],[0, 0, 1])}
    gaiareal_r = {'red': lambda v : np.interp(v, [0,  0.5, 1],[0, 1, 1])[::-1],
            'green': lambda v : np.interp(v, [0, 0.5, 1],[0, 0.5, 1])[::-1],
            'blue': lambda v : np.interp(v, [0, 0.5,  1],[0, 0, 1])[::-1]}
    redpt = [0,20,40,155,156,165,183,184,190, 191, 213,216,255]
    red   = [0,0.59,0,0,0.01,.305,0.976,1, 1,.995,.635,.608,1]
    greenpt = [0,41,  43,69.7,131,156,160,172,192, 193,216,221,255]
    green   = [0, 0,.048, 0.4, 1,.77,.812,  1,  1,.975,0.1,  0,  0]
    bluept  = [   0, 42, 69.7,92, 234, 235,240,245,255]
    blue    = [0.17,  1,  0.4, 0,   0,.008,.15,.355,  1]
    gaiarainbow = {'red': lambda v : np.interp(v, redpt, red),
                   'green': lambda v : np.interp(v, greenpt, green),
                   'blue': lambda v : np.interp(v, bluept, blue)}
    gaiarainbow_r = {'red': lambda v : np.interp(v, redpt, red)[::-1],
                   'green': lambda v : np.interp(v, greenpt, green)[::-1],
                   'blue': lambda v : np.interp(v, bluept, blue)[::-1]}
    
    
    # Set aliases, where colormap exists in matplotlib
    cmap_d['ds9bb'] = cmap_d['afmhot']
    cmap_d['ds9grey'] = cmap_d['gray']
    cmap_d['ds9bb_r'] = cmap_d['afmhot_r']
    cmap_d['ds9grey_r'] = cmap_d['gray_r']    
    # Register all other colormaps
    register_cmap('ds9b', data=ds9b)
    register_cmap('ds9cool', data=ds9cool)
    register_cmap('ds9a', data=ds9a)
    register_cmap('ds9i8', data=ds9i8)
    register_cmap('ds9aips0', data=ds9aips0)
    register_cmap('ds9rainbow', data=ds9rainbow)
    register_cmap('ds9he', data=ds9he)
    register_cmap('ds9heat', data=ds9heat)
    # Register reverse colormaps
    register_cmap('ds9b_r', data=ds9b_r)
    register_cmap('ds9cool_r', data=ds9cool_r)
    register_cmap('ds9a_r', data=ds9a_r)
    register_cmap('ds9i8_r', data=ds9i8_r)
    register_cmap('ds9aips0_r', data=ds9aips0_r)
    register_cmap('ds9rainbow_r', data=ds9rainbow_r)
    register_cmap('ds9he_r', data=ds9he_r)
    register_cmap('ds9heat_r', data=ds9heat_r)
    # Additional gaia map
    register_cmap('real', data=gaiareal)
    register_cmap('real_r', data=gaiareal_r)
    register_cmap('rainbow', data=gaiarainbow)
    register_cmap('rainbow_r', data=gaiarainbow_r)
    # The Normalize class is largely based on code provided by Sarah Graves.
    # http://www.ster.kuleuven.be/~pieterd/python/html/plotting/interactive_colorbar.html


class NavigationToolbar(NavigationToolbar2QT):
    def __init__(self,canvas,parent):
        # Select only a few buttons
        self.iconDir = os.path.join(os.path.dirname(os.path.abspath(__file__)),"icons")
        self.toolitems = [
            #('Home','Go back to original limits','home','home'),
            ('Backward','Go back','back','back'),
            ('Forward','Go forward','forward','forward'),
            ('Pan','Pan figure','move','pan'),
            ('Zoom','Zoom in','zoom_to_rect','zoom'),
        ]
        self.parent = parent
        super().__init__(canvas,parent)

        
        
class MplCanvas(FigureCanvas):
    """ Basic matplotlib canvas class """

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self,QSizePolicy.MinimumExpanding,QSizePolicy.MinimumExpanding)
        FigureCanvas.updateGeometry(self)
        self.compute_initial_figure()

    def sizeHint(self):
        w, h = self.get_width_height()
        return QSize(w, h)

    def minimumSizeHint(self):
        return QSize(5,5)
    
    def compute_initial_figure(self):
        pass

    
class ImageCanvas(MplCanvas):
    """Canvas to plot an image."""
    
    def __init__(self, *args, **kwargs):
        MplCanvas.__init__(self, *args, **kwargs)
        # Define color map
            
    def compute_initial_figure(self, image=None, wcs=None, title=None, cMap = 'real',
                               cMapDir = '_r', stretch='linear', instrument=None, aspect=1):
        self.colorMap = cMap
        self.colorMapDirection = cMapDir
        self.stretch = stretch
        self.aspect = aspect
        if wcs is None:
            '''initial definition when images are not yet read'''
            pass
        else:
            # import time
            # t = time.process_time()
            #print('WCS ', wcs)
            self.wcs = wcs
            try:
                self.fig.delaxes(self.axes)
                self.axes = None
                print("Deleting axes")
            except:
                pass
            self.axes = self.fig.add_axes([0.1,0.1,.8,.8], projection = self.wcs)
            self.axes.coords[0].set_major_formatter('hh:mm:ss')
            self.axes.grid(color='black', ls='dashed')
            self.axes.set_xlabel('R.A.')
            self.axes.set_ylabel('Dec')
            # Colorbar
            self.cbaxes = self.fig.add_axes([0.9,0.1,0.02,0.85])
            # Add title
            if title != None:
                self.title = title
                if title == 'uFlux':
                    title = 'Flux [no atm. corr.]'
                elif title == 'Exp':
                    title = 'Coverage map'
                elif title == 'C0':
                    title = 'Continuum at $\lambda_0$'
                elif title == 'M0':
                    title = 'M$_0$'
                elif title == 'M1':
                    title = 'M$_1$'
                elif title == 'M2':
                    title = 'M$_2$'
                elif title == 'M3':
                    title = 'Normalized Skewness'
                elif title == 'M4':
                    title = 'Kurtosis Excess'
                elif title == 'sv':
                    title = '$FWHM$'
                elif title == 'v':
                   title = 'v'
                elif title == 'L0':
                   title = 'I$_0$'
                elif title == 'L1':
                   title = 'I$_1$'
                elif title == 'v0':
                   title = 'v$_0$'
                elif title == 'v1':
                   title = 'v$_1$'
                elif title == 'd0':
                   title = 'FWHM$_0$'
                elif title == 'd1':
                   title = 'FWHM$_1$'
                elif title == 'unknown':
                   title = ''
                self.fig.suptitle(title)
            # Show image
            self.contrast = 1.
            self.bias = 0.5
            if instrument is not None:
                if instrument == 'GREAT':
                    self.cmin = 0  # np.nanmin(image)
                    self.cmax = 1  # np.nanmax(image)
                else:
                    self.cmin = None
                    self.cmax = None
            else:
                self.cmin = None
                self.cmax = None
            # import time
            # t1 = time.process_time()
            # print('Pre-image is ', t1-t, 's')
            if image is not None:
                self.showImage(image)
            # t2 = time.process_time() 
            # print('Image displayed in ',t2-t, 's')
            self.changed = False            
            # Add ellipse centered on source
            self.pixscale = pixscales(self.wcs)[0]*3600. # Scale in arcsec
            # Apertures
            self.photApertures = []
            self.photApertureSignal = []
            # Contours
            self.contour = None            
            # Activate focus
            self.setFocusPolicy(Qt.ClickFocus)
            self.setFocus()
            # t3 = time.process_time() 
            # print('post image is ', t3-t2,'s')

    def stretchFunc(self, newStretch):
        if newStretch == 'linear':
            stretch = LinearStretch()
        elif newStretch == 'sqrt':
            stretch = SqrtStretch()
        elif newStretch == 'square':
            stretch = SquaredStretch()
        elif newStretch == 'log':
            stretch = LogStretch()
        elif newStretch == 'sinh':
            stretch = SinhStretch()
        elif newStretch == 'asinh':
            stretch = AsinhStretch()
        elif newStretch == 'power':
            stretch = PowerStretch(np.e)
        else:
            print('unknown stretch ...',newStretch)
            stretch = LinearStretch()
        return stretch

    def showImage(self, image):        
        self.oimage = image.copy()
        # Intensity limits
        if self.cmin is None:
            mask =  np.isfinite(image)
            if np.sum(mask) == 0:
                cmin=-10
                cmax=+10
            else:       
                # Check if too many zeros
                vmed0=np.nanmedian(image)
                d0 = np.nanstd(image)
                cmin = vmed0-1*d0
                if cmin < 0:
                    cmin = 0
                cmax = vmed0+5*d0
            self.cmin = cmin
            self.cmax = cmax        
        self.norm = ImageNormalize(vmin=None, vmax=None, stretch=self.stretchFunc(self.stretch))
        # Remove pre-existing image
        try:
            self.image.remove()
        except:
            pass
        #print('aspect is ', self.aspect)
        self.image = self.axes.imshow(image, origin='lower',
                                      cmap=self.colorMap+self.colorMapDirection,
                                      interpolation='nearest',
                                      norm=self.norm, aspect=self.aspect)
        self.fig.colorbar(self.image, cax=self.cbaxes)
        if self.cmin != 0:
            self.image.set_clim([self.cmin,self.cmax])
        else:
            self.cmin, self.cmax = self.image.get_clim()

        # Cursor data format
        def format_coord(x,y):
            """ Redefine how to show the coordinates """
            pixel = np.array([[x, y]], np.float_)
            world = self.wcs.all_pix2world(pixel, 0)                    
            xx = world[0][0]
            yy = world[0][1]
            " Transform coordinates in string "
            radec = SkyCoord(xx*u.deg, yy*u.deg, frame='icrs')
            xx = radec.ra.to_string(u.hour,sep=':',precision=2)
            yy = radec.dec.to_string(sep=':',precision=1)
            return '{:s} {:s} ({:4.0f},{:4.0f})'.format(xx,yy,x,y)
         
        self.axes.format_coord = format_coord
        
    def updateImage(self, image):
        """Update only the image data."""
        self.image.set_data(image)
        self.oimage = image
    
    def updateScale(self,_cmin,_cmax):
        self.image.set_clim([_cmin, _cmax])
        self.cmin = _cmin
        self.cmax = _cmax
        self.fig.canvas.draw_idle()

    def updateNorm(self):
        self.image.set_clim([self.cmin,self.cmax])
        self.fig.canvas.draw_idle()
        
    def refreshImage(self):
        self.fig.canvas.draw_idle()
        self.changed = False

from sospex.moments import histoImage
        
class ImageHistoCanvas(MplCanvas):
    """ Canvas to plot the histogram of image intensity """
    limSignal = pyqtSignal(str)
    levSignal = pyqtSignal(int)
    showLevels = False
    
    def __init__(self, *args, **kwargs):
        MplCanvas.__init__(self, *args, **kwargs)
        self.axes = self.fig.add_axes([0.0,0.4,1.,1.])
        self.axes.yaxis.set_major_formatter(plt.NullFormatter())
        self.axes.spines['top'].set_visible(False)
        self.axes.spines['right'].set_visible(False)
        self.axes.spines['left'].set_visible(False)
        # Start a span selector
        self.span = SpanSelector(self.axes, self.onSelect, 'horizontal', useblit=True,
                                 rectprops=dict(alpha=0.5, facecolor='LightSalmon'),button=1)
        # Initialize contour level vertical lines
        self.lev = []
        self.levels = []
        self._ind = None
        self.changed = True
        
    def compute_initial_figure(self, image=None,xmin=None,xmax=None):
        if image is None:
            ''' initial definition when images are not yet read '''
            pass
        else:
            # Delay the computation of histogram to fist call
            if xmin is not None and xmax is not None:
                self.limits = [xmin, xmax]
            self.median = None
            self.sdev = None
            # Activate focus
            self.setFocusPolicy(Qt.ClickFocus)
            self.setFocus()
        
    def update_figure(self, image=None, percent=None):
        try:
            ima, nbins, xmin, xmax, hmin, hmax, imedian, imin, imax, sdev, epsilon, nh = histoImage(image, percent)
            if hmin != hmax:
                n, self.bins, patches = self.axes.hist(ima, bins=nbins,
                                                       range=(hmin,hmax), fc='k', ec='k')
                # Draw grid (median, median+n*sigma)
                # This should be a collection which could be removed once redrawing the image
                x = imedian
                for i in range(10):
                    if x < imax:
                        self.axes.axvline(x=x,color='black',alpha=0.5)
                        x += sdev                
                x = imedian-sdev
                for i in range(3):
                    if x > imin:
                        self.axes.axvline(x=x,color='black',alpha=0.5)
                        x -= sdev                
                self.min = imin
                self.max = imax
                self.median = imedian
                self.epsilon = epsilon
                self.sdev = sdev
                self.nh = nh
                self.onSelect(xmin,xmax)
            else:
                self.limits = [0, 0]
        except:
            print('Problems with the image')

    def drawLevels(self):
        """ Draw levels as defined in levels"""
        if len(self.levels) > 0:
            for l in self.levels:
                lev = self.axes.axvline(x=l, animated=True, color='cyan')
                self.lev.append(lev)
            #self.span.set_visible(False)
            self.span.active = False
            self.connect()
            self.showLevels = True

    def removeLevels(self):
        """ Delete levels and disconnect """
        if len(self.lev) > 0:
            for lev in self.lev:
                try:
                    lev.remove()
                except:
                    print('level does not exist')
            self.lev = []
            self.levels = []
            self.disconnect()
            #self.span.set_visible(True)
            self.span.active = True
            self.showLevels = False
            print("All levels have been removed")
            
    def onSelect(self,xmin, xmax):
        self.limits = [xmin,xmax]
        try:
            self.shade.remove()
        except:
            pass
        self.limSignal.emit('limits changed')
        self.shade = self.axes.axvspan(self.limits[0], self.limits[1],
                                       facecolor='Lavender', alpha=0.5, linewidth=0)
        # Redefine limits
        x1, x2 = self.axes.get_xlim()
        x2 = self.limits[1] + 6 * self.sdev
        x1 = self.limits[0] - 2 * self.sdev
        self.axes.set_xlim((x1,x2))
        self.fig.canvas.draw_idle()

    def connect(self):
        self.cid_draw = self.fig.canvas.mpl_connect('draw_event', self.draw_callback)
        self.cid_press = self.fig.canvas.mpl_connect('button_press_event',
                                                     self.button_press_callback)
        self.cid_release = self.fig.canvas.mpl_connect('button_release_event',
                                                       self.button_release_callback)
        self.cid_motion = self.fig.canvas.mpl_connect('motion_notify_event',
                                                      self.motion_notify_callback)
        self.cid_key = self.fig.canvas.mpl_connect('key_press_event', self.key_press_callback)
        self.fig.canvas.draw_idle()

    def disconnect(self):
        self.fig.canvas.mpl_disconnect(self.cid_draw)
        self.fig.canvas.mpl_disconnect(self.cid_press)
        self.fig.canvas.mpl_disconnect(self.cid_release)
        self.fig.canvas.mpl_disconnect(self.cid_motion)
        self.fig.canvas.mpl_disconnect(self.cid_key)
        self.fig.canvas.draw_idle()

    def draw_callback(self, event):
        self.background = self.fig.canvas.copy_from_bbox(self.axes.bbox)
        if len(self.lev) > 0:
            for lev in self.lev:
                self.axes.draw_artist(lev)

    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        if not self.showLevels:
            return
        if not event.inaxes:
            return
        if event.button != 1:
            return
        self._ind = self.get_ind_under_point(event)

    def button_release_callback(self, event):
        'whenever a mouse button is released'
        if not self.showLevels:
            return
        if not event.inaxes:
            return
        if event.button != 1:
            return
        # Emit a signal to communicate change of contour (for large images)
        if self.nh > 100000 and self._ind is not None:
            self.levSignal.emit(self._ind)
        self._ind = None
        
    def get_ind_under_point(self, event):
        """get the index of the level if within epsilon tolerance"""
        if not event.inaxes:
            return
        levels = np.array(self.levels)
        d = np.abs(levels - event.xdata)
        indseq, = np.nonzero(d == d.min())
        ind = indseq[0]
        if d[ind] >= self.epsilon:
            ind = None
        # Check if ind is greater than self.levels
        if ind is not None:
            if ind >= np.size(levels):
                ind = np.size(levels)-1
        return ind

    def key_press_callback(self, event):
        'whenever a key is pressed'
        if not event.inaxes:
            return
        if event.key == 't':
            # Toggle between seeing and hide contours
            self.showLevels = not self.showLevels
            for lev in self.lev:
                lev.set_visible(self.showLevels)
            if self.showLevels:
                self.span.active = False
                self.cid_draw = self.fig.canvas.mpl_connect('draw_event', self.draw_callback)
                self.cid_press = self.fig.canvas.mpl_connect('button_press_event',
                                                             self.button_press_callback)
                self.cid_release = self.fig.canvas.mpl_connect('button_release_event',
                                                               self.button_release_callback)
                self.cid_motion = self.fig.canvas.mpl_connect('motion_notify_event',
                                                              self.motion_notify_callback)
            else:
                self._ind = None
                # If levels are not visible, show the span selector
                #self.span.set_visible(True)
                self.span.active = True
                self.fig.canvas.mpl_disconnect(self.cid_draw)
                self.fig.canvas.mpl_disconnect(self.cid_press)
                self.fig.canvas.mpl_disconnect(self.cid_release)
                self.fig.canvas.mpl_disconnect(self.cid_motion)
        elif event.key == 'd':
            ind = self.get_ind_under_point(event)
            if ind is not None:
                # Delete contour level
                self.lev[ind].remove()
                # Remove level from lists
                del self.lev[ind]
                del self.levels[ind]
                # Emit signal to communicate it to images
                self.levSignal.emit(-1000-ind)
                # If there are no more levels, show the span selector
                if len(self.levels) == 0:
                    self.span.active = True
        elif event.key == 'i':
            x = event.xdata
            # Insert new level
            n = 0
            for level in self.levels:
                if x > level:
                    n += 1
            self.levels.insert(n, x)
            lev = self.axes.axvline(x=x,color='cyan',animated='True')
            self.lev.insert(n, lev)
            # Emit signal to communicate it to images (add 1000 to tell that this is a new level)
            self.levSignal.emit(1000+n)
        # Update image
        self.fig.canvas.draw_idle()

    def sort_levels(self):
        # Sort the levels in increasing order
        levels = np.array(self.levels)
        idx = np.argsort(levels)
        self.lev = [self.lev[i] for i in idx]
        self.levels = list(levels[idx])
        print('sorted levels: ', self.levels)
        
    def motion_notify_callback(self, event):
        'on mouse movement - this should be allowed only if the image size is reasonably small'
        if not self.showLevels:
            return
        if self._ind is None:
            return
        if event.button != 1:
            return
        if not event.inaxes:
            return
        x = event.xdata
        # Check to avoid overlapping
        if len(self.levels) > 1:
            if self._ind < len(self.levels) - 1:
                if x >= self.levels[self._ind + 1]:
                    return
            if self._ind > 0:
                if x <= self.levels[self._ind - 1]:
                    return
        self.levels[self._ind] = x
        lev = self.lev[self._ind]
        xl,yl = lev.get_data()
        lev.set_data([x,x],yl)
        self.fig.canvas.restore_region(self.background)
        for lev in self.lev:
            self.axes.draw_artist(lev)
        self.fig.canvas.update()
        #self.fig.canvas.repaint()
        self.fig.canvas.flush_events()
        # Emit a signal to communicate change of contour
        if self.nh <= 100000 and self._ind is not None:
            self.levSignal.emit(self._ind)


class SpectrumCanvas(MplCanvas):
    """ Canvas to plot spectra """
    
    switchSignal = pyqtSignal(str)
    modifyAperture = pyqtSignal(str)
    modifyContinuum = pyqtSignal(str)

    def __init__(self, *args, **kwargs):
        MplCanvas.__init__(self, *args, **kwargs)
        from sospex.lines import define_lines
        self.Lines = define_lines()
        self.fig.set_edgecolor('none')
        self.axes = self.fig.add_axes([0.12,0.17,.8,.70])
        self.axes.format_coord = lambda x, y: "{:8.4f} um  {:10.4f} Jy".format(x,y)        
        # Checks
        self.displayFlux = True
        self.displayUFlux = True
        self.displayAtran = True
        self.displayExposure = True
        self.displayLines = True
        self.displayAuxFlux1 = False # Auxiliary spectral cube
        self.shade = False
        self.regionlimits = None
        self.xunit = 'um'  # Alternatives are THz or km/s
        self.yunit = 'Jy/pix'
        self.xlimits = None
        self.ylimits = None        
        self.axes.spines['top'].set_visible(False)
        self.axes.spines['right'].set_visible(False)
        # Use legend to hide/show lines
        self.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.fig.canvas.mpl_connect('button_release_event', self.onrelease)
        self.dred = None
        self.region = None
        self.guess = None
        self.lguess = None
        self.lines = None
        self.dragged = None
        self.auxiliary1 = None
        self.moments = False
        self.fittedlines = False
        self.fittedaplines = False
        self.rannotation = None
        self.cannotation = None
        self.r = None # radius of aperture
        
    def compute_initial_spectrum(self, name=None, spectrum=None, xmin=None, xmax=None):
        if spectrum is None:
            ''' initial definition when spectrum not yet available '''
        else:
            # Spectrum
            self.name = name
            self.spectrum = spectrum
            self.yunit = self.spectrum.yunit
            self.instrument = spectrum.instrument
            self.drawSpectrum()
            # Activate focus
            self.setFocusPolicy(Qt.ClickFocus)
            self.setFocus()
            
    def drawSpectrum(self):        
        # Initialize
        self.axes.clear()
        self.axes.grid(True, which='both')
        self.axes.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        self.axes.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        s = self.spectrum        
        # Write labels
        #if s.instrument == 'GREAT':
        #    if self.yunit == 'Jy':
        #        self.axes.set_ylabel('Flux [Jy]', picker=True)
        #    else:
        #        self.axes.set_ylabel('Temperature [K]', picker=True)
        #else:
        self.axes.set_ylabel('Flux ['+self.yunit+']')
        ckms = 299792.458  # speed of light in km/s
        if self.xunit == 'um':
            x0 = s.l0 * (1 + s.redshift)
            sy = self.yunit
            self.axes.format_coord = lambda x, y:u"{:8.4f} \u03bcm ({:5.0f} km/s)  {:10.4f} ".format(x, (x / x0 - 1) * ckms, y) + sy 
            self.axes.fmt_xdata = lambda x: "{:.4f}".format(x)        
            self.axes.set_xlabel('$\lambda$ [$\\mu$m]', picker=True)
            self.x = s.wave
            if s.watran is not None:
                self.xa = s.watran
            if s.instrument == 'FIFI-LS':
                self.xr = self.x * (1+s.baryshift)
                if s.watran is not None:
                    self.xar = self.xa * (1+s.baryshift)
            elif s.instrument == 'FORCAST':
                self.xar = self.xa * (1+s.baryshift)
        elif self.xunit == 'THz':
            x0 = s.l0 * (1 + s.redshift)
            sy = self.yunit
            self.axes.format_coord = lambda x, y: u"{:6.4f} THz ({:5.0f} km/s)  {:10.4f} ".format(x, (ckms/x*1.e-3/x0 - 1)*ckms, y)+sy
            self.axes.set_xlabel('$\\nu$ [THz]', picker=True)
            self.x = ckms/s.wave * 1.e-3
            if s.watran is not None:
                self.xa = ckms/s.watran * 1.e-3
            if s.instrument in ['FIFI-LS','FORCAST']:
                self.xr = self.x / (1 + s.baryshift)             
                if s.watran is not None:
                    self.xar = self.xa / (1 + s.baryshift)
                     
        if s.instrument in ['FIFI-LS', 'PACS', 'FORCAST','SPIRE']:
            hiflux = s.flux + 1. * s.eflux  # 1 sigma
            loflux = s.flux - 1. * s.eflux
            # Try to repair NaNs
            try:
                idx = np.isnan(hiflux)
                hiflux[idx] = np.interp(self.x[idx], self.x[~idx], hiflux[~idx])
                idx = np.isnan(loflux)
                loflux[idx] = np.interp(self.x[idx], self.x[~idx], loflux[~idx])
            except:
                print('Everything is NaN')
            self.efluxLine = self.axes.fill_between(self.x, loflux, hiflux,
                                                    color='blue', alpha=0.2)
        self.fluxLine = self.axes.step(self.x, s.flux, color='blue', label='F', 
                                       where='mid', zorder=10)
        self.fluxLayer, = self.fluxLine
        self.contLine = self.axes.plot(self.x, s.continuum, color='skyblue', 
                                       label='Cont', zorder=9)
        self.contLayer, = self.contLine
        # Define limits or adjust to previous limits
        if self.xlimits is None:
            xlim0 = np.min(s.wave)
            xlim1 = np.max(s.wave)
            self.xlimits = (xlim0, xlim1)
            ylim0, ylim1 = self.axes.get_ylim()
            ylim1 *= 1.2
            self.ylimits = (ylim0, ylim1)
        xlim0, xlim1 = self.xlimits
        if self.xunit == 'THz':
            xlim1, xlim0 = ckms / xlim0 * 1.e-3, ckms / xlim1 * 1.e-3
        self.axes.set_xlim(xlim0, xlim1)
        self.axes.set_ylim(self.ylimits)
        # Fake line to have the lines in the legend
        self.linesLine = self.axes.plot([0,0.1], [0,0], color='purple',
                                        alpha=0.4, label='Lines', zorder=11)
        self.linesLayer, = self.linesLine
        # Add redshift value on the plot
        self.zannotation = self.axes.annotate(" cz = {:.1f} km/s".format(ckms*s.redshift),
                                              xy=(-0.15,-0.12), picker=5, xycoords='axes fraction')
        # Add reference wavelength value on the plot
        self.lannotation = self.axes.annotate(" $\\lambda_0$ = {:.4f} $\\mu$m".format(s.l0),
                                              xy=(-0.15,-0.17), picker=5, xycoords='axes fraction')
        if self.auxiliary1:
            print('auxiliary l0 ', self.aux1l0)
            print('auxiliary scale ', self.aux1pixscale)
            self.x1lannotation = self.axes.annotate(" $\\lambda_x$ = {:.4f} $\\mu$m".format(self.aux1l0),
                                              xy=(-0.15,-0.22), picker=5, xycoords='axes fraction', color='cyan')            
        # Check if vel is defined and draw velocity axis  
        try:
            vlims = self.computeVelLimits()         
            try:
                self.fig.delaxes(self.vaxes)
            except:
                pass
            self.vaxes = self.axes.twiny()
            self.vaxes.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
            self.vaxes.set_xlim(vlims)
            self.vaxes.set_xlabel("Velocity [km/s]")
            # Elevate zorder of first axes (to guarantee axes gets the events)
            self.axes.set_zorder(self.vaxes.get_zorder()+1) # put axes in front of vaxes
            self.axes.patch.set_visible(False) # hide the 'canvas' 
        except:
            print('Reference wavelength not defined')

        if s.instrument == 'FIFI-LS':
            try:
                self.fig.delaxes(self.ax2)
                self.fig.delaxes(self.ax3)
                self.fig.delaxes(self.ax4)
            except:
                pass
            self.ax2 = self.axes.twinx()
            self.ax3 = self.axes.twinx()
            self.ax4 = self.axes.twinx()
            self.ax2.set_ylim([0.01,1.1])
            self.ax4.tick_params(labelright='off',right='off')
            self.ax4.set_yticklabels([]) # remove tick labels on the left
            self.atranLine = self.ax2.step(self.xr, s.atran, color='red',
                                           label='Atm', where='mid', zorder=12)
            if s.uatran is not None:
                self.atranLine2 = self.ax2.plot(self.xar, s.uatran, color='red', zorder=14, alpha=0.2)
                self.atranLayer2, = self.atranLine2
            self.exposureLine = self.ax3.step(self.x, s.exposure, where='mid',
                                              color='orange', label='E', zorder=13)
            ymax = np.nanmax(s.flux); ymin = np.nanmin(s.flux)
            yumax = np.nanmax(s.uflux); yumin = np.nanmin(s.uflux)
            if yumax > ymax: ymax=yumax
            if yumin < ymin: ymin=yumin
            self.ax3.set_ylim([0.5,np.nanmax(s.exposure)*1.54])
            self.ufluxLine = self.ax4.step(self.xr, s.uflux, color='green', where='mid',
                                           label='F$_{u}$', zorder=14)
            self.ax4.set_ylim(self.axes.get_ylim())
            self.ufluxLayer, = self.ufluxLine
            self.atranLayer, = self.atranLine
            self.exposureLayer, = self.exposureLine
            lns = self.fluxLine \
                   +self.ufluxLine \
                   +self.atranLine \
                   +self.exposureLine \
                   +self.linesLine
            lines = [self.fluxLayer, 
                     self.ufluxLayer, self.atranLayer,
                     self.exposureLayer, self.linesLayer]
            visibility = [self.displayFlux,
                          self.displayUFlux, 
                          self.displayAtran, self.displayExposure, self.displayLines]            
            # Add axes
            if self.displayExposure:
                self.ax3.get_yaxis().set_tick_params(labelright='on',right='on',
                                  direction='out',pad=5,colors='orange')
            else:
                self.ax3.get_yaxis().set_tick_params(labelright='off',right='off')
            if self.displayAtran:
                self.ax2.get_yaxis().set_tick_params(labelright='on',right='on',
                                  direction='in', pad = -25, colors='red')
            else:
                self.ax2.get_yaxis().set_tick_params(labelright='off',right='off')
        elif s.instrument in ['GREAT','HI','HALPHA','VLA','ALMA','IRAM','CARMA','MMA','PCWI','MUSE','SITELLE']:
            self.displayUFlux = False
            self.displayAtran = False
            self.displayExposure = False
            lns = self.fluxLine + self.linesLine
            lines = [self.fluxLayer,self.linesLayer]
            visibility = [self.displayFlux,self.displayLines]
            if s.instrument == 'GREAT':
                try:
                    self.fig.delaxes(self.Tax)
                except:
                    pass
                self.Taxes = self.axes.twinx()
                #self.Taxes.tick_params(labelright='off',right='off')
                if s.bunit == 'K (Tmb)':
                    if s.yunit == 'Jy':
                        self.Taxes.set_ylabel('T$_{mb}$ [K]')
                    else:
                        self.Taxes.set_ylabel('T$_{mb}$ [K/pix]')
                else:
                    if s.yunit == 'Jy/pix':
                        self.Taxes.set_ylabel('T$^{*}_{A}$ [K]')
                    else:
                        self.Taxes.set_ylabel('T$^{*}_{A}$ [K/pix]')
                ylims = self.axes.get_ylim()
                print('limits in Flux ', ylims)
                print('Tb 2 Jy ', self.spectrum.Tb2Jy)
                print('limits in Tb ', ylims[0]/ self.spectrum.Tb2Jy)
                self.Taxes.set_ylim(ylims[0]/self.spectrum.Tb2Jy, ylims[1]/self.spectrum.Tb2Jy)
        elif s.instrument in ['PACS', 'FORCAST','SPIRE']:
            try:
                self.fig.delaxes(self.ax3)
            except:
                pass
            self.ax3 = self.axes.twinx()
            self.exposureLine = self.ax3.step(self.x, s.exposure, where='mid',
                                              color='orange',label='E',zorder=7)
            self.ax3.set_ylim([np.nanmin(s.exposure)*0.1,np.nanmax(s.exposure)*1.54])
            self.exposureLayer, = self.exposureLine
            self.displayUFlux = False
            self.displayExposure = True
            # Add axes
            if self.displayExposure:
                self.ax3.get_yaxis().set_tick_params(labelright='on', right='on',
                                      direction='out', pad=5, colors='orange')
            else:
                self.ax3.get_yaxis().set_tick_params(labelright='off', right='off')
            if s.instrument == 'FORCAST':
                try:
                    self.fig.delaxes(self.ax2)
                except:
                    pass
                self.ax2 = self.axes.twinx()
                self.atranLine = self.ax2.step(self.xar, s.uatran, color='red', 
                                               label='Atm', where='mid', zorder=8)
                self.ax2.set_ylim([0.01,1.1])
                self.atranLayer, = self.atranLine
                print('Atran layer ', self.atranLayer)
                self.displayAtran = True
                print('Displayed atran ')
                lns = self.fluxLine \
                   +self.atranLine \
                   +self.exposureLine \
                   +self.linesLine
                lines = [self.fluxLayer, self.atranLayer,
                     self.exposureLayer, self.linesLayer]
                visibility = [self.displayFlux, self.displayAtran, 
                          self.displayExposure, self.displayLines]            
                print('Defined visibility')
            else:
                self.displayAtran = False
                lns = self.fluxLine + self.exposureLine + self.linesLine
                lines = [self.fluxLayer, self.exposureLayer, self.linesLayer]
                visibility = [self.displayFlux, self.displayExposure, self.displayLines]
        # It is possible to draw an external cube only if the vaxis is defined
        try:
            try:
                self.fig.delaxes(self.axu1)
            except:
                pass
            if self.auxiliary1:
                # Draw auxiliary flux
                self.axu1 = self.vaxes.twinx()
                #c =  299792.458 # km/s
                v1 = (self.aux1w / self.aux1l0 - 1. - s.redshift) * ckms
                print('aux1l0 ', self.aux1l0)
                self.afluxLine1 = self.axu1.step(v1, self.aflux1, color='cyan', 
                                                 label='F$_{x}$', 
                                                 where='mid', zorder=12)
                # Add scale on the right ?
                self.axu1.get_yaxis().set_tick_params(labelright='on', right='on',
                                          direction='out', pad=5, colors='cyan')
                self.afluxLayer1, = self.afluxLine1
                lns += self.afluxLine1
                lines.append(self.afluxLayer1)
                visibility.append(self.displayAuxFlux1)
        except:
            print('Auxiliary ref. wavelength not defined') 

        # Prepare legend
        print('Prepare legend ')              
        self.labs = [l.get_label() for l in lns]
        print('labels ', self.labs)
        leg = self.axes.legend(lns, self.labs, 
                               bbox_to_anchor=(0.95, -0.20),
                               frameon=False,
                               loc='lower right',
                               ncol=6, borderaxespad=0.,
                               handletextpad=-1, handlelength=0.0)
        # handlelength = 0 makes lines disappear
        # Labels same color as line
        self.textlines = dict()
        self.plottedlines = dict()
        for l, text, line in zip( lns, leg.texts, lines ):
            self.textlines[l.get_label()] = text
            self.plottedlines[l.get_label()] =  line
            text.set_color( l.get_color() )
            text.set_picker(5) # Makes the text of labels pickable
        leg.set_draggable(True)
        self.lined = dict()
        self.labed = dict()
        for legline, origline, txt in zip(leg.get_lines(), lines, leg.texts):
            legline.set_picker(5) # 5pts tolerance
            self.lined[legline] = origline
            self.labed[legline] = txt            

        # Hide lines
        for line, legline, vis in zip(lines, leg.get_lines(), visibility):
            line.set_visible(vis)
            txt = self.labed[legline]
            if vis:
                alpha=1.0
            else:
                alpha=0.2
            legline.set_alpha(alpha)
            txt.set_alpha(alpha)       

        # Shade region considered for the images
        if self.shade == True:
            self.shadeRegion()
        # Add spectral lines
        self.spectralLines()
                
    def spectralLines(self):
        """Draw spectral lines."""
        ckms = 299792.458  # speed of light in km/s
        s = self.spectrum   
        # Add spectral lines
        self.annotations = []
        font = FontProperties(family='DejaVu Sans', size=12)
        xlim0, xlim1 = self.axes.get_xlim()
        ylim0, ylim1 = self.axes.get_ylim()
        #c = 299792458.0  # speed of light in m/s
        if self.xunit == 'THz':
            xlim1, xlim0 = ckms / xlim0 * 1.e-3, ckms / xlim1 * 1.e-3
        dy = ylim1 - ylim0
        for line in self.Lines.keys():
            nline = self.Lines[line][0]
            wline = self.Lines[line][1]*(1.+s.redshift)
            if (wline > xlim0 and wline < xlim1):
                wdiff = abs(s.wave - wline)
                imin = np.argmin(wdiff)
                y = s.flux[imin]
                y1 = y
                if (ylim1-(y+0.1*dy)) > ((y-0.1*dy)-ylim0):
                    y2 = y+0.1*dy
                else:
                    y2 = y-0.1*dy
                if self.xunit == 'um':
                    xline = wline
                elif self.xunit == 'THz':
                    xline = ckms/wline * 1.e-3
                y2 = 0.88  # Axes coordinates
                trans = transforms.blended_transform_factory(self.axes.transData, 
                                                             self.axes.transAxes)
                annotation = self.axes.annotate(nline, xy=(xline,y1),  xytext=(xline, y2),
                                                textcoords = trans,
                                                picker=5,
                                                color='purple', alpha=0.4, ha='center',
                                                arrowprops=dict(edgecolor='purple',facecolor='y', 
                                                                arrowstyle='-',alpha=0.4,
                                                connectionstyle="angle,angleA=0,angleB=90,rad=10"),
                                                rotation = 90, fontstyle = 'italic',
                                                fontproperties=font, visible=self.displayLines)
                annotation.draggable()  # Make annotation draggable
                self.annotations.append(annotation)  

    def computeVelLimits(self):
        """Compute velocity limits."""
        x1, x2 = self.xlimits
        c = 299792.458  # speed of light in km/s
        s = self.spectrum
        vx1 = (x1 / (1 + s.redshift) / s.l0 - 1.) * c  # Velocity at rest
        vx2 = (x2 / (1 + s.redshift) / s.l0 - 1.) * c
        if self.xunit == 'THz':
            return (vx2, vx1)
        else:
            return (vx1, vx2)

    def updateXlim(self):
        """Update xlimits."""
        xlim0, xlim1 = self.xlimits
        vlim = self.computeVelLimits()
        if self.xunit == 'THz':
            c = 299792458.0  # speed of light in m/s
            xlim1, xlim0 = c / xlim0 * 1.e-6, c / xlim1 * 1.e-6
        self.axes.set_xlim(xlim0, xlim1)
        self.vaxes.set_xlim(vlim)
        self.fig.canvas.draw_idle()

    def updateYlim(self,f=None):
        """Update ylimits."""
        # Grab new limits and update flux and 
        ylims = self.ylimits
        self.axes.set_ylim(ylims)
        s = self.spectrum            
        if self.instrument == 'FIFI-LS':
            self.ax4.set_ylim(ylims)
            self.ax3.set_ylim([0.5, np.nanmax(s.exposure) * 1.54])
            self.ax2.set_ylim([0.01, 1.1])
        if self.instrument == 'GREAT':
            self.Taxes.set_ylim(ylims[0]/s.Tb2Jy, ylims[1]/s.Tb2Jy)
        # Adjust line labels if they exist
        try:
           ylim0, ylim1 = ylims
           dy = np.abs(ylim0 - ylim1)
           for a in self.annotations:
              xl, yl = a.xy
              xt, yt = a.get_position()
              # Check unit (wavelength/frequency)
              if self.xunit == 'um':
                 xline = xl
              elif self.xunit == 'THz':
                 c = 299792458.0   
                 xline = c / xl * 1.e-6
              wdiff = abs(s.wave - xline)
              imin = np.argmin(wdiff)
              if f is None:
                 y = s.flux[imin]
              else:
                 y = f[imin]
              y1 = y
              if (ylim1 - (y + 0.2 * dy)) > ((y - 0.2 * dy) - ylim0):
                 y2 = y + 0.2 * dy
              else:
                 y2 = y - 0.2 * dy
              if y1 > ylim1:
                 y1 = ylim1
              if y2 > ylim1:
                 y2 = ylim1 - dy / 20.
              a.xy = (xl, y1)
        except:
           pass
        self.fig.canvas.draw_idle()
    
    def repairNaN(f):
        x = np.arange(len(f))
        idx = np.isnan(f)
        if np.sum(idx) > 0:
            f[idx] = np.interp(x[idx], x[~idx], f[~idx])
        return f

    def updateSpectrum(self, f=None, ef=None, uf=None, exp=None, cont=None, cslope=None, af=None,
                       moments=None, noise=None, atran=None, lines=None, aplines=None, ncell=0):
        try:
            # Remove moments
            try:
                self.arrow1.remove()
                self.arrow2.remove()
                self.gauss.remove()
            except:
                pass
            # Remove fitted lines
            try:
                for voigt in self.voigt:
                    voigt.remove()
            except:
                pass
            try:
                for gauss in self.apfit:
                    gauss.remove()
            except:
                pass
            if atran is not None:
                self.atranLine[0].set_ydata(atran)
                self.axes.draw_artist(self.atranLine[0])
            if cont is not None:
                self.spectrum.continuum = cont
                self.spectrum.cslope = cslope
                self.contLine[0].set_ydata(cont)
                self.axes.draw_artist(self.contLine[0])
            if uf is not None:
                self.ufluxLine[0].set_ydata(uf)
                self.ax4.draw_artist(self.ufluxLine[0])
            if exp is not None:
                self.exposureLine[0].set_ydata(exp)
                self.ax3.draw_artist(self.exposureLine[0])
            if f is not None:
                self.fluxLine[0].set_ydata(f)
                self.axes.draw_artist(self.fluxLine[0])
                ylim0, ylim1 = self.axes.get_ylim()
                xlim0, xlim1 = self.axes.get_xlim()
                n = len(f)
                if self.spectrum.instrument == 'FIFI-LS':
                    # Ignore regions with atm transmission < 0.2
                    atm = self.spectrum.atran
                    mask = (atm > 0.2) & np.isfinite(f)
                    n5 = n // 10
                else:
                    mask = np.ones(n, dtype=bool)
                    n5 = n // 20
                mask[0:n5] = 0
                mask[-n5:] = 0
                if np.sum(mask) > 0: 
                    idx = (mask == 1) & (self.x > xlim0) & (self.x < xlim1)
                    maxf = np.nanmax(f[idx])
                    minf = np.nanmin(f[idx])
                else:
                    maxf = ylim1 / 1.1
                    minf = ylim0 
                if self.spectrum.instrument == 'FIFI-LS':
                    if uf is not None:
                        n = len(uf) // 5
                        umaxf = np.nanmax(uf[n:-n])
                        umaxf = np.nanmax(uf[idx])
                        if umaxf > maxf:
                            maxf = umaxf
                self.ylimits = (minf, maxf*1.2)
                self.axes.set_ylim(self.ylimits)
                self.updateYlim(f=f)
                self.updateGuess(f, ncell)
            if ef is not None:
                path = self.efluxLine.get_paths()[0]
                v = path.vertices
                xf = self.fluxLine[0].get_xdata()
                n1 = (np.where(xf == np.min(v[:,0])))[0][0]
                n2 = (np.where(xf == np.max(v[:,0])))[0][0] 
                if self.xunit == 'THz':
                    y1 = (f - ef)[n2:n1+1]
                    y2 = (f + ef)[n2:n1+1]
                else:
                    y1 = (f - ef)[n1:n2+1]
                    y2 = (f + ef)[n1:n2+1]
                v1 = np.concatenate(([y2[0]],y1,[y2[-1]],y2[::-1],[y2[0]]))
                v1[np.isnan(v1)] = np.nanmedian(v1) # Put medians on NaNs
                v[:,1] = v1
                self.axes.draw_artist(self.efluxLine)
            if af is not None:
                self.afluxLine1[0].set_ydata(af)
                self.vaxes.draw_artist(self.afluxLine1[0])
                # If similar wavelength, it should have same limits as main cube
                dlam = np.abs(self.spectrum.l0 - self.aux1l0)
                print('Update spectrum: dlam ', dlam)
                print('aux1 pixel scale ', self.aux1pixscale)
                dscale = 1 - self.aux1pixscale/self.pixscale
                print('Update spectrum: dscale ', dscale)
                if (dlam < 1):
                    if dscale > 0.1:
                        rescale = (self.aux1pixscale/self.pixscale) ** 2
                        print('rescaling ', rescale)
                        self.afluxLine1[0].set_ydata(af * rescale)
                    self.axu1.set_ylim(self.ylimits)
                else:
                    minf = np.nanmin(af)
                    maxf = np.nanmax(af)
                    self.axu1.set_ylim(minf, maxf*1.1)
            if moments is not None:
                self.moments = True
                # Update position, size, and dispersion from moments
                x = moments[1]
                y = np.nanmedian(cont) # self.spectrum.continuum
                # FWHM of the distribution (assuming Gaussian)
                dx = np.sqrt(2.*np.log(2))* np.sqrt(moments[2])
                # Amplitude of a Gaussian ...
                c = 299792458. # m/s
                A = moments[0] * 1.e26 / np.sqrt(2. * np.pi*moments[2]) * x * x / c * 1.e-6
                style = 'wedge'
                if self.xunit == 'THz':
                    dx = (c*dx)/(x*x-dx*dx) * 1.e-6
                    x = c/x * 1.e-6
                self.arrow1 = FancyArrowPatch((x,y),(x,y+A),arrowstyle=style,
                                              mutation_scale=1.0,color='skyblue')
                self.axes.add_patch(self.arrow1)
                self.arrow2 = FancyArrowPatch((x-dx,y+0.5*A),(x+dx,y+0.5*A),arrowstyle=style,
                                              mutation_scale=1.0,color='skyblue')
                self.axes.add_patch(self.arrow2)
                # Gauss patch
                sig = dx / (np.sqrt(2.*np.log(2)))
                xx = np.arange(x-4*sig,x+4*sig,dx/10.)    
                # Alternatively we can compute the continuum on the range
                y = np.interp(xx, self.contLine[0].get_xdata(), cont)
                gauss = y + A * np.exp(-(xx-x)**2/(2*sig*sig))
                verts = list(zip(xx,gauss))
                self.gauss = Polygon(verts, fill=False, closed=False,color='skyblue')
                self.axes.add_patch(self.gauss)
            if lines is not None:
                self.fittedlines = True
                self.voigt = []
                for iline, line in enumerate(lines):
                    #print('line is ', line)
                    x, sigma, A, alpha, xerr, sigmaerr, Aerr = line
                    c = 299792458. # m/s
                    A = A * 1.e20  * x * x / c  # Retransform in Jy
                    #print('Plot iline, A ', iline, A)
                    # Add redshift
                    #z = self.spectrum.redshift
                    #x /= (1 + z)
                    #sigma /= (1 + z)
                    # Voigt patch 
                    xx = np.arange(x-7*sigma,x+7*sigma,sigma/10.) 
                    if alpha == 0:
                        dx = xx - x
                        A /= np.sqrt(2*np.pi) * sigma
                        y = np.interp(xx, self.contLine[0].get_xdata(), cont)
                        model = y + A * np.exp(-0.5 * (dx / sigma)**2)
                    else:                    
                        sigmag = sigma/np.sqrt(2*np.log(2.))
                        xx2 = (xx-x)**2
                        gauss = np.exp(-xx2/(2*sigmag*sigmag))/(np.sqrt(2*np.pi) * sigmag)
                        cauchy = sigma/np.pi/(xx2+sigma*sigma)
                        y = np.interp(xx, self.contLine[0].get_xdata(), cont)
                        model = y + A * ((1-alpha)* gauss + alpha*cauchy)
                    # Case of Frequency
                    if self.xunit == 'THz':
                        xx = c/xx * 1.e-6
                    verts = list(zip(xx,model))
                    voigt = Polygon(verts, fill=False, closed=False,color='purple',linewidth=3)
                    self.axes.add_patch(voigt)
                    self.voigt.append(voigt)
            if aplines is not None:
                self.aplines = aplines # Conserve the parameters of fitted lines
                self.fittedaplines = True
                self.apfit = []
                for line in aplines:
                    if self.function == 'Voigt':
                        c0, ec0, slope, x, ex, A, eA, sigma, esigma, alpha = line                        
                        c = 299792458. # m/s
                        A = A * 1.e20  * x * x / c  # Retransform in Jy
                        xx = np.arange(x-7*sigma, x+7*sigma, sigma/10.) 
                        dx = xx - x
                        model = c0 + slope * dx
                        s2 = sigma * sigma
                        sg2 = s2/(2*np.log(2.))
                        sg = np.sqrt(2 * np.pi * sg2)
                        dx2 = dx * dx
                        gauss = np.exp(-0.5 * dx2 / sg2) / sg
                        cauchy = sigma / ( np.pi * (dx2 + s2))
                        #factor = (1-alpha)/np.sqrt(np.pi/np.log(2)) + alpha/np.pi
                        #A *= sigma / factor
                        model += A * ((1 - alpha) * gauss + alpha * cauchy) 
                    else:
                        c0, ec0, slope, x, ex, A, eA, sigma, esigma = line
                        c = 299792458. # m/s
                        A = A * 1.e20  * x * x / c  # Retransform in Jy
                        xx = np.arange(x-7*sigma, x+7*sigma, sigma/10.) 
                        dx = xx - x
                        A /= np.sqrt(2*np.pi) * sigma
                        model = c0 + slope * dx + A * np.exp(-0.5 * (dx / sigma)**2)
                    # Case of Frequency
                    if self.xunit == 'THz':
                        c = 299792458. # m/s
                        xx = c / xx * 1.e-6
                    verts = list(zip(xx, model))
                    gauss = Polygon(verts, fill=False, closed=False,
                                    color='violet', linewidth=2)
                    self.axes.add_patch(gauss)
                    self.apfit.append(gauss)
            else:
                self.aplines = None
                self.apfit = None
        except:
            print('Failed to update spectrum')
            pass
        
    def updateGuess(self, f, ncell):
        if self.guess is not None:
            # Reposition the guess to see the markers
            g = self.guess
            x, y = zip(*g.xy)
            x = self.xguess[ncell]
            if self.xunit == 'THz':
                mask = ((self.x < x[0]) & (self.x > x[1])) | ((self.x < x[2]) & (self.x > x[3]))
            else:   
                mask = ((self.x > x[0]) & (self.x < x[1])) | ((self.x > x[2]) & (self.x < x[3]))                        
            med = np.nanmedian(f[mask])
            ymed = np.nanmedian(y)
            dy = med - ymed
            y += dy
            yc = np.nanmedian(y)
            if np.isfinite(yc):
                g.xy = [(i,j) for (i,j) in zip(x,y)]
                g.updateLinesMarkers()
            # Update line guesses
            try:
                if len(self.lines) > 0:
                    for i, line in enumerate(self.lines):
                        if line is None:
                            pass
                        else:
                            gg = self.lguess[line.n][ncell]
                            line.x0 = gg[0]
                            line.fwhm = gg[1]
                            sigma = line.fwhm / 2.355
                            idx = np.abs(self.x - line.x0) < sigma
                            indmin = np.nanargmin(np.abs(self.x-line.x0))
                            yc = self.spectrum.continuum[indmin]
                            if np.isfinite(yc):
                                line.c0 = yc
                                line.cs = self.spectrum.cslope
                                fc = f[idx] - self.spectrum.continuum[idx]
                            else:
                                # compute continuum from median of regions
                                line.c0 = med
                                line.cs = 0.
                                fc = f[idx] - med
                            #A = f[indmin]-line.c0
                            if line.A > 0:
                                A = np.nanmax(fc)
                                if A > 0:
                                    line.A = A
                                else:
                                    pass
                                    #line.A = 0.01
                            else:
                                A = np.nanmin(fc)
                                if A < 0:
                                    line.A = A
                                else:
                                    pass
                                    #line.A = -0.01
                            line.updateCurves()
            except BaseException:
                pass

            
    def shadeRegion(self, limits = None, color=None):
        if limits == None:
            wmin,wmax = self.regionlimits
        else:
            wmin,wmax = limits
        if color == None:
            color = 'Lavender'            
        if self.xunit == 'um':
            xmin = wmin
            xmax = wmax
        elif self.xunit == 'THz':
            c = 299792458.0  # speed of light in m/s
            xmax = c/wmin * 1.e-6
            xmin = c/wmax * 1.e-6
        # Select axes to appear in the background
        ax = self.axes
        zord = self.axes.get_zorder()
        try:
           zord2 = self.ax2.get_zorder()
           if zord2 < zord:
              zord = zord2
              ax = self.ax2
        except:
           pass
        try:
            zord3 = self.ax3.get_zorder()
            if zord3 < zord:
               zord = zord3
               ax = self.ax3
        except:
           pass
        try:
           zord4 = self.ax4.get_zorder()
           if zord4 < zord:
               ax = self.ax4
        except:
           pass
        if color == 'Lavender':
            self.region = ax.axvspan(xmin,xmax,facecolor=color,alpha=1,linewidth=0,zorder=1)
        else:
            self.tmpRegion = ax.axvspan(xmin,xmax,facecolor=color,alpha=1,linewidth=0,zorder=1)

    def shadeSpectrum(self):
        """Shade part of the spectrum used for different operations."""
        # Clean previous region
        try:
            self.region.remove()
        except:
            pass
        # Shade new region
        self.shadeRegion()
        self.shade = True

    def onpick(self, event):
        """React to onpick events."""
        # print('picked event ', event.artist)
        if isinstance(event.artist, Line2D):
            pass
        elif isinstance(event.artist, Text):
            text = event.artist.get_text()
            if event.artist == self.zannotation:
                c = 299792.458 #km/s
                znew = self.getDouble(self.spectrum.redshift * c)
                if znew is not None:
                    znew /= c
                    if znew != self.spectrum.redshift:
                        self.spectrum.redshift = znew
                        for annotation in self.annotations:
                            annotation.remove()
                        self.zannotation.remove()
                        self.lannotation.remove()
                        if self.auxiliary1:
                            self.x1lannotation.remove()
                        self.drawSpectrum()
                        self.fig.canvas.draw_idle()
                        print('Updated spectrum after z change ')
                        # Simulate a release to activate the update of redshift in main program
                        QTest.mouseRelease(self, Qt.LeftButton)
            elif event.artist == self.lannotation:
                lnew = self.getlDouble(self.spectrum.l0)
                if lnew is not None:
                    if lnew != self.spectrum.l0:
                        self.spectrum.l0 = lnew
                        for annotation in self.annotations:
                            annotation.remove()
                        self.lannotation.remove()
                        self.zannotation.remove()
                        self.drawSpectrum()
                        self.fig.canvas.draw_idle()
                        # Simulate a release to activate the update of ref.wavelength in main
                        QTest.mouseRelease(self, Qt.LeftButton)
            elif event.artist == self.rannotation:
                rnew = self.getrDouble(self.r)
                if rnew is not None:
                    if rnew != self.r:
                        self.r = rnew
                        self.modifyAperture.emit(str(self.r))
            elif event.artist == self.cannotation:
                contChoice = self.getContChoice()
                if contChoice is not None:
                    self.guess.fixedContinuum = contChoice
                    if contChoice:
                        self.cannotation.set_text('C Fixed')
                    else:
                        self.cannotation.set_text('C Variable')
                    # If fixed is chosen, ask for value
                    if self.guess.fixedContinuum:
                        self.contLev, okPressed = QInputDialog.getDouble(self, "Continuum level ", "Continuum [Jy]", 0, -10000., 50000., 2)
                        self.guess.contLevel = self.contLev
                        xy = self.guess.xy
                        xg,yg = zip(*xy)
                        xg = np.array(xg); yg = np.array(yg)
                        self.guess.xy = [(x,self.contLev) for x in xg]
                        self.guess.updateLinesMarkers()
                    # Redo fit
                    self.fig.canvas.draw_idle()
                    self.modifyContinuum.emit(str(self.guess.fixedContinuum))
            elif text in ['F', 'F$_{u}$', 'Atm', 'E', 'Lines','F$_{x}$']:
                if text == 'F':
                    self.displayFlux =  not self.displayFlux
                    state = self.displayFlux
                    #label = 'Flux'
                    if self.instrument in ['FIFI-LS', 'PACS','SPIRE']:
                        self.efluxLine.set_visible(self.displayFlux)
                elif text == 'F$_{u}$':
                    self.displayUFlux = not self.displayUFlux
                    state = self.displayUFlux
                    #label = 'Uflux'
                elif text == 'Atm':
                    self.displayAtran = not self.displayAtran
                    try:
                        self.atranLayer2.set_visible(self.displayAtran)
                    except:
                        pass
                    state = self.displayAtran
                    #label = 'Atm'
                    self.ax2.get_yaxis().set_visible(self.displayAtran)
                elif text == 'E':
                    self.displayExposure = not self.displayExposure
                    state = self.displayExposure
                    #label = 'Exp'
                    self.ax3.get_yaxis().set_visible(self.displayExposure)
                elif text == 'F$_{x}$':
                    self.displayAuxFlux1 = not self.displayAuxFlux1
                    #label = 'F$_{x}$'
                    state = self.displayAuxFlux1
                elif text == 'Lines':
                    #label = 'Lines'
                    self.displayLines = not self.displayLines
                    state = self.displayLines
                    self.setLinesVisibility(self.displayLines)
                # Now redisplay
                if self.shade == True:
                    self.shadeRegion()
                txt = self.textlines[text]
                if state:
                    txt.set_alpha(1.0)
                else:
                    txt.set_alpha(0.2)
                self.plottedlines[text].set_visible(state)
                self.fig.canvas.draw_idle()
            else:
                goNext = True
                if self.auxiliary1 is not None:
                    if event.artist == self.x1lannotation:
                        lnew = self.getlDouble(self.aux1l0)
                        if lnew is not None:
                            if lnew != self.aux1l0:
                                self.aux1l0 = lnew
                                for annotation in self.annotations:
                                    annotation.remove()
                                self.x1lannotation.remove()
                                self.lannotation.remove()
                                self.zannotation.remove()
                                self.drawSpectrum()
                                self.fig.canvas.draw_idle()
                                # Simulate a release to activate the update of ref.wavelength in main
                                QTest.mouseRelease(self, Qt.LeftButton)
                        goNext = False
                if goNext:
                    if text == '$\lambda$ [$\mu$m]' or text == '$\\nu$ [THz]':
                        self.switchUnits()
                        self.switchSignal.emit('switched x unit')
                    elif text == 'Flux ['+self.yunit+']' or text == 'Temperature [K]':
                        self.switchFluxUnits()
                        self.switchSignal.emit('switched y unit')
                    else:
                        self.dragged = event.artist
                        self.pick_pos = event.mouseevent.xdata
        else:
            pass
        return True
    
    def pickSwitch(self, artist, mouseevent):
        self.switchUnits()
        self.switchSignal.emit('switched x unit')
        
    def switchUnits(self):
        c = 299792458.0  # speed of light in m/s
        if self.xunit == 'um':
            self.xunit = 'THz'
            if self.yunit == 'Jy':
                self.axes.format_coord = lambda x, y: "{:6.4f} THz {:10.4f} Jy".format(x, y)
            #else:
            #    self.axes.format_coord = lambda x, y: "{:6.4f} THz {:10.4f} K".format(x, y)
        else:
            self.xunit = 'um'
            if self.yunit == 'Jy':
                self.axes.format_coord = lambda x, y: "{:8.4f} um {:10.4f} Jy".format(x, y)
            #else:
            #    self.axes.format_coord = lambda x, y: "{:8.4f} um {:10.4f} K".format(x, y)
        if self.guess is not None:
            # Switch xguess units
            for x in self.xguess:
                x = c / np.array(x) * 1.e-6  # um to THz or viceversa
            # Switch guess
            self.guess.switchUnits()
        # Lines
        if self.lguess is not None:
            for line in self.lines:
                if line is not None:
                    line.switchUnits()
        # Moments computed
        if self.moments:
            xy = self.gauss.get_xy()
            x, y = zip(*xy)
            x = np.asarray(x)
            y = np.asarray(y)
            x = c / x * 1.e-6
            self.gauss.set_xy(list(zip(x, y)))
        # Line fits computed
        if self.fittedlines:
            for voigt in self.voigt:
                xy = voigt.get_xy()
                x, y = zip(*xy)
                x = np.asarray(x)
                y = np.asarray(y)
                x = c / x * 1.e-6
                voigt.set_xy(list(zip(x, y)))
        if self.fittedaplines:
            for gauss in self.apfit:
                xy = gauss.get_xy()
                x, y = zip(*xy)
                x = np.asarray(x)
                y = np.asarray(y)
                x = c / x * 1.e-6
                gauss.set_xy(list(zip(x, y)))
                
        # Redrawing is maybe excessive
        self.drawSpectrum()
        self.fig.canvas.draw_idle()
        
    def switchFluxUnits(self):
        if self.yunit == 'Jy':
            self.yunit = 'K'
            if self.xunit == 'THz':
                self.axes.format_coord = lambda x, y: "{:6.4f} THz {:10.4f} K".format(x, y)
            else:
                self.axes.format_coord = lambda x, y: "{:8.4f} um {:10.4f} K".format(x, y)
            self.spectrum.flux /= self.spectrum.Tb2Jy 
            y0, y1 = self.ylimits
            y0 /= self.spectrum.Tb2Jy
            y1 /= self.spectrum.Tb2Jy
            self.ylimits = (y0, y1)
        else:
            self.yunit = 'Jy'
            if self.xunit == 'THz':
                self.axes.format_coord = lambda x, y: "{:6.4f} THz {:10.4f} Jy".format(x, y)
            else:
                self.axes.format_coord = lambda x, y: "{:8.4f} um {:10.4f} Jy".format(x, y)
            self.spectrum.flux *= self.spectrum.Tb2Jy 
            y0, y1 = self.ylimits
            y0 *= self.spectrum.Tb2Jy
            y1 *= self.spectrum.Tb2Jy
            self.ylimits = (y0, y1)
        # Redraw
        self.drawSpectrum()
        self.fig.canvas.draw_idle()

    def setLinesVisibility(self, visibility=True):
        for annotation in self.annotations:
                annotation.set_visible(visibility)

    def getDouble(self, z):
        znew, okPressed = QInputDialog.getDouble(self, "Redshift", "cz", z, -10000., 50000., 2)
        if okPressed:
            return znew
        else:
            return None

    def getlDouble(self, l):
        lnew, okPressed = QInputDialog.getDouble(self, "Reference wavelength",
                                                 "Ref. wavelength [um]", l, 0., 1000000., 4)
        if okPressed:
            return lnew
        else:
            return None

    def getrDouble(self, r):
        rnew, okPressed = QInputDialog.getDouble(self, "Aperture radius",
                                                 "r [arcsec] ", r, 0.1, 100., 4)
        if okPressed:
            return rnew
        else:
            return None
        
    def getContChoice(self):

        msgBox = QMessageBox()
        msgBox.setText('Perform the fit with:')
        msgBox.addButton('Variable continuum', QMessageBox.ActionRole)
        msgBox.addButton('Fixed continuum', QMessageBox.ActionRole)
        msgBox.addButton('Cancel choice', QMessageBox.ActionRole)
        self.result = msgBox.exec()
        if self.result == 0:
            return False
        elif self.result == 1:
            return True
        else:
            return None

    def onrelease(self, event):
        if self.dragged is not None and self.pick_pos is not None:
            x1 = event.xdata
            x0 = self.pick_pos
            w0 = np.array([self.Lines[line][1]  for  line in self.Lines.keys()])
            wz = w0 * (1. + self.spectrum.redshift)
            if self.xunit == 'um':
                z = (x1 - x0) / x0
                l0 = x0
            elif self.xunit == 'THz':
                z = (x0 - x1) / x1
                c = 299792458.0  # speed of light in m/s
                l0 = c / x0 * 1.e-6
            wdiff = abs(l0 - wz)
            self.spectrum.l0 = (w0[(wdiff == wdiff.min())])[0]
            self.spectrum.redshift = (1. + self.spectrum.redshift) * (1 + z) - 1.
            for annotation in self.annotations:
                annotation.remove()
            self.zannotation.remove()
            self.lannotation.remove()
            self.drawSpectrum()
            self.fig.canvas.draw_idle()
            self.dragged = None
        return True


from lmfit import Parameters, minimize
def residualsPsf(p, dis, data=None, err=None):
    v = p.valuesdict()
    s = v['s']
    A = v['A']
    d2 = (dis/s)**2
    model = A *  np.exp(-0.5 * d2)
    
    if data is None:
        return model
    else:
        if err is None:
            return (model - data.flatten())
        else:
            return (model - data.flatten())/err.flatten()
             
def residualsMoffat(p, dis, data=None, err=None):
    v = p.valuesdict()
    Io = v['Io']
    alpha = v['alpha']
    beta = v['beta']
    model = Io * (1 + (dis/alpha)**2)**(-beta)
    if data is None:
        return model
    else:
        if err is None:
            return (model - data.flatten())
        else:
            return (model - data.flatten())/err.flatten()
    
class PsfCanvas(MplCanvas):
    """ Canvas to plot a PSF """
    
    switchSignal = pyqtSignal(str)
    def __init__(self, *args, **kwargs):
        MplCanvas.__init__(self, *args, **kwargs)
        self.fig.set_edgecolor('none')
        self.axes = self.fig.add_axes([0.12,0.15,.8,.78])
        self.axes.format_coord = lambda x, y: "{:8.4f} arcsec  {:10.4f} ".format(x,y)        
        # Checks
        self.xlimits = None
        self.ylimits = None        
        self.axes.spines['top'].set_visible(False)
        self.axes.spines['right'].set_visible(False)
        
    def compute_initial_psf(self, distance=None, flux=None, xmax=None, 
                            instrument=None, channel=None, order=None,
                            w=None, pix=None):
        self.name = 'PSF'
        self.w = w
        self.instrument = instrument
        self.channel = channel
        self.order = order
        self.pix = pix
        if flux is None:
            ''' initial definition when psf not yet available '''
        else:
            # Spectrum
            self.flux = flux
            self.distance = distance
            self.drawPSF()
            # Activate focus
            self.setFocusPolicy(Qt.ClickFocus)
            self.setFocus()
            
    def drawPSF(self):        
        # Initialize
        self.axes.clear()
        self.axes.set_ylabel('Flux')
        self.axes.set_xlabel('Distance [arcsec]')
        self.axes.grid(True, which='both')
        self.axes.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        self.axes.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        self.pointScatter, = self.axes.plot(self.distance, self.flux, '.', color='blue', label='Data')
        
        # Overplot profile of unresolved line
        if (self.w is None) or (self.instrument is None):
            pass
        else:
            self.unresolvedPsf()
            if self.sigma is None:
                pass
                self.unresolved = None
            else:
                xg = np.arange(0, np.nanmax(self.distance), 0.1)
                yg = np.nanmax(self.flux) * np.exp(- 0.5 * (xg/self.sigma)**2)
                self.unresolved, = self.axes.plot(xg, yg, color = 'red',
                                                  label='FWHM [PSF] = {:2.1f}'.format(self.sigma*2.355))
                        
        # Fit a Moffat function on data
        Io, alpha, beta, fwhm = self.fitPsf()
        print('Moffat with I0 ',Io,' alpha ',alpha,' beta ',beta)
        xf = np.arange(0, np.nanmax(self.distance), 0.1)
        yf = Io * (1 + (xf/alpha)**2)**(-beta)
        self.resolved, = self.axes.plot(xf, yf, color = 'green',
                                        label='FWHM [Fit] = {:2.1f}'.format(fwhm))
        
        self.axes.legend()

    def updatePSF(self, distance, flux):
        # Modify data
        self.distance = distance
        self.flux = flux
        maxflux = np.nanmax(flux)
        self.pointScatter.set_data(distance, flux)
        # Plot new data
        if self.unresolved is not None:
            yg = self.unresolved.get_ydata()
            self.unresolved.set_ydata(yg/np.nanmax(yg) * maxflux)
        
        Io, alpha, beta, fwhm = self.fitPsf()
        xf = np.arange(0, np.nanmax(self.distance), 0.1)
        yf = Io * (1 + (xf/alpha)**2)**(-beta)
        self.resolved.set_data(xf, yf)
        self.resolved.set_label('FWHM [Fit] = {:2.1f}'.format(fwhm))
        self.axes.legend()  # Regenerate legend
        # Reset y limits if needed
        ylims = self.axes.get_ylim()
        if ylims[1] < maxflux:
            self.axes.set_ylim(ylims[0], maxflux*1.1)
        if ylims[1] > maxflux/2.:
            self.axes.set_ylim(ylims[0], maxflux*1.1)
        self.fig.canvas.draw_idle()
        return Io, alpha, beta, fwhm

    def unresolvedPsf(self):
        if self.instrument in ['FIFI-LS', 'GREAT']:
            m1diam = 2.500 # Effective diameter of SOFIA
            #m2diam = 0.350
        elif self.instrument in ['PACS','SPIRE']:
            m1diam = 3.500 # Diameter of Herschel mirror
            #m2diam = 0.308
        else:
            m1diam = None
            print('unknown instrument')            
        if m1diam is None:
            self.sigma = None
        else:
            if self.instrument == 'FIFI-LS':
                tFWHM = 1.013 * self.w *1.e-6/m1diam * 3600. * 180/np.pi
                if self.channel == 'RED':
                    iFWHM =  3.55*(np.sqrt(0.0156*0.0116)*self.w + np.sqrt(1.3214*1.6466)) - 4.5
                else:
                    if self.order == '1':
                        iFWHM = 3.55 * (np.sqrt(0.0151*0.0179)*self.w + np.sqrt(0.8621*0.2169)) - 2.
                    else:
                        iFWHM = 3.55 * (np.sqrt(0.0057*0.0064)*self.w + np.sqrt(1.5699*1.0353)) - 1.0
                FWHM = np.sqrt(tFWHM * tFWHM + iFWHM * iFWHM)
                self.sigma = FWHM / 2.355
            else:
                self.sigma = 1.22*self.w*1e-6/m1diam*180/np.pi*3600 / 2.355
        
    def fitPsfGauss(self):
        '''Fit a Gaussian on data'''
        fit_params = Parameters()
        fit_params.add('A', value=np.nanmax(self.flux))
        if self.sigma is None:
            sigmaguess = 2
        else:
            sigmaguess = self.sigma
        fit_params.add('s', value=sigmaguess, min=.1, max=2*sigmaguess)
        idx = np.isfinite(self.flux)
        if np.sum(idx) > 10:
            if self.pix is None:
                out = minimize(residualsPsf, fit_params, args=(self.distance[idx],), kws={'data': self.flux[idx],})
            else:
                weight =  (self.distance[idx] * self.pix)
                out = minimize(residualsPsf, fit_params, args=(self.distance[idx],),
                               kws={'data': self.flux[idx], 'err':weight})

            A = out.params['A'].value
            sigma = out.params['s'].value
        else:
            A = np.nan
            sigma = np.nan
            print('Not enough good points')
        return A, sigma*2.355

    def fitPsf(self):
        '''Fit a Moffat function on data'''
        fit_params = Parameters()
        fit_params.add('Io', value=np.nanmax(self.flux))
        if self.sigma is None:
            sigmaguess = 2
        else:
            sigmaguess = self.sigma
        fwhm = sigmaguess * 2.355
        beta = 4.7
        alpha = fwhm * 0.5 / np.sqrt(2**(1 / beta) - 1)
        alphamax = 100 * 0.5 / np.sqrt(2**(1 / beta) - 1)
        alphamin = 0.1 * 0.5 / np.sqrt(2**(1 / beta) - 1)
        fit_params.add('alpha', value=alpha, min=alphamin, max=alphamax)        
        fit_params.add('beta', value=4.7, min=3, max=5)
        idx = np.isfinite(self.flux)
        if np.sum(idx) > 10:
            if self.pix is None:
                out = minimize(residualsMoffat, fit_params, args=(self.distance[idx],), kws={'data': self.flux[idx],})
            else:
                weight =  (self.distance[idx] * self.pix)
                out = minimize(residualsMoffat, fit_params, args=(self.distance[idx],),
                               kws={'data': self.flux[idx], 'err':weight})

            Io = out.params['Io'].value
            alpha = out.params['alpha'].value
            beta = out.params['beta'].value
            fwhm = 2*alpha*np.sqrt(2**(1/beta)-1)
        else:
            Io = np.nan
            alpha = np.nan
            beta = np.nan
            fwhm = np.nan
            print('Not enough good points')
        return Io, alpha, beta, fwhm
