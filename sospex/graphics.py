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

from PyQt5.QtWidgets import (QVBoxLayout, QSizePolicy, QInputDialog, QDialog, QListWidget,
                             QListWidgetItem,QPushButton,QLabel,QMessageBox,QScrollArea,QWidget)
from PyQt5.QtGui import QIcon,QFont
from PyQt5.QtCore import Qt, QSize, pyqtSignal
from PyQt5.QtTest import QTest

# https://het.as.utexas.edu/HET/Software/Astropy-1.0/_modules/astropy/visualization/stretch.html
from astropy.visualization import (LinearStretch, SqrtStretch, SquaredStretch, SinhStretch,AsinhStretch,
                                   LogStretch, ImageNormalize, PowerStretch)

class ScrollMessageBox(QMessageBox):
   def __init__(self, l, *args, **kwargs):
      super().__init__(*args, **kwargs)
      font = QFont()
      font.setFamily("Monaco")
      font.setPointSize(10)
      self.setWindowTitle("Header")
      scroll = QScrollArea(self)
      scroll.setWidgetResizable(True)
      self.content = QWidget()
      scroll.setWidget(self.content)
      lay = QVBoxLayout(self.content)
      for item in l:
          label = QLabel(item, self)
          label.setFont(font)
          lay.addWidget(label)
      self.layout().addWidget(scroll, 0, 0, 1, self.layout().columnCount())
      self.setStyleSheet("QScrollArea{min-width:600 px; min-height: 400px}")

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
    # The Normalize class is largely based on code provided by Sarah Graves.
    # http://www.ster.kuleuven.be/~pieterd/python/html/plotting/interactive_colorbar.html


class NavigationToolbar(NavigationToolbar2QT):
    def __init__(self,canvas,parent):
        # Select only a few buttons
        self.iconDir = os.path.join(os.path.dirname(os.path.abspath(__file__)),"icons")
        self.toolitems = [
            ('Home','Go back to original limits','home','home'),
            ('Pan','Pan figure','move','pan'),
            ('Zoom','Zoom in','zoom_to_rect','zoom'),
        ]
        self.parent = parent
        super().__init__(canvas,parent)

        
class cmDialog(QDialog):
    dirSignal = pyqtSignal(str)    

    def __init__(self, cmlist, stlist, clist, currentCM, currentST, currentCC, parent=None):
        super().__init__()
        path0, file0 = os.path.split(__file__)
        self.setWindowTitle('Colors & Stretch')
        layout = QVBoxLayout()
        label1 = QLabel("Color maps")
        self.list = QListWidget(self)
        iconSize = QSize(144,10)
        self.list.setIconSize(iconSize)
        self.list.setSizePolicy(QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding))
        self.list.setMaximumSize(QSize(160,150))
        self.cmlist = cmlist
        for cm in cmlist:
            #item = QListWidgetItem(self.list)
            #item.setText(cm)
            #item.setIcon(QIcon(path0+"/icons/"+cm+".png"))
            #item = QListWidgetItem(QIcon(path0+"/icons/"+cm+".png"),'',self.list)
            #item.setSizeHint(iconSize)
            QListWidgetItem(QIcon(os.path.join(path0,"icons",cm+".png")),'',self.list)
        # For some reason this does not work when in the stylesheet
        stylesheet = "QListWidget::item {"\
                     +"border-style: solid;"\
                     +"border-width:1px;" \
                     +"border-color:transparent;"\
                     +"background-color: transparent;"\
                     +"color: white;"\
                     +"}"\
                     +"QListWidget::item:selected {"\
                     +"border-style: solid;" \
                     +"border-width:1px;" \
                     +"border-color:black;" \
                     +"background-color: transparent;"\
                     +"color: white;"\
                     +"}"            
        self.list.setStyleSheet(stylesheet)
        n = cmlist.index(currentCM)
        self.list.setCurrentRow(n)
        # Button to reverse color map direction
        b1 = QPushButton("Reverse", self)
        b1.clicked.connect(self.reverse)        
        label2 = QLabel("Stretches")        
        self.slist = QListWidget(self)
        self.slist.setSizePolicy(QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding))
        self.slist.setMaximumSize(QSize(160,150))
        for st in stlist:
            QListWidgetItem(QIcon(os.path.join(path0,"icons",st+"_.png")),st,self.slist)
        n = stlist.index(currentST)
        self.slist.setCurrentRow(n)
        label3 = QLabel("Contour color")        
        self.clist = QListWidget(self)
        iconSize = QSize(144,10)
        self.clist.setIconSize(iconSize)
        self.clist.setSizePolicy(QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding))
        self.clist.setMaximumSize(QSize(160,150))
        for col in clist:
            QListWidgetItem(QIcon(os.path.join(path0,"icons",col+".png")),'',self.clist)
        n = clist.index(currentCC)
        self.clist.setCurrentRow(n)
        self.clist.setStyleSheet(stylesheet)        
        # Button with OK to close dialog
        b2 = QPushButton("OK",self)
        b2.clicked.connect(self.end)
        # Layout
        layout.addWidget(label1)
        layout.addWidget(self.list)
        layout.addWidget(b1)
        layout.addWidget(label2)
        layout.addWidget(self.slist)
        layout.addWidget(label3)
        layout.addWidget(self.clist)
        layout.addWidget(b2)
        self.setLayout(layout)

    def end(self):
        self.close()
        
    def reverse(self):
        self.dirSignal.emit('color map reversed')
      
        
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
            
    def compute_initial_figure(self, image=None, wcs=None, title=None, cMap = 'gist_heat',
                               cMapDir = '_r', stretch='linear', instrument=None):
        self.colorMap = cMap
        self.colorMapDirection = cMapDir
        self.stretch = stretch
        if wcs is None:
            '''initial definition when images are not yet read'''
            pass
        else:
            # import time
            # t = time.process_time()
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
                    title = 'Continuum at ref. wavelength $\lambda_0$'
                elif title == 'M0':
                    title = 'M$_0$'
                elif title == 'M1':
                    title = 'M$_1$'
                elif title == 'M2':
                    title = 'M$_2$'
                elif title == 'M3':
                    title = 'M$_3$'
                elif title == 'M4':
                    title = 'M$_4$'
                elif title == 'sv':
                    title = '$\sigma_v$ [km/s]'
                elif title == 'v':
                   title = 'v [km/s]'
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
                cmin = vmed0-2*d0
                cmax = vmed0+5*d0
            self.cmin = cmin
            self.cmax = cmax        
        self.norm = ImageNormalize(vmin=None, vmax=None, stretch=self.stretchFunc(self.stretch))
        # Remove pre-existing image
        try:
            self.image.remove()
        except:
            pass
        self.image = self.axes.imshow(image, origin='lower',
                                      cmap=self.colorMap+self.colorMapDirection,
                                      interpolation='nearest',
                                      norm=self.norm)
        self.fig.colorbar(self.image, cax=self.cbaxes)
        if self.cmin != 0:
            self.image.set_clim([self.cmin,self.cmax])
        else:
            self.cmin, self.cmax = self.image.get_clim()

        # Cursor data format
        def format_coord(x,y):
            """ Redefine how to show the coordinates """
            pixel = np.array([[x, y]], np.float_)
            world = self.wcs.wcs_pix2world(pixel, 1)                    
            xx = world[0][0]
            yy = world[0][1]
            " Transform coordinates in string "
            radec = SkyCoord(xx*u.deg, yy*u.deg, frame='icrs')
            xx = radec.ra.to_string(u.hour,sep=':',precision=1)
            yy = radec.dec.to_string(sep=':',precision=0)
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
            # Print the histogram of finite values
            ima = image.ravel()
            mask = np.isfinite(ima)
            ima = ima[mask]
            # print('image has size', len(ima))
            self.nh = len(ima)
            ima = np.sort(ima)
            s = np.size(ima)
            smin = int(s*0.005)
            smax = min(int(s*0.9995),s-1)
            nbins=256
            self.median = np.median(ima)
            self.sdev   = np.std(ima[smin:smax])
            self.min    = np.min(ima)
            self.max    = np.max(ima)
            self.epsilon = self.sdev/3.            
            # Define the interval containing 99% of the values
            if xmin == None:
                xmin = ima[int(s*0.01)]
            if xmax == None:
                xmax = ima[int(s*0.99)-1]
            # Avoid excessively lower flux
            if xmin < (self.median - 3 * self.sdev):
                xmin = self.median - 3 * self.sdev
            hmin = ima[smin]
            hmax = ima[smax]
            if (hmax - self.median) > 10 *self.sdev:
                hmax = self.median+10*self.sdev
            n, self.bins, patches = self.axes.hist(ima, bins=nbins,
                                                   range=(hmin,hmax), fc='k', ec='k')            
            if np.isfinite(xmin) and np.isfinite(xmax):
                self.onSelect(xmin,xmax)
            else:
                print('Problems with the image')
            # Draw grid (median, median+n*sigma)
            x = self.median
            for i in range(10):
                if x < self.max:
                    self.axes.axvline(x=x,color='black',alpha=0.5)
                    x += self.sdev                
            x = self.median-self.sdev
            for i in range(3):
                if x > self.min:
                    self.axes.axvline(x=x,color='black',alpha=0.5)
                    x -= self.sdev                
            # Activate focus
            self.setFocusPolicy(Qt.ClickFocus)
            self.setFocus()

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
        x1,x2 = self.axes.get_xlim()
        x2 = self.limits[1] + 5 * self.sdev
        x1 = self.limits[0] - 3 * self.sdev
        self.axes.set_xlim((x1,x2))
        self.fig.canvas.draw_idle()

    def connect(self):
        """ When contours are present """
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
        """ When no contours are present """
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
            n = len(self.lev)
            # Add to contour list
            self.levels.append(x)
            lev = self.axes.axvline(x=x,color='cyan',animated='True')
            self.lev.append(lev)
            # Sort the levels in increasing order
            levels = np.array(self.levels)
            idx = np.argsort(levels)
            self.lev = [self.lev[i] for i in idx]
            self.levels = list(levels[idx])
            # Emit signal to communicate it to images (add 1000 to tell that this is a new level)
            n = self.levels.index(x)
            self.levSignal.emit(1000+n)
        # Update image
        self.fig.canvas.draw_idle()

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
    def __init__(self, *args, **kwargs):
        MplCanvas.__init__(self, *args, **kwargs)
        from sospex.lines import define_lines
        self.Lines = define_lines()
        self.fig.set_edgecolor('none')
        self.axes = self.fig.add_axes([0.12,0.15,.8,.78])
        self.axes.format_coord = lambda x, y: "{:8.4f} um  {:10.4f} Jy".format(x,y)        
        # Checks
        self.displayFlux = True
        self.displayUFlux = True
        self.displayAtran = True
        self.displayExposure = True
        self.displayLines = True
        self.shade = False
        self.regionlimits = None
        self.xunit = 'um'  # Alternatives are THz or km/s
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
        self.dragged = None
        
    def compute_initial_spectrum(self, spectrum=None,xmin=None,xmax=None):
        if spectrum is None:
            ''' initial definition when spectrum not yet available '''
        else:
            # Spectrum
            self.spectrum = spectrum
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
        # Write labels
        self.axes.set_ylabel('Flux [Jy]')
        s = self.spectrum        
        if self.xunit == 'um':
            self.axes.fmt_xdata = lambda x: "{:.4f}".format(x)        
            self.axes.set_xlabel('Wavelength [$\\mu$m]', picker=True)
            self.x = s.wave
            if s.instrument == 'FIFI-LS':
                self.xr = self.x * (1+s.baryshift)
        elif self.xunit == 'THz':
            self.axes.format_coord = lambda x, y: "{:6.4f} THz {:10.4f} Jy".format(x, y)
            c = 299792458.0  # speed of light in m/s
            self.axes.set_xlabel('Frequency [THz]', picker=True)
            self.x = c/s.wave * 1.e-6
            if s.instrument == 'FIFI-LS':
                self.xr = self.x / (1 + s.baryshift)             
        #self.fluxLine = self.axes.step(self.x,s.flux,color='blue',label='Flux',zorder=10)
        self.fluxLine = self.axes.step(self.x,s.flux,color='blue',label='F',zorder=10)
        self.fluxLayer, = self.fluxLine
        self.contLine = self.axes.plot(self.x, s.continuum, color='skyblue',label='Cont',zorder=9)
        self.contLayer, = self.contLine
        # Define limits or adjust to previous limits
        if self.xlimits is None:
            xlim0 = np.min(s.wave)
            xlim1 = np.max(s.wave)
            self.xlimits = (xlim0, xlim1)
            self.ylimits = self.axes.get_ylim()
        xlim0,xlim1 = self.xlimits
        if self.xunit == 'THz':
            c = 299792458.0  # speed of light in m/s
            xlim1, xlim0 = c /xlim0 * 1.e-6, c / xlim1 * 1.e-6
        self.axes.set_xlim(xlim0, xlim1)
        self.axes.set_ylim(self.ylimits)
        # Fake line to have the lines in the legend
        self.linesLine = self.axes.plot([0,0.1], [0,0], color='purple',
                                        alpha=0.4, label='Lines', zorder=11)
        self.linesLayer, = self.linesLine
        # Add redshift value on the plot
        c = 299792.458  #km/s
        self.zannotation = self.axes.annotate(" cz = {:.1f} km/s".format(c*s.redshift),
                                              xy=(-0.15,-0.07), picker=5, xycoords='axes fraction')
        # Add reference wavelength value on the plot
        self.lannotation = self.axes.annotate(" $\\lambda_0$ = {:.4f} $\\mu$m".format(s.l0),
                                              xy=(-0.15,-0.12), picker=5, xycoords='axes fraction')        
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
            self.atranLine = self.ax2.step(self.xr, s.atran, color='red', label='Atm', zorder=12)
            self.exposureLine = self.ax3.step(self.x, s.exposure,
                                              color='orange', label='E', zorder=13)
            ymax = np.nanmax(s.flux); ymin = np.nanmin(s.flux)
            yumax = np.nanmax(s.uflux); yumin = np.nanmin(s.uflux)
            if yumax > ymax: ymax=yumax
            if yumin < ymin: ymin=yumin
            self.ax3.set_ylim([0.5,np.nanmax(s.exposure)*1.54])
            self.ufluxLine = self.ax4.step(self.xr, s.uflux, color='green',
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
            lines = [self.fluxLayer, self.ufluxLayer, self.atranLayer,
                     self.exposureLayer, self.linesLayer]
            visibility = [self.displayFlux, self.displayUFlux, self.displayAtran,
                          self.displayExposure, self.displayLines]            
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
        elif s.instrument == 'GREAT':
            self.displayUFlux = False
            self.displayAtran = False
            self.displayExposure = False
            lns = self.fluxLine + self.linesLine
            lines = [self.fluxLayer,self.linesLayer]
            visibility = [self.displayFlux,self.displayLines]
        elif s.instrument == 'PACS':
            try:
                self.fig.delaxes(self.ax3)
            except:
                pass
            self.ax3 = self.axes.twinx()
            self.exposureLine = self.ax3.step(self.x, s.exposure,
                                              color='orange',label='Exp',zorder=7)
            self.ax3.set_ylim([0.5,np.nanmax(s.exposure)*1.54])
            self.exposureLayer, = self.exposureLine
            self.displayUFlux = False
            self.displayAtran = False
            self.displayExposure = True
            lns = self.fluxLine + self.exposureLine + self.linesLine
            lines = [self.fluxLayer, self.exposureLayer, self.linesLayer]
            visibility = [self.displayFlux, self.displayExposure, self.displayLines]
            # Add axes
            if self.displayExposure:
                self.ax3.get_yaxis().set_tick_params(labelright='on', right='on',
                                      direction='out', pad=5, colors='orange')
            else:
                self.ax3.get_yaxis().set_tick_params(labelright='off',right='off')
        # Prepare legend                
        self.labs = [l.get_label() for l in lns]
        leg = self.axes.legend(lns, self.labs, loc='upper center', bbox_to_anchor=(0.5, -0.1),
                               fancybox=True, shadow=True, ncol=5)
        # leg.set_draggable(True)  # works only in matplotlib 3.0.0      
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
                #txt.set_text('')
            legline.set_alpha(alpha)
            txt.set_alpha(alpha)            
        # Shade region considered for the images
        if self.shade == True:
            self.shadeRegion()
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
            print('l0 is not defined')        
        # Add spectral lines
        self.annotations = []
        font = FontProperties(family='DejaVu Sans', size=12)
        xlim0, xlim1 = self.axes.get_xlim()
        ylim0, ylim1 = self.axes.get_ylim()
        c = 299792458.0  # speed of light in m/s
        if self.xunit == 'THz':
            xlim1, xlim0 = c / xlim0 * 1.e-6, c / xlim1 * 1.e-6
        dy = ylim1 - ylim0
        for line in self.Lines.keys():
            nline = self.Lines[line][0]
            wline = self.Lines[line][1]*(1.+s.redshift)
            if (wline > xlim0 and wline < xlim1):
                wdiff = abs(s.wave - wline)
                imin = np.argmin(wdiff)
                y = s.flux[imin]
                y1 = y
                if (ylim1-(y+0.2*dy)) > ((y-0.2*dy)-ylim0):
                    y2 = y+0.2*dy
                else:
                    y2 = y-0.2*dy
                if self.xunit == 'um':
                    xline = wline
                elif self.xunit == 'THz':
                    xline = c/wline * 1.e-6
                y2 = 0.95  # Axes coordinates
                trans = transforms.blended_transform_factory(self.axes.transData, 
                                                             self.axes.transAxes)
                annotation = self.axes.annotate(nline, xy=(xline,y1),  xytext=(xline, y2),
                                                textcoords = trans,
                                                #draggable=True,
                                                picker=5,
                                                color='purple', alpha=0.4, ha='center',
                                                arrowprops=dict(edgecolor='purple',facecolor='y', 
                                                                arrowstyle='-',alpha=0.4,
                                                connectionstyle="angle,angleA=0,angleB=90,rad=10"),
                                                rotation = 90, fontstyle = 'italic',
                                                fontproperties=font, visible=self.displayLines)
                annotation.draggable()
                self.annotations.append(annotation)    

    def computeVelLimits(self):
        """Compute velocity limits."""
        x1, x2 = self.xlimits
        c = 299792.458  # speed of light in km/s
        s = self.spectrum
        vx1 = (x1 / (1 + s.redshift) / s.l0 - 1.) * c
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
              # a.set_position((xl, y2))
        except:
           pass
        self.fig.canvas.draw_idle()

    def updateSpectrum(self, f=None, uf=None, ef=None, exp=None, cont=None,
                       moments=None, noise=None, atran=None):
        try:
            try:
                self.arrow1.remove()
                self.arrow2.remove()
                self.gauss.remove()
            except:
                pass
            if atran is not None:
                self.atranLine[0].set_ydata(atran)
                self.axes.draw_artist(self.atranLine[0])
            if cont is not None:
                self.spectrum.continuum = cont
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
                ylim0,ylim1 = self.axes.get_ylim()
                maxf = np.nanmax(f)
                minf = np.nanmin(f)
                if uf is not None:
                    umaxf = np.nanmax(uf)
                    if umaxf > maxf: maxf = umaxf
                self.axes.set_ylim(minf, maxf*1.1)
                self.ylimits = (minf,maxf*1.1)
                self.updateYlim(f=f)
            if moments is not None:
                # Update position, size, and dispersion from moments
                x = moments[1]; y = np.nanmedian(cont)
                # FWHM of the distribution (assuming Gaussian)
                dx = np.sqrt(2.*np.log(2))* np.sqrt(moments[2])
                # Amplitude of a Gaussian ...
                c = 299792458. # m/s
                A = moments[0] * 1.e26 / np.sqrt(2.*np.pi*moments[2]) * x*x/c*1.e-6
                style = 'wedge'
                if self.xunit == 'THz':
                    dx = (c*dx)/(x*x-dx*dx) * 1.e-6
                    x = c/x * 1.e-6
                #verts = np.array([(x_,g_) for (x_,g_) in zip(xx,gauss)])
                self.arrow1 = FancyArrowPatch((x,y),(x,y+A),arrowstyle=style,
                                              mutation_scale=1.0,color='skyblue')
                self.axes.add_patch(self.arrow1)
                self.arrow2 = FancyArrowPatch((x-dx,y+0.5*A),(x+dx,y+0.5*A),arrowstyle=style,
                                              mutation_scale=1.0,color='skyblue')
                self.axes.add_patch(self.arrow2)
                # Gauss patch
                sig = dx / (np.sqrt(2.*np.log(2)))
                xx = np.arange(x-4*sig,x+4*sig,dx/10.)    
                gauss = y + A * np.exp(-(xx-x)**2/(2*sig*sig))
                verts = list(zip(xx,gauss))
                self.gauss = Polygon(verts, fill=False, closed=False,color='skyblue')
                self.axes.add_patch(self.gauss)
        except:
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
            # print('Line ', event.artist)
            legline = event.artist
            label = legline.get_label()
            origline = self.lined[legline]
            txt = self.labed[legline]
            vis = not origline.get_visible()
            origline.set_visible(vis)
            # Change the alpha on the line in the legend so we can see what lines
            # have been toggled
            if vis:
                legline.set_alpha(1.0)
                txt.set_alpha(1.0)
                # print('label is ',label)
                if label == 'E':
                    self.displayExposure = True
                    self.ax3.tick_params(labelright='on', right='on',
                                         direction='in', pad=-30, colors='orange')
                    txt.set_text('E')
                elif label == 'Atm':
                    self.displayAtran = True
                    self.ax2.get_yaxis().set_tick_params(labelright='on', right='on')            
                    self.ax2.get_yaxis().set_tick_params(which='both', direction='out',colors='red')
                    txt.set_text('Atm')
                elif label == 'F$_{u}$':
                    self.displayUFlux = True
                    txt.set_text('F$_{u}$')
                elif label == 'F':
                    self.displayFlux = True
                    txt.set_text('F')
                elif label == 'Lines':
                    self.displayLines = True
                    txt.set_text('Lines')
            else:
                legline.set_alpha(0.2)
                txt.set_alpha(0.2)
                txt.set_text('')
                if label == 'Exp':
                    self.displayExposure = False
                    self.ax3.get_yaxis().set_tick_params(labelright='off', right='off')
                elif label == 'Atm':
                    self.displayAtran = False
                    self.ax2.get_yaxis().set_tick_params(labelright='off', right='off')            
                elif label == 'Uflux':
                    self.displayUFlux = False
                elif label == 'Flux':
                    self.displayFlux = False
                elif label == 'Lines':
                    self.displayLines = False
            if self.shade == True:
                self.shadeRegion()
            if label == 'Lines':
                self.setLinesVisibility(self.displayLines)
            self.fig.canvas.draw_idle()
        elif isinstance(event.artist, Text):
            text = event.artist.get_text()
            if event.artist == self.zannotation:
                c = 299792.458 #km/s
                znew = self.getDouble(self.spectrum.redshift*c)
                if znew is not None:
                    znew /= c
                    if znew != self.spectrum.redshift:
                        self.spectrum.redshift = znew
                        for annotation in self.annotations:
                            annotation.remove()
                        self.zannotation.remove()
                        self.lannotation.remove()
                        self.drawSpectrum()
                        self.fig.canvas.draw_idle()
                        # Simulate a release to activate the update of redshift in main program
                        QTest.mouseRelease(self, Qt.LeftButton)
            if event.artist == self.lannotation:
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
            else:
                if text == 'Wavelength [$\mu$m]' or text == 'Frequency [THz]':
                    self.switchUnits()
                    self.switchSignal.emit('switched x unit')
                else:
                    self.dragged = event.artist
                    self.pick_pos = event.mouseevent.xdata
        else:
            pass
        return True
    
    def pickSwitch(self, artist, mouseevent):
        # print('artist ', artist)
        # print('switching units ... ', mouseevent)
        self.switchUnits()
        self.switchSignal.emit('switched x unit')
    
    def switchUnits(self):
        if self.xunit == 'um':
            self.xunit = 'THz'
            self.axes.format_coord = lambda x, y: "{:6.4f} THz {:10.4f} Jy".format(x, y)
        else:
            self.xunit = 'um'
            self.axes.format_coord = lambda x, y: "{:8.4f} um {:10.4f} Jy".format(x, y)
        self.axes.clear()
        try:
            self.ax2.clear()
            self.ax3.clear()
            self.ax4.clear()
        except:
            pass
        self.drawSpectrum()
        self.fig.canvas.draw_idle()
        if self.guess is not None:
            self.guess.switchUnits()  
     
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
                                                 "Ref. wavelength [um]", l, 0., 500., 4)
        if okPressed:
            return lnew
        else:
            return None

    def onrelease(self, event):
        if self.dragged is not None and self.pick_pos is not None:
            #print ('old ', self.dragged.get_position(), ' new ', event.xdata)
            x1 = event.xdata
            x0 = self.pick_pos
            #print('pick up position is ',x0)
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
            #print('Reference wavelength is ',self.spectrum.l0)
            self.spectrum.redshift = (1. + self.spectrum.redshift) * (1 + z) - 1.
            for annotation in self.annotations:
                annotation.remove()
            self.zannotation.remove()
            self.lannotation.remove()
            self.drawSpectrum()
            self.fig.canvas.draw_idle()
            self.dragged = None
        return True
