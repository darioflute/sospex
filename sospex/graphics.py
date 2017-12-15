import numpy as np
import os
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT 
from matplotlib.ticker import ScalarFormatter
from matplotlib.font_manager import FontProperties

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

from astropy.wcs.utils import proj_plane_pixel_scales as pixscales


from PyQt5.QtWidgets import (QApplication, QWidget, QMainWindow, QTabWidget, QHBoxLayout,
                             QGroupBox, QVBoxLayout, QSizePolicy, QStatusBar, QSplitter,
                             QToolBar, QAction, QFileDialog, QInputDialog, QTableView, QComboBox, QAbstractItemView,
                             QToolButton)
from PyQt5.QtGui import QIcon, QStandardItem, QStandardItemModel
from PyQt5.QtCore import Qt, QSize, pyqtSignal
from PyQt5.QtTest import QTest

class NavigationToolbar(NavigationToolbar2QT):
    def __init__(self,canvas,parent):
        # Select only a few buttons
        #self.toolitems = [t for t in NavigationToolbar2QT.toolitems if
        #                  t[0] in ('Home', 'Pan', 'Zoom', 'Save')]
        self.iconDir = os.path.join(os.path.dirname(os.path.abspath(__file__)),"icons")

        icon = QIcon(self.iconDir+'exit.png')
        #print('path of icons is ', icon)
        self.toolitems = [
            ('Home','Go back to original limits','home','home'),
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
    """ Canvas to plot an image """
    
    def __init__(self, *args, **kwargs):
        MplCanvas.__init__(self, *args, **kwargs)

            
    def compute_initial_figure(self, image=None, wcs=None, title=None):
        if wcs == None:
            ''' initial definition when images are not yet read '''
            pass
        else:
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
                self.fig.suptitle(title)

            # Show image
            self.showImage(image)
            
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
            

    def showImage(self, image):
        
            self.oimage = image
            self.image = self.axes.imshow(image, cmap='gist_heat_r',interpolation='none')
            self.fig.colorbar(self.image, cax=self.cbaxes)
            # Intensity limits
            vmed0=np.nanmedian(image)
            d0 = np.nanstd(image)
            self.image.set_clim([vmed0-d0,vmed0+4*d0])

            self.changed = False

    def updateScale(self,val):
        _cmin = self.s_cmin.val
        _cmax = self.s_cmax.val
        self.image.set_clim([_cmin, _cmax])
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
        # Test
        xlim0,xlim1 = self.axes.get_xlim()
        ylim0,ylim1 = self.axes.get_ylim()
        bbox_args = dict(boxstyle="round", fc="0.8")
        xl = (xlim0+xlim1)*0.5
        yl = (ylim0+ylim1)*0.5
        a1=self.axes.annotate('Drag me 1', xy=(xl,yl), xycoords='data',
                              #xytext=(.5, .7), textcoords='data',
                              ha="center", va="center",
                              bbox=bbox_args,
                              #arrowprops=arrow_args
                          )
        a1.draggable()

        
    def compute_initial_figure(self, image=None,xmin=None,xmax=None):
        if image is None:
            ''' initial definition when images are not yet read '''
            pass
        else:
            # Print the histogram of finite values
            ima = image.ravel()
            mask = np.isfinite(ima)
            ima = ima[mask]
            print('image has size', len(ima))
            self.nh = len(ima)
            ima = np.sort(ima)
            s = np.size(ima)
            smax = min(int(s*0.9995),s-1)
            nbins=256
            n, self.bins, patches = self.axes.hist(ima, bins=nbins, range=(np.nanmin(ima), ima[smax]), fc='k', ec='k')

            self.median = np.median(ima)
            self.sdev   = np.std(ima-self.median)
            self.min    = np.min(ima)
            self.max    = np.max(ima)
            self.epsilon = self.sdev/3.
            
            # Define the interval containing 99% of the values
            if xmin == None:
                xmin = ima[int(s*0.01)]
            if xmax == None:
                xmax = ima[int(s*0.99)-1]
            self.onSelect(xmin,xmax)

            # Initialize contour level vertical lines
            self.lev = []
            self.levels = []
            self._ind = None

            # Draw grid (median, median+n*sigma)
            x = self.median
            while x < self.max:
                self.axes.axvline(x=x,color='black',alpha=0.5)
                x += self.sdev


    def drawLevels(self):
        """ Draw levels as defined in levels"""

        if len(self.levels) > 0:
            for l in self.levels:
                lev = self.axes.axvline(x=l, animated=True, color='cyan')
                self.lev.append(lev)
            print('there are ',len(self.lev),' contours')
            #self.span.set_visible(False)
            self.span.active = False
            self.connect()
            self.showLevels = True

    def removeLevels(self):
        """ Delete levels and disconnect """

        if len(self.levels) > 0:
            for lev in self.lev:
                lev.remove()
            self.lev = []
            self.levels = []
            self.disconnect()
            #self.span.set_visible(True)
            self.span.active = True
            self.showLevels = False
            print("All levels have been removed")
            
    def onSelect(self,xmin, xmax):
        indmin, indmax = np.searchsorted(self.bins, (xmin, xmax))
        indmax = min(len(self.bins) - 1, indmax)
        self.limits = [self.bins[indmin],self.bins[indmax]]
        try:
            self.shade.remove()
        except:
            pass
        self.limSignal.emit('limits changed')
        self.shade = self.axes.axvspan(self.limits[0],self.limits[1],facecolor='Lavender',alpha=0.5,linewidth=0)
        # Redefine limits
        x1,x2 = self.axes.get_xlim()
        x2 = self.limits[1] + 5 * self.sdev
        self.axes.set_xlim((x1,x2))
        self.fig.canvas.draw_idle()

    def connect(self):
        """ When contours are present """
        self.cid_draw = self.fig.canvas.mpl_connect('draw_event', self.draw_callback)
        self.cid_press = self.fig.canvas.mpl_connect('button_press_event', self.button_press_callback)
        self.cid_release = self.fig.canvas.mpl_connect('button_release_event', self.button_release_callback)
        self.cid_motion = self.fig.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        self.cid_key = self.fig.canvas.mpl_connect('key_press_event', self.key_press_callback)
        self.fig.canvas.draw_idle()
        self.fig.canvas.setFocus()

    def disconnect(self):
        """ When no contours are present """
        self.fig.canvas.mpl_disconnect(self.cid_draw)
        self.fig.canvas.mpl_disconnect(self.cid_press)
        self.fig.canvas.mpl_disconnect(self.cid_release)
        self.fig.canvas.mpl_disconnect(self.cid_motion)
        self.fig.canvas.mpl_disconnect(self.cid_key)
        self.fig.canvas.draw_idle()
        # add something about levels ?


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
        if self.nh > 100000:
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
                #self.span.set_visible(False)
                self.span.active = False
                self.cid_draw = self.fig.canvas.mpl_connect('draw_event', self.draw_callback)
                self.cid_press = self.fig.canvas.mpl_connect('button_press_event', self.button_press_callback)
                self.cid_release = self.fig.canvas.mpl_connect('button_release_event', self.button_release_callback)
                self.cid_motion = self.fig.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
                #self.connect()
            else:
                self._ind = None
                # If levels are not visible, show the span selector
                #self.span.set_visible(True)
                self.span.active = True
                self.fig.canvas.mpl_disconnect(self.cid_draw)
                self.fig.canvas.mpl_disconnect(self.cid_press)
                self.fig.canvas.mpl_disconnect(self.cid_release)
                self.fig.canvas.mpl_disconnect(self.cid_motion)
                #self.disconnect()
        elif event.key == 'd':
            ind = self.get_ind_under_point(event)
            if ind is not None:
                # Delete contour level
                self.lev[ind].remove()
                # Remove level from lists
                del self.lev[ind]
                del self.levels[ind]
                # Emit signal to communicate it to images
                self.levSignal.emit(-ind)
                # If there are no more levels, show the span selector
                if len(self.levels) == 0:
                    #self.span.set_visible(True)
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
            self.levSignal.emit(1000+idx[n])

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
        x, y = event.xdata, event.ydata
        #print(x,y)

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
        if self.nh <= 100000:
            self.levSignal.emit(self._ind)



class SpectrumCanvas(MplCanvas):
    """ Canvas to plot spectra """
    def __init__(self, *args, **kwargs):
        MplCanvas.__init__(self, *args, **kwargs)
        from sospex.lines import define_lines

        self.Lines = define_lines()
        self.fig.set_edgecolor('none')
        self.axes = self.fig.add_axes([0.12,0.15,.8,.78])

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
        self.dragged = None
        self.region = None
        
    def compute_initial_spectrum(self, spectrum=None,xmin=None,xmax=None):
        if spectrum is None:
            ''' initial definition when spectrum not yet available '''
        else:
            # Spectrum
            self.spectrum = spectrum
            self.instrument = spectrum.instrument
            self.drawSpectrum()
        
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
            self.axes.set_xlabel('Wavelength [$\\mu$m]', picker=True)
            self.x = s.wave
            if s.instrument == 'FIFI-LS':
                self.xr = self.x * (1+s.baryshift)
        elif self.xunit == 'THz':
            c = 299792458.0  # speed of light in m/s
            self.axes.set_xlabel('Frequency [THz]', picker=True)
            self.x = c/s.wave * 1.e-6
            if s.instrument == 'FIFI-LS':
                self.xr = self.x / (1+s.baryshift)

             
        self.fluxLine = self.axes.step(self.x,s.flux,color='blue',label='Flux')
        self.fluxLayer, = self.fluxLine

        # Define limits or adjust to previous limits
        if self.xlimits is not None:
            xlim0,xlim1 = self.xlimits
            if self.xunit == 'THz':
                c = 299792458.0  # speed of light in m/s
                xlim1,xlim0 = c/xlim0*1.e-6,c/xlim1*1.e-6
            self.axes.set_xlim(xlim0,xlim1)
            self.axes.set_ylim(self.ylimits)
        else:
            xlim0 = np.min(s.wave)
            xlim1 = np.max(s.wave)
            self.xlimits=(xlim0,xlim1)
            self.ylimits=self.axes.get_ylim()
            self.axes.set_xlim(xlim0,xlim1)
            


        # Fake line to have the lines in the legend
        self.linesLine = self.axes.plot([0,0.1],[0,0],color='purple',alpha=0.4,label='Lines')
        self.linesLayer, = self.linesLine




        # Add redshift value on the plot
        c = 299792.458  #km/s
        self.zannotation = self.axes.annotate(" cz = {:.1f} km/s".format(c*s.redshift), xy=(-0.05,-0.1), picker=5, xycoords='axes fraction')
        #self.zannotation.draggable()
                
        
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
            self.atranLine = self.ax2.step(self.xr, s.atran,color='red',label='Atm')
            self.exposureLine = self.ax3.step(self.x, s.exposure, color='orange',label='Exp')
            ymax = np.nanmax(s.flux); ymin = np.nanmin(s.flux)
            yumax = np.nanmax(s.uflux); yumin = np.nanmin(s.uflux)
            if yumax > ymax: ymax=yumax
            if yumin < ymin: ymin=yumin
            self.ax3.set_ylim([0.5,np.nanmax(s.exposure)*1.54])
            self.ufluxLine = self.ax4.step(self.xr,s.uflux,color='green',label='Uflux')
            self.ax4.set_ylim(self.axes.get_ylim())
            #self.ax1.set_title(spectrum.objname+" ["+spectrum.filegpid+"] @ "+spectrum.obsdate)
            self.ufluxLayer, = self.ufluxLine
            self.atranLayer, = self.atranLine
            self.exposureLayer, = self.exposureLine
            
            
            lns = self.fluxLine \
                  +self.ufluxLine \
                  +self.atranLine \
                  +self.exposureLine \
                  +self.linesLine
            lines = [self.fluxLayer, self.ufluxLayer, self.atranLayer, self.exposureLayer, self.linesLayer]
            visibility = [self.displayFlux, self.displayUFlux, self.displayAtran, self.displayExposure, self.displayLines]
            
            # Add axes
            if self.displayExposure:
                self.ax3.get_yaxis().set_tick_params(labelright='on',right='on',direction='out',pad=5,colors='orange')
            else:
                self.ax3.get_yaxis().set_tick_params(labelright='off',right='off')
            if self.displayAtran:
                self.ax2.get_yaxis().set_tick_params(labelright='on',right='on', direction='in', pad = -25, colors='red')
            else:
                self.ax2.get_yaxis().set_tick_params(labelright='off',right='off')            
        elif s.instrument == 'GREAT':
            self.displayUFlux = False
            self.displayAtran = False
            self.displayExposure = False
            lns = self.fluxLine + self.linesLine
            lines = [self.fluxLayer,self.linesLayer]
            visibility = [self.displayFlux,self.displayLines]


                
        # Prepare legend                
        self.labs = [l.get_label() for l in lns]
        leg = self.axes.legend(lns, self.labs, loc='upper center', bbox_to_anchor=(0.5, -0.1),
                               fancybox=True, shadow=True, ncol=5)
        leg.draggable()
        
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
        xlim0,xlim1 = self.axes.get_xlim()
        ylim0,ylim1 = self.axes.get_ylim()
        if self.xunit == 'THz':
            c = 299792458.0  # speed of light in m/s
            xlim1,xlim0 = c/xlim0*1.e-6,c/xlim1*1.e-6
        dy = ylim1-ylim0

        for line in self.Lines.keys():
            nline = self.Lines[line][0]
            wline = self.Lines[line][1]*(1.+s.redshift)
            if (wline > xlim0 and wline < xlim1):
                wdiff = abs(s.wave - wline)
                y = s.flux[(wdiff == wdiff.min())]
                y1 = y
                if (ylim1-(y+0.2*dy)) > ((y-0.2*dy)-ylim0):
                    y2 = y+0.2*dy
                else:
                    y2 = y-0.2*dy
                c = 299792.458  # speed of light in km/s
                if self.xunit == 'um':
                    xline = wline
                    vline = (xline/(1+s.redshift)/s.l0-1.)*c
                elif self.xunit == 'THz':
                    xline = c/wline * 1.e-9
                    vline = (c/xline/(1+s.redshift)/s.l0*1.e-3-1.)*c
                annotation = self.axes.annotate(nline, xy=(xline,y1),  xytext=(xline, y2), color='purple', alpha=0.4,
                                                arrowprops=dict(edgecolor='purple',facecolor='y', arrowstyle='-',alpha=0.4,
                                                connectionstyle="angle,angleA=0,angleB=90,rad=10"),
                                                rotation = 90, fontstyle = 'italic', fontproperties=font,
                                                visible=self.displayLines,)
                annotation.draggable()
                self.annotations.append(annotation)     

    def computeVelLimits(self):
        x1,x2 = self.axes.get_xlim()
        c = 299792.458  # speed of light in km/s
        if self.xunit == 'um':
            vx1 = (x1/(1+s.redshift)/s.l0-1.)*c
            vx2 = (x2/(1+s.redshift)/s.l0-1.)*c
        elif self.xunit == 'THz':
            vx1 = (c/x1/(1+s.redshift)/s.l0*1.e-3-1.)*c
            vx2 = (c/x2/(1+s.redshift)/s.l0*1.e-3-1.)*c

        return (vx1,vx2)
            

    def updateXlim(self):
        """ update xlimits """
        xlim0,xlim1 = sc.xlimits
        vlim = sc.computeVelLimits(xlim0,xlim1)
        if self.xunit == 'THz':
            c = 299792458.0  # speed of light in m/s
            xlim1,xlim0 = c/xlim0*1.e-6,c/xlim1*1.e-6
        sc.axes.set_xlim(xlim0,xlim1)
        sc.vaxes.set_xlim(vlim)
        sc.fig.canvas.draw_idle()

        
    def updateSpectrum(self,f,uf=None,exp=None):

        try:
            self.fluxLine[0].set_ydata(f)
            self.axes.draw_artist(self.fluxLine[0])
            if uf is not None:
                self.ufluxLine[0].set_ydata(uf)
                self.ax4.draw_artist(self.ufluxLine[0])
            if exp is not None:
                self.exposureLine[0].set_ydata(exp)
                self.ax3.draw_artist(self.exposureLine[0])
            ylim0,ylim1 = self.axes.get_ylim()
            if np.nanmax(f) > ylim1:
                self.axes.set_ylim(ylim0, np.nanmax(f)*1.1)
                self.axes.draw_artist(self.axes.yaxis)
            self.fig.canvas.draw_idle()
            #self.fig.canvas.update()
            #self.fig.canvas.flush_events()
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

        if color == 'Lavender':
            self.region = self.axes.axvspan(xmin,xmax,facecolor=color,alpha=0.5,linewidth=0)
        else:
            self.tmpRegion = self.axes.axvspan(xmin,xmax,facecolor=color,alpha=0.5,linewidth=0)
        
                
    def shadeSpectrum(self):
        """ shade part of the spectrum used for different operations """

        # Clean previous region
        try:
            print("previous shade removed")
            self.region.remove()
        except:
            pass

        # Shade new region
        self.shadeRegion()
        self.shade = True
        
                
    def onpick(self, event):
        """ React to onpick events """

        if isinstance(event.artist, Line2D):
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
                if label == 'Exp':
                    self.displayExposure = True
                    self.ax3.tick_params(labelright='on',right='on',direction='in',pad=-30,colors='orange')
                    txt.set_text('Exp')
                elif label == 'Atm':
                    self.displayAtran = True
                    self.ax2.get_yaxis().set_tick_params(labelright='on',right='on')            
                    self.ax2.get_yaxis().set_tick_params(which='both', direction='out',colors='red')
                    txt.set_text('Atm')
                elif label == 'Uflux':
                    self.displayUFlux = True
                    txt.set_text('Uflux')
                elif label == 'Flux':
                    self.displayFlux = True
                    txt.set_text('Flux')
                elif label == 'Lines':
                    self.displayLines = True
                    txt.set_text('Lines')
            else:
                legline.set_alpha(0.2)
                txt.set_alpha(0.2)
                txt.set_text('')
                if label == 'Exp':
                    self.displayExposure = False
                    self.ax3.get_yaxis().set_tick_params(labelright='off',right='off')
                elif label == 'Atm':
                    self.displayAtran = False
                    self.ax2.get_yaxis().set_tick_params(labelright='off',right='off')            
                elif label == 'Uflux':
                    self.displayUFlux = False
                elif label == 'Flux':
                    self.displayFlux = False
                elif label == 'Lines':
                    self.displayLines = False
            if self.shade == True:
                self.shadeRegion()
            for annotation in self.annotations:
                annotation.set_visible(self.displayLines)
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
                        self.drawSpectrum()
                        self.fig.canvas.draw_idle()
                        # Simulate a release to activate the update or redshift in main program
                        QTest.mouseRelease(self, Qt.LeftButton)
            else:
                if text == 'Wavelength [$\mu$m]' or text == 'Frequency [THz]':
                    if self.xunit == 'um':
                        self.xunit = 'THz'
                    else:
                        self.xunit = 'um'

                    self.axes.clear()
                    try:
                        self.ax2.clear()
                        self.ax3.clear()
                        self.ax4.clear()
                    except:
                        pass
                    self.drawSpectrum()
                    self.fig.canvas.draw_idle()
                else:
                    print(text, event.mouseevent.xdata)
                    self.dragged = event.artist
                    self.pick_pos = event.mouseevent.xdata
                
        else:
            pass
        return True

    def getDouble(self,z):
        znew, okPressed = QInputDialog.getDouble(self, "Redshift","cz:", z, -1000., 30000., 2)
        if okPressed:
            return znew
        else:
            return None
        

    def onrelease(self, event):
        if self.dragged is not None :
            print ('old ', self.dragged.get_position(), ' new ', event.xdata)
            x1 = event.xdata
            x0 = self.pick_pos
            print('pick up position is ',x0)
            w0 = np.array([self.Lines[line][1]  for  line in self.Lines.keys()])
            wz = w0*(1.+self.spectrum.redshift)
            if self.xunit == 'um':
                z = (x1-x0)/x0
                l0 = x0
            elif self.xunit == 'THz':
                z = (x0-x1)/x1
                c = 299792458.0  # speed of light in m/s
                l0 = c/x0 * 1.e-6
            wdiff = abs(l0 - wz)
            self.spectrum.l0 = w0[(wdiff == wdiff.min())]
            print('Reference wavelength is ',self.spectrum.l0)
            self.spectrum.redshift = (1.+self.spectrum.redshift)*(1+z)-1.
            for annotation in self.annotations:
                annotation.remove()
            self.zannotation.remove()
            self.drawSpectrum()
            self.fig.canvas.draw_idle()
            self.dragged = None
        return True
