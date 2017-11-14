import numpy as np

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT 
from matplotlib.ticker import ScalarFormatter

# Matplotlib parameters
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family']='STIXGeneral'
rcParams['font.size']=13
rcParams['mathtext.fontset']='stix'
rcParams['legend.numpoints']=1

from matplotlib.widgets import SpanSelector
from astropy.wcs.utils import proj_plane_pixel_scales as pixscales


from PyQt5.QtWidgets import (QApplication, QWidget, QMainWindow, QTabWidget, QHBoxLayout,
                             QGroupBox, QVBoxLayout, QSizePolicy, QStatusBar, QSplitter,
                             QToolBar, QAction, QFileDialog, QTableView, QComboBox, QAbstractItemView,
                             QToolButton)
from PyQt5.QtGui import QIcon, QStandardItem, QStandardItemModel
from PyQt5.QtCore import Qt, QSize, pyqtSignal


class NavigationToolbar(NavigationToolbar2QT):
    # Select only a few buttons
    toolitems = [t for t in NavigationToolbar2QT.toolitems if
                 t[0] in ('Home', 'Pan', 'Zoom', 'Save')]


class imageNavigationToolbar(NavigationToolbar2QT):
    # Select only a few buttons
    toolitems = [t for t in NavigationToolbar2QT.toolitems if
                 t[0] in ('Home', 'Pan', 'Zoom', 'Save')]
    #toolitems.append(('Levels','Set levels of the image',self.path0+'/icons/levels.png,))


    #levelsAction = self.createAction(self.path0+'/icons/levels.png','Adjust image levels','Ctrl+L',self.changeVisibility)


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
            self.axes = self.fig.add_axes([0.1,0.1,.8,.8], projection = self.wcs)
            self.oimage = image
            self.image = self.axes.imshow(image, cmap='gist_heat_r',interpolation='none')
            self.axes.coords[0].set_major_formatter('hh:mm:ss')
            self.axes.grid(color='black', ls='dashed')
            self.axes.set_xlabel('R.A.')
            self.axes.set_ylabel('Dec')

            # Plot with North up (corners are clockwise from left-bottom)
            # corners = self.wcsn.calc_footprint()
            # self.flip = False
            # if corners[0,1] < 0:
            #     if corners[0,1] > corners[1,1]:
            #         ylim = self.axes.get_ylim()
            #         self.axes.set_ylim([ylim[1],ylim[0]])
            #         self.flip = True

            # Add title
            if title != None:
                self.fig.suptitle(title)

            # Colorbar
            cbaxes = self.fig.add_axes([0.9,0.1,0.02,0.85])
            self.fig.colorbar(self.image, cax=cbaxes)

            # Sliders to adjust intensity
            #self.ax_cmin = self.figure.add_axes([0.1, 0.01, 0.8, 0.01])
            #self.ax_cmax = self.figure.add_axes([0.1, 0.04, 0.8, 0.01])
            #self.ax_cmin.clear()
            #self.ax_cmax.clear()
            vmed0=np.nanmedian(image)
            d0 = np.nanstd(image)
            #self.s_cmin = Slider(self.ax_cmin, 'low', vmed0-2*d0, vmed0+5*d0, valinit=vmed0-d0, facecolor='goldenrod')
            #self.s_cmax = Slider(self.ax_cmax, 'high', vmed0-2*d0, vmed0+5*d0, valinit=vmed0+4*d0, facecolor='goldenrod')
            self.image.set_clim([vmed0-d0,vmed0+4*d0])
            #self.s_cmin.valtext.set_visible(False)
            #self.s_cmax.valtext.set_visible(False)
            #self.slider1=self.s_cmin.on_changed(self.updateScale)
            #self.slider2=self.s_cmax.on_changed(self.updateScale)


            # Add ellipse centered on source
            pixscale = pixscales(self.wcs)*3600. # Scale in arcsec
            #if self.flip:
            #    theta2= 0
            #    theta1 = 100
            #else:
            #    theta1=0
            #    theta2=100


            # Ellipse
            # self.arcell = self.ArcEll((xc,yc), 5/pixscale[0], 5/pixscale[1], 'Lime', 30)
            # for a in self.arcell:
            #     self.axes.add_patch(a)
            #     self.drrEllipse = DragResizeRotateEllipse(self.arcell)

            self.changed = False



    def updateScale(self,val):
        _cmin = self.s_cmin.val
        _cmax = self.s_cmax.val
        self.image.set_clim([_cmin, _cmax])
        self.fig.canvas.draw_idle()

    def ArcEll(self,pos,w,h,color,theta):

        arcell = []
        arcell.append(Ellipse(pos,w,h,edgecolor=color,facecolor='none'))
        arcell.append(Arc(pos,w,h, edgecolor=color, facecolor='none',theta1=0 -theta,theta2=0 +theta,lw=4))
        arcell.append(Arc(pos,w,h, edgecolor=color, facecolor='none',theta1=90 -theta,theta2=90 +theta,lw=4))
        arcell.append(Arc(pos,w,h, edgecolor=color, facecolor='none',theta1=180 -theta,theta2=180 +theta,lw=4))
        arcell.append(Arc(pos,w,h, edgecolor=color, facecolor='none',theta1=270 -theta,theta2=270 +theta,lw=4))
        return arcell

class ImageHistoCanvas(MplCanvas):
    """ Canvas to plot the histogram of image intensity """
    def __init__(self, *args, **kwargs):
        MplCanvas.__init__(self, *args, **kwargs)

    mySignal = pyqtSignal(str)
        
    def compute_initial_figure(self, image=None,xmin=None,xmax=None):
        if image is None:
            ''' initial definition when images are not yet read '''
            pass
        else:
            self.axes = self.fig.add_axes([0.0,0.4,1.,1.])
            self.axes.yaxis.set_major_formatter(plt.NullFormatter())
            self.axes.spines['top'].set_visible(False)
            self.axes.spines['right'].set_visible(False)
            self.axes.spines['left'].set_visible(False)
            # Print the histogram of finite values
            ima = image.ravel()
            mask = np.isfinite(ima)
            ima = ima[mask]
            ima = np.sort(ima)
            s = np.size(ima)
            smax = min(int(s*0.9995),s-1)
            nbins=256
            n, self.bins, patches = self.axes.hist(ima, bins=nbins, range=(np.nanmin(ima), ima[smax]), fc='k', ec='k')
            # Define the interval containing 99% of the values

            self.x = np.arange(s)
            if xmin == None:
                xmin = ima[int(s*0.01)]
            if xmax == None:
                xmax = ima[int(s*0.99)-1]
            self.onSelect(xmin,xmax)
            # Start a span selector
            self.span = SpanSelector(self.axes, self.onSelect, 'horizontal', useblit=True,
                                     rectprops=dict(alpha=0.5, facecolor='LightSalmon'),button=1)
            # Define the way to draw/read a shaded interval (check from sospex)
            # Communicate values to the image (from the tab we know the canvas, but I have to communicate changes to the main window)


    def onSelect(self,xmin, xmax):
        indmin, indmax = np.searchsorted(self.bins, (xmin, xmax))
        indmax = min(len(self.bins) - 1, indmax)
        self.limits = [self.bins[indmin],self.bins[indmax]]
        try:
            self.shade.remove()
        except:
            pass
        self.mySignal.emit('limits changed')
        self.shade = self.axes.axvspan(self.limits[0],self.limits[1],facecolor='Lavender',alpha=0.5,linewidth=0)
        self.fig.canvas.draw_idle()


class SpectrumCanvas(MplCanvas):
    """ Canvas to plot spectra """
    def __init__(self, *args, **kwargs):
        MplCanvas.__init__(self, *args, **kwargs)
        #self.fig.set_tight_layout(True)
        self.fig.set_edgecolor('none')
        
        self.axes = self.fig.add_axes([0.1,0.1,.85,.85])
        self.axes.grid(True, which='both')
        self.axes.set_xlabel('Wavelength [$\\mu$m]')
        self.axes.set_ylabel('Flux [Jy]')
        
        # Checks
        self.displayFlux = True
        self.displayUFlux = True
        self.displayAtran = True
        self.displayExposure = True
        self.displayLines = False

        self.axes.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        self.axes.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        #self.axes.xaxis.set_major_formatter(plt.NullFormatter())
        #self.axes.yaxis.set_major_formatter(plt.NullFormatter())
        self.axes.spines['top'].set_visible(False)
        self.axes.spines['right'].set_visible(False)
        #self.axes.spines['left'].set_visible(False)

        # Use legend to hide/show lines
        self.fig.canvas.mpl_connect('pick_event', self.onpick)
        #mySignal = pyqtSignal(str)
        
    def compute_initial_spectrum(self, spectrum=None,xmin=None,xmax=None):
        if spectrum is None:
            ''' initial definition when spectrum not yet available '''
        else:
            s = spectrum
            self.instrument = s.instrument
            self.fluxLine = self.axes.step(s.wave,s.flux,color='blue',label='Flux')
            self.fluxLayer, = self.fluxLine
            if s.instrument == 'FIFI-LS':
                try:
                    self.ax2.remove()
                    self.ax3.remove()
                    self.ax4.remove()
                except:
                    pass
                self.ax2 = self.axes.twinx()
                self.ax3 = self.axes.twinx()
                self.ax4 = self.axes.twinx()
                self.ax2.set_ylim([0,1.1])
                self.ax4.tick_params(labelright='off',right='off')
                self.atranLine = self.ax2.step(s.wave*(1.+s.baryshift), s.atran,color='red',label='Atm Trans')
                self.exposureLine = self.ax3.step(s.wave, s.exposure, color='orange',label='Exposure')
                ymax = np.nanmax(s.flux); ymin = np.nanmin(s.flux)
                yumax = np.nanmax(s.uflux); yumin = np.nanmin(s.uflux)
                if yumax > ymax: ymax=yumax
                if yumin < ymin: ymin=yumin
                self.ax3.set_ylim([0.5,np.nanmax(s.exposure)*1.04])
                self.ufluxLine = self.ax4.step(s.wave*(1.+s.baryshift),s.uflux,color='green',label='Unc Flux')
                self.ax4.set_ylim(self.axes.get_ylim())
                self.axes.set_xlim([np.min(s.wave),np.max(s.wave)])
                #self.ax1.set_title(spectrum.objname+" ["+spectrum.filegpid+"] @ "+spectrum.obsdate)
                self.ufluxLayer, = self.ufluxLine
                self.atranLayer, = self.atranLine
                self.exposureLayer, = self.exposureLine

                
                lns = self.fluxLine \
                      +self.ufluxLine \
                      +self.atranLine \
                      +self.exposureLine
                lines = [self.fluxLayer, self.ufluxLayer, self.atranLayer, self.exposureLayer]
                visibility = [self.displayFlux, self.displayUFlux, self.displayAtran, self.displayExposure]

                
            elif s.instrument == 'GREAT':
                self.displayUFlux = False
                self.displayAtran = False
                self.displayExposure = False
                lns = self.fluxLine
                lines = [self.fluxLayer]
                visibility = [self.displayFlux]

            # Add axes
            if self.displayExposure:
                self.ax3.tick_params(labelright='on',right='on',direction='in',pad=-30,colors='orange')
            if self.displayAtran:
                self.ax2.get_yaxis().set_tick_params(labelright='on',right='on')            
                self.ax2.get_yaxis().set_tick_params(which='both', direction='out',colors='red')

            # Prepare legend                
            self.labs = [l.get_label() for l in lns]
            leg = self.axes.legend(lns, self.labs, loc='upper left',frameon=False,framealpha=0.0)
            self.lined = dict()
            self.labed = dict()
            for legline, origline, txt in zip(leg.get_lines(), lines, leg.texts):
                legline.set_picker(5) # 5pts tolerance
                self.lined[legline] = origline
                self.labed[legline] = txt

            # Hide lines
            for line, legline, vis in zip(lines, leg.get_lines(), visibility):
                line.set_visible(vis)
                if vis:
                    alpha=1.0
                else:
                    alpha=0.2
                legline.set_alpha(alpha)
                txt = self.labed[legline]
                txt.set_alpha(alpha)


                
    def onpick(self, event):
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
            if label == 'Exposure':
                self.displayExposure = True
                self.ax3.tick_params(labelright='on',right='on',direction='in',pad=-30,colors='orange')
            elif label == 'Atm Trans':
                self.displayAtran = True
                self.ax2.get_yaxis().set_tick_params(labelright='on',right='on')            
                self.ax2.get_yaxis().set_tick_params(which='both', direction='out',colors='red')
            elif label == 'Unc Flux':
                self.displayUFlux = True
            elif label == 'Flux':
                self.displayFlux = True
        else:
            legline.set_alpha(0.2)
            txt.set_alpha(0.2)
            if label == 'Exposure':
                self.displayExposure = False
                self.ax3.get_yaxis().set_tick_params(labelright='off',right='off')
            elif label == 'Atm Trans':
                self.displayAtran = False
                self.ax2.get_yaxis().set_tick_params(labelright='off',right='off')            
            elif label == 'Unc Flux':
                self.displayUFlux = False
            elif label == 'Flux':
                self.displayFlux = False
        self.fig.canvas.draw_idle()

            
    # def update_spectrum(self):
        
        #self.fluxLayer.set_visible(self.displayFlux)
        #self.ufluxLayer.set_visible(self.displayUFlux)
        #self.atranLayer.set_visible(self.displayAtran)
        #self.exposureLayer.set_visible(self.displayExposure)
        #if self.instrument == 'FIFI-LS':
        #    if self.displayExposure:
        #        self.ax3.tick_params(labelright='on',right='on',direction='in',pad=-30,colors='orange')
        #    else:
        #        self.ax3.get_yaxis().set_tick_params(labelright='off',right='off')
        #    if self.displayAtran:
        #        self.ax2.get_yaxis().set_tick_params(labelright='on',right='on')            
        #        self.ax2.get_yaxis().set_tick_params(which='both', direction='out',colors='red')
        #    else:
        #        self.ax2.get_yaxis().set_tick_params(labelright='off',right='off')            


