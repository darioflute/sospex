import numpy as np

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
                             QToolBar, QAction, QFileDialog, QTableView, QComboBox, QAbstractItemView,
                             QToolButton)
from PyQt5.QtGui import QIcon, QStandardItem, QStandardItemModel
from PyQt5.QtCore import Qt, QSize, pyqtSignal


class NavigationToolbar(NavigationToolbar2QT):
    # Select only a few buttons
    toolitems = [t for t in NavigationToolbar2QT.toolitems if
                 t[0] in ('Home', 'Pan', 'Zoom', 'Save')]



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
            

            # Plot with North up (corners are clockwise from left-bottom)
            # corners = self.wcsn.calc_footprint()
            # self.flip = False
            # if corners[0,1] < 0:
            #     if corners[0,1] > corners[1,1]:
            #         ylim = self.axes.get_ylim()
            #         self.axes.set_ylim([ylim[1],ylim[0]])
            #         self.flip = True




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
        self.axes = self.fig.add_axes([0.0,0.4,1.,1.])
        self.axes.yaxis.set_major_formatter(plt.NullFormatter())
        self.axes.spines['top'].set_visible(False)
        self.axes.spines['right'].set_visible(False)
        self.axes.spines['left'].set_visible(False)
        # Start a span selector
        self.span = SpanSelector(self.axes, self.onSelect, 'horizontal', useblit=True,
                                 rectprops=dict(alpha=0.5, facecolor='LightSalmon'),button=1)

    mySignal = pyqtSignal(str)
        
    def compute_initial_figure(self, image=None,xmin=None,xmax=None):
        if image is None:
            ''' initial definition when images are not yet read '''
            pass
        else:
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
            if xmin == None:
                xmin = ima[int(s*0.01)]
            if xmax == None:
                xmax = ima[int(s*0.99)-1]
            self.onSelect(xmin,xmax)


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
        from lines import define_lines

        self.Lines = define_lines()
        self.fig.set_edgecolor('none')
        self.axes = self.fig.add_axes([0.12,0.12,.82,.8])
        
        # Checks
        self.displayFlux = True
        self.displayUFlux = True
        self.displayAtran = True
        self.displayExposure = True
        self.displayLines = False
        self.shade = False
        self.regionlimits = None
        self.xunit = 'um'  # Alternatives are THz or km/s
        
        self.axes.spines['top'].set_visible(False)
        self.axes.spines['right'].set_visible(False)

        # Use legend to hide/show lines
        self.fig.canvas.mpl_connect('pick_event', self.onpick)
        
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
        self.axes.grid(True, which='both')
        self.axes.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        self.axes.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        
        # Write labels
        self.axes.set_ylabel('Flux [Jy]')
        s = self.spectrum
        
        if self.xunit == 'um':
            self.axes.set_xlabel('Wavelength [$\\mu$m]', picker=True)
            x = s.wave
            if s.instrument == 'FIFI-LS':
                xr = x * (1+s.baryshift)
        elif self.xunit == 'THz':
            c = 299792458.0  # speed of light in m/s
            self.axes.set_xlabel('Frequency [THz]', picker=True)
            x = c/s.wave * 1.e-6
            if s.instrument == 'FIFI-LS':
                xr = x / (1+s.baryshift)

        x0 = np.min(x); x1 = np.max(x)
        self.axes.set_xlim([x0,x1])
                
        self.fluxLine = self.axes.step(x,s.flux,color='blue',label='Flux')
        self.fluxLayer, = self.fluxLine
                
        if s.instrument == 'FIFI-LS':
            try:
                self.fig.delaxes(self.ax2)
                self.fig.delaxes(self.ax3)
                self.fig.delaxes(self.ax4)
                #self.ax2.remove()
                #self.ax3.remove()
                #self.ax4.remove()
            except:
                pass
            self.ax2 = self.axes.twinx()
            self.ax3 = self.axes.twinx()
            self.ax4 = self.axes.twinx()
            self.ax2.set_ylim([0.01,1.1])
            self.ax4.tick_params(labelright='off',right='off')
            self.atranLine = self.ax2.step(xr, s.atran,color='red',label='Atm Trans')
            self.exposureLine = self.ax3.step(x, s.exposure, color='orange',label='Exposure')
            ymax = np.nanmax(s.flux); ymin = np.nanmin(s.flux)
            yumax = np.nanmax(s.uflux); yumin = np.nanmin(s.uflux)
            if yumax > ymax: ymax=yumax
            if yumin < ymin: ymin=yumin
            self.ax3.set_ylim([0.5,np.nanmax(s.exposure)*1.54])
            self.ufluxLine = self.ax4.step(xr,s.uflux,color='green',label='Unc Flux')
            self.ax4.set_ylim(self.axes.get_ylim())
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


        # Check if vel is defined and draw velocity axis
        #if s.l0 is not None:
        try:
            #tl = self.axes.get_xticks()
            #tl = tl[1:-1]
            #print('ticks positions ', tl)
            x1,x2 = self.axes.get_xlim()
            c = 299792.458  # speed of light in km/s
            if self.xunit == 'um':
                #v = (tl/s.l0-1.)*c
                #vtl = ["%.1f" % z for z in v]
                vx1 = (x1/s.l0-1.)*c
                vx2 = (x2/s.l0-1.)*c
            elif self.xunit == 'THz':
                #v = (c/s.l0/tl*1.e-3-1.)*c
                #vtl = ["%.1f" % z for z in v]
                vx1 = (c/s.l0/x1*1.e-3-1.)*c
                vx2 = (c/s.l0/x2*1.e-3-1.)*c
            #print('ticks labels ',vtl)
            try:
                self.fig.delaxes(self.vaxes)
                #self.vaxes.remove()
            except:
                pass
            self.vaxes = self.axes.twiny()
            self.vaxes.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
            self.vaxes.set_xlim((vx1,vx2))
            #self.vaxes.set_xticks(v)
            #self.vaxes.set_xticklabels(vtl)
            self.vaxes.set_xlabel("Velocity [km/s]")
        except:
            print('l0 is not defined')
        # Add axes
        if self.displayExposure:
            self.ax3.tick_params(labelright='on',right='on',direction='out',pad=0,colors='orange')
        if self.displayAtran:
            self.ax2.get_yaxis().set_tick_params(labelright='on',right='on')            
            self.ax2.get_yaxis().set_tick_params(which='both', direction='in', pad = -20, colors='red')

        # Add spectral lines
        self.annotations = []
        font = FontProperties(family='DejaVu Sans', size=12)
        xlim0,xlim1 = self.axes.get_xlim()
        ylim0,ylim1 = self.axes.get_ylim()
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
                annotation = self.axes.annotate(nline, xy=(wline,y1),  xytext=(wline, y2), color='blue', alpha=0.4,
                                                arrowprops=dict(color='blue',facecolor='y', arrowstyle='-',alpha=0.4),
                                                rotation = 90, fontstyle = 'italic', fontproperties=font, visible=self.displayLines,)
                self.annotations.append(annotation)     

        for annotation in self.annotations:
            annotation.set_visible(True)
                
        # Prepare legend                
        self.labs = [l.get_label() for l in lns]
        leg = self.axes.legend(lns, self.labs, loc='best',frameon=False,framealpha=0.0)
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
            
        # Shade region considered for the images
        if self.shade == True:
            self.shadeRegion()

    def shadeRegion(self):

        wmin,wmax = self.regionlimits
        if self.xunit == 'um':
            self.region = self.axes.axvspan(wmin,wmax,facecolor='Lavender',alpha=0.5,linewidth=0)
        elif self.xunit == 'THz':
            c = 299792458.0  # speed of light in m/s
            fmax = c/wmin * 1.e-6
            fmin = c/wmax * 1.e-6
            self.region = self.axes.axvspan(fmin,fmax,facecolor='Lavender',alpha=0.5,linewidth=0)            
        
                
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
            if self.shade == True:
                self.shadeRegion()
            self.fig.canvas.draw_idle()
        elif isinstance(event.artist, Text):            
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
            pass
