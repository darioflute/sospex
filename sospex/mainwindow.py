#!/usr/bin/env python
import sys,os
import numpy as np
from PyQt5.QtWidgets import (QApplication, QWidget, QMainWindow, QTabWidget, QTabBar,QHBoxLayout,
                             QVBoxLayout, QSizePolicy, QStatusBar, QSplitter,
                             QToolBar, QAction, QFileDialog,  QTableView, QComboBox, QAbstractItemView,
                             QMessageBox, QInputDialog, QDialog, QLabel, QProxyStyle,QStyle)
from PyQt5.QtGui import QIcon, QStandardItem, QStandardItemModel, QPixmap, QMovie
from PyQt5.QtCore import Qt, QSize, QTimer, QThread, pyqtSignal, pyqtSlot

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.widgets import SpanSelector, PolygonSelector, RectangleSelector, EllipseSelector
from matplotlib.patches import Polygon
from scipy.spatial import cKDTree as KDTree  # cKDTree implementation in C++ (KDTree in pure Python)
import scipy.ndimage as ndimage

# To avoid excessive warning messages
import warnings

from sospex.moments import ( multiFitContinuum, multiComputeMoments, ContParams, ContFitParams, 
                            SlicerDialog, FitCubeDialog, multiFitLines, residualsPsf, guessParams)
from sospex.graphics import  (NavigationToolbar, ImageCanvas, ImageHistoCanvas, SpectrumCanvas,
                       cmDialog, ds9cmap, ScrollMessageBox, PsfCanvas)
from sospex.apertures import (photoAperture, PolygonInteractor, EllipseInteractor,
                              RectangleInteractor, PixelInteractor)
from sospex.specobj import specCube, specCubeAstro, Spectrum, ExtSpectrum
from sospex.cloud import cloudImage
from sospex.interactors import (SliderInteractor, SliceInteractor, DistanceSelector,
                                VoronoiInteractor, LineInteractor, PsfInteractor,
                                InteractorManager, SegmentsSelector, SegmentsInteractor)
from sospex.inout import exportAperture, importAperture, exportGuesses, importGuesses

class MyProxyStyle(QProxyStyle):
    pass
    def pixelMetric(self, QStyle_PixelMetric, option=None, widget=None):

        if QStyle_PixelMetric == QStyle.PM_SmallIconSize:
            return 40
        else:
            return QProxyStyle.pixelMetric(self, QStyle_PixelMetric, option, widget)


#class UpdateTabs(QObject):
#    newImage = pyqtSignal([cloudImage])

class DownloadThread(QThread):
    """Thread to download images from web archives."""
    updateTabs = pyqtSignal([cloudImage])
    sendMessage = pyqtSignal([str])

    def __init__(self,lon,lat,xsize,ysize,band, parent = None):
        super().__init__(parent)
        self.lon = lon
        self.lat = lat
        self.xsize = xsize
        self.ysize = ysize
        self.band = band
        self.parent = parent
        
    def __del__(self):
        self.wait()

    @pyqtSlot()
    def run(self):
        downloadedImage = cloudImage(self.lon,self.lat,self.xsize,self.ysize,self.band)
        if downloadedImage.data is not None:
            #self.updateTabs.newImage.emit(downloadedImage)
            self.updateTabs.emit(downloadedImage)
            message = 'New image downloaded'
        else:
            message = 'The selected survey does not cover the displayed image'
        print(message)
        self.sendMessage.emit(message)
        # Disconnect signal at the end of the thread
        #self.updateTabs.disconnect()
        #self.updateTabs.newImage.disconnect()
        
    def stop(self):
        self.terminate()
        
class UpdateHistogram(QThread):
    #updateHisto = updateHisto()
    sendMessage = pyqtSignal([str])
                
    def __init__(self, image, clim, parent=None):
        super().__init__(parent)
        self.image = image
        self.clim = clim
        
    def run(self):
        pass
        
class GUI (QMainWindow):
    """Main GUI window."""
    
    def __init__(self):
        super().__init__()
        self.title = 'SOSPEX: SOFIA Spectral Explorer'
        self.left = 10
        self.top = 10
        self.width = 640
        self.height = 480
        # Avoid Python deleting object in unpredictable order on exit
        self.setAttribute(Qt.WA_DeleteOnClose)
        # Default color map for images
        ds9cmap()
        self.colorMap = 'real'
        self.colorMapDirection = '_r'
        self.stretchMap = 'linear'
        self.colorContour = ['cyan', 'yellow']
        self.tabContour = [None, None]
        # Set initial state of all spectral tab (all true is computed)
        self.all = False
        # Set status of continuum fit (all if continuum defined for the whole cube)
        self.fitcont = False
        # Set default number of lines to fit across the cube
        self.abslines = 0
        self.emslines = 0
        # Default kernel
        self.kernel = 1
        # Default number of cells
        self.ncells = 1
        # Initial press setting 
        self.press = None
        self.press2 = None
        # Get the path of the package
        self.path0, file0 = os.path.split(__file__)
        # Define style
        with open(os.path.join(self.path0,'yellow-stylesheet.css'),"r") as fh:
            self.setStyleSheet(fh.read())
        self.initUI()
 
    def initUI(self):
        """Define the user interface."""
        self.setWindowTitle(self.title)
        #self.setWindowIcon(QIcon(self.path0+'/icons/sospex.png'))
        self.setGeometry(self.left, self.top, self.width, self.height)
        # Create main widget
        wid = QWidget()
        self.setCentralWidget(wid)
        # Main layout is horizontal
        mainLayout = QHBoxLayout(wid)
        # Horizontal splitter
        self.hsplitter = QSplitter(Qt.Horizontal)        
        # Create main panels
        self.createImagePanel()
        self.createSpectralPanel()
        # Add panels to splitter
        self.hsplitter.addWidget(self.imagePanel)
        self.hsplitter.addWidget(self.spectralPanel)
        # Add panels to main layout
        mainLayout.addWidget(self.hsplitter)
        wid.setLayout(mainLayout)
        self.show()
        # Timer for periodical events
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.blinkTab)
        # Load lines
        from sospex.lines import define_lines
        self.Lines = define_lines()
        # Welcome message
        self.welcomeMessage()
        # Menu
        self.createMenu()

    def welcomeMessage(self):
        self.wbox = QMessageBox()
        pixmap = QPixmap(os.path.join(self.path0,'icons','sospex.png'))
        self.wbox.setIconPixmap(pixmap)
        self.wbox.setText("Welcome to SOSPEX")
        self.wbox.setInformativeText('SOFIA SPectral EXplorer\n\n '+\
                                     '* Click on folder icon to load spectra\n\n'+\
                                     '* Click on running men icon to exit\n\n'+\
                                     '* Click on question mark for further help')
        self.wbox.show()

    def createMenu(self):
        """Menu with all the actions (more complete than the icons)."""
        # File import/save
        bar = self.menuBar()
        file = bar.addMenu("File")
        file.addAction(QAction("Quit",self,shortcut='Ctrl+q',triggered=self.fileQuit))
        file.addAction(QAction("Open cube",self,shortcut='Ctrl+n',triggered=self.newFile))
        file.addAction(QAction("Reload cube",self,shortcut='Ctrl+n',triggered=self.reloadFile))
        file.addAction(QAction("Import image",self,shortcut='Ctrl+d',
                               triggered=self.selectDownloadImage))
        cube = file.addMenu("Save cube")
        cube.addAction(QAction('Trimmed', self, shortcut='',triggered=self.trimCube))
        cube.addAction(QAction('Cropped', self, shortcut='',triggered=self.cropCube))
        cube.addAction(QAction('Continuum subtracted', self, shortcut='',triggered=self.savelCube))
        cube.addAction(QAction('Current status', self, shortcut='',triggered=self.saveMaskedCube))
        file.addAction(QAction('Save image', self, shortcut='',triggered=self.saveFits))
        file.addAction(QAction('Save spectrum', self, shortcut='',triggered=self.saveSpectrum))
        aperture = file.addMenu("Aperture I/O")
        aperture.addAction(QAction('Export',self,shortcut='',triggered=self.exportApertureAction))
        aperture.addAction(QAction('Import',self,shortcut='',triggered=self.importApertureAction))
        guesses = file.addMenu("Guesses I/O")
        guesses.addAction(QAction('Export', self, shortcut='', triggered=self.exportGuessesAction))
        guesses.addAction(QAction('Import', self, shortcut='', triggered=self.importGuessesAction))        
        # View
        view = bar.addMenu("View")
        slice = view.addMenu("Show slider")
        slice.addAction(QAction('on channel', self, shortcut='', triggered=self.initializeSlider))
        slice.addAction(QAction('on slice', self, shortcut='', triggered=self.initializeSlicer))
        slice.addAction(QAction('no', self, shortcut='', triggered=self.removeSliders))
        level = view.addMenu("Image levels")
        level.addAction(QAction('100.0%',self,shortcut='', triggered=self.changeVisibility100))
        level.addAction(QAction('99.5%',self,shortcut='', triggered=self.changeVisibility995))
        level.addAction(QAction('99.0%',self,shortcut='', triggered=self.changeVisibility990))
        level.addAction(QAction('98.0%',self,shortcut='', triggered=self.changeVisibility980))
        level.addAction(QAction('95.0%',self,shortcut='', triggered=self.changeVisibility950))
        level.addAction(QAction('90.0%',self,shortcut='', triggered=self.changeVisibility900))
        level.addAction(QAction('80.0%',self,shortcut='', triggered=self.changeVisibility800))
        #self.menuHisto = QAction('Image levels',self,shortcut='',
        #                         checkable = True, triggered=self.changeVisibility)
        #view.addAction(self.menuHisto)
        view.addAction(QAction('Colors and stretch',self,shortcut='',triggered=self.changeColorMap))
        view.addAction(QAction('Show header',self,shortcut='',triggered=self.showHeader))
        magnify = view.addMenu("Magnify image")
        magnify.addAction(QAction('+10%',self,shortcut='',triggered=self.zoomUp))
        magnify.addAction(QAction('-10%',self,shortcut='',triggered=self.zoomDown))
        kernel = view.addMenu("Choose kernel for spectrum")
        self.kernel1 = QAction('1 pixel',self,shortcut='',checkable=True,
                               triggered=self.kernel1pixel)
        self.kernel5 = QAction('5 pixels',self,shortcut='',checkable=True,
                               triggered=self.kernel5pixel)
        self.kernel9 = QAction('9 pixels',self,shortcut='',checkable=True,
                               triggered=self.kernel9pixel)
        kernel.addAction(self.kernel1)
        kernel.addAction(self.kernel5)
        kernel.addAction(self.kernel9)
        self.kernel1.setChecked(True)
        contours = view.addMenu("Contours")
        # Checkable is for checking if contour is 'on', valid only for
        # one set of contours. In the future, more than one contour will be available
        contours.addAction(QAction('Add contour',self,shortcut='',
                                    checkable=True,triggered=self.overlapContours))
        contours.addAction(QAction('Clear contours',self,shortcut='',
                                    checkable=True,triggered=self.overlapContours))
        # Tools
        tools = bar.addMenu("Tools")
        flux = tools.addMenu("Recompute flux")
        flux.addAction(QAction('.. with atm. transmission at ref. wavelength',self,shortcut='',
                               triggered=self.fluxRefWavAT))        
        flux.addAction(QAction('.. with median atm. transmission',self,shortcut='',
                               triggered=self.fluxMedianAT))        
        flux.addAction(QAction('.. with new atm. transmission',self,shortcut='',
                               triggered=self.fluxNewAT))        
        erase = tools.addMenu("Mask part of cube")
        erase.addAction(QAction('.. lower than min contour level',self,shortcut='',
                                triggered=self.maskCubeContour))
        erase.addAction(QAction('.. inside a polygon',self,shortcut='',
                                triggered=self.maskCubeInsidePolygon))
        erase.addAction(QAction('.. outside a polygon',self,shortcut='',
                                triggered=self.maskCubeOutsidePolygon))
        continuum = tools.addMenu("Continuum and moments")
        continuum.addAction(QAction('Define guesses',self,shortcut='',
                                    triggered=self.guessContinuum))
        continuum.addAction(QAction('Compute',self,shortcut='',
                                    triggered=self.ContMomLines))
        #continuum.addAction(QAction('Fit all cube',self,shortcut='',triggered=self.fitContAll))
        #continuum.addAction(QAction('Fit inside region',self,shortcut='',
        #                            triggered=self.fitContRegion))
        #continuum.addAction(QAction('Set continuum to zero ',self,shortcut='',
        #                            triggered=self.setContinuumZero))        
        #continuum.addAction(QAction('Set continuum to medians ',self,shortcut='',
        #                            triggered=self.setContinuumMedian))        
        #moments = tools.addMenu("Compute moments")
        # moments.addAction(QAction('Define slice',self,shortcut='',triggered=self.sliceCube))
        #moments.addAction(QAction('Compute all cube',self,shortcut='',
         #                         triggered=self.computeMomentsAll))
        #moments.addAction(QAction('Compute inside region',self,shortcut='',
         #                         triggered=self.computeRegion))
        tools.addAction(QAction(u'Recompute C\u2080, v, \u03c3\u1d65',self,shortcut='',
                                triggered=self.computeVelocities))
        apertures = tools.addMenu("Select aperture")
        apertures.addAction(QAction('Square',self,shortcut='',triggered=self.selectSquareAperture))
        apertures.addAction(QAction('Rectangle',self,shortcut='',
                                    triggered=self.selectRectangleAperture))
        apertures.addAction(QAction('Circle',self,shortcut='',triggered=self.selectCircleAperture))
        apertures.addAction(QAction('Ellipse',self,shortcut='',
                                    triggered=self.selectEllipseAperture))
        apertures.addAction(QAction('Polygon',self,shortcut='',
                                    triggered=self.selectPolygonAperture))
        # Help 
        help = bar.addMenu("Help")
        help.addAction(QAction('About', self, shortcut='Ctrl+a',triggered=self.about))
        help.addAction(QAction('Tutorials', self, shortcut='Ctrl+h',triggered=self.onHelp))
        help.addAction(QAction('Issues', self, shortcut='Ctrl+i',triggered=self.onIssue))
        bar.setNativeMenuBar(False)
        
    def exportApertureAction(self):
        exportAperture(self)
        
    def importApertureAction(self):
        importAperture(self)
        
    def exportGuessesAction(self):
        exportGuesses(self)
        
    def importGuessesAction(self):
        importGuesses(self)
        print('Imported ', self.ncells, ' cells.')

    def showHeader(self):
        """ Show header of the spectral cube """
        # NOTE: the function repr can be used to format the header printout
        try:
            header = self.specCube.header
            hv = header.values()
            hk = header.keys()
            hc = header.comments
            h = []
            cc = False
            for k,v,c in zip(hk,hv,hc):
                if k == 'HISTORY':
                    pass
                elif k == 'COMMENT':
                    if cc:
                        h.append('        {0:20}'.format(v)+''+c)
                    else:
                        h.append('\n        {0:20}'.format(v)+''+c)
                        cc = True
                else:
                    if cc:
                        h.append('\n')
                        cc = False
                    if isinstance(v, str):
                        vs = v.rjust(20)
                    else:
                        vs = str(v)
                        vs=vs.rjust(20)
                    k = k.ljust(8)
                    c = c.ljust(30)
                    h.append('{k:8s} = {v:20s} {c:30s}'.format(k=k,v=vs,c=c))    
            #message = '\n'.join(h)
            #QMessageBox.about(self, "Header", message)
            msgBox = ScrollMessageBox(h, None)
            msgBox.setWindowTitle("Header")
            msgBox.exec_()
        except:
            self.sb.showMessage("First load a spectral cube ", 1000)

    def about(self):
        from sospex import __version__
        # Get path of the package
        file=open(os.path.join(self.path0,"copyright.txt"),"r")
        message=file.read()
        message = 'SOSPEX - version: ' + __version__ + '\n' + message
        QMessageBox.about(self, "About", message)

    def createImagePanel(self):
        """Panel to display images."""
        #self.imagePanel = QGroupBox("")
        self.imagePanel = QWidget()
        layout = QVBoxLayout(self.imagePanel)
        # Tabs with images        
        self.itabs = QTabWidget()
        self.itabs.setTabsClosable(True)
        self.itabs.tabCloseRequested.connect(self.removeTab)
        #self.itabs.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding)
        self.itabs.setSizePolicy(QSizePolicy.Minimum,QSizePolicy.Minimum)
        self.itabs.currentChanged.connect(self.onITabChange)  # things to do when changing tab
        self.tabi = []
        self.ici  = []
        self.ihi  = []
        self.ihcid = []
        self.icid1 = []
        self.icid2 = []
        self.icid3 = []
        self.icid4 = []
        # Add widgets to panel
        layout.addWidget(self.itabs)

    def addSpectrum(self, b):
        ''' Add a tab with a spectrum '''
        t = QWidget()
        t.layout = QHBoxLayout(t)
        t.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Ignored) # Avoid expansion
        self.stabs.addTab(t, b)
        sc = SpectrumCanvas(t, width=5.5, height=5.25, dpi=100)
        sc.switchSignal.connect(self.switchUnits)
        # Toolbar
        toolbar = QToolBar()
        toolbar.setOrientation(Qt.Vertical)
        # Check if new aperture
        # Add actions to toolbar
        toolbar.addAction(self.cutAction)
        toolbar.addAction(self.specAction)
        print('Aperture is ', b)
        if b in ['All', 'Pix']:        
            toolbar.addAction(self.guessAction)
        else:
            print('Aperture special actions ')
            toolbar.addAction(self.guessLinesAction)
            toolbar.addAction(self.fitLinesAction)            
        # toolbar.addAction(self.sliceAction)
        toolbar.addAction(self.slideAction)
        toolbar.addSeparator()
        toolbar.addSeparator()
        toolbar.addAction(self.hresizeAction)
        toolbar.addAction(self.vresizeAction)
        # Navigation toolbar
        sc.toolbar = NavigationToolbar(sc, self)
        #sc.toolbar.resize(100,10) 
        # Panel
        panel = QWidget()
        panel.layout = QVBoxLayout(panel)
        panel.layout.addWidget(sc)
        panel.layout.addWidget(sc.toolbar)
        t.layout.addWidget(panel)
        t.layout.addWidget(toolbar)
        self.stabs.resize(self.stabs.minimumSizeHint())  # Avoid expansion
        # connect image and histogram to  events
        scid1 = sc.mpl_connect('button_release_event', self.onDraw2)
        scid2 = sc.mpl_connect('scroll_event', self.onWheel2)
        scid3 = sc.mpl_connect('key_press_event', self.onKeyPress2)
        scid4 = sc.mpl_connect('key_release_event', self.onKeyRelease2)
        scid5 = sc.mpl_connect('motion_notify_event', self.onMotion2)
        scid6 = sc.mpl_connect('button_press_event', self.onPress2)
        self.ctrlIsHeld = False
        return t, sc, scid1, scid2, scid3, scid4, scid5, scid6

    def addPsfPlot(self):
        '''Add a tab with the plot of the PSF'''
        t = QWidget()
        t.layout = QVBoxLayout(t)
        t.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Ignored) # Avoid expansion
        sc = PsfCanvas()
        # Toolbar
        toolbar = QToolBar()
        # Navigation toolbar
        # Add actions to toolbar
        toolbar.addAction(self.centroidAction)
        toolbar.addAction(self.centerAction)
        sc.toolbar = NavigationToolbar(sc, self)
        foot = QWidget()
        foot.layout = QHBoxLayout(foot)
        foot.layout.addWidget(toolbar)
        foot.layout.addWidget(sc.toolbar)
        # 
        self.stabs.addTab(t, 'PSF')
        t.layout.addWidget(sc)
        t.layout.addWidget(foot)
        self.stabs.resize(self.stabs.minimumSizeHint())  # Avoid expansion
        self.stabi.append(t)
        self.sci.append(sc)
        self.scid1.append(None)
        self.scid2.append(None)
        self.scid3.append(None)
        self.scid4.append(None)
        self.scid5.append(None)
        self.scid6.append(None)
        # Initialize plot
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        ny, nx = np.shape(ic.oimage)
        xi = np.arange(nx); yi = np.arange(ny)
        xi, yi = np.meshgrid(xi, yi)
        self.imagePoints = (xi, yi)
        # Compute distances and fluxes
        dist, flux = self.computePsfData()
        sc.compute_initial_psf(distance=dist, flux=flux,
                               instrument=self.specCube.instrument,
                               w=self.specCube.l0,
                               pix=ic.pixscale)

    def addImage(self, b):
        '''Add a tab with an image.'''
        t = QWidget()
        t.layout = QHBoxLayout(t)
        t.setSizePolicy(QSizePolicy.Ignored,QSizePolicy.Ignored) # Avoid expansion
        if b == 'sv':
            self.itabs.addTab(t, u'\u03c3\u1d65')  # unicode for sigma
        elif b == 'Flux':
            self.itabs.addTab(t, u'F')  # unicode
        elif b == 'uFlux':
            self.itabs.addTab(t, u'F\u1d64')  # unicode
            #self.itabs.addTab(t, u'F\u1d64\u2099')  # unicode (the n does not work on mac os-x)
            #self.itabs.addTab(t, "F<sub>unc</sub>")  # markup
        elif b == 'Exp':
            self.itabs.addTab(t, u'E')  # unicode
        elif b == 'C0':
            self.itabs.addTab(t, u'C\u2080')  # unicode for 0
        elif b == 'M0':
            self.itabs.addTab(t, u'M\u2080')  # unicode for 0
        elif b == 'M1':
            self.itabs.addTab(t, u'M\u2081')  # unicode for 1
        elif b == 'M2':
            self.itabs.addTab(t, u'M\u2082')  # unicode for 2
        elif b == 'M3':
            #self.itabs.addTab(t, u'M\u2083')  # unicode for 3
            self.itabs.addTab(t, u'Sk')  # Skewness
        elif b == 'M4':
            #self.itabs.addTab(t, u'M\u2084')  # unicode for 4
            self.itabs.addTab(t, u'K')  # Kurtosis
        elif b == 'L0':
            self.itabs.addTab(t, u'I\u2080')  # unicode 0
        elif b == 'L1':
            self.itabs.addTab(t, u'I\u2081')  # unicode 1      
        else:
            self.itabs.addTab(t, b)
        ic = ImageCanvas(t, width=11, height=10.5, dpi=100)
        if b in ['Flux','uFlux','Exp','C0','M0','M1','M2','M3','M4','L0','L1']:
            ic.crota2 = self.specCube.crota2
        # No contours available
        ic.contours = None
        ic.contour0 = None
        ih = ImageHistoCanvas(t, width=11, height=0.5, dpi=100)
        ih.setVisible(False)
        ic.toolbar = NavigationToolbar(ic, self)
        # Toolbar
        toolbar = QToolBar()
        toolbar.setOrientation(Qt.Vertical)
        toolbar.addAction(self.levelsAction)
        toolbar.addAction(self.cmapAction)
        toolbar.addAction(self.blinkAction)
        toolbar.addAction(self.contoursAction)
        #toolbar.addAction(self.momentAction)
        toolbar.addAction(self.cropAction)
        toolbar.addAction(self.cloudAction)
        toolbar.addAction(self.fitsAction)
        #toolbar.addAction(self.fitregionAction)
        toolbar.addAction(self.fitContAction)
        # toolbar.addAction(self.compMomAction)
        toolbar.addAction(self.maskAction)
        toolbar.addAction(self.distanceAction)
        toolbar.addAction(self.psfAction)
        toolbar.addSeparator()
        #toolbar.addWidget(self.apertureAction)        
        # Panel
        panel = QWidget()
        panel.layout = QVBoxLayout(panel)
        panel.layout.addWidget(ic)
        panel.layout.addWidget(ih)
        panel.layout.addWidget(ic.toolbar)        
        t.layout.addWidget(toolbar)
        t.layout.addWidget(panel)
        self.itabs.resize(self.itabs.minimumSizeHint())  # Avoid expansion
        # connect image and histogram to  events
        cidh = ih.limSignal.connect(self.onChangeIntensity)
        cid1 = ic.mpl_connect('button_release_event', self.onDraw)
        cid2 = ic.mpl_connect('scroll_event', self.onWheel)
        cid3 = ic.mpl_connect('motion_notify_event', self.onMotion)
        cid4 = ic.mpl_connect('button_press_event', self.onPress)        
        return t, ic, ih, cidh, cid1, cid2, cid3, cid4

    def removeTab(self, itab, confirm=True):
        if confirm:
            close_msg = "Are you sure you want to close this image tab ?"
            reply = QMessageBox.question(self, 'Message', close_msg, QMessageBox.Yes, QMessageBox.No)
            if reply == QMessageBox.No:
                return
        widget = self.itabs.widget(itab)
        if widget is not None:
            widget.deleteLater()
        self.itabs.removeTab(itab)
        # Disconnect and remove canvases
        tab = self.tabi[itab]
        ima = self.ici[itab]
        his = self.ihi[itab]
        hcid = self.ihcid[itab]
        c1 = self.icid1[itab]
        c2 = self.icid2[itab]
        c3 = self.icid3[itab]
        c4 = self.icid4[itab]
        his.mpl_disconnect(hcid)
        ima.mpl_disconnect(c1)
        ima.mpl_disconnect(c2)
        ima.mpl_disconnect(c3)
        ima.mpl_disconnect(c4)
        self.tabi.remove(tab)
        self.ici.remove(ima)
        self.ihi.remove(his)
        self.ihcid.remove(hcid)
        self.icid1.remove(c1)
        self.icid2.remove(c2)
        self.icid3.remove(c3)
        self.icid4.remove(c4)
        ima = None
        his = None
        # Put back to None the band
        b = self.bands[itab]
        if b == 'C0': self.C0 = None
        elif b == 'M0': self.M0 = None
        elif b == 'M1': self.M1 = None
        elif b == 'M2': self.M2 = None
        elif b == 'M3': self.M3 = None
        elif b == 'M4': self.M4 = None
        elif b == 'v' : self.v  = None
        elif b == 'sv': self.sv = None
        elif b == 'L0': self.L0 = None
        elif b == 'L1': self.L1 = None
        # Remove band from band list
        del self.bands[itab]

    def removeSpecTab(self, stab, confirm=True):
        if confirm:
            close_msg = "Are you sure you want to close this tab ?"
            reply = QMessageBox.question(self, 'Message', close_msg, QMessageBox.Yes, QMessageBox.No)
            if reply == QMessageBox.No:
                return
        tabname = self.stabs.tabText(stab)
        print('Removing spectral tab no ',stab, ' name ', tabname,' canvas name ',self.sci[stab].name)
        itab = self.itabs.currentIndex()
        ic0 = self.ici[itab]
        if tabname in ['All','Pix']:
            pass
        elif tabname == 'PSF':
            # Remove PSF apertures
            self.PsfI.mySignal.disconnect()
            self.PsfI.disconnect()
            del self.PsfI
            self.PsfI = None
            ic0.fig.canvas.draw_idle()
        else:
            # Once the tab is removed, also the relative aperture should be removed
            # n = stab-1
            scname = self.sci[stab].name
            n = int(scname)
            print('removing aperture ',n,' type: ',ic0.photApertures[n].type)
            for ic in self.ici:
                ap = ic.photApertures[n]
                #aps = 
                ic.photApertureSignal[n]
                #print('removing the aperture: ',ap.type)
                ap.mySignal.disconnect()
                ap.disconnect()
                del ic.photApertureSignal[n]
                del ic.photApertures[n]
            # Remove photoAperture
            del self.photoApertures[n]
            # Remove from self.spectra
            self.spectra.remove(scname)
            # Redraw apertures
            ic0.fig.canvas.draw_idle()
        
        # Remove tab
        widget = self.stabs.widget(stab)
        if widget is not None:
            widget.deleteLater()
        # Disconnect and remove canvases
        tab = self.stabi[stab]
        spec = self.sci[stab]
        c1 = self.scid1[stab]
        c2 = self.scid2[stab]
        c3 = self.scid3[stab]
        c4 = self.scid4[stab]
        c5 = self.scid5[stab]
        c6 = self.scid6[stab]
        if tabname != 'PSF':
            spec.mpl_disconnect(c1)
            spec.mpl_disconnect(c2)
            spec.mpl_disconnect(c3)
            spec.mpl_disconnect(c4)
            spec.mpl_disconnect(c5)
            spec.mpl_disconnect(c6)
        self.stabi.remove(tab)
        self.sci.remove(spec)
        self.scid1.remove(c1)
        self.scid2.remove(c2)
        self.scid3.remove(c3)
        self.scid4.remove(c4)
        self.scid5.remove(c5)
        self.scid6.remove(c6)
        spec = None
        # Remove stab from list
        self.stabs.removeTab(stab)
        #print('New number of spectral tabs ', len(self.stabs))
        # Rename aperture tabs
        i = 1
        for stab in range(len(self.stabs)):
            tabname = self.stabs.tabText(stab)
            #print('tab name ', tabname)
            if tabname in ['All','Pix','PSF']:
                pass
            else:
                apname = "{:d}".format(i)
                #print('Change ', tabname,' into ', apname)
                self.stabs.setTabText(stab, apname)
                self.sci[stab].name = apname
                i += 1
        #if len(self.stabs) > 2:
        #    for i in range(2, len(self.stabs)):
        #        apname = "{:d}".format(i-1)
        #        self.stabs.setTabText(i,apname)
        
    def onITabChange(self, itab):
        ''' When tab changes check if latest update of ellipse are implemented '''
        if len(self.ici) == 0:
            return
        if itab < len(self.ici):
            ima = self.ici[itab]
            if len(self.stabs) > 1:
                # Check if vertices are activated correctly
                # istab = self.stabs.currentIndex()
                #nap = len(self.stabs)-1
                n = self.nAper()
                # n = istab-1  # aperture number
                # Activate interactor (toogle on) and disactivate
                #nap = len(ima.photApertures)
                for iap, ap in enumerate(ima.photApertures):
                    #ap = ima.photApertures[iap]
                    if iap == n:
                        ap.showverts = True
                    else:
                        ap.showverts = False
                    ap.updateMarkers()
                    ap.line.set_visible(ap.showverts)
                ima.changed = True
            if ima.changed:
                # Compute contours the 1st time
                if self.contours == 'off':
                    ima.contour0 = None
                if ima.contour0 is not None:
                    if isinstance(ima.contour0, int):
                        # print('Computing contours ...')
                        from astropy.nddata import Cutout2D
                        itab = ima.contour0
                        ic0 = self.ici[itab]
                        ih0 = self.ihi[itab]
                        # Consider only the part visible in the image
                        x = ima.axes.get_xlim()
                        y = ima.axes.get_ylim()
                        verts = [(i, j) for i, j in zip(x, y)]
                        adverts = np.array([(ima.wcs.all_pix2world(x, y, 0)) for (x,y) in verts])                
                        verts = [(ic0.wcs.all_world2pix(ra, dec, 0)) for (ra,dec) in adverts]
                        x, y = zip(*verts)
                        center =  ((x[0]+x[1])*0.5,(y[0]+y[1])*0.5)
                        size = (np.abs((y[1]-y[0]).astype(int)),np.abs((x[1]-x[0]).astype(int)))
                        co = Cutout2D(ic0.oimage, center, size, wcs=ic0.wcs, mode='partial')
                        # If the contour image is much large, one should lower the resolution
                        # to speed-up the computation
                        ny0, nx0 = np.shape(co.data)
                        ny ,nx = np.shape(ima.oimage)
                        levs =  sorted(ih0.levels)
                        # Smooth data
                        #ismo = ndimage.gaussian_filter(co.data, sigma=0.1, order=0)
                        ismo = co.data
                        if nx0 > 2 * nx:
                            # print('Reprojecting image for contours')
                            from reproject import reproject_interp
                            from astropy.io import fits
                            hdu = fits.PrimaryHDU(ima.oimage)
                            hdu.header.extend(ima.wcs.to_header())
                            # print("header", hdu.header)
                            hdu0 = fits.PrimaryHDU(ismo)
                            hdu0.header.extend(co.wcs.to_header())
                            array, footprint = reproject_interp(hdu0, hdu.header)
                            ima.contour = ima.axes.contour(array, levs, colors=self.colorContour[0])
                        else:
                            ima.contour = ima.axes.contour(ismo, levs, colors=self.colorContour[0],
                                                           transform=ima.axes.get_transform(co.wcs))
                        ima.contour0 = None
                ima.fig.canvas.draw_idle()
                ima.changed = False
            if self.blink == 'select':
                # Select 2nd tab and start blinking until blink status changes ...
                self.btab[1] = itab
                self.blink = 'on'
                self.timer.start(1000)
        
    def onChangeIntensity(self, event):
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        ih = self.ihi[itab]
        # apply intensity limits to the relative figure
        ic.image.set_clim(ih.limits)
        ic.fig.canvas.draw_idle()

    def onHelp(self, event):
        import webbrowser
        webbrowser.open('file://'+os.path.abspath(os.path.join(self.path0,'help','Help.html')))

    def onIssue(self, event):
        import webbrowser
        webbrowser.open('https://github.com/darioflute/sospex/issues')

    def onPress(self, event):
        if event.inaxes:
            self.press = event.xdata, event.ydata
        else:
            self.press = None
            return

    def onMotion(self, event):
        """ Update spectrum when moving an aperture on the image """
        # Update bias and contrast with third mouse button
        if self.press is None:
            return
        if event.inaxes:
            dx = event.xdata - self.press[0]
            dy = event.ydata - self.press[1]
            # self.press = (event.x, event.y)
        else:
            return
        if event.button == 3:
            # Get current bias and contrast
            itab = self.itabs.currentIndex()
            ih = self.ihi[itab]
            cmin = ih.limits[0]
            cmax = ih.limits[1]
            if np.abs(dx) >= np.abs(dy):
                percent = 0.01
                diff = (cmax-cmin)*percent*np.sign(dx)
                cmin += diff
                cmax += diff
            else:
                percent = 0.04
                mid  = (cmin+cmax)*0.5
                diff = np.abs(cmax-cmin)*0.5*(1.+percent*np.sign(dy))           
                cmin = mid - diff
                cmax = mid + diff
            ih.onSelect(cmin,cmax)
        elif event.button == 2:
            # Pan using the mouse middle button
            itab = self.itabs.currentIndex()
            ic = self.ici[itab]
            xlim = np.asarray(ic.axes.get_xlim())
            ylim = np.asarray(ic.axes.get_ylim())
            x = xlim - dx
            y = ylim - dy
            ic.axes.set_xlim(x)
            ic.axes.set_ylim(y)
            ic.fig.canvas.draw_idle()

    def updateAperture(self):
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        ic0 = self.ici[0]
        #nap = self.stabs.currentIndex()-1
        nap = self.nAper()
        aper = ic.photApertures[nap]
        aper0 = ic0.photApertures[nap]
        if aper.type == 'Polygon':
            verts = aper.poly.get_xy()
            adverts = np.array([(ic.wcs.all_pix2world(x, y, 0)) for (x,y) in verts])                
            verts = [(ic0.wcs.all_world2pix(ra, dec, 0)) for (ra,dec) in adverts]
            aper0.poly.set_xy(verts)
        elif aper.type == 'Ellipse' or aper.type == 'Circle':
            x0,y0 = aper.ellipse.center
            w0    = aper.ellipse.width
            h0    = aper.ellipse.height
            angle = aper.ellipse.angle
            ra0,dec0 = ic.wcs.all_pix2world(x0, y0, 0)
            ws = w0 * ic.pixscale; hs = h0 * ic.pixscale
            x0,y0 = ic0.wcs.all_world2pix(ra0, dec0, 0)
            w0 = ws/ic0.pixscale; h0 = hs/ic0.pixscale
            aper0.ellipse.center = x0, y0
            aper0.ellipse.width = w0
            aper0.ellipse.height = h0
            aper0.ellipse.angle = angle
        elif aper.type == 'Rectangle' or aper.type == 'Square' or aper.type == 'Pixel':
            x0,y0 = aper.rect.get_xy()
            w0    = aper.rect.get_width()
            h0    = aper.rect.get_height()
            angle = aper.rect.angle
            ra0,dec0 = ic.wcs.all_pix2world(x0, y0, 0)
            ws = w0 * ic.pixscale; hs = h0 * ic.pixscale
            x0,y0 = ic0.wcs.all_world2pix(ra0, dec0, 0)
            w0 = ws/ic0.pixscale; h0 = hs/ic0.pixscale
            aper0.rect.set_xy((x0,y0))
            aper0.rect.set_width(w0)
            aper0.rect.set_height(h0)
            aper0.rect.angle = angle

    def onRemoveContinuum(self, event):
        istab = self.stabs.currentIndex()
        sc = self.sci[istab]
        # print('remove continuum event: ', event)
        if event[:-2] == 'line deleted':
            n = int(event[-2:])
            print('Disconnect line ', n)
            sc.lines[n].disconnect
            sc.lines[n] = None
            for line in sc.lines:
                print('line ', line)
            # No need to update guesses since each line has an ID (line.n)
        elif event == 'segments deleted':
            sc.guess.disconnect()
            sc.guess = None
            if sc.lines is not None:
                for line in sc.lines:
                    line.disconnect()
        sc.fig.canvas.draw_idle()

    def onRemoveAperture(self,event):
        """Interpret signal from apertures."""        
        itab = self.itabs.currentIndex()
        istab = self.stabs.currentIndex()
        # n = istab-1
        n = self.nAper()
        ap = self.ici[itab].photApertures[n]
        apertype = ap.__class__.__name__
        if (event == 'rectangle deleted' and apertype == 'RectangleInteractor') or \
           (event == 'ellipse deleted' and apertype == 'EllipseInteractor') \
           or (event == 'polygon deleted' and apertype == 'PolygonInteractor'):
            self.stabs.currentChanged.disconnect()
            self.removeSpecTab(istab, False)
            self.stabs.setCurrentIndex(1)  # Pixel tab
            self.stabs.currentChanged.connect(self.onSTabChange)
        else:
            print(event)

    def onModifiedAperture(self, event):
        """Update spectrum when aperture is modified."""
        itab = self.itabs.currentIndex()
        # Grab aperture in the flux image to compute the new fluxes
        istab = self.stabs.currentIndex()
        sc = self.sci[istab]
        # print('modify aperture on ', sc.name)
        if sc.name == 'PSF': # Ideally it should go to the stab which correspond to the interactor
            # Check activated aperture
            ic = self.ici[itab]
            for i, ap in enumerate(ic.photApertures):
                if ap.showverts:
                    n = i
            if n == 0:       
                istab = self.spectra.index('Pix')
            else:
                istab = n
            self.stabs.setCurrentIndex(istab)
            sc = self.sci[istab]
        else:
            n = self.nAper()
        if istab > 0:
            #print('changing apertures')
            #self.onDraw('changed aperture')
            s = self.specCube
            # I should find a way to know if the aperture has changed
            if itab != 0:
                self.updateAperture()
            #aperture = self.ici[0].photApertures[istab-1].aperture
            aperture = self.ici[0].photApertures[n].aperture
            path = aperture.get_path()
            transform = aperture.get_patch_transform()
            npath = transform.transform_path(path)
            inpoints = s.points[npath.contains_points(s.points)]
            xx,yy = inpoints.T
            # If tab is pix, then compute the center to decide which Voronoi cells it belongs
            if istab == 1:
                if self.ncells > 1:
                    xc = np.median(xx)
                    yc = np.median(yy)
                    if np.isfinite(xc) & np.isfinite(yc):
                        yc = int(yc)
                        xc = int(xc)
                        ny, nx = np.shape(self.regions)
                        if yc >= ny:
                            yc = ny - 1
                        if xc >= nx:
                            xc = nx - 1
                        ncell = self.regions[int(yc), int(xc)]
                    else:
                        ncell = 0
                else:
                    ncell = 0
            else:
                ncell = 0
            if s.instrument == 'GREAT':
                t2j = self.specCube.Tb2Jy
            else:
                t2j = 1.
            if istab == 1 and self.continuum is not None:
                xc=np.median(xx); yc = np.median(yy)
                if np.isfinite(xc) & np.isfinite(yc):
                    i = int(np.rint(xc)); j = int(np.rint(yc))
                    try:
                        cont = self.continuum[:, j, i]*t2j
                        cslope=self.Cs[j,i]*t2j
                    except:
                        cont = None
                        cslope=None
                    try:
                        moments = [self.M0[j, i]*t2j, self.M1[j, i], self.M2[j, i],
                                   self.M3[j, i], self.M4[j, i]]
                        noise = self.noise[j, i]
                    except:
                        moments = None
                        noise = None
                    try:
                        lines = []
                        for line in self.lines:
                            lines.append([line[0][j, i], line[1][j, i], 
                                          line[2][j, i], line[3][j, i]])
                    except:
                        lines = None
                else:
                    cont = None
                    cslope=None
                    moments = None
                    noise = None
                    lines = None                    
            else:
                cont = None
                cslope=None
                moments = None
                noise = None
                lines = None
            if istab == 1: # case of pixel (with different kernels)
                fluxAll = np.nanmean(s.flux[:, yy, xx], axis=1)
                if sc.auxiliary:
                    # Get the pixel of the auxiliary cube from aperture pixel in an image
                    x0, y0 = aperture.get_xy()
                    ic0 = self.ici[0]
                    ra0, dec0 = ic0.wcs.all_pix2world(x0, y0, 0) 
                    xxa, yya = self.auxSpecCube.wcs.all_world2pix(ra0, dec0, 0)
                    xxa = int (xxa // 1)
                    yya = int (yya // 1)
                    try:
                        afluxAll = self.auxSpecCube.flux[:, yya, xxa]
                        # Normalization
                        # I should pass the cube and show the scale on the right
                        #nmax = len(fluxAll) // 3
                        #fmax = np.nanmax(fluxAll[nmax:-nmax]) # Avoid borders
                        #fmin = np.nanmedian(fluxAll)
                        #anmax = len(afluxAll) // 3
                        #afmax = np.nanmax(afluxAll[anmax:-anmax]) # Avoid borders
                        #afmin = np.nanmedian(afluxAll)
                        #afluxAll = (afluxAll - afmin) / (afmax - afmin) * (fmax - fmin) + fmin
                        # afluxAll *= np.nanmax(fluxAll)/np.nanmax(afluxAll)
                    except:
                        afluxAll = self.auxSpecCube.flux[:,0,0] * np.nan
                        #print('Pixel out of map')
                    # Normalization
                    sc.aflux = afluxAll
            else:
                fluxAll = np.nansum(s.flux[:, yy, xx], axis=1)
            sc.spectrum.flux = fluxAll
            if s.instrument in ['GREAT','HI','VLA','ALMA','MUSE','IRAM','CARMA','PCWI']:
                if sc.auxiliary:
                    sc.updateSpectrum(f=fluxAll*t2j, af=afluxAll, cont=cont, cslope=cslope, moments=moments, 
                                      lines=lines, noise=noise, ncell=ncell)
                else:
                    sc.updateSpectrum(f=fluxAll*t2j, cont=cont,cslope=cslope, moments=moments, 
                                      lines=lines, noise=noise, ncell=ncell)
            elif s.instrument in ['PACS','FORCAST']:
                if istab == 1:
                    expAll = np.nanmean(s.exposure[:, yy, xx], axis=1)
                else:
                    expAll = np.nansum(s.exposure[:, yy, xx], axis=1)
                if sc.auxiliary:
                    sc.updateSpectrum(f=fluxAll, af=afluxAll, exp=expAll, cont=cont, cslope=cslope,
                                      moments=moments, lines=lines,
                                      noise=noise, ncell=ncell)
                else:
                    sc.updateSpectrum(f=fluxAll, exp=expAll, cont=cont, cslope=cslope,
                                      moments=moments, lines=lines,
                                      noise=noise, ncell=ncell)
            elif s.instrument == 'FIFI-LS':
                if istab == 1:
                    ufluxAll = np.nanmean(s.uflux[:, yy, xx], axis=1)
                    efluxAll = np.sqrt(np.nanmean(s.eflux[:,yy,xx]**2, axis=1))
                    expAll = np.nanmean(s.exposure[:, yy, xx], axis=1)
                    #print('modified flux ',np.size(xx))
                else:
                    ufluxAll = np.nansum(s.uflux[:, yy, xx], axis=1)
                    efluxAll = np.sqrt(np.nansum(s.eflux[:, yy, xx]**2, axis=1))
                    expAll = np.nansum(s.exposure[:, yy, xx], axis=1)                    
                sc.spectrum.uflux = ufluxAll
                sc.spectrum.eflux = efluxAll
                sc.spectrum.exposure = expAll
                #print('updating aperture ...')
                if sc.auxiliary:
                    sc.updateSpectrum(f=fluxAll, af=afluxAll, uf=ufluxAll, exp=expAll,
                                      cont=cont,cslope=cslope, moments=moments, lines=lines, 
                                      noise=noise, ncell=ncell)
                else:
                    sc.updateSpectrum(f=fluxAll, uf=ufluxAll, exp=expAll,
                                      cont=cont,cslope=cslope, moments=moments, lines=lines, 
                                      noise=noise, ncell=ncell)

           
    def onDraw(self,event):
        if len(self.ici) <= 1:
            return
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        # Deselect pan option on release of mouse
        if ic.toolbar._active == "PAN":
            ic.toolbar.pan()
        # Update patch in all the images
        # a status should be added to the apertures to avoid unnecessary redrawings
        istab = self.stabs.currentIndex()
        ici = self.ici.copy()
        ici.remove(ic)
        tabname = self.stabs.tabText(istab)
        if tabname in ['All', 'PSF']:
            return
        else:
            if tabname == 'Pix':
                ntab = 0
            else:
                ntab = int(tabname)
            aper = ic.photApertures[ntab]
            #apertype = aper.__class__.__name__
            if aper.type == 'Ellipse' or aper.type == 'Circle':
                x0,y0 = aper.ellipse.center
                w0    = aper.ellipse.width
                h0    = aper.ellipse.height
                angle = aper.ellipse.angle
                ra0, dec0 = ic.wcs.all_pix2world(x0, y0, 0)
                ws = w0 * ic.pixscale; hs = h0 * ic.pixscale
                for ima in ici:
                    x0,y0 = ima.wcs.all_world2pix(ra0, dec0, 0)
                    w0 = ws/ima.pixscale; h0 = hs / ima.pixscale
                    ap = ima.photApertures[ntab]
                    ap.ellipse.center = x0, y0
                    ap.ellipse.width = w0
                    ap.ellipse.height = h0
                    ap.ellipse.angle = angle - ic.crota2 + ima.crota2
                    ap.updateMarkers()
                    ima.changed = True
            if aper.type == 'Rectangle' or aper.type == 'Square' or aper.type == 'Pixel':
                x0,y0 = aper.rect.get_xy()
                w0    = aper.rect.get_width()
                h0    = aper.rect.get_height()
                angle = aper.rect.angle
                ra0,dec0 = ic.wcs.all_pix2world(x0, y0, 0)
                ws = w0 * ic.pixscale; hs = h0 * ic.pixscale
                for ima in ici:
                    x0,y0 = ima.wcs.all_world2pix(ra0, dec0, 0)
                    w0 = ws / ima.pixscale
                    h0 = hs / ima.pixscale
                    ap = ima.photApertures[ntab]
                    ap.rect.set_xy((x0,y0))
                    ap.rect.set_width(w0)
                    ap.rect.set_height(h0)
                    ap.rect.angle = angle  - ic.crota2 + ima.crota2
                    ap.updateMarkers()
                    ima.changed = True
            elif aper.type == 'Polygon':
                verts = aper.poly.get_xy()
                adverts = np.array([(ic.wcs.all_pix2world(x, y, 0)) for (x, y) in verts])                
                for ima in ici:
                    verts = [(ima.wcs.all_world2pix(ra, dec, 0)) for (ra, dec) in adverts]
                    ap = ima.photApertures[ntab]
                    ap.poly.set_xy(verts)
                    ap.updateMarkers()
                    ima.changed = True

    def onWheel(self,event):
        """Enable zoom with mouse wheel and propagate changes to other tabs."""
        eb = event.button
        self.zoomImage(eb)

    def zoomUp(self):
        self.zoomImage('up')
        
    def zoomDown(self):
        self.zoomImage('down')
        
    def zoomImage(self, eb):
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        curr_xlim = ic.axes.get_xlim()
        curr_ylim = ic.axes.get_ylim()
        curr_x0 = (curr_xlim[0]+curr_xlim[1])*0.5
        curr_y0 = (curr_ylim[0]+curr_ylim[1])*0.5
        if eb == 'up':
            factor=0.9
        elif eb == 'down':
            factor=1.1
        new_width = (curr_xlim[1]-curr_xlim[0])*factor*0.5
        new_height= (curr_ylim[1]-curr_ylim[0])*factor*0.5
        x = [curr_x0-new_width,curr_x0+new_width]
        y = [curr_y0-new_height,curr_y0+new_height]
        ic.axes.set_xlim(x)
        ic.axes.set_ylim(y)
        ic.fig.canvas.draw_idle()
        ici = self.ici.copy()
        ici.remove(ic)
        ra, dec = ic.wcs.all_pix2world(x, y, 0)
        for ima in ici:
            x,y = ima.wcs.all_world2pix(ra, dec, 0)
            ima.axes.set_xlim(x)
            ima.axes.set_ylim(y)
            ima.changed = True

    def onDraw2(self,event):
        stab = self.stabs.currentIndex()        
        sc = self.sci[stab]
        if sc.spectrum.redshift != self.specCube.redshift:
            flags = QMessageBox.Yes 
            flags |= QMessageBox.No
            question = "Do you want to update the redshift ?"
            response = QMessageBox.question(self, "Question",
                                                  question,
                                                  flags)            
            if response == QMessageBox.Yes:
                self.sb.showMessage("Updating the redshift ", 2000)
                self.specCube.redshift = sc.spectrum.redshift
                # Check if all has been already computed 
                if self.all == False:
                    try:
                        self.computeAll()
                    except BaseException:
                        pass
                # Propagate to other tabs
                for sc_ in self.sci:
                    # Check if tab is psf
                    if sc_.name != 'PSF':
                        sc_.spectrum.redshift = self.specCube.redshift
                        for annotation in sc_.annotations:
                            annotation.remove()
                        sc_.lannotation.remove()
                        sc_.zannotation.remove()
                        sc_.drawSpectrum()
                        sc_.fig.canvas.draw_idle()
                    else:
                        pass
            elif QMessageBox.No:
                self.sb.showMessage("Redshift value unchanged ", 2000)
                sc.spectrum.redshift = self.specCube.redshift
                for annotation in sc.annotations:
                    annotation.remove()
                sc.zannotation.remove()
                sc.lannotation.remove()
                sc.drawSpectrum()
                sc.fig.canvas.draw_idle()
            else:
                pass           
        if sc.spectrum.l0 != self.specCube.l0:
            flags = QMessageBox.Yes 
            flags |= QMessageBox.No
            question = "Do you want to update the reference wavelength ?"
            response = QMessageBox.question(self, "Question",
                                                  question,
                                                  flags)            
            if response == QMessageBox.Yes:
                self.sb.showMessage("Updating the reference wavelength ", 2000)
                self.specCube.l0 = sc.spectrum.l0
                # Check if all has been already computed 
                if self.all == False:
                    try:
                        self.computeAll()
                    except BaseException:
                        pass
                # Propagate to other tabs
                for sc_ in self.sci:
                    if sc_.name != 'PSF':
                        sc_.spectrum.l0 = self.specCube.l0
                        for annotation in sc_.annotations:
                            annotation.remove()
                        sc_.lannotation.remove()
                        sc_.zannotation.remove()
                        sc_.drawSpectrum()
                        sc_.fig.canvas.draw_idle()
                    else:
                        pass
            elif QMessageBox.No:
                self.sb.showMessage("Redshift value unchanged ", 2000)
                sc.spectrum.l0 = self.specCube.l0
                for annotation in sc.annotations:
                    annotation.remove()
                sc.zannotation.remove()
                sc.lannotation.remove()
                sc.drawSpectrum()
                sc.fig.canvas.draw_idle()
            else:
                pass                       
        # Deselect pan & zoom options on mouse release
        if sc.toolbar._active == "PAN":
            sc.toolbar.pan()
        if sc.toolbar._active == "ZOOM":
            sc.toolbar.zoom()

    def onKeyPress2(self, event):
        if event.key in ['control', 'cmd', 'shift', 'alt']:
            self.ctrlIsHeld = True
            
    def onPress2(self, event):
        if event.inaxes:
            self.press2 = event.xdata, event.ydata
        else:
            self.press2 = None

    def onKeyRelease2(self, event):
        if event.key in ['control', 'cmd', 'shift', 'alt']:
            self.ctrlIsHeld = False

    def onWheel2(self, event):
        """Wheel zooms/unzooms spectrum."""
        itab = self.stabs.currentIndex()
        sc = self.sci[itab]
        if event.inaxes:
            # zoom/unzoom 
            eb = event.button
            curr_xlim = sc.xlimits
            curr_ylim = sc.ylimits
            curr_x0 = (curr_xlim[0]+curr_xlim[1])*0.5
            curr_y0 = (curr_ylim[0]+curr_ylim[1])*0.5
            if eb == 'up':
                factor=0.9
            elif eb == 'down':
                factor=1.1
            new_width = (curr_xlim[1]-curr_xlim[0])*factor*0.5
            new_height= (curr_ylim[1]-curr_ylim[0])*factor*0.5
            sc.xlimits = (curr_x0-new_width,curr_x0+new_width)
            sc.updateXlim()
            sc.ylimits = (curr_y0-new_height,curr_y0+new_height)
            sc.updateYlim()
               
    def onMotion2(self, event):
        """Reacts to movements of the mouse."""
        if self.press2 is None:
            return
        if event.inaxes:
            dx = event.xdata - self.press2[0]
            dy = event.ydata - self.press2[1]
        else:
            return
        if event.button == 2:
            itab = self.stabs.currentIndex()
            sc = self.sci[itab]
            if np.abs(dx) > 0:
                sc.xlimits = (sc.xlimits[0] - dx, sc.xlimits[1] - dx)
                sc.updateXlim()
            if np.abs(dy) > 0:
                sc.ylimits = (sc.ylimits[0] - dy, sc.ylimits[1] - dy)
                sc.updateYlim()
        
    def createSpectralPanel(self):
        """Panel to plot spectra."""
        #self.spectralPanel = QGroupBox("")
        self.spectralPanel = QWidget()
        layout = QVBoxLayout(self.spectralPanel)
        # Toolbar
        self.createToolbar()
        # Tabs with plots        
        self.stabs = QTabWidget()
        self.stabs.setTabsClosable(True)
        self.stabs.tabCloseRequested.connect(self.removeSpecTab)
        self.stabs.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding)
        self.stabs.currentChanged.connect(self.onSTabChange)  # things to do when changing tab
        self.stabi = []
        self.sci  = []
        self.scid1 = []
        self.scid2 = []
        self.scid3 = []
        self.scid4 = []
        self.scid5 = []
        self.scid6 = []
        # Status bar
        self.sb = QStatusBar()
        self.sb.showMessage("Click the folder icon to load a cube !", 10000)
        # Add widgets to panel
        banner = QWidget()
        banner.layout = QHBoxLayout(banner)
        banner.layout.addWidget(self.tb)
        banner.layout.addWidget(self.sb)
        layout.addWidget(self.stabs)
        layout.addWidget(banner)

    def createToolbar(self):
        """Toolbar with main commands."""
        # Toolbar definition
        self.tb = QToolBar()
        self.tb.setMovable(True)
        self.tb.setObjectName('toolbar')
        # Actions
        self.helpAction = self.createAction(os.path.join(self.path0,'icons','help.png'),
                                            'Help','Ctrl+q',self.onHelp)
        self.issueAction = self.createAction(self.path0+'/icons/issue.png',
                                             'Report an issue','Ctrl+q',self.onIssue)
        self.quitAction = self.createAction(os.path.join(self.path0,'icons','exit.png'),
                                            'Quit program', 'Ctrl+q',self.fileQuit)
        self.startAction = self.createAction(os.path.join(self.path0,'icons','open.png'),
                                             'Load new observation', 'Ctrl+n',self.newFile)
        self.reloadAction = self.createAction(os.path.join(self.path0,'icons','reload.png'),
                                              'Reload observation', 'Ctrl+R',self.reloadFile)
        self.levelsAction = self.createAction(os.path.join(self.path0,'icons','levels.png'),
                                              'Adjust image levels',
                                              'Ctrl+L',self.changeVisibility)
        self.cmapAction = self.createAction(os.path.join(self.path0,'icons','rainbow.png'),
                                            'Choose color map', 'Ctrl+m',self.changeColorMap)
        self.blink = 'off'
        self.blinkAction = self.createAction(os.path.join(self.path0,'icons','blink.png'),
                                             'Blink between 2 images', 'Ctrl+B',self.blinkImages)
        self.momentAction = self.createAction(os.path.join(self.path0,'icons','map.png'),
                                              'Compute moment 0', 'Ctrl+m',self.zeroMoment)
        self.contours = 'off'
        self.contoursAction = self.createAction(os.path.join(self.path0,'icons','contours.png'),
                                                'Overlap contours', 'Ctrl+c',self.overlapContours)
        self.apertureAction = self.createApertureAction()
        self.fitAction = self.createFitAction()
        self.cutAction = self.createAction(os.path.join(self.path0,'icons','cut.png'),
                                           'Trim cube', 'Ctrl+k',self.trimCube)
        self.cropAction = self.createAction(os.path.join(self.path0,'icons','crop.png'),
                                            'Crop the cube', 'Ctrl+K',self.cropCube)
        self.sliceAction = self.createAction(os.path.join(self.path0,'icons','slice.png'),
                                             'Define a slice to compute moments and/or display',
                                             'Ctrl+K',self.sliceCube)
        self.slideAction = self.createAction(os.path.join(self.path0,'icons','slidecube.png'),
                                             'Display a slice of the cube and show a slider',
                                             'Ctrl+S',self.selectSlider)
        self.maskAction =  self.createAction(os.path.join(self.path0,'icons','eraser.png'),
                                             'Erase a region', '', self.maskCube)
        self.cloudAction = self.createAction(os.path.join(self.path0,'icons','cloud.png'),
                                             'Download image from cloud',
                                             'Ctrl+D', self.selectDownloadImage)
        self.fitsAction =  self.createAction(os.path.join(self.path0,'icons','download.png'),
                                             'Save the image as a FITS/PNG/JPG/PDF file',
                                             'Ctrl+S',self.saveFits)
        self.specAction = self.createAction(os.path.join(self.path0,'icons','download.png'),
                                            'Save the spectrum as a ASCII/FITS/PNG/JPG/PDF file',
                                            'Ctrl+S',self.saveSpectrum)
        self.vresizeAction = self.createAction(os.path.join(self.path0,'icons','vresize.png'),
                                               'Resize image vertically',
                                               'Ctrl+V',self.vresizeSpectrum)
        self.hresizeAction = self.createAction(os.path.join(self.path0,'icons','hresize.png'),
                                               'Resize image horizontally',
                                               'Ctrl+H',self.hresizeSpectrum)
        self.guessAction = self.createAction(os.path.join(self.path0,'icons','guessCont.png'),
                                             'Define cube fitting parameters',
                                             'Ctrl+g',self.guessContinuum)
        self.guessLinesAction = self.createAction(os.path.join(self.path0,'icons','guess.png'),
                                             'Define lines fitting parameters',
                                             'Ctrl+G',self.guessApLines)
        self.fitLinesAction = self.createAction(os.path.join(self.path0,'icons','fitline.png'),
                                             'Define lines fitting parameters',
                                             'Ctrl+F',self.fitApLines)
        self.fitContAction =self.createAction(os.path.join(self.path0,'icons','fit.png'),
                                              'Fit continuum/moments/lines',
                                              'Ctrl+F',self.ContMomLines)
        self.compMomAction = self.createAction(os.path.join(self.path0,'icons','computeMoments.png'),
                                               'Compute moments','Ctrl+g',self.chooseComputeMoments)
        self.fitregionAction = self.createAction(os.path.join(self.path0,'icons','fitregion.png'),
                                                 'Fit baseline and compute moments inside region',
                                                 'Ctrl+f',self.fitRegion)
        self.distanceAction = self.createAction(os.path.join(self.path0,'icons','squareset.png'),
                                                'Measure distance', '', self.computeDistance)
        self.psfAction = self.createAction(os.path.join(self.path0,'icons','psf.png'),
                                                'Estimate FWHM of PSF', '', self.estimatePSF)
        self.centroidAction = self.createAction(os.path.join(self.path0,'icons','centroid.png'),
                                                'Centroid of the PSF', '', self.centroidPSF)
        self.centerAction = self.createAction(os.path.join(self.path0,'icons','center.png'),
                                                'Fit center of the PSF', '', self.centerPSF)
        self.spacer = QWidget()
        self.spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.tb.addAction(self.startAction)
        self.tb.addAction(self.reloadAction)
        self.tb.addAction(self.helpAction)
        self.tb.addAction(self.issueAction)
        self.tb.addWidget(self.apertureAction)        
        #self.tb.addWidget(self.fitAction)        
        self.tb.addAction(self.quitAction)

    def openContinuumTab(self):        
        """Clear previous continuum estimate and open new tab."""
        s = self.specCube
        self.continuum = np.full((s.nz,s.ny,s.nx), np.nan) # Fit of continuum
        self.Cmask = np.zeros((s.nz,s.ny,s.nx), dtype=bool) # Spectral cube mask (for fitting the continuum)
        self.C0 = np.full((s.ny,s.nx), np.nan) # Continuum at ref. wavelength
        # Open tabs if they do not exist
        if 'C0' not in self.bands:
            self.addBand('C0')
        else:
            itab = self.bands.index('C0')
            ic = self.ici[itab]
            ic.showImage(self.C0)
            ic.fig.canvas.draw_idle()

    def kernel1pixel(self):
        # Set kernel to 1 pixel
        self.setKernel(1)

    def kernel5pixel(self):
        # Set kernel to 1 pixel
        self.setKernel(5)

    def kernel9pixel(self):
        # Set kernel to 1 pixel
        self.setKernel(9)

    def setKernel(self, size):
        """Set the kernel to compute the spectrum in the Pixel tab (and the continuum)."""
        # Change to pixel tab
        istab = self.spectra.index('Pix')
        if self.stabs.currentIndex() != istab:
            self.stabs.setCurrentIndex(istab)
        # Get the reference scale
        ic0 = self.ici[0]
        w0 = ic0.pixscale
        k1=False
        k5=False
        k9=False
        if size == 1:
            self.kernel = 1
            theta = 0
            k1 = True
        elif size == 5:
            self.kernel = 5
            theta = 45
            w0 *= np.sqrt(2.)*1.5
            k5 = True
        elif size == 9:
            self.kernel = 9
            w0 *= 2.6
            theta = 0.
            k9 = True
        else:
            self.kernel = 1
            theta = 0.
            k1 = True
        self.kernel1.setChecked(k1)
        self.kernel5.setChecked(k5)
        self.kernel9.setChecked(k9)
        # Update the pixel marker with new kernel
        for ic in self.ici:
            w = w0/ic.pixscale;
            pixel = ic.photApertures[0]
            x,y = pixel.rect.get_xy()
            theta_ = pixel.rect.angle
            w_  = pixel.rect.get_width()
            xc = x + w_/np.sqrt(2.) * np.sin((45.-theta_)*np.pi/180.)
            yc = y + w_/np.sqrt(2.) * np.cos((45.-theta_)*np.pi/180.)
            pixel.rect.set_width(w)
            pixel.rect.set_height(w)
            pixel.rect.angle = theta
            x = xc - w/np.sqrt(2.) * np.sin((45.-theta)*np.pi/180.)
            y = yc - w/np.sqrt(2.) * np.cos((45.-theta)*np.pi/180.)
            pixel.rect.set_xy((x,y))
            ic.changed = True
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        ic.fig.canvas.draw_idle()
        # Compute the new aperture flux
        self.onModifiedAperture(itab)

    def guessContinuum(self):
        """Create a first guess for fitting the continuum."""
        # Change to pixel tab
        istab = self.spectra.index('Pix')
        if self.stabs.currentIndex() != istab:
            self.stabs.setCurrentIndex(istab)
        sc = self.sci[istab]
        # Change to flux map
        imtab = self.bands.index('Flux')
        ic0 = self.ici[imtab]
        if self.itabs.currentIndex() != imtab:
            self.itabs.setCurrentIndex(imtab)
        w0 = ic0.pixscale
        # Dialog to select continuum fit paramenters
        self.CP = ContParams(self.kernel)
        if self.CP.exec_() == QDialog.Accepted:
            function, boundary, kernel, regions, emlines, ablines = self.CP.save()
            if function == 'Constant':
                self.zeroDeg = True
            else:
                self.zeroDeg = False
            if boundary == 'None':
                self.positiveContinuum = False
            else:
                self.positiveContinuum = True
            k1=k5=k9=False
            self.kernel = int(kernel)
            if kernel == 1:
                theta = 0.
                k1= True
            elif kernel == 5:
                w0 *= np.sqrt(2.)*1.5 
                theta = 45.
                k5 = True
            elif kernel == 9:
                w0 *= 2.6
                theta = 0.
                k9 = True
            else:
                self.kernel = 1
                theta = 0.
                k1 = True
            self.kernel1.setChecked(k1)
            self.kernel5.setChecked(k5)
            self.kernel9.setChecked(k9)
            self.ncells = int(regions)   # Number of Voronoi cells
            print('selected ', self.ncells, ' regions')
            sc.abslines = int(ablines)  # Number of absorption lines
            sc.emslines = int(emlines)  # Number of emission lines
            # Create tessellation
            nx = self.specCube.nx
            ny = self.specCube.ny
            if self.ncells == 1:
                self.sites = np.array([[nx // 2, ny // 2]])
            elif self.ncells == 4:
                dx = nx // 3
                dy = ny // 3
                self.sites = np.array([[dx, dy], [dx * 2.02, dy],
                                       [dx * 2, 2 * dy], [dx * 1.01, 2 * dy]
                                       ])
            else:
                self.computeSites(self.ncells)
            if sc.displayLines == True:
                sc.displayLines = False
                sc.setLinesVisibility(sc.displayLines)
                sc.fig.canvas.draw_idle()
        else:
            return
        # Update the pixel marker with new kernel
        for ic in self.ici:
            w = w0/ic.pixscale;
            pixel = ic.photApertures[0]
            x,y = pixel.rect.get_xy()
            theta_ = pixel.rect.angle
            w_  = pixel.rect.get_width()
            xc = x + w_/np.sqrt(2.) * np.sin((45.-theta_)*np.pi/180.)
            yc = y + w_/np.sqrt(2.) * np.cos((45.-theta_)*np.pi/180.)
            pixel.rect.set_width(w)
            pixel.rect.set_height(w)
            pixel.rect.angle = theta
            x = xc - w/np.sqrt(2.) * np.sin((45.-theta)*np.pi/180.)
            y = yc - w/np.sqrt(2.) * np.cos((45.-theta)*np.pi/180.)
            pixel.rect.set_xy((x,y))
            ic.changed = True
        # Help on status bar
        self.sb.showMessage("Click and drag the mouse over the spectrum to select two continuum regions ",
                            2000)        
        # Delete previous guess and select new one
        if sc.guess is not None:
            self.onRemoveContinuum('segments deleted')
        # I don't know how to disactivate the draggable feature, so I hide the annotations
        # when selecting the continuum
        try:
            sc.lguess = []
            sc.xguess = []
        except:
            pass
        self.CS = SegmentsSelector(sc.axes, sc.fig, self.onContinuumSelect, zD=self.zeroDeg)
        self.openContinuumTab()
        # Create Voronoi sites, KDTree, plot Voronoi ridges on image
        self.removeVI()
        if self.ncells > 1:
            self.createVI()
        else:
            ic0.fig.canvas.draw_idle()

    def guessApLines(self):
        """Create a guess of continuum and lines for an aperture."""
        # We should be able to mask part of spectrum (like in showspectra ... )
        self.GP = guessParams()
        if self.GP.exec_() == QDialog.Accepted:
            cont, em, ab = self.GP.save()
            if cont == 'Constant':
                self.zeroDeg = True
            else:
                self.zeroDeg = False
            try:
                self.onRemoveContinuum('all')
            except BaseException:
                pass
            istab = self.stabs.currentIndex()
            sc = self.sci[istab]
            sc.emslines = int(em)
            sc.abslines = int(ab)
            self.CS = SegmentsSelector(sc.axes, sc.fig, self.onContinuumSelect, zD=self.zeroDeg)
        else:
            return     
        
    def fitApLines(self):
        """Fit lines after defining guesses."""
        istab = self.stabs.currentIndex()
        sc = self.sci[istab]
        if sc.emslines + sc.abslines == 0:
            return
        else:
            print('I can fit lines')
 
    def removeVI(self):
        """Remove Voronoi tessellation if there."""
        try:
            self.VI.removeRidges()
            self.VI.disconnect
            self.VI.modSignal.disconnect()
            self.VI = None
        except BaseException:
            pass
        
    def createVI(self):
        """Create new Voronoi tessellation."""
        nx = self.specCube.nx
        ny = self.specCube.ny
        imtab = self.bands.index('Flux')
        ic0 = self.ici[imtab]
        tree = KDTree(self.sites)
        tq = tree.query(self.specCube.points)
        self.regions = tq[1].reshape(ny, nx)
        self.VI = VoronoiInteractor(ic0.axes, self.sites)
        self.VI.modSignal.connect(self.updateKDTree)
            
    def computeSites(self, n):
        nx = self.specCube.nx
        ny = self.specCube.ny
        area = nx * ny
        l = np.sqrt(area  * 2 / (3 * np.sqrt(3) * self.ncells))
        l2 = l * np.sqrt(3) / 2.
        nxc = int((nx - 2 * l)// (3 * l))
        nyc = int((ny - 2 * l2) // l2)
        # Second estimate from sides
        ly = ny / np.sqrt(3.) / (nyc / 2. + 1)
        lx = nx / (3  * nxc + 2)
        l = np.max(np.array([lx,ly]))
        l2 = l * np.sqrt(3) / 2.
        nxc = int(nx // (3 * l) + 1)
        nyc = int(ny // l2)
        self.sites = []
        for j in range(nyc):
            y = l2 * (j + 1)
            for i in range(nxc):
                if j % 2 == 1:
                    x = l * (3 * i )
                else:
                    x = l * (1.5 + (nxc - 1 - i) * 3)
                self.sites.append([x, y])   
        self.sites = np.array(self.sites)
        self.ncells = len(self.sites)
        print('Computed ', self.ncells, ' regions')

    def updateKDTree(self, event):
        """React to modification of Voronoi cells."""
        sc = self.sci[self.spectra.index('Pix')]
        nx = self.specCube.nx
        ny = self.specCube.ny
        self.sites = self.VI.sites
        tree = KDTree(self.sites)
        tq = tree.query(self.specCube.points)
        self.regions = tq[1].reshape(ny, nx)
        self.ncells = len(self.sites)
        if event == 'voronoi modified':
            print('number of cells is now ', self.ncells)
            pass
        elif event == 'one voronoi site added':
            try:
                # Add to the list the last value
                newguess = sc.xguess[-1]
                print('Adding guess ', newguess)
                sc.xguess.append(newguess)
                if len(sc.lines) > 0:
                    for line in sc.lines:
                        if line is None:
                            pass
                        else:
                            newguess = sc.lguess[line.n][-1]
                            sc.lguess[line.n].append(newguess)
            except BaseException:
                message = 'Please, define the continuum guess !'
                self.sb.showMessage(message, 4000)
                print(message)
        else:
            print('Site '+event+' removed')
            try:
                ind = int(event)
                del sc.xguess[ind]
                if len(sc.lines) > 0:
                    for line in sc.lines:
                        if line is None:
                            pass
                        else:
                            del sc.lguess[line.n][ind]
            except BaseException:
                message = 'Please, define the continuum guess !'
                self.sb.showMessage(message, 4000)
                print(message)

    def addBand(self, band):
        """Add a band and display it in a new tab."""
        # Open tab and display the image
        self.bands.append(band)
        t, ic, ih, h, c1, c2, c3, c4 = self.addImage(band)
        self.tabi.append(t)
        self.ici.append(ic)
        self.ihi.append(ih)
        self.ihcid.append(h)
        self.icid1.append(c1)
        self.icid2.append(c2)
        self.icid3.append(c3)
        self.icid4.append(c4)
        # image
        s = self.specCube
        if band == 'M0':
            image = self.M0
            if s.instrument == 'GREAT':
                image *= self.specCube.Tb2Jy
        elif band == 'M1':
            image = self.M1
        elif band == 'M2':
            image = self.M2
        elif band == 'M3':
            image = self.M3
        elif band == 'M4':
            image = self.M4
        elif band == 'C0':
            image = self.C0
        elif band == 'v':
            image = self.v
        elif band == 'sv':
            image = self.sv   
        elif band == 'L0':
            image = self.L0
        elif band == 'L1':
            image = self.L1
        if s.instrument == 'PCWI':
            aspect = s.ypixscale/s.pixscale
        else:
            aspect = 1.

        ic.compute_initial_figure(image=image, wcs=self.specCube.wcs, title=band,
                                  cMap=self.colorMap, cMapDir=self.colorMapDirection,
                                  stretch=self.stretchMap, aspect=aspect)
        if band in ['M0','M1','M2','M3','M4','C0','v','sv','L0','L1']:
            ic.crota2 = self.specCube.crota2
        ih = self.ihi[self.bands.index(band)]
        ih.compute_initial_figure(image=None)
        # Add apertures
        self.addApertures(ic)
        # Add contours
        self.addContours(ic) 
        # Align with flux cube
        ic0 = self.ici[0]
        x = ic0.axes.get_xlim()
        y = ic0.axes.get_ylim()
        #ra, dec = ic0.wcs.wcs_pix2world(x, y, 0)
        #x, y = ic.wcs.wcs_world2pix(ra, dec, 0)            
        ic.axes.set_xlim(x)
        ic.axes.set_ylim(y)
        ic.zoomlimits = (x,y)
        # Callback to propagate axes limit changes among images
        #ic.cid = ic.axes.callbacks.connect('xlim_changed' and 'ylim_changed', self.doZoomAll)
        ic.cid = ic.axes.callbacks.connect('ylim_changed', self.doZoomAll)

    def onContinuumSelect(self, verts):
        istab = self.stabs.currentIndex()
        sc = self.sci[istab]
        # Order the x coordinates of the verts
        x, y = zip(*verts)
        x = np.array(x); y = np.array(y)
        # Order increasing if wavelength, decreasing if frequency
        idx = np.argsort(x)
        if sc.xunit == 'THz':
            idx = idx[::-1]
        x = x[idx]
        y = y[idx]
        verts = [(i,j) for (i,j) in zip(x,y)]
        sc.guess = SegmentsInteractor(sc.axes, verts, self.zeroDeg)
        sc.guess.modSignal.connect(self.onModifiedGuess)
        sc.guess.mySignal.connect(self.onRemoveContinuum)
        interactors = [sc.guess]
        # Add lines
        if sc.emslines > 0:
            sc.lines = self.addLines(sc.emslines, x, 'emission')
        else:
            sc.lines = []
        if sc.abslines > 0:
            sc.lines.extend(self.addLines(sc.abslines, x, 'absorption', sc.emslines))
        interactors.extend(sc.lines)
        sc.interactorManager = InteractorManager(sc.axes, interactors)
        
        if istab == self.spectra.index('Pix'):
            if sc.emslines + sc.abslines > 0:
                sc.lguess = []
                for line in sc.lines:
                    sc.lguess.append([[line.x0, line.fwhm, line.A]] * self.ncells) 
            # Generate a list of limits connected to each Voronoi cell
            xg,yg = zip(*sc.guess.xy)
            print('x limits of guess ', xg)
            sc.xguess = [xg] * self.ncells
        else:
            if sc.emslines + sc.abslines > 0:
                sc.lguess = []
                for line in sc.lines:
                    sc.lguess.append([line.x0, line.fwhm, line.A]) 
        #
        if sc.displayLines == False:
            sc.displayLines = True
            sc.setLinesVisibility(sc.displayLines)
            sc.fig.canvas.draw_idle()
            
    def getLineVelocity(self):
        czline, okPressed = QInputDialog.getDouble(self, "Velocity of line ", "cz", 0, -10000., 50000., 2)
        if okPressed:
            return czline
        else:
            return None
            
    def addLines(self, n, x, type, nstart=0, x0s=None, fwhms=None, As=None):
        lines = []
        #istab = self.spectra.index('Pix')
        istab = self.stabs.currentIndex()
        sc = self.sci[istab]
        nid = nstart
        colors = np.array(['orange','forestgreen','dodgerblue','violet','salmon','peru'])
        for i in range(n):
            if fwhms is None:
                dx = (x[2] - x[1]) / (2 * n)
                fwhm = dx * 0.5
            else:
                fwhm = fwhms
            #x0 = x[1] + dx + i * 2 * dx
            #print('x0 precomputed is ', x0)
            # Ask for cz 
            if x0s is None:
                #czline = self.getLineVelocity()
                #w0 = self.specCube.l0
                #c =  299792.458 # km/s
                #x0 = w0 * (czline / c + self.specCube.redshift + 1)
                x0 = x[1] + dx + i * 2 * dx
            else:
                x0 = x0s * ( 1 + self.specCube.redshift)
            c0 = sc.guess.intcpt + x0 * sc.guess.slope  # continuum at line center
            if As is None:
                idx = (sc.x > (x0 - dx)) & (sc.x < (x0 + dx))
                if type == 'emission':
                    A = np.nanmax(sc.spectrum.flux[idx]) - c0
                else:
                    A = np.nanmin(sc.spectrum.flux[idx]) - c0
            else:
                A = As
            LI = LineInteractor(sc.axes, c0, sc.guess.slope, x0, A, fwhm, nid,color=colors[i])
            LI.modSignal.connect(self.onModifiedGuess)
            LI.mySignal.connect(self.onRemoveContinuum)
            lines.append(LI)
            nid += 1
        return lines

    def onModifiedGuess(self, event):
        """Pass modification to the xguess limits of the Voronoi cell."""
        istab = self.stabs.currentIndex()
        if istab == 1:
            if event == 'continuum guess modified':
                # Grab new values 
                sc = self.sci[istab]
                x, y = zip(*sc.guess.xy)
                # Identify cursor position of pixel-aperture on image
                aperture = self.ici[0].photApertures[0].aperture
                xc, yc = aperture.xy
                if self.ncells > 1:
                    ncell = self.regions[int(yc), int(xc)]
                else:
                    ncell = 0
                # Update the guess limits for cell
                sc.xguess[ncell] = x
            elif event[:-2] == 'line guess modified':
                # Grab new values 
                sc = self.sci[istab]
                # Identify cursor position of pixel-aperture on image
                aperture = self.ici[0].photApertures[0].aperture
                xc, yc = aperture.xy
                if self.ncells > 1:
                    ncell = self.regions[int(yc), int(xc)]
                else:
                    ncell = 0
                # N line
                nline = int(event[-2:])
                line = sc.lines[nline]
                sc.lguess[nline][ncell] = [line.x0, line.fwhm, line.A]
                sc.fig.canvas.draw_idle()

    def ContMomLines(self):
        """Dialog to select fit options for the cube."""
        #if self.continuum is not None:
        sc = self.sci[self.spectra.index('Pix')]
        # print('fitcont: ', self.fitcont)
        if sc.guess is not None:
            if self.fitcont:
                moments = True
                options = []
            else:
                moments = False
                options = []
            if sc.abslines + sc.emslines > 0:
                if self.fitcont:
                    lines = True
                else:
                    lines = False
            else:
                lines = False
            if self.ncells > 1:
                options.extend(['Set to medians', 'Set to lowest 25% medians', 
                                'Set to lowest 50% medians', 'Offset of lowest 25% uncorrected medians',
                                'Fit region', 'Fit all cube', 'Set to zero'])
            else:
                options.extend(['Set to medians', 'Set to lowest 25% medians', 
                                'Set to lowest 50% medians', 'Offset of lowest 25% uncorrected medians',
                                'Fit all cube', 'Set to zero'])
        else:
            options = ['Set to medians', 'Set to lowest 25% medians', 
                       'Set to lowest 50% medians','Offset of lowest 25% uncorrected medians','Set to zero']
            moments = False
            lines = False
        FCD = FitCubeDialog(options, moments, lines)
        if FCD.exec_() == QDialog.Accepted:
            coption, moption, loption = FCD.save()
            if coption is None:
                pass
            else:
                if coption == 'Fit all cube':
                    self.fitContAll()
                elif coption == 'Fit region':
                    self.fitContRegion()
                elif coption == 'Set to medians':
                    self.setContinuumMedian()
                elif coption == 'Set to lowest 25% medians':
                    self.setContinuumMedianPercent(25)
                elif coption == 'Set to lowest 50% medians':
                    self.setContinuumMedianPercent(50)
                elif coption == 'Offset of lowest 25% uncorrected medians':
                    self.setContinuumMedianPercent(33, True)
                elif coption == 'Set to zero':
                    self.setContinuumZero()
            if moption is None:
                pass
            else:
                if moption == 'Region':
                    print('Set continuum to median and compute region moments')
                    self.setContinuumMedian()
                    self.computeMomentsRegion()
                elif moption == 'All':
                    print('Compute all cube moments')
                    self.computeMomentsAll()
            if loption is None:
                pass
            else:
                if loption == 'Region':
                    self.fitLinesRegion()
                elif loption == 'All':
                    self.fitLinesAll()
        else:
            message = 'Define a guess for the continuum on the spectrum panel'
            self.sb.showMessage(message, 4000)

    def fitCont(self):
        """Options to fit the continuum."""
        if self.continuum is not None:
            if self.ncells > 1:
                options = ['Fit all cube', 'Fit region', 'Set to zero', 'Set to medians']
            else:
                options = ['Fit all cube', 'Set to zero', 'Set to medians']
        else:
            options = ['Set to zero', 'Set to medians']
        CFP = ContFitParams(options)
        if CFP.exec_() == QDialog.Accepted:

            option = CFP.save()
            if option == 'Fit all cube':
                self.fitContAll()
            elif option == 'Fit region':
                self.fitContRegion()
            elif option == 'Set to medians':
                self.setContinuumMedian()
            elif option == 'Set to zero':
                self.setContinuumZero()

            else:
                pass
        else:
            if self.continuum is None:
                message = 'Define a guess for the continuum on the spectrum panel'
                self.sb.showMessage(message, 4000)
            else:
                pass

    def setContinuumZero(self):
        """Set continuum to zero."""
        self.openContinuumTab()        

        nz,ny,nx = np.shape(self.specCube.flux)
        self.continuum = np.zeros((nz,ny,nx))
        self.C0 = np.zeros((ny,nx))
        self.Cs = np.zeros((ny,nx))
        self.refreshContinuum()
        self.fitcont = True

    def getContinuumGuess(self, ncell=None):
        sc = self.sci[self.spectra.index('Pix')]
        if sc.guess is None:
            return None, None, None, None  
        elif self.ncells == 1:
            # guess values for the continuum in wavelength
            xy = sc.guess.xy
            xg,yg = zip(*xy)
            xg = np.array(xg); yg = np.array(yg)
            if sc.xunit == 'THz':

                c = 299792458.0  # speed of light in m/s
                xg = c/xg * 1.e-6  # THz to um
        else:
            xg = sc.xguess[ncell]
            xg = np.array(xg)
            if sc.xunit == 'THz':
                c = 299792458.0  # speed of light in m/s
                xg = c/xg * 1.e-6  # THz to um
        # Compute i0,i1,i2,i3 from xy
        i0 = np.argmin(np.abs(self.specCube.wave-xg[0]))
        i1 = np.argmin(np.abs(self.specCube.wave-xg[1]))
        i2 = np.argmin(np.abs(self.specCube.wave-xg[2]))
        i3 = np.argmin(np.abs(self.specCube.wave-xg[3]))
        #print('cell ', ncell, ' xguess ', xg)
        #print('indeces ', i0, i1, i2, i3)
        return i0, i1, i2, i3

    def setContinuumMedian(self):
        """Compute continuum by using the median signal per pixel."""
        self.openContinuumTab()

        sc = self.sci[self.spectra.index('Pix')]
        #print('ncells ', self.ncells)
        if sc.guess is None:
            self.C0 = np.nanmedian(self.specCube.flux, axis=0)
        else:
            if self.ncells == 1:
                i0, i1, i2, i3 = self.getContinuumGuess()
                mask = np.zeros(len(self.specCube.flux), dtype=bool)
                mask[i0:i1] = True
                mask[i2:i3] = True
                self.C0 = np.nanmedian(self.specCube.flux[mask,:,:], axis=0)
            else:
                # Otherwise, find the regions
                for ncell in range(self.ncells):
                    i0, i1, i2, i3 = self.getContinuumGuess(ncell)
                    #print('cell ',ncell, 'is', i0,i1,i2,i3)
                    mask = np.zeros(len(self.specCube.flux), dtype=bool)
                    mask[i0:i1] = True
                    mask[i2:i3] = True
                    j, i = np.where(self.regions == ncell)
                    self.C0[j, i] = np.nanmedian((self.specCube.flux[:, j, i])[mask,:], axis=0)
        self.continuum = np.broadcast_to(self.C0, np.shape(self.specCube.flux))
        self.Cs = self.C0.copy() * 0. # Set all slopes to 0
        self.refreshContinuum()
        self.fitcont = True
        
    def setContinuumMedianPercent(self, percent, uncorrected=False):
        """Compute continuum by using the median signal per pixel."""
        self.openContinuumTab()

        if self.specCube.instrument != 'FIFI-LS':
            uncorrected = False

        sc = self.sci[self.spectra.index('Pix')]
        n = 100 // percent 
        #print('ncells ', self.ncells)
        if sc.guess is None:
            if uncorrected:
                sflux = np.sort(self.specCube.uflux, axis=0)
            else:
                sflux = np.sort(self.specCube.flux, axis=0)
            print('shape of flux is ',np.shape(sflux))
            nz, ny, nx = np.shape(sflux)
            npc = nz // n
            self.C0 = np.nanmedian(sflux[0: npc, :, :], axis=0)
        else:
            if self.ncells == 1:
                i0, i1, i2, i3 = self.getContinuumGuess()
                npc = (i3-i0) // n
                if uncorrected:
                    sflux = np.sort(self.specCube.uflux[i0:i3, :, :], axis=0)
                else:
                    sflux = np.sort(self.specCube.flux[i0:i3, :, :], axis=0)
                #print('n10 ', n50, 'sflux ', np.shape(sflux))
                self.C0 = np.nanmedian(sflux[0:npc, :, :], axis=0)
            else:
                # Otherwise, find the regions
                for ncell in range(self.ncells):
                    i0, i1, i2, i3 = self.getContinuumGuess(ncell)
                    j, i = np.where(self.regions == ncell)
                    npc = (i3-i0) // n 
                    if uncorrected:
                        sflux = np.sort(self.specCube.uflux[i0:i3,j,i], axis=0)
                    else:
                        sflux = np.sort(self.specCube.flux[i0:i3,j,i], axis=0)
                    self.C0[j, i] = np.nanmedian(sflux[0:npc,:], axis=0)
        # If uncorrected, put to zero positive offsets
        if uncorrected:
            idx = self.C0 > 0
            self.C0[idx] = 0
        self.continuum = np.broadcast_to(self.C0, np.shape(self.specCube.flux))
        self.Cs = self.C0.copy() * 0. # Set all slopes to 0
        self.refreshContinuum()
        self.fitcont = True
    

    def refreshContinuum(self):
        """Refresh the plotted image of the continuuum."""
        itab = self.bands.index('C0')
        ic = self.ici[itab]
        ic.showImage(image=self.C0)
        ic.image.format_cursor_data = lambda z: "{:.2e} Jy".format(float(z))        
        ih = self.ihi[itab]
        ih.compute_initial_figure(image = self.C0)
        ih.update_figure(image = self.C0)
        self.addContours(ic)
        # Update limits of image
        ic.image.set_clim(ih.limits)
        ic.changed = True
        # Update continuum on pixel tab
        sc = self.sci[self.spectra.index('Pix')]
        ic = self.ici[0]
        xc,yc = ic.photApertures[0].xy[0]  # marker coordinates (center of rectangle)
        i = int(np.rint(xc)); j = int(np.rint(yc))
        if self.specCube.instrument == 'GREAT':
            t2j = self.specCube.Tb2Jy
        else:
            t2j = 1.
        sc.updateSpectrum(cont=self.continuum[:,j,i]*t2j)
        sc.fig.canvas.draw_idle()

    def chooseComputeMoments(self):
        """Options to compute the moments."""
        if self.continuum is not None:
            cmask = np.isfinite(self.continuum)
            cmask0 = np.sum(cmask)
            if cmask0 > 0:
                sc = self.sci[self.spectra.index('Pix')]
                # Check if there is a region defined to compute moments
                if sc.regionlimits is not None:
                    # Dialog to choose between fitting the entire cube or only a region of it
                    msgBox = QMessageBox()
                    msgBox.setText('Compute the moments over:')
                    msgBox.addButton('All cube', QMessageBox.ActionRole)
                    msgBox.addButton('Region', QMessageBox.ActionRole)
                    msgBox.addButton('Cancel', QMessageBox.ActionRole)
                    self.result = msgBox.exec()
                    if self.result == 0:
                        self.computeMomentsAll()
                    elif self.result == 1:
                        self.computeRegion()
                    else:
                        pass
                else:
                    message = 'First define a region where to compute the moments'
                    self.sb.showMessage(message, 4000)
            else:
                message = 'First fit the continuum'
                self.sb.showMessage(message, 4000)                
        else:
            message = 'First define a guess for the continuum'
            self.sb.showMessage(message, 4000)
            
    def moveIndex(self, oldpos, newpos):
        
        bigList = [
                self.bands,
                self.tabi,
                self.ici,
                self.ihi,
                self.ihcid,
                self.icid1,
                self.icid2,
                self.icid3,
                self.icid4,
        ]

        for elist in bigList:
            item = elist[oldpos]
            elist.remove(item)
            elist.insert(newpos, item)

    def computeVelocities(self):
        """ Compute velocity and vel. dispersion from the moments """
        if self.M0 is not None:
            itab = self.itabs.currentIndex()
            origband = self.bands[itab]
            # Update tabs
            newbands = ['v','sv']
            sbands = [self.v,self.sv]
            for new, sb in zip(newbands, sbands):
                if new not in self.bands:
                    self.addBand(new)
                else:
                    itab = self.bands.index(new)
                    self.removeTab(itab, False)
                    self.addBand(new)
            # Compute velocity fields
            c = 299792.458 # km/s
            w0 = self.specCube.l0
            z = self.specCube.redshift
            self.v = (self.M1/w0 -1. - z)* c # km/s
            self.sv = np.sqrt(self.M2) * c/w0  # km/s
            # Refresh the plotted images
            bands = ['v', 'sv']
            sbands = [self.v, self.sv]
            for b,sb in zip(bands,sbands):
                itab = self.bands.index(b)
                ic = self.ici[itab]
                ic.showImage(image=sb)
                if (b == 'sv') & (self.specCube.instrument == 'FIFI-LS'):
                    R = self.specCube.getResolutionFIFI()
                    c = 299792.458 # km/s
                    sR = c/R/2.355
                    ic.fig.suptitle('$\sigma_v$ [km/s]  (R = '+str(int(R)) + 
                                    ", $\sigma_R$ = {:.1f} km/s)".format(sR))
                ic.image.format_cursor_data = lambda z: "{:.1f} km/s".format(float(z))
                ih = self.ihi[itab]
                ih.compute_initial_figure(image = sb)
                ih.update_figure(image = sb)
                ic.image.set_clim(ih.limits)
                # Adopt same limits as flux image
                ic0=self.ici[0]
                ic.image.axes.set_xlim(ic0.image.axes.get_xlim())
                ic.image.axes.set_ylim(ic0.image.axes.get_ylim())
                ic.changed = True
            # Reorder tabs
            for new in ['v','sv']:  
                itab = self.bands.index(new)
                itab0 = self.bands.index('M0')
                if new == 'v':
                    newpos = itab0 + 1
                else:
                    newpos = itab0 + 2
                self.itabs.tabBar().moveTab(itab, newpos) # move v, sv after M0
                # Move elements in associated lists
                self.moveIndex(itab, newpos)
            # Go back to original tab
            itab = self.bands.index(origband)
            self.itabs.setCurrentIndex(itab)
        else:
            message = 'First compute the moments'
            self.sb.showMessage(message, 4000)

    def fitRegion(self):
        """Fit the guess over a square region."""
        # 0) Check if guess is defined
        sc = self.sci[1]
        if sc.guess is None:
            self.sb.showMessage("Please, define a guess for the continuum", 4000)
        else:
            # 0) Toggle pixel marker
            itab = self.itabs.currentIndex()
            ic = self.ici[itab]
            pixel = ic.photApertures[0]
            pixel.showverts = False
            pixel.line.set_visible(pixel.showverts)
            # 1) Define region
            self.RS = RectangleSelector(ic.axes, self.fitContRegion,
                                        drawtype='box', useblit=False,
                                        button=[1, 3],  # don't use middle button
                                        minspanx=5, minspany=5,
                                        spancoords='pixels',
                                        rectprops = dict(facecolor='g', edgecolor = 'g',
                                                         alpha=0.8, fill=False),
                                        lineprops = dict(color='g', linestyle='-',linewidth = 2,
                                                         alpha=0.8),
                                        interactive=False)
            self.RS.to_draw.set_visible(False)
            self.RS.set_visible(True)
            #self.RS.state.add('center')

    def computeRegion(self):
        """Fit the guess over a square region."""
        # 0) Check if guess is defined
        sc = self.sci[1]
        if sc.regionlimits is None:
            self.sb.showMessage("Please, define the region to compute the moments ", 4000)
        else:
            # 0) Toggle pixel marker
            itab = self.itabs.currentIndex()
            ic = self.ici[itab]
            pixel = ic.photApertures[0]
            pixel.showverts = False
            pixel.line.set_visible(pixel.showverts)
            # 1) Define region
            self.RS = RectangleSelector(ic.axes, self.computeMomentsRegion,
                                        drawtype='box', useblit=False,
                                        button=[1, 3],  # don't use middle button
                                        minspanx=5, minspany=5,
                                        spancoords='pixels',
                                        rectprops = dict(facecolor='g', edgecolor = 'g',
                                                         alpha=0.8, fill=False),
                                        lineprops = dict(color='g', linestyle='-',
                                                         linewidth = 2, alpha=0.8),
                                        interactive=False)
            self.RS.to_draw.set_visible(False)
            self.RS.set_visible(True)

    def fitContRegionOld(self, eclick, erelease):
        """Fit the continuum in a selected region."""
        # Select points in the region
        points = self.onFitRect(eclick,erelease)
        # Create the mask for the continuum
        self.continuumMask(points)
        # Execute the fit
        self.fitContinuum(points)

    def onFitRect(self, eclick, erelease):
        """eclick and erelease are the press and release events."""        
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
        # if tab is not flux, uflux, exp, or m0 recompute the values
        itab = self.itabs.currentIndex()
        band = self.bands[itab]
        if (band == 'Flux') | (band == 'uFlux') | (band =='Exp') | (band == 'M0'):
            pass
        else:
            # Transform into original cube coordinates
            ic = self.ici[itab]
            ic0 = self.ici[0]
            ra1,dec1 = ic.wcs.all_pix2world(x1, y1, 0)
            ra2,dec2 = ic.wcs.all_pix2world(x2, y2, 0)
            x1,y1 = ic0.wcs.all_world2pix(ra1, dec1, 0)
            x2,y2 = ic0.wcs.all_world2pix(ra2, dec2, 0)
        # Disactive the selector
        self.disactiveSelectors()
        # Untoggle pixel marker
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        pixel = ic.photApertures[0]
        pixel.showverts = True
        pixel.line.set_visible(pixel.showverts)
        # Find the points  (create a meshgrid of the selected rectangle)
        x1 = int(np.rint(x1))
        x2 = int(np.rint(x2))
        y1 = int(np.rint(y1))
        y2 = int(np.rint(y2))
        xi = np.arange(x1,x2)
        yi = np.arange(y1,y2)
        xi,yi = np.meshgrid(xi,yi)
        points = np.array([np.ravel(xi),np.ravel(yi)]).transpose()
        return points
    
    def fitContRegion(self):
        """Fit continuum inside region occupied by the cursor."""
        # position of the cursor
        if self.ncells > 1:
            aperture = self.ici[0].photApertures[0].aperture
            xc, yc = aperture.xy
            ncell = self.regions[int(yc), int(xc)]
            exp = np.nansum(self.specCube.exposure, axis=0)
            mask = (self.regions == ncell) & (exp > 0)
            mask[0, :] = False
            mask[-1, :] = False
            mask[:, 0] = False
            mask[:, -1] = False
            yi, xi = np.where(mask == True)
            points = np.c_[xi, yi]
            # Create the mask for the continuum
            self.continuumMask(points)
            # Execute the fit
            self.fitContinuum(points)
        else:
            print('There are no defined regions')
            return

    def fitContAll(self):
        """Fit continuum all over the cube."""
        # For efficiency, points with no exposure should be masked
        exp = np.nansum(self.specCube.exposure, axis=0)
        mask = exp > 0
        mask[0, :] = False
        mask[-1, :] = False
        mask[:, 0] = False
        mask[:, -1] = False
        yi, xi = np.where(mask == True)
        points = np.c_[xi, yi]
        # Continuum mask
        self.continuumMask(points, True)
        # Fit
        self.fitContinuum(points)
        # Set flag
        self.fitcont = True

    def continuumMask(self, points, smooth=False):
        from scipy.signal import convolve2d as convolve
        # Update masks
        if self.ncells <= 1:
            i0, i1, i2, i3 = self.getContinuumGuess()
            for p in points:
                i, j = p
                self.Cmask[:,j,i] = 0
                self.Cmask[i0:i1,j,i] = 1
                self.Cmask[i2:i3,j,i] = 1
        else:
            # Do some smoothing before computing the Cmask
            if smooth:
                # A) Define the limits for each region
                nx = self.specCube.nx
                ny = self.specCube.ny
                i0map = np.ones((ny,nx))
                i1map = np.ones((ny,nx))
                i2map = np.ones((ny,nx))
                i3map = np.ones((ny,nx))
                for cell in range(self.ncells):
                    i0, i1, i2, i3 = self.getContinuumGuess(cell)
                    idx = np.where(self.regions == cell)
                    i0map[idx] = i0
                    i1map[idx] = i1
                    i2map[idx] = i2
                    i3map[idx] = i3 
                # B) Convolve each limit map with a kernel
                kernel = np.array([[1/16., 1/8., 1/16.], 
                                   [1/8., 1/4., 1/8.],
                                   [1/16.,1/8.,1/16.]])
                i0map = convolve(i0map, kernel, boundary='symm', mode='same')
                i1map = convolve(i1map, kernel, boundary='symm', mode='same')
                i2map = convolve(i2map, kernel, boundary='symm', mode='same')
                i3map = convolve(i3map, kernel, boundary='symm', mode='same')
                # C) Define the continuum mask
                for p in points:
                    i, j = p
                    i0 = int(i0map[j, i])
                    i1 = int(i1map[j, i])
                    i2 = int(i2map[j, i])
                    i3 = int(i3map[j, i])                   
                    self.Cmask[:,j,i] = 0
                    self.Cmask[i0:i1,j,i] = 1
                    self.Cmask[i2:i3,j,i] = 1                
            else:
                for p in points:
                    i, j = p
                    ncell = self.regions[j, i]
                    i0, i1, i2, i3 = self.getContinuumGuess(ncell)
                    self.Cmask[:,j,i] = 0
                    self.Cmask[i0:i1,j,i] = 1
                    self.Cmask[i2:i3,j,i] = 1
        print('Continuum mask computed')
        # Eventually apply smoothing on the borders of the xguess map to avoid sudden change of continuum

    def fitContinuum(self, points):
        """Fit the continuum on a selected set of points."""
        sc = self.sci[self.spectra.index('Pix')]
        intcp = sc.guess.intcpt
        slope = sc.guess.slope
        c, c0, cs = multiFitContinuum(self.Cmask, self.specCube.wave, self.specCube.flux,
                                  self.continuum, self.C0, self.specCube.l0,
                                  points, slope, intcp, self.positiveContinuum,
                                  self.kernel, exp=self.specCube.exposure)
        self.continuum = c
        print('max continuum is ',np.nanmax(self.continuum))
        self.C0 = c0
        self.Cs = cs
        # Refresh the plotted image
        itab = self.bands.index('C0')
        ic = self.ici[itab]
        if self.specCube.instrument == 'GREAT':
            t2j = self.specCube.Tb2Jy
        else:
            t2j = 1.
        ic.showImage(image=self.C0*t2j)
        ic.image.format_cursor_data = lambda z: "{:.2e} Jy".format(float(z))        
        ih = self.ihi[itab]
        ih.compute_initial_figure(image=self.C0*t2j)
        ih.update_figure(image=self.C0*t2j)
        self.addContours(ic)
        # Update limits of image
        ic.image.set_clim(ih.limits)
        # Adopt same limits as flux image
        ic0 = self.ici[0]
        ic.image.axes.set_xlim(ic0.image.axes.get_xlim())
        ic.image.axes.set_ylim(ic0.image.axes.get_ylim())
        # Differ the change
        ic.changed = True
        # Update continuum on pixel tab
        sc = self.sci[self.spectra.index('Pix')]
        ic = self.ici[0]
        xc,yc = ic.photApertures[0].xy[0]  # marker coordinates (center of rectangle)
        i = int(np.rint(xc)); j = int(np.rint(yc))
        sc.updateSpectrum(cont=self.continuum[:,j,i]*t2j, cslope = self.Cs[j,i]*t2j)
        sc.fig.canvas.draw_idle()

    def computeMomentsAll(self):
        """Compute moments all over the cube."""
        # Remove contours
        try:
            self.removeContours()
        except:
            pass
        # Define region (whole image)
        nx = self.specCube.nx
        ny = self.specCube.ny
        xi = np.arange(nx)
        yi = np.arange(ny)
        xi,yi = np.meshgrid(xi,yi)
        points = np.array([np.ravel(xi),np.ravel(yi)]).transpose()
        # Define moments
        self.defineMoments()
        # Compute moment mask
        self.momentsMask(points, True)
        # Compute moments
        self.computeMoments(points)


    def defineMoments(self):
        """Define mask, moments, and velocities."""
        # In the case they are not already defined ...
        if self.M0 is None:
            s = self.specCube
            self.Mmask = np.zeros((s.nz,s.ny,s.nx), dtype=bool) # Spectral cube mask (for computing the moments)
            self.M0 = np.full((s.ny,s.nx), np.nan) # 0th moment
            self.M1 = np.full((s.ny,s.nx), np.nan) # 1st moment
            self.M2 = np.full((s.ny,s.nx), np.nan) # 2nd moment
            self.M3 = np.full((s.ny,s.nx), np.nan) # 3rd moment
            self.M4 = np.full((s.ny,s.nx), np.nan) # 4th moment
            self.v = np.full((s.ny,s.nx), np.nan) # velocity field
            self.sv = np.full((s.ny,s.nx), np.nan) # vel. disp. field
            # Current canvas
            itab = self.itabs.currentIndex()
            ic0 = self.ici[itab]
            # Create/update moment tabs
            #newbands = ['M0','M1','M2','M3','M4']
            #sbands = [self.M0, self.M1, self.M2, self.M3,self.M4]
            newbands = ['M0', 'M3', 'M4']
            if self.specCube.instrument == 'GREAT':
                t2j = self.specCube.Tb2Jy
            else:
                t2j = 1.
            sbands = [self.M0*t2j, self.M3, self.M4]
            for new,sb in zip(newbands,sbands):
                if new not in self.bands:
                    self.addBand(new)
                else:
                    itab = self.bands.index(new)
                    ic = self.ici[itab]
                    ic.showImage(sb)
                    if ic == ic0:
                        ic.fig.canvas.draw_idle()
                    else:
                        ic.changed = True
                        
    def defineLines(self):
        """Define line structure."""
        if self.L0 is None:
            s = self.specCube
            self.Mmask = np.zeros((s.nz,s.ny,s.nx), dtype=bool) # Spectral cube mask (for computing the moments)
            self.L0 = np.full((s.ny,s.nx), np.nan) #  1st line integral
            self.L1 = np.full((s.ny,s.nx), np.nan) #  2nd line integral
            sc = self.sci[self.spectra.index('Pix')]
            nlines = len(sc.lguess)
            self.lines = np.full((nlines, 4, s.ny, s.nx), np.nan) # Fitted parameters of lines
            # Current canvas
            itab = self.itabs.currentIndex()
            ic0 = self.ici[itab]
            # Create/update moment tabs
            if nlines == 1:
                newbands = ['L0']
                sbands = [self.L0]                
            elif nlines == 2:
                newbands = ['L0','L1']
                sbands = [self.L0, self.L1]
            for new,sb in zip(newbands,sbands):
                if new not in self.bands:
                    self.addBand(new)
                else:
                    itab = self.bands.index(new)
                    ic = self.ici[itab]
                    ic.showImage(sb)
                    if ic == ic0:
                        ic.fig.canvas.draw_idle()
                    else:
                        ic.changed = True            

    def computeMomentsRegionOld(self, eclick, erelease):
        """Compute moments in a defined region."""
        # Select points in the region
        points = self.onFitRect(eclick,erelease)
        # Define moments
        self.defineMoments()
        # Compute moments
        self.computeMoments(points)
        
    def computeMomentsRegion(self):
        """Fit continuum inside region occupied by the cursor."""
        # position of the cursor
        if self.ncells > 1:
            # Remove contours
            try:
                self.removeContours()
            except:
                pass

            aperture = self.ici[0].photApertures[0].aperture
            xc, yc = aperture.xy
            ncell = self.regions[int(yc), int(xc)]
            print('region is ', ncell)
            exp = np.nansum(self.specCube.exposure, axis=0)
            mask = (self.regions == ncell) & (exp > 0)
            mask[0, :] = False
            mask[-1, :] = False
            mask[:, 0] = False
            mask[:, -1] = False
            print('region has ',np.sum(mask),' points')
            yi, xi = np.where(mask == True)
            points = np.c_[xi, yi]
            #print('Points defined ... for computing moments, points: ', np.shape(points))
            # Define moments
            self.defineMoments()
            # Compute moment mask
            self.momentsMask(points)
            # Compute moments inside the region
            self.computeMoments(points)
        else:
            print('There are no defined regions')
            return
        
    def fitLinesAll(self):
        """Fit defined lines in the cell occupied by the cursor."""
            
        for ncell in range(self.ncells):
            exp = np.nansum(self.specCube.exposure, axis=0)
            mask = (self.regions == ncell) & (exp > 0)
            mask[0, :] = False
            mask[-1, :] = False
            mask[:, 0] = False
            mask[:, -1] = False
            print('Cell ', ncell, ' has ',np.sum(mask),' points')
            if np.sum(mask) > 0:
                yi, xi = np.where(mask == True)
                points = np.c_[xi, yi]
                print('Points defined ... for computing moments, points: ', np.shape(points))
                # Define line fit results
                self.defineLines()
                # Compute the mask
                self.momentsMask(points)
                # Fit the lines inside the cell
                self.fitLinesOnly(points, ncell)  
        self.fitLinesDisplay()
        
        
    def fitLinesRegion(self):
        """Fit defined lines in the cell occupied by the cursor."""
        if self.ncells > 1:
            aperture = self.ici[0].photApertures[0].aperture
            xc, yc = aperture.xy
            ncell = self.regions[int(yc), int(xc)]
            print('region is ', ncell)
            exp = np.nansum(self.specCube.exposure, axis=0)
            mask = (self.regions == ncell) & (exp > 0)
            mask[0, :] = False
            mask[-1, :] = False
            mask[:, 0] = False
            mask[:, -1] = False
            print('region has ',np.sum(mask),' points')
            if np.sum(mask) > 0:
                yi, xi = np.where(mask == True)
                points = np.c_[xi, yi]
                print('Points defined ... for computing moments, points: ', np.shape(points))
                # Define line fit results
                self.defineLines()
                # Compute the mask
                self.momentsMask(points)
                # Fit the lines inside the cell
                self.fitLines(points)
        else:
            print('There are no defined regions')
            return
           

    def momentsMask(self, points, smooth=False):
        from scipy.signal import convolve2d as convolve
        # Update masks
        if self.ncells < 1:
            i0, i1, i2, i3 = self.getContinuumGuess()
            for p in points:
                i, j = p
                self.Mmask[:,j,i] = 0
                #self.Mmask[i1:i2,j,i] = 1
                self.Mmask[i0:i3,j,i] = 1
        else:
            # Do some smoothing before computing the Cmask
            if smooth:
                # A) Define the limits for each region
                nx = self.specCube.nx
                ny = self.specCube.ny
                i1map = np.ones((ny,nx))
                i2map = np.ones((ny,nx))
                for cell in range(self.ncells):
                    i0, i1, i2, i3 = self.getContinuumGuess(cell)
                    idx = np.where(self.regions == cell)
                    #i1map[idx] = i1
                    #i2map[idx] = i2
                    i1map[idx] = i0
                    i2map[idx] = i3
                # B) Convolve each limit map with a kernel
                kernel = np.array([[1/16., 1/8., 1/16.], 
                                   [1/8., 1/4., 1/8.],
                                   [1/16.,1/8.,1/16.]])
                i1map = convolve(i1map, kernel, boundary='symm', mode='same')
                i2map = convolve(i2map, kernel, boundary='symm', mode='same')
                # C) Define the continuum mask
                for p in points:
                    i, j = p
                    i1 = int(i1map[j, i])
                    i2 = int(i2map[j, i])
                    self.Mmask[:,j,i] = 0
                    self.Mmask[i1:i2,j,i] = 1
            else:
                for p in points:
                    i, j = p
                    ncell = self.regions[j, i]
                    i0, i1, i2, i3 = self.getContinuumGuess(ncell)
                    self.Mmask[:,j,i] = 0
                    self.Mmask[i1:i2,j,i] = 1
        print('Moment mask computed')
        # Eventually apply smoothing on the borders of the xguess map to avoid sudden change of continuum

    def computeMoments(self, points):
        """Compute moments and velocities."""
        m = self.Mmask
        moments = [self.M0, self.M1, self.M2, self.M3, self.M4]
        f = self.specCube.flux
        w = self.specCube.wave
        c = self.continuum
        moments, self.noise = multiComputeMoments(m, w, f, c, moments,points)
        self.M0, self.M1, self.M2, self.M3, self.M4 = moments
        # Refresh the plotted images
        #bands = ['M0', 'M1', 'M2', 'M3', 'M4']
        #sbands = [self.M0, self.M1, self.M2, self.M3, self.M4]
        bands = ['M0', 'M3', 'M4']
        if self.specCube.instrument == 'GREAT':
            t2j = self.specCube.Tb2Jy
        else:
            t2j = 1.
        sbands = [self.M0*t2j, self.M3, self.M4]
        for b,sb in zip(bands,sbands):
            itab = self.bands.index(b)
            ic = self.ici[itab]
            ic.showImage(image=sb)
            if  b == 'M1':
                ic.image.format_cursor_data = lambda z: "{:.4f} um".format(float(z))
            elif  b == 'M2':
                ic.image.format_cursor_data = lambda z: "{:.4f} um2".format(float(z))
            elif  b == 'M0':
                ic.image.format_cursor_data = lambda z: "{:.2e} W/m2".format(float(z))
            else:
                ic.image.format_cursor_data = lambda z: "{:.2e} ".format(float(z))
            ih = self.ihi[itab]
            ih.compute_initial_figure(image=sb)
            ih.update_figure(image=sb)
            self.addContours(ic)
            ic.image.set_clim(ih.limits)
            # Adopt same limits as flux image
            ic0 = self.ici[0]
            ic.image.axes.set_xlim(ic0.image.axes.get_xlim())
            ic.image.axes.set_ylim(ic0.image.axes.get_ylim())
            ic.changed = True
        # Refresh current image (if a moment)
        itab = self.itabs.currentIndex()
        if self.bands[itab] in bands:
            ic = self.ici[itab]
            ic.fig.canvas.draw_idle()
        # Update moments on pixel tab
        sc = self.sci[self.spectra.index('Pix')]
        ic = self.ici[0]
        xc,yc = ic.photApertures[0].xy[0]
        i = int(np.rint(xc)); j = int(np.rint(yc))
        if self.specCube.instrument == 'GREAT':
            t2j = self.specCube.Tb2Jy
        else:
            t2j = 1
        moments = [self.M0[j, i] * t2j, self.M1[j, i], self.M2[j, i], self.M3[j, i], self.M4[j, i]]
        sc.updateSpectrum(cont=self.continuum[:, j, i]* t2j, cslope=self.Cs[j,i]* t2j,
                          moments=moments, noise=self.noise[j, i]* t2j)
        sc.fig.canvas.draw_idle()
        # Compute Velocities
        self.computeVelocities()
        
    def fitLines(self, points):
        """Fit lines inside a defined region."""
        m = self.Mmask
        # Find cell and guesses
        sc = self.sci[self.spectra.index('Pix')]
        # Select cell          
        aperture = self.ici[0].photApertures[0].aperture
        xc, yc = aperture.xy
        ncell = self.regions[int(yc), int(xc)]
        lineguesses = []
        for guess in sc.lguess:
            g = guess[ncell]
            lineguesses.append(g)
        f = self.specCube.flux
        w = self.specCube.wave
        c = self.continuum
        multiFitLines(m, w, f, c, lineguesses, self.lines, points)
        print('Number of lines ',len(self.lines))
        # Update L0 and L1 (first two lines)
        self.L0 = self.lines[0][2] # the plane no 2 corresponds to the amplitude
        if len(self.lines) == 2:
            self.L1 = self.lines[1][2]
        # Then display them
        if len(self.lines) == 2:
            bands = ['L0', 'L1']
            sbands = [self.L0, self.L1]
        else:
            bands = ['L0']
            sbands = [self.L0]
        for b, sb in zip(bands, sbands):
            itab = self.bands.index(b)
            ic = self.ici[itab]
            ic.showImage(image=sb)
            ic.image.format_cursor_data = lambda z: "{:.2e} W/m2".format(float(z))
            ih = self.ihi[itab]
            ih.compute_initial_figure(image=sb)
            ih.update_figure(image=sb)
            self.addContours(ic)
            ic.image.set_clim(ih.limits)
            # Adopt same limits as flux image
            ic0 = self.ici[0]
            ic.image.axes.set_xlim(ic0.image.axes.get_xlim())
            ic.image.axes.set_ylim(ic0.image.axes.get_ylim())
            ic.changed = True
        # Refresh current image (if a line)
        itab = self.itabs.currentIndex()
        if self.bands[itab] in bands:
            ic = self.ici[itab]
            ic.fig.canvas.draw_idle()
        # Update lines on pixel tab
        # sc = self.sci[self.spectra.index('Pix')]
        ic = self.ici[0]
        xc,yc = ic.photApertures[0].xy[0]
        i = int(np.rint(xc)); j = int(np.rint(yc))
        lines = []
        for line in self.lines:
            lines.append([line[0][j, i], line[1][j, i], line[2][j, i], line[3][j, i]])
        sc.updateSpectrum(cont=self.continuum[:, j, i], cslope=self.Cs[j,i], lines=lines)
        sc.fig.canvas.draw_idle()
        
    def fitLinesOnly(self, points, ncell):
        """Fit lines inside a defined region."""
        m = self.Mmask
        # Find cell and guesses
        sc = self.sci[self.spectra.index('Pix')]
        lineguesses = []
        for guess in sc.lguess:
            g = guess[ncell]
            lineguesses.append(g)
        f = self.specCube.flux
        w = self.specCube.wave
        c = self.continuum
        multiFitLines(m, w, f, c, lineguesses, self.lines, points)
        # Update L0 and L1 (first two lines)
        self.L0 = self.lines[0][2] # the plane no 2 corresponds to the amplitude
        if len(self.lines) == 2:
            self.L1 = self.lines[1][2]
        
    def fitLinesDisplay(self):
        # Display images after fitting
        if len(self.lines) == 2:
            bands = ['L0', 'L1']
            sbands = [self.L0, self.L1]
        else:
            bands = ['L0']
            sbands = [self.L0]
        for b, sb in zip(bands, sbands):
            itab = self.bands.index(b)
            ic = self.ici[itab]
            ic.showImage(image=sb)
            ic.image.format_cursor_data = lambda z: "{:.2e} W/m2".format(float(z))
            ih = self.ihi[itab]
            ih.compute_initial_figure(image=sb)
            ih.update_figure(image=sb)
            self.addContours(ic)
            ic.image.set_clim(ih.limits)
            # Adopt same limits as flux image
            ic0 = self.ici[0]
            ic.image.axes.set_xlim(ic0.image.axes.get_xlim())
            ic.image.axes.set_ylim(ic0.image.axes.get_ylim())
            ic.changed = True
        # Refresh current image (if a line)
        itab = self.itabs.currentIndex()
        if self.bands[itab] in bands:
            ic = self.ici[itab]
            ic.fig.canvas.draw_idle()
        # Update lines on pixel tab
        ic = self.ici[0]
        xc,yc = ic.photApertures[0].xy[0]
        i = int(np.rint(xc)); j = int(np.rint(yc))
        lines = []
        for line in self.lines:
            lines.append([line[0][j, i], line[1][j, i], line[2][j, i], line[3][j, i]])
        sc = self.sci[self.spectra.index('Pix')]
        sc.updateSpectrum(cont=self.continuum[:, j, i], cslope=self.Cs[j,i],lines=lines)
        sc.fig.canvas.draw_idle()
 


    def createAction(self,icon,text,shortcut,action):
        act = QAction(QIcon(icon), text, self)
        act.setShortcut(shortcut)
        act.triggered.connect(action)
        return act

    def createApertureAction(self):
        """Create combo box for choosing an aperture."""
        self.apertures = [['apertures','Square','Rectangle'],
                     ['Circle','Ellipse','Polygon']]
        self.model = QStandardItemModel()
        for d in self.apertures:                
            row = []
            for text in d:
                item = QStandardItem(QIcon(os.path.join(self.path0,'icons',text.lower()+'.png')),"")
                item.setTextAlignment(Qt.AlignCenter)
                if text != 'apertures':
                    item.setToolTip("Choose a "+text)
                else:
                    item.setToolTip("Choose an aperture ")
                #item.setCheckable(True)
                #item.setCheckState(False)
                row.append(item)
            self.model.appendRow(row)
        self.apView = QTableView()
        # Remove headers and grid
        self.apView.verticalHeader().setVisible(False) 
        self.apView.horizontalHeader().setVisible(False)
        self.apView.setShowGrid(False)
        self.apView.setIconSize(QSize(20,20))
        apertureAction = QComboBox()
        apertureAction.setToolTip("Choose an aperture\n")
        #apertureAction.SizeAdjustPolicy(QComboBox.AdjustToContentsOnFirstShow)
        apertureAction.SizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
        apertureAction.setIconSize(QSize(20, 20))
        apertureAction.setView(self.apView)
        apertureAction.setModel(self.model)
        self.apView.setModel(self.model)
        self.apView.resizeColumnsToContents()  # Once defined the model, resize the column width
        self.apView.setSelectionMode(QAbstractItemView.SingleSelection)
        self.apView.setSelectionBehavior(QAbstractItemView.SelectItems)
        #sp = self.apView.sizePolicy()
        #sp.setHorizontalPolicy(QSizePolicy.MinimumExpanding)
        self.apView.setMinimumWidth(120)  
        self.apView.setMinimumHeight(70)  
        #self.apView.setSizePolicy(sp)
        apertureAction.activated.connect(self.chooseAperture)
        return apertureAction

    def createFitAction(self):
        """Create combo box for choosing an aperture."""
        self.fitoptions = [['gaussfit','continuum','list'],
                           ['location','dispersion','check']]
        self.fmodel = QStandardItemModel()
        for d in self.fitoptions:                
            row = []
            for text in d:
                item = QStandardItem(QIcon(os.path.join(self.path0,'icons',text+'.png')),"")
                item.setTextAlignment(Qt.AlignCenter)
                if text == 'gaussfit':
                    item.setToolTip("Choose an option")
                elif text == 'list':
                    item.setToolTip("Select manually all the constrains")
                elif text == 'check':
                    item.setToolTip("Fit all over the cube")
                else:
                    item.setToolTip("Constrain the "+text)
                row.append(item)
            self.fmodel.appendRow(row)
        self.fitView = QTableView()
        # Remove headers and grid
        self.fitView.verticalHeader().setVisible(False) 
        self.fitView.horizontalHeader().setVisible(False)
        self.fitView.setShowGrid(False)
        self.fitView.setIconSize(QSize(20, 20))
        fitAction = QComboBox()
        fitAction.setToolTip("Fit line and continuum (under construction)\n")
        #apertureAction.SizeAdjustPolicy(QComboBox.AdjustToContentsOnFirstShow)
        fitAction.SizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
        fitAction.setIconSize(QSize(24,24))
        fitAction.setView(self.fitView)
        fitAction.setModel(self.fmodel)
        self.fitView.setModel(self.fmodel)
        self.fitView.resizeColumnsToContents()  # Once defined the model, resize the column width
        self.fitView.setSelectionMode(QAbstractItemView.SingleSelection)
        self.fitView.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.fitView.setMinimumWidth(120)  
        self.fitView.setMinimumHeight(70)  
        fitAction.activated.connect(self.chooseFitOption)
        return fitAction

    def onPolySelect(self, verts):
        self.disactiveSelectors()
        # 1 vertices in RA,Dec coords
        itab = self.itabs.currentIndex()
        ic0 = self.ici[itab]
        adverts = np.array([(ic0.wcs.all_pix2world(x, y, 0)) for (x,y) in verts])
        self.drawNewPolygonAperture(adverts)
        
    def drawNewPolygonAperture(self, adverts):
        """Save aperture with vertices in ra,dec coordinates."""
        n = len(self.photoApertures)
        self.photoApertures.append(photoAperture(n,'polygon',adverts))
        for ic in self.ici:
            # First adjust vertices to astrometry (they are in xy coords)
            verts = [(ic.wcs.all_world2pix(ra, dec, 0)) for (ra,dec) in adverts]
            poly  = PolygonInteractor(ic.axes, verts)
            ic.photApertures.append(poly)
            cidap=poly.mySignal.connect(self.onRemoveAperture)
            ic.photApertureSignal.append(cidap)
            poly.modSignal.connect(self.onModifiedAperture)
        self.drawNewSpectrum(n)

    def drawNewSpectrum(self, n):        
        """Add tab with the flux inside the aperture."""
        apname = "{:d}".format(n)
        self.spectra.append(apname)
        t, sc, scid1, scid2, scid3, scid4, scid5, scid6 = self.addSpectrum(apname)
        self.stabi.append(t)
        self.sci.append(sc)
        self.scid1.append(scid1)
        self.scid2.append(scid2)
        self.scid3.append(scid3)
        self.scid4.append(scid4)
        self.scid5.append(scid5)
        self.scid6.append(scid6)
        # Draw spectrum from polygon
        aperture = self.ici[0].photApertures[n].aperture
        path = aperture.get_path()
        transform = aperture.get_patch_transform()
        npath = transform.transform_path(path)
        s = self.specCube
        inpoints = s.points[npath.contains_points(s.points)]
        xx,yy = inpoints.T        
        fluxAll = np.nansum(s.flux[:,yy,xx], axis=1)
        if s.instrument == 'GREAT':
            spec = Spectrum(s.wave, fluxAll*s.Tb2Jy, instrument=s.instrument, 
                            redshift=s.redshift, l0=s.l0, Tb2Jy=s.Tb2Jy, bunit=s.bunit)
        elif s.instrument in ['HI','VLA','ALMA','MUSE','IRAM','CARMA','PCWI']:
            spec = Spectrum(s.wave, fluxAll, instrument=s.instrument, redshift=s.redshift, l0=s.l0)
        elif s.instrument == 'FIFI-LS':
            ufluxAll = np.nansum(s.uflux[:,yy,xx], axis=1)
            expAll = np.nansum(s.exposure[:,yy,xx], axis=1)
            efluxAll = np.sqrt(np.nansum(s.eflux[:,yy,xx]**2, axis=1))
            spec = Spectrum(s.wave, fluxAll, eflux=efluxAll, uflux=ufluxAll,
                            exposure=expAll, atran = s.atran, instrument=s.instrument,
                            redshift=s.redshift, baryshift=s.baryshift, l0=s.l0, 
                            watran=s.watran, uatran=s.uatran)
        elif s.instrument in ['PACS','FORCAST']:
            expAll = np.nansum(s.exposure[:,yy,xx], axis=1)
            spec = Spectrum(s.wave, fluxAll, exposure=expAll, instrument=s.instrument,
                            redshift=s.redshift, l0=s.l0 )
        # Inherit the x-units of pix 
        istab = self.spectra.index('Pix')
        sc.xunit = self.sci[istab].xunit
        sc.compute_initial_spectrum(name=apname, spectrum=spec)
        self.specZoomlimits = [sc.xlimits,sc.ylimits]
        sc.cid = sc.axes.callbacks.connect('xlim_changed' and 'ylim_changed', self.doZoomSpec)
        # Start the span selector to show only part of the cube
        sc.span = SpanSelector(sc.axes, self.onSelect, 'horizontal', useblit=True,
                               rectprops=dict(alpha=0.3, facecolor='LightGreen'))
        sc.span.active = False
        # Select new tab
        self.stabs.setCurrentIndex(len(self.stabs)-1)
                   
    def onRectSelect(self, eclick, erelease):
        """eclick and erelease are the press and release events."""       
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
        if self.selAp == 'Square' or self.selAp == 'Rectangle':
            x0=x1;y0=y1
            w  = np.abs(x2-x1)
            h  = np.abs(y2-y1)
        else:
            x0 = (x1+x2)*0.5
            y0 = (y1+y2)*0.5
            w  = np.abs(x2-x1)
            h  = np.abs(y2-y1)
        self.newSelectedAperture(x0,y0,w,h,self.selAp)

    def newSelectedAperture(self, x0, y0, w, h, selAp):       
        itab = self.itabs.currentIndex()
        ic0 = self.ici[itab]
        r0, d0 = ic0.wcs.all_pix2world(x0, y0, 0)
        ws = w * ic0.pixscale
        hs = h * ic0.pixscale
        self.drawNewAperture(selAp, r0, d0, ws, hs, 0.)

    def drawNewAperture(self, selAp, r0, d0, ws, hs, angle):
        """Draw new selected aperture."""
        n = len(self.photoApertures)
        itab = self.itabs.currentIndex()
        ic0 = self.ici[itab]
        if selAp == 'Square':
            self.disactiveSelectors()
            # Define square
            data = [r0, d0, ws]
            self.photoApertures.append(photoAperture(n, 'square', data))
            for ic in self.ici:
                x0, y0 = ic.wcs.all_world2pix(r0, d0, 0)
                w = ws / ic.pixscale
                h = hs / ic.pixscale
                square = RectangleInteractor(ic.axes, (x0, y0), w, 
                                             angle=angle - ic0.crota2 + ic.crota2)
                ic.photApertures.append(square)
                cidap = square.mySignal.connect(self.onRemoveAperture)
                ic.photApertureSignal.append(cidap)
                square.modSignal.connect(self.onModifiedAperture)
        elif selAp == 'Rectangle':
            self.disactiveSelectors()
            # Define rectangle
            data = [r0, d0, ws, hs]
            self.photoApertures.append(photoAperture(n, 'rectangle', data))
            for ic in self.ici:
                x0, y0 = ic.wcs.all_world2pix(r0, d0, 0)
                w = ws / ic.pixscale
                h = hs / ic.pixscale
                rectangle = RectangleInteractor(ic.axes, (x0, y0), w, h, 
                                                angle=angle - ic0.crota2 + ic.crota2)
                ic.photApertures.append(rectangle)
                cidap=rectangle.mySignal.connect(self.onRemoveAperture)
                ic.photApertureSignal.append(cidap)
                #cidapm=
                rectangle.modSignal.connect(self.onModifiedAperture)
        elif selAp == 'Circle':
            self.disactiveSelectors()
            # Define circle
            data = [r0, d0, ws]
            self.photoApertures.append(photoAperture(n, 'circle', data))
            for ic in self.ici:
                x0, y0 = ic.wcs.all_world2pix(r0, d0, 0)
                w = ws / ic.pixscale
                h = hs / ic.pixscale
                circle = EllipseInteractor(ic.axes, (x0, y0), w)
                ic.photApertures.append(circle)
                cidap=circle.mySignal.connect(self.onRemoveAperture)
                ic.photApertureSignal.append(cidap)
                circle.modSignal.connect(self.onModifiedAperture)
        elif selAp == 'Ellipse':
            self.disactiveSelectors()
            # Define ellipse
            data = [r0, d0, ws, hs]
            self.photoApertures.append(photoAperture(n, 'ellipse', data))
            for ic in self.ici:
                x0, y0 = ic.wcs.all_world2pix(r0, d0, 0)
                w = ws / ic.pixscale
                h = hs / ic.pixscale
                ellipse = EllipseInteractor(ic.axes, (x0, y0), w, h,
                                            angle=angle - ic0.crota2 + ic.crota2)
                ic.photApertures.append(ellipse)
                cidap=ellipse.mySignal.connect(self.onRemoveAperture)
                ic.photApertureSignal.append(cidap)
                ellipse.modSignal.connect(self.onModifiedAperture)
        self.drawNewSpectrum(n)

    def disactiveSelectors(self):
        """Disactive all selectors, in the case more than one is selected."""
        if self.RS is not None:
            self.RS.set_active(False)
            for artist in self.RS.artists:
                artist.remove()
            self.RS = None
        if self.ES is not None:
            self.ES.set_active(False)
            for artist in self.ES.artists:
                artist.remove()
            self.ES = None
        if self.PS is not None:
            self.PS.set_active(False)
            for artist in self.PS.artists:
                artist.remove()
            self.PS = None
        if self.LS is not None:
            self.LS.set_active(False)
            for artist in self.LS.artists:
                artist.remove()
            self.LS = None

    def chooseFitOption(self, i):
        """Choosing a fit option."""
        index  = self.fitView.selectionModel().currentIndex()
        i = index.row()
        j = index.column()
        self.selFit = self.fitoptions[i][j]
        if self.selFit == 'list':
            self.sb.showMessage("You chose to open the "+self.selFit, 1000)
        if self.selFit == 'check':
            self.sb.showMessage("You chose to compute the fit for all the pixels ", 1000)
        elif self.selFit == 'gaussfit':
            self.sb.showMessage("Choose a fitting option ", 1000)
        else:
            self.sb.showMessage("You chose the action "+self.selFit, 1000)
        #put back to the 0-th item
        self.fitAction.setCurrentIndex(0)

    def chooseAperture(self, i):
        """Choosing an aperture."""
        index  = self.apView.selectionModel().currentIndex()
        i = index.row()
        j = index.column()
        self.selAp = self.apertures[i][j]
        if self.selAp == 'Ellipse':
            self.sb.showMessage("You chose an "+self.selAp, 1000)
        elif self.selAp == 'apertures':
            self.sb.showMessage("Choose an aperture shape ", 1000)
        else:
            self.sb.showMessage("You chose a "+self.selAp, 1000)
        if len(self.ici) == 0:
            self.sb.showMessage("Start by opening a new image ", 1000)
            self.apertureAction.setCurrentIndex(0)
            return
        self.activateAperture()

    def selectSquareAperture(self):
        self.selAp = 'Square'
        self.sb.showMessage("You chose a "+self.selAp, 1000)
        self.activateAperture()
        
    def selectRectangleAperture(self):
        self.selAp = 'Rectangle'
        self.sb.showMessage("You chose a "+self.selAp, 1000)
        self.activateAperture()

    def selectCircleAperture(self):
        self.selAp = 'Circle'
        self.sb.showMessage("You chose a "+self.selAp, 1000)
        self.activateAperture()
        
    def selectEllipseAperture(self):
        self.selAp = 'Ellipse'
        self.sb.showMessage("You chose an "+self.selAp, 1000)
        self.activateAperture()
        
    def selectPolygonAperture(self):
        self.selAp = 'Polygon'
        self.sb.showMessage("You chose a "+self.selAp, 1000)
        self.activateAperture()

    def activateAperture(self):
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        print('aperture is ',self.selAp)
        # Deactivate current aperture
        n = self.nAper()
        if n >= 0:
            print('deactivate aperture ', n)
            ap = ic.photApertures[n]
            ap.showverts = False
            ap.line.set_visible(ap.showverts)
        if self.selAp == 'Polygon':
            self.PS = PolygonSelector(ic.axes, self.onPolySelect,
                                      lineprops=dict(linestyle='-',color='g'),
                                      useblit=True,markerprops=dict(marker='o',mec='g'),
                                      vertex_select_radius=15)
        elif self.selAp == 'Rectangle':
            self.RS = RectangleSelector(ic.axes, self.onRectSelect,
                                        drawtype='box', useblit=True,
                                        button=[1, 3],  # don't use middle button
                                        minspanx=5, minspany=5,
                                        spancoords='pixels',
                                        rectprops = dict(facecolor='g', edgecolor = 'g',
                                                         alpha=0.8, fill=False),
                                        lineprops = dict(color='g', linestyle='-',linewidth = 2,
                                                         alpha=0.8),
                                        interactive=False)
            self.RS.state.add('center')
        elif self.selAp == 'Square':
            self.RS = RectangleSelector(ic.axes, self.onRectSelect,
                                        drawtype='box', useblit=True,
                                        button=[1, 3],  # don't use middle button
                                        minspanx=5, minspany=5,
                                        spancoords='pixels',
                                        rectprops = dict(facecolor='g', edgecolor = 'g',
                                                         alpha=0.8, fill=False),
                                        lineprops = dict(color='g', linestyle='-',linewidth = 2,
                                                         alpha=0.8),
                                        interactive=False)
            self.RS.state.add('square')
            self.RS.state.add('center')
        elif self.selAp == 'Ellipse':
            self.ES = EllipseSelector(ic.axes, self.onRectSelect,
                                      drawtype='line', useblit=True,
                                      button=[1, 3],  # don't use middle button
                                      minspanx=5, minspany=5,
                                      spancoords='pixels',
                                      rectprops = dict(facecolor='g', edgecolor = 'g',
                                                       alpha=0.8, fill=False),
                                      lineprops = dict(color='g', linestyle='-',linewidth = 2,
                                                       alpha=0.8),
                                      interactive=False)
            self.ES.state.add('center')
        elif self.selAp == 'Circle':
            self.ES = EllipseSelector(ic.axes, self.onRectSelect,
                                      drawtype='line', useblit=True,
                                      button=[1, 3],  # don't use middle button
                                      minspanx=5, minspany=5,
                                      spancoords='pixels',
                                      rectprops = dict(facecolor='g', edgecolor = 'g',
                                                       alpha=0.8, fill=False),
                                      lineprops = dict(color='g', linestyle='-',linewidth = 2,
                                                       alpha=0.8),
                                      interactive=False)
            self.ES.state.add('center')        
            self.ES.state.add('square')        
        if self.selAp != 'apertures':
            ic.fig.canvas.draw_idle()
        #put back to the 0-th item
        self.apertureAction.setCurrentIndex(0)

    def selectDownloadImage(self):
        """Select the image to download."""
        selectDI = QInputDialog()
        selectDI.setStyleSheet("* { font-size: 14pt; }")
        selectDI.setOption(QInputDialog.UseListViewForComboBoxItems)
        selectDI.setWindowTitle("Select image to download")
        selectDI.setLabelText("Selection")
        imagelist = ['local image', 'local spectral cube',
                     'sdss-u','sdss-g','sdss-r','sdss-i','sdss-z',
                     'panstarrs-g','panstarrs-r','panstarrs-i','panstarrs-z','panstarrs-y',
                     '2mass-j','2mass-h','2mass-k',
                     'wise1','wise2','wise3','wise4',
                     'first','nvss','sumss']
        selectDI.setComboBoxItems(imagelist)
        select = selectDI.exec_()
        if select == QDialog.Accepted:
            self.downloadImage(selectDI.textValue())

    def downloadImage(self, band):
        """Download an image covering the cube."""
        # Compute center and size of image (in arcmin)
        nz,ny,nx = np.shape(self.specCube.flux)
        lon, lat = self.specCube.wcs.celestial.all_pix2world(ny//2, nx//2, 0)
        xsize = nx * self.specCube.pixscale / 60. #size in arcmin
        ysize = ny * self.specCube.pixscale / 60. #size in arcmin
        # Compute center and size (arcmin) of the displayed image 
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        x = ic.axes.get_xlim()
        y = ic.axes.get_ylim()
        ra,dec = ic.wcs.all_pix2world(x, y, 0)
        lon = np.mean(ra)
        lat = np.mean(dec)
        xsize = np.abs(ra[0] - ra[1]) * np.cos(lat * np.pi / 180.) * 60.
        ysize = np.abs(dec[0] - dec[1]) * 60.
        if band == 'local image':
            # Download the local fits
            downloadedImage = cloudImage(lon, lat, xsize, ysize, band)
            if downloadedImage.data is not None:
                self.newImageTab(downloadedImage)
                message = 'New image downloaded'
            else:
                message = 'The selected survey does not cover the displayed image'
            self.newImageMessage(message)
        elif band == 'local spectral cube':
            print('Uploading auxiliary spectral cube')
            downloadedCube = cloudImage(lon, lat, xsize, ysize, band)
            self.auxSpecCube = downloadedCube.data
            if self.auxSpecCube is not None:
                message = 'Auxiliary spectral cube downloaded'
                # Add the spectrum to the Pix tab
                istab = self.spectra.index('Pix')
                sc = self.sci[istab]
                sc.displayAuxFlux = True
                sc.auxiliary = True
                sc.auxl0 = self.auxSpecCube.l0
                sc.auxw = self.auxSpecCube.wave
                #sc.drawSpectrum()
                self.onModifiedAperture('Auxiliary spectrum')  # Update the spectrum plot
                sc.spectrum.l0 = self.specCube.l0
                sc.drawSpectrum()
                sc.fig.canvas.draw_idle()
            else:
                message = 'The selected survey does not cover the displayed image'
            self.newImageMessage(message)
        else:
            # Here call the thread
            self.downloadThread = DownloadThread(lon, lat, xsize, ysize, band, parent=self)
            #self.downloadThread.updateTabs.newImage.connect(self.newImage)
            self.downloadThread.updateTabs.connect(self.newImage)
            self.downloadThread.sendMessage.connect(self.newImageMessage)
            self.downloadThread.start()
            # and start the spinning messagebox
            self.msgbox = QMessageBox()
            label = QLabel(self.msgbox)
            pixmap = QPixmap(os.path.join(self.path0, 'icons', 'niet.png'))
            label.setPixmap(pixmap)
            movie = QMovie(os.path.join(self.path0, 'icons', 'spinplane.gif'))
            label.setMovie(movie)
            movie.jumpToFrame(0)
            movie.start()
            label.resize(QSize(200,200))
            self.msgbox.setIconPixmap(pixmap)
            self.msgbox.setText("Querying " + band + " ... ")
            self.msgbox.exec_()
       
    def newImageMessage(self, message):
        """Message sent from download thread.""" 
        print('new message called')
        self.sb.showMessage(message, 5000)
        try:
            self.msgbox.done(1)
        except:
            pass

    def newImage(self, downloadedImage):
        "Save and display."
        print('Save and display called')
        print('image ', np.shape(downloadedImage.data))
        self.saveDownloadedFits(downloadedImage)
        self.newImageTab(downloadedImage)

    def newImageTab(self, downloadedImage):
        """Open  a tab and display the new image."""
        image = downloadedImage.data
        if np.issubdtype(image[0,0], int):
            print('Image contains integers')
        else:
            pass
            # Put to NaN values of mode if this is more common than 10%
            #ima = image.ravel()
            #(_, idx, counts) = np.unique(ima, return_index=True, return_counts=True)
            #index = idx[np.argmax(counts)]
            #mode = ima[index]
            #m0 = image == mode
            #n0 = np.sum(m0)
            #nx, ny = np.shape(image)
            #p0 = n0 / (nx * ny)
            #print('Mode of image is ', mode)
            #image = np.array(image, dtype=float)
            #if (p0 > 0.1):
            #    image[m0] = np.nan
        mask = np.isfinite(image)
        if np.sum(mask) == 0:
            self.sb.showMessage("The selected survey does not cover the displayed image", 2000)
        else:
            self.sb.showMessage("Image downloaded", 2000)
            band = downloadedImage.source
            self.bands.append(band)
            t,ic,ih,h,c1,c2,c3,c4 = self.addImage(band)
            self.tabi.append(t)
            self.ici.append(ic)
            self.ihi.append(ih)
            self.ihcid.append(h)
            self.icid1.append(c1)
            self.icid2.append(c2)
            self.icid3.append(c3)
            self.icid4.append(c4)           
            image = downloadedImage.data
            wcs = downloadedImage.wcs
            ic.compute_initial_figure(image=image,wcs=wcs,title=band,
                                      cMap=self.colorMap,cMapDir=self.colorMapDirection,
                                      stretch=self.stretchMap)
            try:
                ic.crota2 = downloadedImage.crota2
            except:
                ic.crota2 = 0.
            ih.compute_initial_figure(image=image)
            ih.update_figure(image=image)
            ic.image.set_clim(ih.limits)
            ic.changed = True
            # Check rotation angle
            # Align with spectral cube
            ic0 = self.ici[0]
            x = ic0.axes.get_xlim()
            y = ic0.axes.get_ylim()
            xmin = np.min(x)
            xmax = np.max(x)
            ymin = np.min(y)
            ymax = np.max(y)
            x = (xmin, xmax, xmin, xmax)
            y = (ymin, ymin, ymax, ymax)
            ra,dec = ic0.wcs.all_pix2world(x, y, 0)
            x,y = ic.wcs.all_world2pix(ra, dec, 0)            
            xmin = np.min(x)
            xmax = np.max(x)
            ymin = np.min(y)
            ymax = np.max(y)
            ic.axes.set_xlim(xmin, xmax)
            ic.axes.set_ylim(ymin, ymax)
            ic.zoomlimits = [(xmin,xmax),(ymin,ymax)]
            ic.changed = True
            # Add existing apertures
            self.addApertures(ic)
            # Add existing contours
            self.addContours(ic)
            # Callback to propagate axes limit changes to other images
            # ic.cid = ic.axes.callbacks.connect('xlim_changed' and 'ylim_changed', self.doZoomAll)
            ic.cid = ic.axes.callbacks.connect( 'ylim_changed', self.doZoomAll)

    def saveDownloadedFits(self, downloadedImage):
        """Save the downloaded FITS image."""
        from astropy.io import fits        
        # Dialog to save file
        fd = QFileDialog()
        fd.setLabelText(QFileDialog.Accept, "Save as")
        fd.setNameFilters(["Fits Files (*.fits)","All Files (*)"])
        fd.setOptions(QFileDialog.DontUseNativeDialog)
        fd.setViewMode(QFileDialog.List)
        fd.selectFile(downloadedImage.source+'.fits')        
        if (fd.exec()):
            fileName = fd.selectedFiles()
            outfile = fileName[0]
            # Check the 
            filename, file_extension = os.path.splitext(outfile)
            basename = os.path.basename(filename)            
            # Primary header
            image = downloadedImage.data
            wcs   = downloadedImage.wcs
            header = wcs.to_header()
            header.remove('WCSAXES')
            header['OBJECT'] = (self.specCube.objname, 'Object Name')
            header['INSTRUME'] = (basename, 'Instrument')
            hdu = fits.PrimaryHDU(image)
            hdu.header.extend(header)
            hdul = fits.HDUList([hdu])
            hdul.writeto(outfile,overwrite=True) # clobber true  allows rewriting
            hdul.close()
            
    def uploadSpectrum(self, event):
        """Upload existing spectrum."""        
        fd = QFileDialog()
        fd.setLabelText(QFileDialog.Accept, "Import")
        fd.setNameFilters(["Fits Files (*.fits, *.fits.gz)","WXY fits files (*WXY*.fits*)", "All Files (*)"])
        fd.setOptions(QFileDialog.DontUseNativeDialog)
        fd.setViewMode(QFileDialog.List)
        fd.setFileMode(QFileDialog.ExistingFile)
        if (fd.exec()):
            fileName= fd.selectedFiles()
            # Read external spectrum
            self.extSpectrum = ExtSpectrum(fileName[0])            
            # Plot over selected tab
            istab = self.stabs.currentIndex()
            sc = self.sci[istab]
            sc.extspecLayer, = sc.axes.plot(self.extSpectrum.wave,self.extSpectrum.flux, color='orange')
            sc.displayExtSpec = True

    def addApertures(self, ic):
        """ Add apertures already defined on new image """
        ic0 = self.ici[0]
        for aper in ic0.photApertures:
            #apertype = aper.__class__.__name__
            if aper.type == 'Ellipse' or aper.type == 'Circle':
                x0,y0 = aper.ellipse.center
                w0    = aper.ellipse.width
                h0    = aper.ellipse.height
                angle = aper.ellipse.angle - ic0.crota2 + ic.crota2
                ra0,dec0 = ic0.wcs.all_pix2world(x0, y0, 0)
                ws = w0 * ic0.pixscale; hs = h0 * ic0.pixscale
                # Add ellipse
                x0,y0 = ic.wcs.all_world2pix(ra0, dec0, 0)
                w0 = ws/ic.pixscale; h0 = hs/ic.pixscale
                ellipse = EllipseInteractor(ic.axes, (x0,y0),w0,h0,angle)
                ellipse.type = aper.type
                ellipse.showverts = aper.showverts
                ic.photApertures.append(ellipse)
                cidap=ellipse.mySignal.connect(self.onRemoveAperture)
                ic.photApertureSignal.append(cidap)
                ellipse.modSignal.connect(self.onModifiedAperture)
            elif aper.type == 'Rectangle' or aper.type == 'Square':
                x0,y0 = aper.rect.get_xy()
                w0    = aper.rect.get_width()
                h0    = aper.rect.get_height()
                #print(type(h0))
                angle = aper.rect.angle - ic0.crota2 + ic.crota2
                ra0,dec0 = ic0.wcs.all_pix2world(x0, y0, 0)
                ws = w0 * ic0.pixscale; hs = h0 * ic0.pixscale
                # Add rectangle
                x0,y0 = ic.wcs.all_world2pix(ra0, dec0, 0)
                w0 = ws/ic.pixscale; h0 = hs/ic.pixscale
                rectangle = RectangleInteractor(ic.axes, (x0,y0),w0,h0,angle)
                rectangle.type = aper.type
                rectangle.showverts = aper.showverts
                #rectangle.rect.set_xy((x0,y0))
                rectangle.updateMarkers()
                ic.photApertures.append(rectangle)
                cidap=rectangle.mySignal.connect(self.onRemoveAperture)
                ic.photApertureSignal.append(cidap)
                rectangle.modSignal.connect(self.onModifiedAperture)
            elif aper.type == 'Polygon':
                verts = aper.poly.get_xy()
                adverts = np.array([(ic0.wcs.all_pix2world(x, y, 0)) for (x,y) in verts])
                verts = [(ic.wcs.all_world2pix(ra, dec, 0)) for (ra,dec) in adverts]
                # Add polygon
                poly = PolygonInteractor(ic.axes,verts)                
                poly.showverts = aper.showverts
                ic.photApertures.append(poly)
                cidap=poly.mySignal.connect(self.onRemoveAperture)
                ic.photApertureSignal.append(cidap)
                poly.modSignal.connect(self.onModifiedAperture)
            elif aper.type == 'Pixel':
                x0,y0 = aper.rect.get_xy()
                w0    = aper.rect.get_width()
                h0    = aper.rect.get_height()
                angle = aper.rect.angle - ic0.crota2 + ic.crota2
                ra0,dec0 = ic0.wcs.all_pix2world(x0, y0, 0)
                ws = w0 * ic0.pixscale; hs = h0 * ic0.pixscale
                # Add rectangle
                x0,y0 = ic.wcs.all_world2pix(ra0, dec0, 0)
                w0 = ws/ic.pixscale; h0 = hs/ic.pixscale
                pixel = PixelInteractor(ic.axes, (x0,y0), w0)
                pixel.type = aper.type
                pixel.showverts = aper.showverts
                pixel.rect.set_xy((x0,y0))
                pixel.updateMarkers()
                ic.photApertures.append(pixel)
                cidap=pixel.mySignal.connect(self.onRemoveAperture)
                ic.photApertureSignal.append(cidap)
                #cidapm=
                pixel.modSignal.connect(self.onModifiedAperture)

    def blinkTab(self):
        ''' keep switching between two tabs until blink changes state '''
        itab = self.itabs.currentIndex()
        if itab == self.btab[0]:
            i = 1
        else:
            i = 0
        self.itabs.setCurrentIndex(self.btab[i])
        
    def blinkImages(self, event):
        ''' Blink between two images in different tabs or stop blinking'''
        if self.blink == 'off':
            self.btab = [self.itabs.currentIndex(),0]
            self.sb.showMessage("Select another tab to blink / click again to stop blinking", 2000)
            self.blink = 'select'
        else:
            self.blink = 'off'
            self.timer.stop()
        
    def fileQuit(self):
        """ Quitting the program """
        self.close()

    def trimCubeOld(self):  
        """ Trim the cube """
        self.sb.showMessage("Drag the mouse over the slice of the cube to trim ", 2000)
        self.trimcube = 'on'
        istab = self.spectra.index('All')
        self.stabs.setCurrentIndex(istab)
        sc = self.sci[istab]
        #sc.span.set_visible(True)
        sc.span.set_active(True)
    
    def trimMessage(self):
        self.wbox = QMessageBox()
        self.wbox.setText("How to trim a cube")
        self.wbox.setInformativeText('* Open slicer in the spectral panel\n\n'+\
                                '* Adjust the limits to the region of interest\n\n'+\
                                '* Trim the rest using this command again.')
        self.wbox.show()
        
    def trimCube(self):
        """Trimming the cube."""
        if self.slicer is None:
            # Message to define a slice first
            self.trimMessage()
        else:
            # Find indices of the slice
            xmin = self.slicer.xl
            xmax = self.slicer.xr
            istab = self.stabs.currentIndex()
            sc = self.sci[istab]
            if sc.xunit == 'THz':
                c = 299792458.0  # speed of light in m/s
                xmin, xmax = c/xmax*1.e-6, c/xmin*1.e-6
            indmin, indmax = np.searchsorted(self.specCube.wave, (xmin, xmax))
            indmax = min(len(self.specCube.wave) - 1, indmax)
            size = indmax-indmin
            nz, nx, ny = np.shape(self.specCube.flux)
            if size == nx:
                self.sb.showMessage("No cutting needed ", 2000)
            else:
                flags = QMessageBox.Yes | QMessageBox.No
                question = "Do you want to trim the cube to the slice selected on the spectrum ?"
                response = QMessageBox.question(self, "Question", question, flags)            
            if response == QMessageBox.Yes:
                self.sb.showMessage("Trimming the cube ", 2000)
                self.trimCube1D(indmin,indmax)
                self.saveCube()
                # Load trimmed cube
                self.reloadFile()
            elif QMessageBox.No:
                self.sb.showMessage("Trimming aborted ", 2000)
            else:
                pass
        
    def cropCube(self):
        """ Crop part of the cube """
        self.sb.showMessage("Crop the cube using the zoomed image shown ", 2000)
        # Get limits and center
        ic0 = self.ici[0]
        xlimits = ic0.axes.get_xlim()
        ylimits = ic0.axes.get_ylim()
        center =  ((xlimits[0]+xlimits[1])*0.5,(ylimits[0]+ylimits[1])*0.5)
        size = (np.abs((ylimits[1]-ylimits[0]).astype(int)),np.abs((xlimits[1]-xlimits[0]).astype(int)))
        nz,nx,ny = np.shape(self.specCube.flux)
        if size[0] == nx and size[1] == ny:
            self.sb.showMessage("No cropping needed ", 2000)
        else:
            flags = QMessageBox.Yes 
            flags |= QMessageBox.No
            question = "Do you want to crop the part of the cube shown on the image ?"
            response = QMessageBox.question(self, "Question",
                                                  question,
                                                  flags)            
            if response == QMessageBox.Yes:
                self.sb.showMessage("Cropping the cube ", 2000)
                self.cropCube2D(center,size)
                self.saveCube()
                # Load cropped cube
                self.reloadFile()
            elif response == QMessageBox.No:
                self.sb.showMessage("Cropping aborted ", 2000)
            else:
                pass

    def cropCube2D(self, center, size):
        """ Generate cropped cube """
        from astropy.nddata import Cutout2D
        ic0 = self.ici[0]
        co = Cutout2D(ic0.oimage,center,size,wcs=ic0.wcs)
        bb = co.bbox_original
        self.specCube.flux = self.specCube.flux[:,bb[0][0]:bb[0][1]+1,bb[1][0]:bb[1][1]+1]
        self.specCube.wcs = co.wcs
        if self.specCube.instrument == 'FIFI-LS':
            self.specCube.eflux = self.specCube.eflux[:,bb[0][0]:bb[0][1]+1,bb[1][0]:bb[1][1]+1]
            self.specCube.uflux = self.specCube.uflux[:,bb[0][0]:bb[0][1]+1,bb[1][0]:bb[1][1]+1]
            self.specCube.euflux = self.specCube.euflux[:,bb[0][0]:bb[0][1]+1,bb[1][0]:bb[1][1]+1]
            self.specCube.x = self.specCube.x[bb[1][0]:bb[1][1]+1]
            self.specCube.y = self.specCube.y[bb[0][0]:bb[0][1]+1]
            self.specCube.exposure = self.specCube.exposure[:,bb[0][0]:bb[0][1]+1,bb[1][0]:bb[1][1]+1]
        elif self.specCube.instrument == 'PACS':
            self.specCube.exposure = self.specCube.exposure[:,bb[0][0]:bb[0][1]+1,bb[1][0]:bb[1][1]+1]
        elif self.specCube.instrument == 'FORCAST':
            self.specCube.eflux = self.specCube.eflux[:,bb[0][0]:bb[0][1]+1,bb[1][0]:bb[1][1]+1]
            self.specCube.exposure = self.specCube.exposure[:,bb[0][0]:bb[0][1]+1,bb[1][0]:bb[1][1]+1]
        # Create a grid of points
        nz,ny,nx = np.shape(self.specCube.flux)
        xi = np.arange(nx); yi = np.arange(ny)
        xi,yi = np.meshgrid(xi,yi)
        self.points = np.array([np.ravel(xi),np.ravel(yi)]).transpose()

    def trimCube1D(self,xmin,xmax):
        """ Generate trimmed cube """
        self.specCube.flux = self.specCube.flux[xmin:xmax,:,:]
        self.specCube.wave = self.specCube.wave[xmin:xmax]
        nz,ny,nx = np.shape(self.specCube.flux)
        self.specCube.n = nz
        # Cut the cubes
        if self.specCube.instrument == 'FIFI-LS':
            self.specCube.eflux = self.specCube.eflux[xmin:xmax,:,:]
            self.specCube.uflux = self.specCube.uflux[xmin:xmax,:,:]
            self.specCube.euflux = self.specCube.euflux[xmin:xmax,:,:]
            self.specCube.exposure = self.specCube.exposure[xmin:xmax,:,:]
            self.specCube.atran = self.specCube.atran[xmin:xmax]
            self.specCube.response = self.specCube.response[xmin:xmax]
        if self.specCube.instrument == 'PACS':
            self.specCube.exposure = self.specCube.exposure[xmin:xmax,:,:]
        if self.specCube.instrument == 'FORCAST':
            self.specCube.eflux = self.specCube.eflux[xmin:xmax,:,:]
            self.specCube.exposure = self.specCube.exposure[xmin:xmax,:,:]

    def saveFits(self):
        """ Save the displayed image as a FITS file """
        from astropy.io import fits        
        # Dialog to save file
        fd = QFileDialog()
        fd.setLabelText(QFileDialog.Accept, "Save as")
        fd.setNameFilters(["Fits Files (*.fits)","PNG Files (*.png)",
                           "JPG Files (*.jpg)","PDF Files (*.pdf)","All Files (*)"])
        fd.setOptions(QFileDialog.DontUseNativeDialog)
        fd.setViewMode(QFileDialog.List)
        if (fd.exec()):
            fileName = fd.selectedFiles()
            outfile = fileName[0]
            itab = self.itabs.currentIndex()
            ic = self.ici[itab]
            band = self.bands[itab]
            if band == 'Flux' or band == 'Uflux' or band == 'Exp':
                instrument = self.specCube.instrument
            else:
                instrument = band
            # Check the 
            filename, file_extension = os.path.splitext(outfile)
            if file_extension == '.fits':
                # Primary header
                header = ic.wcs.to_header()
                header.remove('WCSAXES')
                header['INSTRUME'] = instrument
                header['OBJECT'] = (self.specCube.objname, 'Object Name')
                hdu = fits.PrimaryHDU(ic.oimage)
                hdu.header.extend(header)
                hdul = fits.HDUList([hdu])
                hdul.writeto(outfile,overwrite=True) # clobber true  allows rewriting
                hdul.close()
            elif file_extension == '.png' or file_extension == '.pdf' or file_extension == '.jpg':
                ic.fig.savefig(outfile)
            else:
                message = 'extension has to be *.fits, *.png, *.jpg or *.pdf' 
                print(message)
                self.sb.showMessage(message, 2000)

    def saveSpectrum(self):
        """ Save the displayed spectrum as a FITS/ASCII file or as PNG/PDF image """
        from astropy.io import fits        
        # Dialog to save file
        fd = QFileDialog()
        fd.setLabelText(QFileDialog.Accept, "Save as")
        fd.setNameFilters(["Fits Files (*.fits)","PNG Files (*.png)","JPG Files (*.jpg)",
                           "PDF Files (*.pdf)","ASCII Files (*.txt)", "CSV Files (*.csv)","All Files (*)"])
        fd.setOptions(QFileDialog.DontUseNativeDialog)
        fd.setViewMode(QFileDialog.List)
        if (fd.exec()):
            fileName = fd.selectedFiles()
            outfile = fileName[0]
            filename, file_extension = os.path.splitext(outfile)
            # Tabs
            istab = self.stabs.currentIndex()
            itab = self.itabs.currentIndex()
            sc = self.sci[istab]
            ic = self.ici[itab]
            n = istab-1
            # Compute area of the aperture
            if n >= 0:
                aperture = self.ici[0].photApertures[n].aperture
                path = aperture.get_path()
                transform = aperture.get_patch_transform()
                npath = transform.transform_path(path)
                s = self.specCube
                inpoints = s.points[npath.contains_points(s.points)]
                npoints = np.size(inpoints)/2
                ps = s.pixscale/3600.
                area = npoints*ps*ps            
            if file_extension == '.fits':
                # Primary header
                hdu = fits.PrimaryHDU()
                hdu.header['OBJ_NAME'] = (self.specCube.objname, 'Object Name')
                hdu.header['INSTRUME'] = (self.specCube.instrument, 'SOFIA instrument')
                hdu.header['REDSHIFT'] = (self.specCube.redshift, 'Object Redshift')
                if self.specCube.instrument == 'FIFI-LS':
                    hdu.header['BARYSHFT'] = (self.specCube.baryshift, 'Barycentric shift')
                if n >= 0:
                    aper = ic.photApertures[n]
                    if aper.type == 'Ellipse':
                        x0,y0 = aper.ellipse.center
                        pixel = np.array([[x0, y0]], np.float_)
                        world = ic.wcs.all_pix2world(pixel, 0)
                        hdu.header['APERTURE']=('Ellipse','Type of photometric aperture')
                        hdu.header['RA'] = (world[0][0], 'RA of aperture center')
                        hdu.header['DEC'] = (world[0][1], 'Dec of aperture center')
                        hdu.header['ANGLE'] = (aper.ellipse.angle, 'Angle of elliptical aperture [degs]')
                        hdu.header['MAJAX'] = (aper.ellipse.width*ic.pixscale*0.5, 'Major semi-axis of elliptical aperture [arcsec]')
                        hdu.header['MINAX'] = (aper.ellipse.height*ic.pixscale*0.5, 'Minor semi-axis of elliptical aperture [arcsec]')
                    elif aper.type == 'Circle':
                        x0,y0 = aper.ellipse.center
                        pixel = np.array([[x0, y0]], np.float_)
                        world = ic.wcs.all_pix2world(pixel, 0)
                        hdu.header['APERTURE']=('Circle','Type of photometric aperture')
                        hdu.header['RA'] = (world[0][0]/15., 'RA of aperture center [hours]')
                        hdu.header['DEC'] = (world[0][1], 'Dec of aperture center [degs]')
                        hdu.header['RADIUS'] = (aper.ellipse.height*ic.pixscale*0.5, 'Radius of circular aperture [arcsec]')
                    elif aper.type == 'Square':
                        x0,y0 = aper.xy[0]
                        pixel = np.array([[x0, y0]], np.float_)
                        world = ic.wcs.all_pix2world(pixel, 0)
                        hdu.header['APERTURE']=('Square','Type of photometric aperture')
                        hdu.header['RA'] = (world[0][0]/15., 'RA of aperture center [hours]')
                        hdu.header['DEC'] = (world[0][1], 'Dec of aperture center [degs]')
                        hdu.header['ANGLE'] = (aper.rect.angle, 'Angle of square aperture')
                        hdu.header['SIDE'] = (aper.rect.get_height()*ic.pixscale, 'Side of square aperture [arcsec]')
                    elif aper.type == 'Rectangle':
                        x0,y0 = aper.xy[0]
                        pixel = np.array([[x0, y0]], np.float_)
                        world = ic.wcs.all_pix2world(pixel, 0)
                        hdu.header['APERTURE']=('Rectangle','Type of photometric aperture')
                        hdu.header['RA'] = (world[0][0]/15., 'RA of aperture center [hours]')
                        hdu.header['DEC'] = (world[0][1], 'Dec of aperture center [degs]')
                        hdu.header['WIDTH'] = (aper.rect.get_width()*ic.pixscale, 'Width of rectangle aperture [arcsec]')
                        hdu.header['HEIGHT'] = (aper.rect.get_height()*ic.pixscale, 'Height of rectangle aperture [arcsec]')
                        hdu.header['ANGLE'] = (aper.rect.angle, 'Angle of rectangle aperture [degs]')
                    elif aper.type == 'Polygon':
                        hdu.header['APERTURE']=('Polygon','Type of photometric aperture')
                        xy = np.asarray(aper.poly.xy)
                        world = ic.wcs.all_pix2world(xy, 0)
                        i = 0
                        for w in world:
                            hdu.header['RA_PT'+"{:03d}".format(i)] = (w[0]/15.,'RA [hours] of polygon aperture point no {:d}'.format(i))
                            hdu.header['DECPT'+"{:03d}".format(i)] = (w[1],'Dec [degs] of polygon aperture point no {:d}'.format(i))
                            i += 1
                    hdu.header['AREA'] = (area,'Area aperture in sq. degs.')
                # Add extensions
                hdu1 = self.addExtension(sc.spectrum.wave,'WAVELENGTH','um',None)
                hdu2 = self.addExtension(sc.spectrum.flux,'FLUX','Jy',None)
                hdlist = [hdu,hdu1,hdu2]
                if self.specCube.instrument == 'FIFI-LS':
                    hdu3 = self.addExtension(sc.spectrum.eflux,'FLUX_ERROR','Jy',None)
                    hdu4 = self.addExtension(sc.spectrum.uflux,'UNCORR_FLUX','Jy',None)
                    hdu5 = self.addExtension(sc.spectrum.exposure,'EXPOSURE','s',None)
                    hdu6 = self.addExtension(self.specCube.atran,'ATM_TRANS','Norm',None)
                    hdlist.append(hdu3)
                    hdlist.append(hdu4)
                    hdlist.append(hdu5)
                    hdlist.append(hdu6)
                if self.specCube.instrument == 'FORCAST':
                    hdu3 = self.addExtension(sc.spectrum.eflux,'FLUX_ERROR','Jy',None)
                    hdu4 = self.addExtension(sc.spectrum.exposure,'EXPOSURE','s',None)
                    hdlist.append(hdu3)
                    hdlist.append(hdu4)
                # Save file
                hdul = fits.HDUList(hdlist)
                hdul.writeto(outfile,overwrite=True) # clobber true  allows rewriting
                hdul.close()
            elif file_extension == '.txt' or file_extension == '.csv':
                header = "# Object name: "+self.specCube.objname
                header += "\n# Instrument: "+self.specCube.instrument
                header += "\n# z: {:.8f}".format(self.specCube.redshift)
                header += "\n# Ref. Wav.: {:.8f}".format(self.specCube.l0)
                if n >= 0:
                    aper = ic.photApertures[n]
                    if aper.type == 'Ellipse':
                        x0,y0 = aper.ellipse.center
                        pixel = np.array([[x0, y0]], np.float_)
                        world = ic.wcs.all_pix2world(pixel, 0)                    
                        header += '\n# Aperture: Ellipse'
                        header += '\n# Center: {:.5f} {:.6f}'.format(world[0][0], world[0][1])
                        header += '\n# Angle: {:.1f} degs'.format(aper.ellipse.angle)
                        header += '\n# Axes: {:.1f} {:.1f} [arcsec]'.format(aper.ellipse.width*ic.pixscale,aper.ellipse.height*ic.pixscale)
                    elif aper.type == 'Circle':
                        x0,y0 = aper.ellipse.center
                        pixel = np.array([[x0, y0]], np.float_)
                        world = ic.wcs.all_pix2world(pixel, 0)                    
                        header += '\n# Aperture: Circle'
                        header += '\n# Center: {:.5f} {:.6f}'.format(world[0][0], world[0][1])
                        header += '\n# Radius: {:.1f} [arcsec]'.format(aper.ellipse.height*ic.pixscale)
                    elif aper.type == 'Square':
                        x0,y0 = aper.xy[0]
                        pixel = np.array([[x0, y0]], np.float_)
                        world = ic.wcs.all_pix2world(pixel, 0)
                        header += '\n# Aperture: Square'
                        header += '\n# Center: {:.5f} {:.6f}'.format(world[0][0], world[0][1])
                        header += '\n# Side: {:.1f} [arcsec]'.format(aper.rect.get_height()*ic.pixscale)
                        header += '\n# Angle: {:.1f} degs'.format(aper.rect.angle)
                    elif aper.type == 'Rectangle':
                        x0,y0 = aper.xy[0]
                        pixel = np.array([[x0, y0]], np.float_)
                        world = ic.wcs.all_pix2world(pixel, 0)
                        header += '\n# Aperture: Rectangle'
                        header += '\n# Center: {:.5f} {:.6f}'.format(world[0][0], world[0][1])
                        header += '\n# Height: {:.1f} [arcsec]'.format(aper.rect.get_height()*ic.pixscale)
                        header += '\n# Width: {:.1f} [arcsec]'.format(aper.rect.get_width()*ic.pixscale)
                        header += '\n# Angle: {:.1f} degs'.format(aper.rect.angle)
                    elif aper.type == 'Polygon':
                        header += '\n# Aperture: Polygon'
                        xy = np.asarray(aper.poly.xy)
                        world = ic.wcs.all_pix2world(xy, 0)
                        i = 0
                        for w in world:
                            header += '\n# Point {:03d}: {:.5f}h {:.6f}d'.format(i,w[0]/15.,w[1])
                            i += 1
                    header += '\n# Area aperture: {:.5f} degs'.format(area)
                #
                w = sc.spectrum.wave
                f = sc.spectrum.flux
                if file_extension == '.txt':
                    delimiter = ' '
                else:
                    delimiter = ','
                if self.specCube.instrument == 'FIFI-LS':
                    uf = sc.spectrum.uflux
                    ef = sc.spectrum.eflux
                    e  = sc.spectrum.exposure
                    a  = self.specCube.atran
                    # Normal ASCII file
                    fmt = delimiter.join(["%10.6e"]*6)
                    with open(outfile, 'wb') as file:
                        file.write(header.encode())
                        file.write(b'\n"wavelength","flux","eflux","uflux","exposure","atran"\n')
                        np.savetxt(file, np.column_stack((w,f,ef,uf,e,a)), fmt=fmt, delimiter=delimiter)
                else:
                    fmt = delimiter.join(["%10.6e"]*2)
                    with open(outfile, 'wb') as file:
                        file.write(header.encode())
                        file.write(b'\n"wavelength","flux"\n')
                        np.savetxt(file, np.column_stack((w,f)), fmt=fmt, delimiter=delimiter)                
            elif file_extension == '.png' or file_extension == '.pdf' or file_extension == '.jpg':
                sc.fig.savefig(outfile)
            else:
                message = "Extension has to be *.fits, *.txt, *.csv, *.png, *.jpg, or *.pdf "
                self.sb.showMessage(message, 2000)
                print(message)
    
    def savelCube(self):
        """ Save a cube after subtracting the continuum """
        if self.continuum is not None:
            self.saveCube(cont=True)
        else:
            message = "No continuum has been fit to the cube"
            self.sb.showMessage(message, 2000)
            print(message)

    def saveMaskedCube(self):
        """ Save a cube after masking """
        message = "Saving current version of the cube"
        self.saveCube()
        print(message)
        
    def saveCube(self, cont=False):
        """Save a trimmed/cropped cube."""
        from astropy.io import fits
        # Dialog to save file
        fd = QFileDialog()
        fd.setLabelText(QFileDialog.Accept, "Save as")
        fd.setNameFilters(["Fits Files (*.fits)","All Files (*)"])
        fd.setOptions(QFileDialog.DontUseNativeDialog)
        fd.setViewMode(QFileDialog.List)            
        if (fd.exec()):
            fileName = fd.selectedFiles()
            outfile = fileName[0]
            # Update name in cube to allow automatical reloading
            self.specCube.filename = outfile
            if cont and self.continuum is not None:
                flux = self.specCube.flux - self.continuum
            else:
                flux = self.specCube.flux
            # Reusable header
            header = self.specCube.wcs.to_header()
            header.remove('WCSAXES')
            header['CRPIX3'] = (self.specCube.crpix3,'Reference pixel')
            header['CRVAL3'] = (self.specCube.crval3,'Reference pixel value')
            header['CDELT3'] = (self.specCube.cdelt3,'Increment')
            header['NAXIS3'] = (self.specCube.n,'3rd dimension')
            header['INSTRUME'] = (self.specCube.instrument, 'Instrument')
            header['DATE-OBS'] = (self.specCube.obsdate, 'Date of the observation')  
            try:
                header['BUNIT'] = self.specCube.header['BUNIT']
            except:
                print('No flux unit defined in original cube')
            if self.specCube.instrument == 'FIFI-LS':
                header['CUNIT3'] = ('um','Wavelength unit')
                header['OBJ_NAME'] = (self.specCube.objname, 'Object Name')
                header['REDSHIFT'] = (self.specCube.redshift, 'Redshift')
                header['FILEGPID'] = (self.specCube.filegpid, 'File Group ID')
                header['BARYSHFT'] = (self.specCube.baryshift, 'Barycentric shift')
                header['RESOLUN'] = (self.specCube.resolution, 'Spectral resolution')
                header['ZA_START'] = (self.specCube.za[0],'Zenith angle [degrees]')
                header['ZA_END'] = (self.specCube.za[1],'Zenith angle [degrees]')
                header['ALTI_STA'] = (self.specCube.altitude[0],'Equivalent aircraft pressure altitude [feet]')
                header['ALTI_END'] = (self.specCube.altitude[1],'Equivalent aircraft pressure altitude [feet]')
                header['PIXSCAL'] = (self.specCube.pixscale,'Pixel scale [arcsec]' )
                header['RESOLUN'] = (self.specCube.resolution,'Spectral resolution')
                header['DETCHAN'] = (self.specCube.header['DETCHAN'],'Data comes from ..')
                header['NAXIS'] = (3,'Number of axis')  
                exptime = self.specCube.header['EXPTIME']
                nexp = self.specCube.header['NEXP']
                header['EXPTIME'] = (exptime,'On-source exposure time in seconds')
                header['NEXP'] = (nexp,'Number of exposure in source')
                # Primary header
                hdu = fits.PrimaryHDU()
                hdu.header.extend(header)
                if cont and self.continuum is not None:
                    uflux = self.specCube.uflux - self.continuum
                else:
                    uflux = self.specCube.uflux
                # Extensions
                hdu1 = self.addExtension(flux,'FLUX','Jy',header)
                hdu2 = self.addExtension(self.specCube.eflux,'ERROR','Jy',header)
                hdu3 = self.addExtension(uflux,'UNCORRECTED_FLUX','Jy',header)
                hdu4 = self.addExtension(self.specCube.euflux,'UNCORRECTED_ERROR','Jy',header)
                hdu5 = self.addExtension(self.specCube.wave,'WAVELENGTH','um',None)
                hdu6 = self.addExtension(self.specCube.x,'X',None,None)
                hdu7 = self.addExtension(self.specCube.y,'Y',None,None)
                hdu8 = self.addExtension(self.specCube.atran,'TRANSMISSION',None,None)
                hdu9 = self.addExtension(self.specCube.response,'RESPONSE',None,None)
                hdu10 = self.addExtension(self.specCube.exposure * nexp / exptime,'EXPOSURE_MAP',None,header)
                if self.specCube.watran is not None:
                    uatran = np.array([self.specCube.watran, self.specCube.uatran])
                    hdu11 = self.addExtension(uatran,'UNSMOOTHED_TRANSMISSION',None,None)
                    hdul = fits.HDUList([hdu, hdu1, hdu2, hdu3, hdu4, hdu5, hdu6, hdu7, hdu8, hdu9, hdu10, hdu11])
                else:
                    hdul = fits.HDUList([hdu, hdu1, hdu2, hdu3, hdu4, hdu5, hdu6, hdu7, hdu8, hdu9, hdu10])            
                #hdul.info()    
                hdul.writeto(outfile,overwrite=True) 
                hdul.close()
            elif self.specCube.instrument == 'GREAT':
                header['OBJECT'] = (self.specCube.objname, 'Object Name')
                c = 299792458.0  # speed of light in m/s 
                header['VELO-LSR'] = self.specCube.redshift * c
                header['RESTFREQ'] = self.specCube.header['RESTFREQ']
                header['CUNIT3'] = ('m/s','Velocity unit')
                #eta_fss=0.97
                #eta_mb =0.67
                #calib = 971.
                #factor = calib*eta_fss*eta_mb
                temperature = flux #/ factor  # Flux is already a Tb
                # Primary header
                hdu = fits.PrimaryHDU(temperature)
                header['NAXIS'] = (3,'Number of axis')
                header['BMAJ'] = self.specCube.header['BMAJ']
                header['BMIN'] = self.specCube.header['BMIN']
                header['DATAMAX'] = self.specCube.header['DATAMAX']
                hdu.header.extend(header)
                hdul = fits.HDUList([hdu])
                hdul.writeto(outfile,overwrite=True) 
                hdul.close()
            elif self.specCube.instrument == 'FORCAST':
                header['OBJECT'] = (self.specCube.objname, 'Object Name')
                c = 299792458.0  # speed of light in m/s 
                header['CUNIT3'] = ('um','Units of the spectral axis')
                # Primary header
                hdu = fits.PrimaryHDU()
                hdu.header.extend(header)
                # Extensions
                hdu1 = self.addExtension(flux,'FLUX','Jy',header)
                hdu2 = self.addExtension(self.specCube.eflux**2,'VARIANCE','Jy2',header)
                hdu3 = self.addExtension(self.specCube.exposure * self.specCube.nz,
                                         'EXPOSURE', None, header)
                hdul = fits.HDUList([hdu, hdu1, hdu2, hdu3])
                hdul.writeto(outfile, overwrite=True) 
                hdul.close()
            elif self.specCube.instrument == 'PACS':
                """ Experimental """
                header['NAXIS'] = (3,'Number of axis')
                header['OBJECT'] = (self.specCube.objname, 'Object Name')
                c = 299792.458  # speed of light in km/s 
                header['REDSHFTV'] = self.specCube.redshift * c
                # Primary header
                hdu = fits.PrimaryHDU()
                hdu.header.extend(header)
                # Extensions
                hdu1 = self.addExtension(flux,'image','Jy',header)
                hdu2 = self.addExtension(self.specCube.exposure,'coverage',None,header)
                # Writing the wavelength table
                wave = [np.array([[x] for x in self.specCube.wave])]
                layer=np.array([np.arange(len(self.specCube.wave))])
                n = str(len(self.specCube.wave))
                col1 = fits.Column(name='wavelen', format=n+'D', array=wave)
                col2 = fits.Column(name='layer',format = n+'J', array=layer)
                coldefs = fits.ColDefs([col1, col2])
                hdw = fits.BinTableHDU.from_columns(coldefs)
                hdw.header['EXTNAME'] = 'wcs-tab'
                hdw.header.extend(header)
                hdul = fits.HDUList([hdu, hdu1, hdu2, hdw])            
                hdul.writeto(outfile,overwrite=True) 
                hdul.close()
            elif self.specCube.instrument == 'CARMA':
                header['NAXIS'] = (3,'Number of axis')
                c = 299792458.0  # speed of light in m/s 
                header['TELESCOP'] = self.specCube.header['TELESCOP']
                header['BMAJ'] = self.specCube.header['BMAJ']
                header['BMIN'] = self.specCube.header['BMIN']
                header['BPA'] = self.specCube.header['BPA']
                header['BSCALE'] = self.specCube.header['BSCALE']
                header['BZERO'] = self.specCube.header['BZERO']
                header['RESTFRQ'] = self.specCube.header['RESTFRQ']
                header['ALTRVAL'] = self.specCube.header['ALTRVAL']
                header['ALTRPIX'] = self.specCube.header['ALTRPIX']
                header['CTYPE3'] = self.specCube.header['CTYPE3']
                header['CRPIX3'] = self.specCube.header['CRPIX3']
                header['CRVAL3'] = self.specCube.header['CRVAL3']
                header['CDELT3'] = self.specCube.header['CDELT3']
                header['CUNIT3'] = ('m/s','Velocity unit')
                hdu = fits.PrimaryHDU(flux)
                hdu.header.extend(header)
                hdul = fits.HDUList([hdu])
                hdul.writeto(outfile,overwrite=True) 
                hdul.close()
            else:
                pass  

    def addExtension(self,data, extname, unit, hdr):
        from astropy.io import fits
        hdu = fits.ImageHDU()
        hdu.data = data
        hdu.header['EXTNAME']=(extname)
        if unit !=None: hdu.header['BUNIT']=(unit)
        if hdr != None: hdu.header.extend(hdr)
        return hdu
        
    def sliceCube(self):
        """Select part of the cube."""
        self.sb.showMessage("Define slice of the cube ", 1000)
        self.slice = 'on'
        istab = self.stabs.currentIndex()
        sc = self.sci[istab]
        ## toggle off continuum
        try:
            SI = sc.guess
            SI.showverts = False
            SI.line.set_visible(SI.showverts)
        except:
            pass
        ## Hide lines
        if sc.displayLines == True:
            sc.displayLines = False
            sc.setLinesVisibility(sc.displayLines)
            sc.fig.canvas.draw_idle()
        try:
            self.slider.disconnect()
            self.slider = None
        except BaseException:
            pass
        sc.span.active = True
        
    def nAper(self):
        '''Return the number of the aperture'''
        istab = self.stabs.currentIndex()
        tabname = self.stabs.tabText(istab)
        # print('Name of tab is ', tabname)
        if tabname in ['All','PSF']:
            n = -1
        elif tabname == 'Pix':
            n = 0
        else:
            n = int(tabname)
        return n
        
    def computeDistance(self):
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        # Disactivate aperture
        n = self.nAper()
        if n >= 0:
            ap = ic.photApertures[n]
            ap.showverts = False
            ap.line.set_visible(ap.showverts)
        self.DS = DistanceSelector(ic.axes, ic.fig, ic.wcs, self.printDistance)
        
    def printDistance(self, xy):
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        n = self.nAper()
        if n >= 0:
            ap = ic.photApertures[n]
            ap.showverts = True
            ap.line.set_visible(ap.showverts)        
            ic.fig.canvas.draw_idle()
        
    def estimatePSF(self):
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        # Remove previous PSF
        try:                
            tabnames = [self.stabs.tabText(i) for i in range(len(self.stabs))]
            if 'PSF' in tabnames:
                # Remove PSF tab
                istab = tabnames.index('PSF')
                print('tab ', istab, tabnames[istab])
                self.stabs.currentChanged.disconnect()
                self.removeSpecTab(istab, False)
                self.stabs.setCurrentIndex(1)  # Pixel tab
                self.stabs.currentChanged.connect(self.onSTabChange)
                print('tab removed')
        except:
            pass
        # Disactivate aperture
        n = self.nAper()
        if n >= 0:
            ap = ic.photApertures[n]
            ap.showverts = False
            ap.line.set_visible(ap.showverts)
        # Estimate new PSF
        self.PsfS = EllipseSelector(ic.axes, self.onDrawPSF, drawtype='line', useblit=True,
                                  button=[1, 3], minspanx=5, minspany=5, spancoords='pixels',
                                  rectprops = dict(facecolor='g',
                                                   edgecolor='g',
                                                   alpha=0.8, fill=False),
                                  lineprops = dict(color='g',
                                                   linestyle='-',
                                                   linewidth=2,
                                                   alpha=0.8),
                                  interactive=False)
                                  #state_modifier_keys = dict(square='square', center='center')
                                  #)
        self.PsfS.state.add('center')
        self.PsfS.state.add('square')
        ic.fig.canvas.draw_idle()                                      
        
    def onDrawPSF(self, eclick, erelease):
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
        # print('x1,y1,x2,y2',x1,y1,x2,y2)
        x0 = (x1 + x2) * 0.5
        y0 = (y1 + y2) * 0.5
        radius = np.hypot(x1-x2,y1-y2) * 0.5
        # print('center, radius ', (x0, y0), radius)
        # Deactivate selector
        self.PsfS.set_active(False)
        for artist in self.PsfS.artists:
            artist.remove()
        self.PsfS = None
        self.PsfI = PsfInteractor(ic.axes, (x0, y0), radius)
        # Connect signals
        self.PsfCidr = self.PsfI.mySignal.connect(self.onRemovePsf)
        self.PsfCidm = self.PsfI.modSignal.connect(self.onModifiedPsf)
        # Make apertures active again
        n = self.nAper()
        if n >= 0:
            ap = ic.photApertures[n]
            ap.showverts = True
            ap.line.set_visible(ap.showverts)        
            ic.fig.canvas.draw_idle()
        # print('Call routine to compute PSF and display on right panel')
        self.addPsfPlot()
        # Move to tab with plot
        self.stabs.setCurrentIndex(len(self.stabs)-1)
        
    def onRemovePsf(self, event):
        """Remove the PSF from the image and the tab with PSF profile."""
        #self.PsfI.disconnect(self.PsfCidr)
        #self.PsfI.disconnect(self.PsfCidm)
        #itab = self.itabs.currentIndex()
        istab = self.stabs.currentIndex()
        self.stabs.currentChanged.disconnect()
        self.removeSpecTab(istab, False)
        self.stabs.setCurrentIndex(1)  # Pixel tab
        self.stabs.currentChanged.connect(self.onSTabChange)
    
    def computePsfData(self):
         # Get center and radii and compute distances
        x0, y0 = self.PsfI.innerCircle.center
        xi, yi = self.imagePoints
        distance = np.hypot(xi - x0, yi - y0).ravel()
        r0 = self.PsfI.inRadius
        r1 = self.PsfI.outRadius

        # Continuum
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        image = ic.oimage.ravel()
        mc = (distance > r0) & (distance <= r1)
        annulus = image[mc]
        cont = np.nanmedian(annulus)
        
        # New fluxes
        mi = distance <= r0
        dist = distance[mi] * ic.pixscale  # transformed into arcsec
        flux = image[mi] - cont

        return dist, flux       
    
    def onModifiedPsf(self, event):
        """Update the profile after changing the aperture to compute the PSF."""   
        istab = self.stabs.currentIndex()
        sc = self.sci[istab]
        if sc.name != 'PSF':
            tabnames = [self.stabs.tabText(i) for i in range(len(self.stabs))]
            i = tabnames.index('PSF')
            self.stabs.setCurrentIndex(i)
            sc = self.sci[i]
        # Compute new distance and flux
        dist, flux = self.computePsfData()
        # Update figure of PSF
        sc.updatePSF(dist, flux)    
        
    def centroidPSF(self, event):
        """Recenter the PSF aperture on centroid of intensity."""
        x0, y0 = self.PsfI.innerCircle.center
        xi, yi = self.imagePoints
        distance = np.hypot(xi - x0, yi - y0)
        r0 = self.PsfI.inRadius
        mc = distance <= r0
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        image = ic.oimage
        ima2 = image[mc]
        neg = ima2 < 0
        ima2[neg] = 0
        totima = np.nansum(ima2)
        x0 = np.nansum(ima2 * xi[mc]) / totima
        y0 = np.nansum(ima2 * yi[mc]) / totima
        self.PsfI.innerCircle.center = x0, y0
        self.PsfI.outerCircle.center = x0, y0
        # Update interactor
        self.PsfI.updateInteractor()
        # Update plot
        self.onModifiedPsf('center changed')
        
    def centerPSF(self, event):
        """Recenter the PSF aperture optimizing the center with a fit"""
        from lmfit import Parameters, minimize
        x0, y0 = self.PsfI.innerCircle.center
        xi, yi = self.imagePoints
        distance = np.hypot(xi - x0, yi - y0)
        r0 = self.PsfI.inRadius
        r1 = self.PsfI.outRadius
        mc = distance <= r0
        ma = (distance > r0) & (distance <= r1)
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        stab = self.stabs.currentIndex()
        sc = self.sci[stab]
        image = ic.oimage
        annulus = image[ma]
        cont = np.nanmedian(annulus)
        image -= cont
        fit_params = Parameters()
        fit_params.add('A', value=np.nanmax(image[mc]))
        if sc.sigma is None:
            sigmaguess = 5
        else:
            sigmaguess = sc.sigma
        try:
            fit_params.add('s', value=sigmaguess, min=1, max=2*sigmaguess)
            fit_params.add('x0', value=x0)
            fit_params.add('y0', value=y0)
            idx = np.isfinite(image)
            if ic.pixscale is None:
                out = minimize(residualsPsf, fit_params, args=(xi[idx],yi[idx]),
                               kws={'data': image[idx],})
            else:
                weight = (distance * ic.pixscale)
                weight[~mc] = np.nanmax(weight)
                out = minimize(residualsPsf, fit_params, args=(xi[idx],yi[idx]),
                           kws={'data': image[idx], 'err':weight[idx]})
            x0 = out.params['x0'].value
            y0 = out.params['y0'].value
            sc.sigma = out.params['s'].value
            self.PsfI.innerCircle.center = x0, y0
            self.PsfI.outerCircle.center = x0, y0
            # Update interactor
            self.PsfI.updateInteractor()
            # Update plot
            self.onModifiedPsf('center changed')
        except:
            print('Fit has failed')
        
    def maskCube(self):
        """Mask a slice of the cube."""
        # Dialog to choose between masking with contour level or polygon
        msgBox = QMessageBox()
        msgBox.setText('Mask the region')
        msgBox.addButton('lower than minimum contour', QMessageBox.ActionRole)
        msgBox.addButton('inside a polygon', QMessageBox.ActionRole)
        msgBox.addButton('outside a polygon', QMessageBox.ActionRole)
        msgBox.addButton('cancel', QMessageBox.ActionRole)
        self.result = msgBox.exec()            
        if self.result == 0:
            self.maskCubeContour()
        elif self.result == 1:
            self.maskCubeInsidePolygon()
        elif self.result == 2:
            self.maskCubeOutsidePolygon()
        else:
            pass
            
    def maskCubeContour(self):
        self.sb.showMessage("Masking data", 2000)
        # Find the tab with levels
        itab = None
        for ih in self.ihi:
            if len(ih.lev) > 0:
                # Identify minimum level
                minlev = min(ih.levels)
                itab = self.ihi.index(ih)
                # Remove lowest level contour
                ind = ih.levels.index(minlev)
                ih.lev[ind].remove()
                # Remove level from lists
                del ih.lev[ind]
                del ih.levels[ind]
                ih.fig.canvas.draw_idle()
                # Modify contours
                self.onModifyContours(-1000-ind)
        if itab is not None:
            band = self.bands[itab]
            if band in ['Flux','uFlux','Exp','C0','M0','M1','M2','M3','M4','v','sv']:
                ic = self.ici[itab]
                mask2d = ic.oimage < minlev
                xi = np.arange(self.specCube.nx); yi = np.arange(self.specCube.ny)
                xx,yy = np.meshgrid(xi,yi)
                xx = xx[mask2d]
                yy = yy[mask2d]
                # Mask images and cubes
                self.specCube.flux[:,yy,xx] = np.nan
                icis = [self.ici[0]]
                if self.specCube.instrument == 'FIFI-LS':
                    print('Masking uflux')
                    self.specCube.uflux[:,yy,xx] = np.nan
                    icis.append(self.ici[1])
                for ic0 in icis:
                    image = ic0.oimage
                    image[yy,xx] = np.nan
                    cmin = ic0.cmin; cmax=ic0.cmax
                    ic0.updateImage(image)
                    ic0.updateScale(cmin,cmax)
                # mask C0, Mi, v, sv
                sbands = [self.C0, self.M0, self.M1, self.M2, self.M3, self.M4, 
                          self.v, self.sv,self.L0, self.L1]
                bands = ['C0','M0','M1','M2','M3','M4','v','sv','L0','L1']
                for b,sb in zip(bands,sbands):
                    if sb is not None:
                        itab = self.bands.index(b)
                        sb[yy,xx] = np.nan
                        ic = self.ici[itab]
                        ih = self.ihi[itab]
                        clim = ic.image.get_clim()
                        ic.updateImage(sb)
                        ic.fig.canvas.draw_idle()
                        ih.axes.cla()
                        ih.compute_initial_figure(image=sb, xmin=clim[0],xmax=clim[1])
                        ic.updateScale(clim[0],clim[1])
            else:
                self.sb.showMessage('Contours are considered only in cube or derivated images',4000)
        else:
            self.sb.showMessage("Please, define contours before masking ", 4000)                                

    def maskCubeInsidePolygon(self):
        self.sb.showMessage("Draw the region to mask", 2000)        
        # Start a Selector to define a polygon aperture
        itab = self.itabs.currentIndex()
        band = self.bands[itab]
        if band not in ['Flux','uFlux','Exp','C0','M0','M1','M2','M3','M4','v','sv','L0','L1']:
            itab = 0
            self.itabs.setCurrentIndex(itab)            
        ic = self.ici[itab]
        if ic.toolbar._active == 'ZOOM':
            ic.toolbar.zoom()  # turn off zoom
            #self.doZoomAll()
            x = ic.axes.get_xlim()
            y = ic.axes.get_ylim()
            self.zoomlimits = [x,y]            
        self.LS = PolygonSelector(ic.axes, onselect=self.insideMask)

    def maskCubeOutsidePolygon(self):
        self.sb.showMessage("Draw the region to mask", 2000)        
        # Start a Selector to define a polygon aperture
        itab = self.itabs.currentIndex()
        band = self.bands[itab]
        if band not in ['Flux','uFlux','Exp','C0','M0','M1','M2','M3','M4','v','sv','L0','L1']:
            itab = 0
            self.itabs.setCurrentIndex(itab)            
        ic = self.ici[itab]
        if ic.toolbar._active == 'ZOOM':
            ic.toolbar.zoom()  # turn off zoom
            x = ic.axes.get_xlim()
            y = ic.axes.get_ylim()
            self.zoomlimits = [x,y]            
        self.LS = PolygonSelector(ic.axes, onselect=self.outsideMask)

    def insideMask(self,verts):
        self.onMask(verts,True)

    def outsideMask(self, verts):
        self.onMask(verts,False)

    def onMask(self, verts, inside = True):
        """ Uses the vertices of the mask to mask the cube (and moments) """
        s= self.specCube
        poly = Polygon(list(verts), fill=False, closed=True, color='lime')
        self.disactiveSelectors()        
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        ic.axes.add_patch(poly)
        """ Ask to mask points """
        flags = QMessageBox.Yes 
        flags |= QMessageBox.No
        if inside:
            question = "Do you want to mask the data inside the polygon ?"
        else:
            question = "Do you want to mask the data outside the polygon ?"
        response = QMessageBox.question(self, "Question",
                                        question,
                                        flags)            
        if response == QMessageBox.Yes:
            path = poly.get_path()
            transform = poly.get_patch_transform()
            npath = transform.transform_path(path)
            # We could transform this path into another one with different astrometry here before checking the points inside
            insidePts = npath.contains_points(s.points)
            if inside:
                inpoints = s.points[insidePts]
            else:
                inpoints = s.points[~insidePts]
            xx,yy = inpoints.T
            poly.remove()            
            self.sb.showMessage("Masking data ", 2000)
            self.specCube.flux[:,yy,xx] = np.nan
            icis = [self.ici[0]]
            if self.specCube.instrument == 'FIFI-LS':
                self.specCube.uflux[:,yy,xx] = np.nan
                icis.append(self.ici[1])
            for ic0 in icis:
                image = ic0.oimage
                image[yy,xx] = np.nan
                cmin = ic0.cmin; cmax=ic0.cmax
                ic0.updateImage(image)
                ic0.updateScale(cmin,cmax)
            # mask C0, Mi, v, sv
            sbands = [self.C0, self.M0, self.M1, self.M2, self.M3, self.M4, self.v, self.sv,self.L0,self.L1]
            bands = ['C0','M0','M1','M2','M3','M4','v','sv','L0','L1']
            for b,sb in zip(bands,sbands):
                if sb is not None:
                    itab = self.bands.index(b)
                    sb[yy,xx] = np.nan
                    ic0 = self.ici[itab]
                    ih = self.ihi[itab]
                    clim = ic0.image.get_clim()
                    ic0.updateImage(sb)
                    ic0.updateScale(clim[0],clim[1])
                    ih.axes.cla()
                    ih.compute_initial_figure(image=sb, xmin=clim[0],xmax=clim[1])
            x,y = ic.zoomlimits
            ic.axes.set_xlim(x)
            ic.axes.set_ylim(y)
            ic.fig.canvas.draw_idle()
        elif QMessageBox.No:
            poly.remove()
        
    def computeZeroMoment(self):        
        c = 299792458. # m/s
        w = self.specCube.wave
        dw = [] 
        dw.append([w[1]-w[0]])
        dw.append(list((w[2:]-w[:-2])*0.5))
        dw.append([w[-1]-w[-2]])
        dw = np.concatenate(dw)
        nz,ny,nx = np.shape(self.specCube.flux)
        self.M0 = np.zeros((ny,nx))
        for i in range(nx):
            for j in range(ny):
                Snu = self.specCube.flux[:,j,i]
                Slambda = c*(Snu-np.nanmedian(Snu))/(w*w)*1.e6   # [Jy * Hz / um]
                self.M0[j,i] = np.nansum(Slambda*dw)*1.e-26 # [Jy Hz]  (W/m2 = Jy*Hz*1.e-26)
        
    def zeroMoment(self):
        """ Compute and display zero moment of flux """
        self.computeZeroMoment()        
        band = 'M0'
        # Open tab and display the image
        self.bands.append(band)
        t,ic,ih,h,c1,c2,c3,c4 = self.addImage(band)
        self.tabi.append(t)
        self.ici.append(ic)
        self.ihi.append(ih)
        self.ihcid.append(h)
        self.icid1.append(c1)
        self.icid2.append(c2)
        self.icid3.append(c3)
        self.icid4.append(c4)        
        ic.compute_initial_figure(image=self.M0,wcs=self.specCube.wcs,title=band,
                                  cMap=self.colorMap,cMapDir=self.colorMapDirection,stretch=self.stretchMap)
        # Callback to propagate axes limit changes among images
        #ic.cid = ic.axes.callbacks.connect('xlim_changed' and 'ylim_changed', self.doZoomAll)
        ic.cid = ic.axes.callbacks.connect('ylim_changed', self.doZoomAll)
        ih = self.ihi[self.bands.index(band)]
        clim = ic.image.get_clim()
        ih.compute_initial_figure(image=self.M0,xmin=clim[0],xmax=clim[1])
        # Add apertures
        self.addApertures(ic)
        # Add contours
        self.addContours(ic) 
        # Align with spectral cube
        ic0 = self.ici[0]
        x = ic0.axes.get_xlim()
        y = ic0.axes.get_ylim()
        ra,dec = ic0.wcs.all_pix2world(x, y, 0)
        x,y = ic.wcs.all_world2pix(ra, dec, 0)            
        ic.axes.set_xlim(x)
        ic.axes.set_ylim(y)
        ic.changed = True

    def addContours(self, ic):
        """ Add existing contours to newly added images """
        if self.contours == 'on':
            ih0 = None
            for ih,ic_ in zip(self.ihi, self.ici):
                if len(ih.levels) > 0:
                    ih0 = ih
                    ic0 = ic_
            if ih0 is not None:
                itab = self.ici.index(ic0)
                ic.contour0 = itab
                #ic.contour = ic.axes.contour(ic0.oimage, ih0.levels, colors=self.colorContour,
                #                             transform=ic.axes.get_transform(ic0.wcs))
                ic.changed = False
            else:
                pass
            
    def overlapContours(self):
        """ Compute contours and overlap/remove them on images """
        if self.contours == 'off':
            self.sb.showMessage("Click again to remove contours", 2000)
            self.drawContours()
            self.contours = 'on'
            #self.menuContours.setChecked(True)
        else:
            self.contours = 'off'
            self.tabContour[0] = None
            #self.menuContours.setChecked(False)
            # Remove level lines in histogram 
            for ih in self.ihi:
                if len(ih.lev) > 0:
                    print('There are ',len(ih.lev),' contour levels')
                    ih.levSignal.disconnect()
                    ih.removeLevels()
            # Remove contours
            for ic in self.ici:
                if ic.contour is not None:
                    for coll in ic.contour.collections:
                        coll.remove()
                    ic.contour = None
                    ic.changed = True
            # Update current tab
            itab = self.itabs.currentIndex()
            ic0 = self.ici[itab]
            ic0.fig.canvas.draw_idle()
            ic0.changed = False
            self.sb.showMessage('Contours deleted ', 1000)

    def drawContours(self):
        """ Draw contours of image """
        itab = self.itabs.currentIndex()
        ic0 = self.ici[itab]
        ih0 = self.ihi[itab]
        # Conserve tab number
        self.tabContour[0] = itab
        if ih0.median is None or ih0.sdev is None:
            print('updating histogram')
            ih0.update_figure(image=ic0.oimage)
        if self.bands[itab] == 'Cov':
            ih0.levels = list(np.arange(ih0.min,ih0.max,(ih0.max-ih0.min)/8))
        else:
            levels = ih0.median + np.array([-1,0,1,2,3,5,10]) * ih0.sdev
            mask = levels < ih0.max
            ih0.levels = list(levels[mask])
        #print('Contour levels are: ',ih0.levels)
        if self.specCube.instrument != 'FIFI-LS':
            ismo = ndimage.gaussian_filter(ic0.oimage, sigma=1.0, order=0)
        else:
            ismo = ic0.oimage
        ic0.contour = ic0.axes.contour(ismo, ih0.levels, colors=self.colorContour[0])
        ic0.fig.canvas.draw_idle()
        # Add levels to histogram
        ih0.drawLevels()
        # Connect signal event to action
        #cidh0=
        ih0.levSignal.connect(self.onModifyContours)
        # Update contours on all other images
        ici = self.ici.copy()
        ici.remove(ic0)
        for ic in ici:
            #ic.contour = ic.axes.contour(ic0.oimage,ih0.levels, colors=self.colorContour,
            #                             transform=ic.axes.get_transform(ic0.wcs))
            ic.contour0 = itab
            ic.changed = True

    def onModifyContours(self, n):
        """ Called by mysignal in the histogram canvas if contour levels change """
        # In some cases, computing contours can be computationally intense. So
        # we call threads in the case of new/modified contours
        # Unfortunately this does not work because matplotlib is not thread safe.        
        itab = self.itabs.currentIndex()
        ic0 = self.ici[itab]
        ih0 = self.ihi[itab]
        ncontours = len(ic0.contour.collections)
        if ic0.contour is not None:
            if n >= 1000:
                # Add a brand new level
                n -= 1000
                if self.specCube.instrument != 'FIFI-LS':
                    ismo = ndimage.gaussian_filter(ic0.oimage, sigma=0.5, order=0)
                else:
                    ismo = ic0.oimage
                new = ic0.axes.contour(ismo, [ih0.levels[n]], colors=self.colorContour[0])
                # Insert new contour in the contour collection
                ic0.contour.collections.insert(n, new.collections[0])                
            elif n <= -1000:
                n += 1000
                # Remove contour from image
                if n < ncontours:
                    ic0.axes.collections.remove(ic0.contour.collections[n])
                    # Delete element from contour collection list
                    del ic0.contour.collections[n]
                else:
                    print('contour ',n,' does not exist')
            else:
                # Move level by removing contour and adding it at new level
                if n < ncontours:
                    ic0.axes.collections.remove(ic0.contour.collections[n])
                    new = ic0.axes.contour(ic0.oimage, [ih0.levels[n]], colors=self.colorContour[0])
                    ic0.contour.collections[n] = new.collections[0]
                else:
                    print('contour ',n,' does not exist')
            self.modifyOtherImagesContours(itab)
                
    def modifyOtherImagesContours(self, i0):
        """ Once the new contours are computed, propagate them to other images """        
        ic0 = self.ici[i0]
        # ih0 = self.ihi[i0]
        ic0.fig.canvas.draw_idle()        
        ici = self.ici.copy()
        ici.remove(ic0)
        for ic in ici:
            if ic.contour is not None:
                # Remove previous contours
                for coll in ic.contour.collections:
                    coll.remove()
                ic.contour = None
                # Compute new contours
                ic.contour0 = i0
                #levels =  sorted(ih0.levels)
                #ic.contour = ic.axes.contour(ic0.oimage, levels, colors=self.colorContour,
                #                             transform=ic.axes.get_transform(ic0.wcs))
                # Differ drawing until changing tab
                ic.changed = True

    # Redefine flux for FIFI-LS cubes
    def fluxMedianAT(self):
        try:
            # Check if FIFI-LS cube
            if self.specCube.instrument == 'FIFI-LS':
                # substitute flux with uncorrected flux divided by the median atmospheric transmission
                at = self.specCube.atran
                atmed = np.nanmedian(at)
                atran = self.specCube.atran*0.+atmed
                self.specCube.atran = atran
                #self.specCube.atran[:] = atmed
                # The uncorrected flux should be interpolated over the corrected flux wavelength
                # grid after applying the baryshift correction ...
                self.computeFluxAtm(atran)
                # self.specCube.flux = self.specCube.uflux/atmed
                # Redraw the spectrum
                for istab in range(len(self.stabs)):
                    if self.stabs.tabText(istab) != 'PSF':
                        self.sci[istab].updateSpectrum(atran=atran)
                # tab with total flux
                self.doZoomAll('new AT')
                # tabs with apertures
                self.onModifiedAperture('new AT')
                # Update image
                self.slideCube('new AT')
            else:
                self.sb.showMessage("This operation is possible with FIFI-LS cubes only", 2000)    
        except:
            self.sb.showMessage("First choose a cube ", 1000)
            
            
    def fluxRefWavAT(self):
        try:
            # Check if FIFI-LS cube
            if self.specCube.instrument == 'FIFI-LS':
                # substitute flux with uncorrected flux divided by the atm transmission at ref wav
                at = self.specCube.atran
                idx = np.argmin(np.abs(self.specCube.wave-self.specCube.l0 * (1+self.specCube.redshift)))
                atran = self.specCube.atran*0.+at[idx]
                #self.specCube.atran[:] = atmed
                # The uncorrected flux should be interpolated over the corrected flux wavelength
                # grid after applying the baryshift correction ...
                self.computeFluxAtm(at[idx])
                self.specCube.atran = atran
                # self.specCube.flux = self.specCube.uflux/atmed
                # Redraw the spectrum
                for istab in range(len(self.stabs)):
                    if self.stabs.tabText(istab) != 'PSF':
                        self.sci[istab].updateSpectrum(atran=atran)
                # tab with total flux
                self.doZoomAll('new AT')
                # tabs with apertures
                self.onModifiedAperture('new AT')
                # Update image
                self.slideCube('new AT')
            else:
                self.sb.showMessage("This operation is possible with FIFI-LS cubes only", 2000)    
        except:
            self.sb.showMessage("First choose a cube ", 1000)
        
    def computeFluxAtm(self, atmed):
        nz, ny, nx = np.shape(self.specCube.uflux)
        x = self.specCube.wave
        xr = x * (1+self.specCube.baryshift)
        for i in range(ny):
            for j in range(nx):
                uflux = self.specCube.uflux[:, i, j]
                offset = 0
                #offset = np.nanmedian(uflux)
                self.specCube.flux[:, i, j] = np.interp(x, xr, (uflux - offset)/ atmed + offset)
        
    def readAtran(self, detchan, order):
        import os
        from astropy.io import fits
        path0, file0 = os.path.split(__file__)
        # print('Path for Atran is ', path0)
        if detchan == 'BLUE':
            file = 'AtranBlue'+str(order)+'.fits.gz'
        else:
            file = 'AtranRed.fits.gz'
        path = path0+'/data/'
        hdl = fits.open(path+file)
        wt = hdl['WAVELENGTH'].data
        atran = hdl['ATRAN'].data
        altitudes = hdl['ALTITUDE'].data
        wvzs = hdl['WVZ'].data
        hdl.close()
        return (wt, atran, altitudes, wvzs)
            
    def fluxNewAT(self):
        try:
            # Check if FIFI-LS cube
            print('instrument is ', self.specCube.instrument)
            if self.specCube.instrument == 'FIFI-LS':
                # Call a dialog showing the altitude and elevation angle in the header and asking for zenithal water vapor
                za = 0.5 * (self.specCube.za[0] + self.specCube.za[1])
                altitude = 0.5 * (self.specCube.altitude[0] + self.specCube.altitude[1])
                wvz = self.getWVZ(altitude,za)
                print('Selected Zenithal Water Vapor: ', wvz)
                # Download a new AT curve and apply it to the uflux
                atrandata = self.readAtran(self.specCube.channel, self.specCube.order)
                wt, atran, altitudes, wvs = atrandata
                imin = np.argmin(np.abs(altitudes-altitude))
                at = atran[imin]
                angle = za * np.pi/180.
                cos_angle = np.cos(angle)
                #depth = 1. / cos_angle  # Flat Earth approximation
                r = 6383.5/50.  # assuming r_earth = 6371 km, altitude = 12.5 km, and 50 km of more stratosphere
                rcos = r * cos_angle
                depth = -rcos + np.sqrt(rcos * rcos + 1 + 2 * r) # Taking into account Earth curvature
                imin = np.argmin(np.abs(wvs-wvz))
                at = at[imin]**depth
                atmed = np.interp(self.specCube.wave, wt , at)
                atran = atmed
                atmed[atmed < 0.3] = np.nan   # Do not correct for too low atmospheric transmission
                self.computeFluxAtm(atmed)
                self.specCube.atran = atran # update atran
                for istab in range(len(self.stabs)):
                    if self.stabs.tabText(istab) != 'PSF':
                        self.sci[istab].updateSpectrum(atran=atran)
                # tab with total flux
                self.doZoomAll('new AT')
                # tabs with apertures
                self.onModifiedAperture('new AT')
                # Update image
                self.slideCube('new AT')
            else:
                self.sb.showMessage("This operation is possible with FIFI-LS cubes only", 2000)    
        except:
            self.sb.showMessage("First choose a cube ", 1000)
        
    def getWVZ(self, alt, za):
        wvzline, okPressed = QInputDialog.getDouble(self, "Altitude: " + str(alt) + 
                                                    " Zenithal Angle: " + str(za),
                                                    "Water vapor at zenith", 5, 1.5, 10., 2)
        if okPressed:
            return wvzline
        else:
            return None

    def reloadFile(self):
        """Reload the file."""        
        try:
            filename = self.specCube.filename
            self.loadFile(filename)
            self.initializeImages()
            self.initializeSpectra()
            self.initializeSlider()
            if self.specCube.instrument in ['GREAT','HI','VLA','ALMA','MUSE','IRAM','CARMA','PCWI']:
                self.specCube.computeExpFromNan()
                #idx = np.isfinite(self.specCube.flux)
                #print('No of bad ', np.sum(~idx))
                if self.specCube.instrument  == 'GREAT':
                    self.slideCube('exp computed')
            self.all = False
            self.fitcont = False
            # Set default number of lines to fit across the cube
            self.abslines = 0
            self.emslines = 0
            # Default to one region
            self.ncells = 1
        except:
            self.sb.showMessage("ERROR: You have to load a file first", 2000)
            return
        
    def newFile(self):
        """Display a new image."""
        # First remove contours
        try:
            self.removeContours()
        except:
            pass
        fd = QFileDialog()
        fd.setLabelText(QFileDialog.Accept, "Import")
        fd.setNameFilters(["Fits Files (*.fit*)", "WXY fits files (*WXY*.fits*)", "All Files (*)"])
        fd.setOptions(QFileDialog.DontUseNativeDialog)
        fd.setViewMode(QFileDialog.List)
        fd.setFileMode(QFileDialog.ExistingFile)
        if (fd.exec()):
            fileName= fd.selectedFiles()
            print('Reading file ', fileName[0])
            # Save the file path for future reference
            self.pathFile, file = os.path.split(fileName[0])
            self.loadFile(fileName[0])
            try:
                self.initializeImages()
                print('images initialized ')
                self.initializeSpectra()
                print('spectra initialized ')
                if self.specCube.instrument in ['GREAT','HI','VLA','ALMA','MUSE','IRAM','CARMA','PCWI']:
                    #print('compute Exp from Nan')
                    self.specCube.computeExpFromNan()
                self.all = False
                self.fitcont = False
                # Set default number of lines to fit across the cube
                self.abslines = 0
                self.emslines = 0
                self.initializeSlider()
                if self.specCube.instrument == 'GREAT':
                    self.slideCube('Exp computed')
            except:
                print('No spectral cube is defined')
                pass
            
    def loadFile(self, infile):
        # Read the spectral cube
        try:
            self.specCube = specCube(infile)
        except:
            self.sb.showMessage("fitsio cannot read this file ", 1000)
            try:
                self.specCube = specCubeAstro(infile)
            except:
                self.sb.showMessage("ERROR: The selected file is not a good spectral cube ", 2000)
                return
        # Delete pre-existing spectral tabs
        try:
            for stab in reversed(range(len(self.sci))):
                self.removeSpecTab(stab, False)
            #print('all spectral tabs removed')
        except:
            pass
        # Delete pre-existing image tabs
        try:
            # Remove tabs, image and histo canvases and disconnect them
            # The removal is done in reversed order to get all the tabs
            for itab in reversed(range(len(self.ici))):
                self.removeTab(itab, False)
        except BaseException:
            pass
        # Update window title (to include object name)
        self.setWindowTitle(self.title + " [ "+self.specCube.objname+" - "+
                                            self.specCube.instrument+" ]")
        # Initialize
        self.tabi = []
        self.ici  = []
        self.ihi  = []
        self.ihcid = []
        self.icid1 = []
        self.icid2 = []
        self.icid3 = []
        self.icid4 = []
        self.stabi = []
        self.sci  = []
        self.scid1 = []
        self.scid2 = []
        self.scid3 = []
        self.scid4 = []
        self.scid5 = []
        self.scid6 = []
        # Open new tabs and display it
        print('Instrument ', self.specCube.instrument)
        if self.specCube.instrument == 'FIFI-LS':
            self.bands = ['Flux','uFlux','Exp']
            self.spectra = ['All','Pix']
        elif self.specCube.instrument in ['GREAT','HI','VLA','ALMA','MUSE','IRAM','CARMA','PCWI']: 
            self.bands = ['Flux']
            self.spectra = ['All','Pix']
        elif self.specCube.instrument == 'FORCAST':
            self.bands = ['Flux','Exp']
            self.spectra = ['All','Pix']            
        elif self.specCube.instrument == 'PACS':
            self.bands = ['Flux','Exp']
            self.spectra = ['All','Pix']
        else:
            self.spectra = []
            self.bands = []
        self.photoApertures = []
        for b in self.bands:
            t, ic, ih, h, c1, c2, c3, c4 = self.addImage(b)
            self.tabi.append(t)
            self.ici.append(ic)
            self.ihi.append(ih)
            self.ihcid.append(h)
            self.icid1.append(c1)
            self.icid2.append(c2)
            self.icid3.append(c3)
            self.icid4.append(c4)
        # Make tab 'Flux' unclosable
        self.itabs.tabBar().setTabButton(0,QTabBar.LeftSide,None)
        self.itabs.tabBar().setTabButton(0,QTabBar.RightSide,None)
        for s in self.spectra:
            t, sc, scid1, scid2, scid3, scid4, scid5, scid6 = self.addSpectrum(s)
            self.stabi.append(t)
            self.sci.append(sc)
            self.scid1.append(scid1)
            self.scid2.append(scid2)
            self.scid3.append(scid3)
            self.scid4.append(scid4)
            self.scid5.append(scid5)
            self.scid6.append(scid6)
        self.stabs.tabBar().setTabButton(0,QTabBar.LeftSide,None)
        self.stabs.tabBar().setTabButton(0,QTabBar.RightSide,None)
        # Make tab 'Pix' unclosable
        self.stabs.tabBar().setTabButton(1,QTabBar.LeftSide,None)
        self.stabs.tabBar().setTabButton(1,QTabBar.RightSide,None)
        # Start selecting Pixel tab
        istab = self.spectra.index('Pix')
        if self.stabs.currentIndex() != istab:
            self.stabs.currentChanged.disconnect()
            self.stabs.setCurrentIndex(istab)  # Pixel tab
            self.stabs.currentChanged.connect(self.onSTabChange)
           
    def initializeImages(self):
        import time
        #t = time.process_time()
        s = self.specCube
        # Compute initial images
        # print('Initialize images')
        print('bands ',self.bands)
        for ima in self.bands:
            #print('ima is ', ima)
            ts = time.process_time()
            ic = self.ici[self.bands.index(ima)]
            if ima == 'Flux':
                image = s.flux[s.n0,:,:]
                if s.instrument == 'GREAT':
                    image *= s.Tb2Jy
            elif ima == 'uFlux':
                image = s.uflux[s.n0,:,:]
            elif ima == 'Exp':
                #print('exposure exists ! ', np.shape(s.exposure))
                print('n0 is ',s.n0)
                image = s.exposure[s.n0,:,:]
                print('image is ', np.shape(image), image)
            elif ima == 'M0':
                self.computeZeroMoment()
                image = self.M0
            else:
                pass
            if (ima == 'Flux') & (s.instrument == 'PCWI'):
                aspect = s.ypixscale/s.pixscale
            else:
                aspect = 1.
            print('aspect in mainwindow is ',aspect)
            # t0 = time.process_time()
            # print('Image prepared in ', t0-ts, 's')
            ic.compute_initial_figure(image=image,wcs=s.wcs,title=ima,cMap=self.colorMap,
                                      cMapDir=self.colorMapDirection,stretch=self.stretchMap,
                                      instrument = s.instrument, aspect=aspect)
            t1 = time.process_time() 
            print('Image displayed in ', t1-ts,' s')
            # print('select output format')
            if ima == 'Exp':
                ic.image.format_cursor_data = lambda z: "{:10.3f} s".format(float(z))
            else:
                ic.image.format_cursor_data = lambda z: "{:10.4f} Jy".format(float(z))
            # print('callback for limits')
            # Callback to propagate axes limit changes among images
            #ic.cid = ic.axes.callbacks.connect('xlim_changed' or 'ylim_changed',
            #                                   self.doZoomAll)
            ic.cid = ic.axes.callbacks.connect('ylim_changed',
                                               self.doZoomAll)
            # print('Define histogram')
            ih = self.ihi[self.bands.index(ima)]
            #ih.changed = False
            try:
                clim = ic.image.get_clim()
            except BaseException:
                clim = [0, 0]
            ih.compute_initial_figure(image=image, xmin=clim[0], xmax=clim[1])
            t2 = time.process_time() 
            print('Histogram computed ', t2-t1, ' s')
            # temporary ...
            x = ic.axes.get_xlim()
            y = ic.axes.get_ylim()
            ic.zoomlimits = [x,y]
            print('zoom limits ', x, y)
        # Re-initialize variables
        self.contours = 'off'
        self.blink = 'off'
        self.slice = 'off'
        self.trimcube = 'off'
        print('variables off ..')
        self.continuum = None
        self.L0 = None
        self.L1 = None
        self.M0 = None
        self.M1 = None
        self.M2 = None
        self.M3 = None
        self.M4 = None
        self.C0 = None
        self.v  = None
        self.sv = None
        # Selectors
        self.PS = None
        self.ES = None
        self.RS = None
        self.LS = None
        print('all cleared ..')
        return image, clim
            
    def initializeSpectra(self):
        print('initializing spectrum ...')
        s = self.specCube
        # Compute initial pixel spectrum
        print('spectra ', self.spectra)
        spectrum = self.spectra[1]
        sc = self.sci[self.spectra.index(spectrum)]
        # Add pixel aperture
        ic0 = self.ici[0]
        x0 = s.nx // 2
        y0 = s.ny // 2
        r0,d0 = ic0.wcs.all_pix2world(x0, y0, 0)
        ws = ic0.pixscale       
        n = len(self.photoApertures)
        # Define pixel aperture
        data = [r0,d0,ws]
        print('aperture at ', data)
        self.photoApertures.append(photoAperture(n,'pixel',data))
        for ic in self.ici:
            x0, y0 = ic.wcs.all_world2pix(r0, d0, 0)
            w = ws / ic.pixscale
            pixel = PixelInteractor(ic.axes, (x0, y0), w)
            ic.photApertures.append(pixel)
            cidap=pixel.mySignal.connect(self.onRemoveAperture)
            ic.photApertureSignal.append(cidap)
            pixel.modSignal.connect(self.onModifiedAperture)
        x0 = s.nx // 2
        y0 = s.ny // 2
        fluxAll = s.flux[:,y0,x0]
        if s.instrument == 'GREAT':
            #print('max flux is ', np.nanmax(fluxAll*s.Tb2Jy))
            spec = Spectrum(s.wave, fluxAll*s.Tb2Jy, instrument=s.instrument,
                            redshift=s.redshift, l0=s.l0, Tb2Jy=s.Tb2Jy, bunit=s.bunit)
        elif s.instrument in ['HI','VLA','ALMA','MUSE','IRAM','CARMA','PCWI']:
            spec = Spectrum(s.wave, fluxAll,instrument=s.instrument,
                            redshift=s.redshift, l0=s.l0)
        elif s.instrument == 'PACS':
            expAll = s.exposure[:, y0, x0]
            spec = Spectrum(s.wave, fluxAll, exposure=expAll, instrument=s.instrument,
                            redshift=s.redshift, l0=s.l0)
        elif s.instrument == 'FORCAST':
            expAll = s.exposure[:, y0, x0]
            efluxAll = s.eflux[:, y0, x0]
            spec = Spectrum(s.wave, fluxAll, eflux=efluxAll,
                            exposure=expAll, instrument=s.instrument,
                            redshift=s.redshift, l0=s.l0)
        elif s.instrument == 'FIFI-LS':
            ufluxAll = s.uflux[:, y0, x0]
            expAll = s.exposure[:, y0, x0]
            efluxAll = s.eflux[:, y0, x0]
            spec = Spectrum(s.wave, fluxAll, eflux=efluxAll, uflux= ufluxAll,
                            exposure=expAll, atran = s.atran, instrument=s.instrument,
                            redshift=s.redshift, baryshift=s.baryshift, l0=s.l0,
                            watran=s.watran, uatran=s.uatran)
            #print('spectrum defined')
        print('compute initial spectrum')
        sc.compute_initial_spectrum(name='Pix', spectrum=spec)
        print('initial spectrum computed')
        self.specZoomlimits = [sc.xlimits, sc.ylimits]
        sc.cid = sc.axes.callbacks.connect('xlim_changed' or 'ylim_changed', self.doZoomSpec)
        print('define span selector')
        sc.span = SpanSelector(sc.axes, self.onSelect, 'horizontal', useblit=True,
                               rectprops=dict(alpha=0.3, facecolor='LightGreen'))
        sc.span.active = False
        #print('s.n0 ', s.n0, ' len(wave) ', len(s.wave))
        if (s.n0 <= 0) or (s.n0 >= (len(s.wave)-2)):
            s.n0 = len(s.wave)//2
        #print('s.n0 updated ', s.n0)
        wave0 = s.wave[s.n0]
        dwave = (s.wave[s.n0+1]-wave0)*0.5
        sc.regionlimits = wave0-dwave,wave0+dwave
        self.slider = None
        self.slicer = None

    def selectSlider(self):
        SD = SlicerDialog()
        if SD.exec_() == QDialog.Accepted:
            option = SD.save()
            if option == 'Channel':
                self.initializeSlider()
            elif option == 'Cube slice':
                self.initializeSlicer()
            elif option == 'None':
                self.removeSliders()
            else:
                pass
            
    def removeSliders(self):
        try:
            self.slicer.disconnect()
            self.slicer = None
        except BaseException:
            pass
        try:
            self.slider.disconnect()
            self.slider = None
        except BaseException:
            pass
        
    def initializeSlicer(self):
        if self.slicer is not None:
            return
        # Number or channels on one side
        ndw = 3
        s = self.specCube
        w0 = s.wave[s.n0]
        dw = s.wave[s.n0+1]-w0
        sc = self.sci[self.stabs.currentIndex()]
        try:
            x = self.slider.x
            if sc.xunit == 'THz':
                c = 299792458.0  # speed of light in m/s
                w0 = c / x * 1.e-6
            else:
                w0 = x            
            self.slider.disconnect()
            self.slider = None
        except BaseException:
            pass
        if sc.xunit == 'THz':
            c = 299792458.0  # speed of light in m/s
            r = w0/dw
            w0 = c / w0 * 1.e-6
            dw = c / dw * 1.e-6 / np.abs(r * r - 0.25)
        self.slicer = SliceInteractor(sc.axes, w0 - ndw * dw, w0 + ndw * dw)
        self.slicer.modSignal.connect(self.updateImages)
        sc.fig.canvas.draw_idle()  

    def initializeSlider(self):
        s = self.specCube
        sc = self.sci[self.stabs.currentIndex()]
        print('no is ', s.n0)
        w0 = s.wave[s.n0]
        dw = s.wave[s.n0+1]-w0
        try:
            xm = 0.5 * (self.slicer.xl + self.slicer.xr)
            if sc.xunit == 'THz':
                c = 299792458.0  # speed of light in m/s
                w0 = c / xm * 1.e-6
            else:
                w0 = xm
            self.slicer.disconnect()
            self.slicer = None
        except BaseException:
            pass
        try:
            x = self.slider.x
            if sc.xunit == 'THz':
                c = 299792458.0  # speed of light in m/s
                w0 = c / x * 1.e-6
            else:
                w0 = x            
            self.slider.disconnect()
            self.slider = None
        except BaseException:
            pass
        if sc.xunit == 'THz':
            c = 299792458.0  # speed of light in m/s
            r = w0/dw
            w0 = c / w0 * 1.e-6
            dw = c / dw * 1.e-6 / np.abs(r * r - 0.25)
        #print('slider at ',w0,' with width ',dw)
        self.slider = SliderInteractor(sc.axes, w0, dw)
        self.slider.modSignal.connect(self.slideCube)

    def slideCube(self, event):
        """Slide over the depth of the cube once the slider moves."""
        # Capture the position of the slider and convert it in position in the cube
        # Here we are just using wavelengths, make this more generic with frequency
        w = self.slider.x
        stab = self.stabs.currentIndex()
        sc = self.sci[stab]
        if sc.xunit == 'THz':
            c = 299792458.0  # speed of light in m/s
            w = c / w * 1.e-6
        n = np.argmin(np.abs(self.specCube.wave - w))
       # Display channel n of the spectral cube
        if self.specCube.instrument in ['GREAT','HI','VLA','ALMA','MUSE','IRAM','CARMA','PCWI']:
            imas = ['Flux']
        elif self.specCube.instrument in ['PACS', 'FORCAST']:
            imas = ['Flux','Exp']
        elif self.specCube.instrument == 'FIFI-LS':
            imas = ['Flux','uFlux','Exp']
        # x,y = self.zoomlimits
        itab = self.itabs.currentIndex()
        ic0 = self.ici[itab]
        for ima in imas:
            if ima == 'Flux':
                if self.specCube.instrument == 'GREAT':
                    #print('update image of great cube')
                    #idx = np.isfinite(self.specCube.flux[n,:,:])
                    #print('NaN are ',np.sum(~idx))
                    image = self.specCube.flux[n,:,:] * self.specCube.Tb2Jy
                else:
                    image = self.specCube.flux[n,:,:]
            elif ima == 'uFlux':
                image = self.specCube.uflux[n,:,:]
            elif ima == 'Exp':
                image = self.specCube.exposure[n,:,:]
            else:
                pass
            ic = self.ici[self.bands.index(ima)]
            # Update only data to go faster
            ic.updateImage(image)
            if ic == ic0:
                ic.fig.canvas.draw_idle()
                ic.changed = False
                # draw_idle seems faster than update ...
                # ic.axes.draw_artist(ic.image)
                # ic.fig.canvas.update()
                # ic.fig.canvas.flush_events()
            else:
                ic.changed = True
            ih = self.ihi[self.bands.index(ima)]
            ih.changed = True

    def switchUnits(self, event):
        """React to switch in units of the spectrum canvas."""
        #print('event is ', event)
        if event == 'switched x unit':
            if self.slider is not None:
                self.sliderSwitchUnits()
            if self.slicer is not None:
                self.slicerSwitchUnits()
            # change units in other tabs
            stab = self.stabs.currentIndex()
            sc0 = self.sci[stab]
            sci = self.sci.copy()
            sci.remove(sc0)
            if self.all == False:
                sci.remove(self.sci[0])
                self.sci[0].xunit = 'THz'
            for sc in sci:
                sc.switchUnits()
        elif event == 'switched y unit':
            print('Switch Jy/K')
            # Change units in other tabs
            stab = self.stabs.currentIndex()
            sc0 = self.sci[stab]
            sci = self.sci.copy()
            sci.remove(sc0)
            if self.all == False:
                sci.remove(self.sci[0])
                self.sci[0].yunit = 'K'
            for sc in sci:
                sc.switchFluxUnits()
            

    def sliderSwitchUnits(self):
        c = 299792458.0  # speed of light in m/s
        r = self.slider.x/self.slider.dx
        self.slider.x = c / self.slider.x * 1.e-6
        self.slider.dx = c / self.slider.dx * 1.e-6 / np.abs(r * r - 0.25)
        # print('dx ', self.slider.dx)
        self.slider.redraw(0)
        
    def slicerSwitchUnits(self):
        c = 299792458.0  # speed of light in m/s
        xl = self.slicer.xl
        xr = self.slicer.xr
        self.slicer.xl = c / xr * 1.e-6
        self.slicer.xr = c / xl * 1.e-6
        self.slicer.redraw()

    def computeAll(self):
        """Compute initial total spectrum."""
        print('Computing total spectrum')
        s = self.specCube
        spectrum = self.spectra[0]
        sc = self.sci[self.spectra.index(spectrum)]
        fluxAll = np.nansum(s.flux, axis=(1,2))
        if s.instrument == 'GREAT':
            spec = Spectrum(s.wave, fluxAll*s.Tb2Jy, instrument=s.instrument,
                            redshift=s.redshift, l0=s.l0, Tb2Jy=s.Tb2Jy, bunit=s.bunit)
        elif s.instrument in ['HI','VLA','ALMA','MUSE','IRAM','CARMA','PCWI']:
            spec = Spectrum(s.wave, fluxAll, instrument=s.instrument,
                            redshift=s.redshift, l0=s.l0)
        elif s.instrument in ['PACS', 'FORCAST']:
            expAll = np.nansum(s.exposure, axis=(1,2))
            spec = Spectrum(s.wave, fluxAll, exposure=expAll,instrument=s.instrument,
                            redshift=s.redshift, l0=s.l0 )
        elif s.instrument == 'FIFI-LS':
            ufluxAll = np.nansum(s.uflux, axis=(1,2))
            expAll = np.nansum(s.exposure, axis=(1,2))
            efluxAll = np.sqrt(np.nansum(s.eflux*s.eflux, axis=(1,2)))
            spec = Spectrum(s.wave, fluxAll, eflux=efluxAll, uflux= ufluxAll,
                            exposure=expAll, atran = s.atran, instrument=s.instrument,
                            redshift=s.redshift, baryshift = s.baryshift, l0=s.l0)
        sc.compute_initial_spectrum(name='All', spectrum=spec)
        self.specZoomlimits = [sc.xlimits,sc.ylimits]
        sc.cid = sc.axes.callbacks.connect('xlim_changed' or 'ylim_changed', self.doZoomSpec)
        # Start the span selector to show only part of the cube
        sc.span = SpanSelector(sc.axes, self.onSelect, 'horizontal', useblit=True,
                               rectprops=dict(alpha=0.3, facecolor='LightGreen'))
        sc.span.active = False
        self.all = True

    def onSelect(self, xmin, xmax):
        """ Consider only a slice of the cube when computing the image """
        if self.slice == 'on':
            # Find indices of the shaded region
            sc = self.sci[self.stabs.currentIndex()]
            self.removeContours()
            # Draw region on spectrum
            self.slicer = SliceInteractor(sc.axes, xmin, xmax)
            self.slicer.modSignal.connect(self.updateImages)
            # Hide span selector
            sc.span.active = False
            # Show the lines
            sc.displayLines = True
            sc.setLinesVisibility(sc.displayLines)
            # Redraw the canvas
            sc.fig.canvas.draw_idle()
            # Toggle on the continuum
            try:
                SI = sc.guess
                SI.showverts = True
                SI.line.set_visible(SI.showverts)
            except BaseException:
                pass
            # Update images (flux, uflux, coverage)
            self.updateImages()
            self.slice = 'off'
        elif self.trimcube == 'on':
            # Find indices of the shaded region
            sc = self.sci[self.spectra.index('All')]
            if sc.xunit == 'THz':
                c = 299792458.0  # speed of light in m/s
                xmin, xmax = c/xmax*1.e-6, c/xmin*1.e-6
            sc.shadeRegion([xmin,xmax],'LightYellow')
            sc.fig.canvas.draw_idle()
            sc.span.active = False
            indmin, indmax = np.searchsorted(self.specCube.wave, (xmin, xmax))
            indmax = min(len(self.specCube.wave) - 1, indmax)
            size = indmax-indmin
            nz, nx, ny = np.shape(self.specCube.flux)
            if size == nx:
                self.sb.showMessage("No cutting needed ", 2000)
            else:
                flags = QMessageBox.Yes 
                flags |= QMessageBox.No
                question = "Do you want to trim the cube to the part selected on the image ?"
                response = QMessageBox.question(self, "Question",
                                                question,
                                                flags)            
                if response == QMessageBox.Yes:
                    self.sb.showMessage("Trimming the cube ", 2000)
                    self.trimCube1D(indmin,indmax)
                    self.saveCube()
                    # Load trimmed cube
                    self.reloadFile()
                elif QMessageBox.No:
                    self.sb.showMessage("Trimming aborted ", 2000)
                else:
                    pass
            self.trimcube = 'off'
            sc.tmpRegion.remove()
            sc.fig.canvas.draw_idle()

    def updateImages(self):
        """Update images once the cursors of the slice move."""
        sc = self.sci[self.spectra.index('Pix')]
        xmin = self.slicer.xl
        xmax = self.slicer.xr
        if sc.xunit == 'THz':
            c = 299792458.0  # speed of light in m/s
            xmin, xmax = c/xmax*1.e-6, c/xmin*1.e-6
        indmin, indmax = np.searchsorted(self.specCube.wave, (xmin, xmax))
        indmax = min(len(self.specCube.wave) - 1, indmax)
        sc.regionlimits = [xmin,xmax]
        if self.specCube.instrument in ['GREAT','HI','VLA','ALMA','MUSE','IRAM','CARMA','PCWI']:
            imas = ['Flux']
        elif self.specCube.instrument in ['PACS','FORCAST']:
            imas = ['Flux','Exp']
        elif self.specCube.instrument == 'FIFI-LS':
            imas = ['Flux','uFlux','Exp']
        # x,y = self.zoomlimits
        for ima in imas:
            ic = self.ici[self.bands.index(ima)]
            #ih = self.ihi[self.bands.index(ima)]
            if ima == 'Flux':
                image = np.nanmean(self.specCube.flux[indmin:indmax,:,:], axis=0)
            elif ima == 'uFlux':
                image = np.nanmean(self.specCube.uflux[indmin:indmax,:,:], axis=0)
            elif ima == 'Exp':
                image = np.nansum(self.specCube.exposure[indmin:indmax,:,:], axis=0)
            else:
                pass
            # Update image
            ic.updateImage(image)
            ic.fig.canvas.draw_idle()            
        
    def removeContours(self):
        """Remove previous contours on image and histogram."""
        for ic in self.ici:
            if ic.contour is not None:
                for coll in ic.contour.collections:
                    coll.remove()
                ic.contour = None
                ic.changed = True
        self.contours = 'off'
        # Redraw current image
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        ic.fig.canvas.draw_idle()
        # Remove contour lines in the histogram
        for ih in self.ihi:
            if len(ih.lev) > 0:
                print('There are ',len(ih.lev),len(ih.levels),' contour levels')
                ih.levSignal.disconnect()
                ih.removeLevels()
        # Update tabContour
        self.tabContour[0] = None

    def doZoomAll(self, event):
        """Propagate limit changes to all images."""
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        if ic.axes == event: # only consider axes on screen (not other tabs)
            #print('limits changed ....')
            self.zoomAll(itab)

    def zoomAll(self, itab):
        """Update total spectrum."""
        #print('zoomall called')
        s = self.specCube
        spectrum = self.spectra.index('All')
        sc = self.sci[spectrum]
        ic = self.ici[itab]
        if ic.toolbar._active == 'ZOOM':
            print('called from toolbar', self.bands[itab])
            ic.toolbar.zoom()  # turn off zoom
        ymin, ymax = ic.axes.get_ylim()
        xmin, xmax = ic.axes.get_xlim()
        ic.zoomlimits = [(xmin, xmax), (ymin, ymax)]
        x = [xmin, xmax, xmin, xmax]
        y = [ymin, ymin, ymax, ymax]
        ra, dec = ic.wcs.all_pix2world(x, y, 0)
        # If not in the flux image, compute values for the flux image
        band = self.bands.index('Flux')
        if itab != band:
            ic0 = self.ici[band]
            x, y = ic0.wcs.all_world2pix(ra, dec, 0)
        x = [np.min(x), np.max(x)]
        y = [np.min(y), np.max(y)]
        # Compute limits for new total spectrum
        x0 = int(np.min(x))
        x1 = int(np.max(x))
        y0 = int(np.min(y))
        y1 = int(np.max(y))
        if x0 < 0: x0=0
        if y0 < 0: y0=0
        if x1 >= s.nx: x1 = s.nx-1
        if y1 >= s.ny: y1 = s.ny-1
        # Set new values in other image tabs
        ici = self.ici.copy()
        ici.remove(ic)
        for ima in ici:
            x, y = ima.wcs.all_world2pix(ra, dec, 0)
            x = [np.min(x), np.max(x)]
            y = [np.min(y), np.max(y)]
            ima.axes.callbacks.disconnect(ima.cid)
            ima.axes.set_xlim(x)
            ima.axes.set_ylim(y)
            ima.zoomlimits = (x, y) 
            ima.changed = True
            ima.cid = ima.axes.callbacks.connect('ylim_changed', self.doZoomAll)
        fluxAll = np.nansum(self.specCube.flux[:, y0:y1, x0:x1], axis=(1, 2))
        if s.instrument in ['GREAT']:
            t2j = self.specCube.Tb2Jy
            sc.updateSpectrum(f=fluxAll*t2j)
        elif s.instrument in ['HI','VLA','ALMA','MUSE','IRAM','CARMA','PCWI']:
            sc.updateSpectrum(f=fluxAll)
        elif s.instrument in ['PACS','FORCAST']:
            expAll = np.nansum(s.exposure[:, y0:y1, x0:x1], axis=(1, 2))
            sc.updateSpectrum(f=fluxAll,exp=expAll)            
        elif self.specCube.instrument == 'FIFI-LS':
            ufluxAll = np.nansum(s.uflux[:, y0:y1, x0:x1], axis=(1, 2))
            expAll = np.nansum(s.exposure[:, y0:y1, x0:x1], axis=(1, 2))
            sc.updateSpectrum(f=fluxAll, uf=ufluxAll, exp=expAll)

    def doZoomSpec(self,event):
        """ In the future impose the same limits to all the spectral tabs """
        stab = self.stabs.currentIndex()
        sc = self.sci[stab]
        if sc.toolbar._active == 'ZOOM':
            sc.toolbar.zoom()  # turn off zoom
        xmin,xmax = sc.axes.get_xlim()
        if sc.xunit == 'THz':
            c = 299792458.0  # speed of light in m/s
            xmin, xmax = c/xmax*1.e-6, c/xmin*1.e-6  # Transform in THz as expected by onSelect            
        sc.xlimits = (xmin,xmax)
        sc.ylimits = sc.axes.get_ylim()
        self.specZoomlimits = [sc.xlimits,sc.ylimits]

    def changeVisibility100(self):
        self.changeVisibility(percent=100)

    def changeVisibility995(self):
        self.changeVisibility(percent=99.5)

    def changeVisibility990(self):
        self.changeVisibility(percent=99.0)

    def changeVisibility980(self):
        self.changeVisibility(percent=98.0)

    def changeVisibility950(self):
        self.changeVisibility(percent=95.0)

    def changeVisibility900(self):
        self.changeVisibility(percent=90.0)

    def changeVisibility800(self):
        self.changeVisibility(percent=80.0)

    def changeVisibility(self, dumb=None, percent=None):
        """ Hide/show the histogram of image intensities """
        #print(percent)
        #print(dumb)   understand why percent is not None is first argument
        try:
            itab = self.itabs.currentIndex()
            ih = self.ihi[itab]
            if percent is None:
                state = ih.isVisible()
                ih.setVisible(not state)
            else:
                ih.changed = True
                ih.setVisible(True)
            if ih.isVisible():
                if ih.changed:
                    ih.axes.clear()
                    ih.update_figure(image=self.ici[itab].oimage, percent=percent)
                    ih.fig.canvas.draw_idle()
                    if percent is None:
                        ih.changed = False
                    #if percent is not None:
                    #    ic = self.ici[itab]
                    #    ic.fig.canvas.draw_idle()
        except:
            self.sb.showMessage("First choose a cube ", 1000)

    def changeColorMap(self):
        """ Change a color map for the images """
        # check DS9 colormaps in http://nbviewer.jupyter.org/gist/adonath/c9a97d2f2d964ae7b9eb
        # stretching and normalization are treated at: http://docs.astropy.org/en/stable/visualization/normalization.html
        # normalization as 'sqrt' for instance ... norm = simple_norm(image,'sqrt') with astropy.visualization
        # Great adaption: https://github.com/glue-viz/ds9norm
        # In the colormap dialog I should add a list of possible stretches (sqrt, log , ...)
        if len(self.ihi) > 0:
            self.CMlist = ['real','gist_heat','afmhot','ds9heat','gist_earth','gist_gray','inferno','ocean','plasma','seismic','jet',
                           'ds9a','ds9b','ds9cool','ds9i8','ds9aips0','ds9rainbow','ds9he']
            self.STlist = ['linear','sqrt','square','log','power','sinh','asinh']
            self.CClist = ['cyan','lime','magenta','red','blue','purple','black','white','yellow']
            self.selectCM = cmDialog(self.CMlist,self.STlist, self.CClist, self.colorMap, self.stretchMap, self.colorContour)
            self.selectCM.list.currentRowChanged.connect(self.updateColorMap)
            self.selectCM.slist.currentRowChanged.connect(self.updateStretchMap)
            self.selectCM.clist.currentRowChanged.connect(self.updateColorContour)
            self.selectCM.clist2.currentRowChanged.connect(self.updateColorContour2)
            self.selectCM.dirSignal.connect(self.reverseColorMap)
            self.selectCM.exec_()
        else:
            self.sb.showMessage("First choose a cube ", 1000)            
            
    def updateColorContour(self, newRow):
        """Update the stretch of the color map."""
        newColor = self.CClist[newRow]
        if newColor != self.colorContour[0]:
            self.colorContour[0] = newColor
            if self.contours == 'on':
                for ic in self.ici:
                    try:
                        for c in ic.contour.collections:
                            c.set_color(self.colorContour[0])
                        ic.fig.canvas.draw_idle()
                    except BaseException:
                        pass
                # Update color level ticks in histogram
                ih = self.ihi[self.tabContour[0]]
                for lev in ih.lev:
                    lev.set_color(self.colorContour[0])
                ih.fig.canvas.draw_idle()

    def updateColorContour2(self, newRow):
        """Update the stretch of the color map."""
        newColor = self.CClist[newRow]
        if newColor != self.colorContour[1]:
            self.colorContour[1] = newColor
            if self.contours2 == 'on':
                for ic in self.ici:
                    try:
                        for c in ic.contour2.collections:
                            c.set_color(self.colorContour[1])
                        ic.fig.canvas.draw_idle()
                    except BaseException:
                        pass
                # Update color level ticks in histogram
                ih = self.ihi[self.tabContour[1]]
                for lev in ih.lev:
                    lev.set_color(self.colorContour[1])
                ih.fig.canvas.draw_idle()

    def updateStretchMap(self, newRow):
        """ Update the stretch of the color map """
        from astropy.visualization import ImageNormalize
        newStretch = self.STlist[newRow]
        if newStretch != self.stretchMap:
            print('new stretch is ', newStretch)
            self.stretchMap = newStretch
            for ic in self.ici:
                ic.stretch = self.stretchMap
                ic.norm = ImageNormalize(vmin=None, vmax=None,
                                         stretch=ic.stretchFunc(ic.stretch))
                ic.image.norm = ic.norm
                ic.changed = True
            itab  = self.itabs.currentIndex()
            ic = self.ici[itab]
            ic.image.autoscale()  # trick to update colorbar (https://github.com/matplotlib/matplotlib/issues/5424)
            ic.refreshImage()
            #ic.fig.canvas.draw_idle()
            #ic.changed = False

    def updateColorMap(self, newRow):
        """ Update the color map of the image tabs """
        newCM = self.CMlist[newRow]
        if newCM != self.colorMap:
            self.colorMap = newCM
            for ic in self.ici:
                ic.colorMap = self.colorMap
                ic.image.set_cmap(ic.colorMap+ic.colorMapDirection)
                ic.changed = True
            itab  = self.itabs.currentIndex()
            ic = self.ici[itab]
            ic.refreshImage()
            #ic.fig.canvas.draw_idle()
            #ic.changed = False
                
    def reverseColorMap(self, reverse):
        """ Reverse color map direction """
        if self.colorMapDirection == "":
            self.colorMapDirection = "_r"
        else:
            self.colorMapDirection = ""

        for ic in self.ici:
            ic.colorMapDirection = self.colorMapDirection
            ic.image.set_cmap(ic.colorMap+ic.colorMapDirection)
            ic.changed = True
        itab  = self.itabs.currentIndex()
        ic = self.ici[itab]
        ic.fig.canvas.draw_idle()
        ic.changed = False

    def onSTabChange(self, stab):
        # This should activate an aperture (put dots on it and/or change color)
        if len(self.stabs) > 1:
            itab  = self.itabs.currentIndex()
            ic = self.ici[itab]
            #nap = len(self.stabs)-1
            nap = len(ic.photApertures)
            istab = self.stabs.currentIndex()
            n = self.nAper()
            tabname = self.stabs.tabText(istab)
            #n = istab-1  # aperture number
            # Activate interactor (toogle on) and disactivate
            for iap in range(nap):
                ap = ic.photApertures[iap]
                ic.photApertureSignal[iap]
                if iap == n:
                    ap.showverts = True
                    #ic.photApertureSignal[iap]=ap.mySignal.connect(self.onRemovePolyAperture)
                else:
                    ap.showverts = False
                    #ap.mySignal.disconnect()
                ap.line.set_visible(ap.showverts)
            ic.fig.canvas.draw_idle()
            #
            if tabname != 'PSF':
                if self.slider is not None:
                    x = self.slider.x
                    dx = self.slider.dx
                    self.slider.disconnect()
                    sc = self.sci[istab]
                    self.slider = SliderInteractor(sc.axes, x, dx)
                    self.slider.modSignal.connect(self.slideCube)
                if self.slicer is not None:
                    xl = self.slicer.xl
                    xr = self.slicer.xr
                    self.slicer.disconnect()
                    sc = self.sci[istab]
                    self.slicer = SliceInteractor(sc.axes, xl, xr)
                    self.slicer.modSignal.connect(self.updateImages)
                tabnames = [self.stabs.tabText(i) for i in range(len(self.stabs))]
                if 'PSF' in tabnames:
                    self.PsfI.showverts = False
                    self.PsfI.line.set_visible(self.PsfI.showverts)   
            else:
                self.PsfI.showverts = True
                self.PsfI.line.set_visible(self.PsfI.showverts)
        # Delayed computation of all tab
        if stab == 0:
            if self.all == False:
                try:
                    self.computeAll()
                except BaseException:
                    pass

    def hresizeSpectrum(self):
        """ Expand spectrum to maximum wavelength range """
        istab = self.stabs.currentIndex()
        sc = self.sci[istab]
        xlim0 = np.min(sc.spectrum.wave)
        xlim1 = np.max(sc.spectrum.wave)
        sc.xlimits=(xlim0,xlim1)
        sc.updateXlim()
        
    def vresizeSpectrum(self):
        """ Expand spectrum to maximum flux range """
        istab = self.stabs.currentIndex()
        sc = self.sci[istab]
        spec = sc.spectrum
        mask = np.isfinite(spec.flux)
        if np.sum(mask) > 10:
            ylim0 = np.nanmin(spec.flux)
            ylim1 = np.nanmax(spec.flux)
        else:
            ylim0, ylim1 = sc.axes.get_ylim()
        if sc.instrument == 'FIFI-LS':
            u0 = np.nanmin(spec.uflux)
            u1 = np.nanmax(spec.uflux)
            if u0 < ylim0: ylim0 = u0
            if u1 > ylim1: ylim1 = u1
        if sc.displayAuxFlux:
            sc.vaxes.set_ylim(np.nanmin(sc.aflux), np.nanmax(sc.aflux))
        # Slightly higher maximum
        sc.ylimits = (ylim0,ylim1*1.1)
        sc.updateYlim()     
        
def main():
    from sospex import __version__
    # QApplication.setStyle('Fusion')
    # Ignore warnings
    warnings.filterwarnings('ignore')
    print('sospex version: ', __version__)
    app = QApplication(sys.argv)
    #myStyle = MyProxyStyle('Motif')    # The proxy style should be based on an existing style,
                                        # like 'Windows', 'Motif', 'Plastique', 'Fusion', ...
    #app.setStyle(myStyle)
    gui = GUI()    
    # Adjust geometry to size of the screen
    screen_resolution = app.desktop().screenGeometry()
    width = screen_resolution.width()
    gui.setGeometry(width*0.005, width*0.005, width*0.99, width*0.5)
    gui.hsplitter.setSizes ([width*0.48,width*0.48])
    # Add an icon for the application
    app.setWindowIcon(QIcon(os.path.join(gui.path0,'icons','sospex.png')))
    app.setApplicationName('SOSPEX')
    app.setApplicationVersion(__version__)
    sys.exit(app.exec_())
