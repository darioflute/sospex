#!/usr/bin/env python
import sys,os
import numpy as np
from PyQt5.QtWidgets import (QApplication, QWidget, QMainWindow, QTabWidget, QHBoxLayout,
                             QGroupBox, QVBoxLayout, QSizePolicy, QStatusBar, QSplitter,
                             QToolBar, QAction, QFileDialog,  QTableView, QComboBox, QAbstractItemView,
                             QToolButton, QMessageBox, QPushButton, QInputDialog, QDialog)
from PyQt5.QtGui import QIcon, QStandardItem, QStandardItemModel
from PyQt5.QtCore import Qt, QSize, QTimer, pyqtSignal

import matplotlib
matplotlib.use('Qt5Agg')
from graphics import  NavigationToolbar
from matplotlib.widgets import SpanSelector, PolygonSelector, RectangleSelector, EllipseSelector
from matplotlib.patches import Ellipse, Rectangle, Circle, Ellipse, Polygon
from matplotlib.path import Path

# To find points inside a:
#     Polygon
# #path = Path(polygon)
# path = polygon.get_path()   # path in patch coords
# transform = polygon.get_patch_transform()
# npath = transform.transform_path(path)
# inpoints = spectrum.points[npath.contains_points(spectrum.points)]
# xx,yy = inpoints.T
#
#    Ellipse/Circle, Rectangle/Square
# path = ellipse.get_path()   # path in patch coords
# transform = ellipse.get_patch_transform()   # transform from patch coords to data coords
# npath = transform.transform_path(path)
# inpoints = spectrum.points[npath.contains_points(spectrum.points)]
# xx,yy = inpoints.T
#
# Circle
# path = circle.get_path()
# transform = circle.get_transform()   # Apply radius
# npath = transform.transform_path(path)
# inpoints = spectrum.points[npath.contains_points(spectrum.points)]
# xx,yy = inpoints.T

from apertures import photoAperture

class GUI (QMainWindow):
 
    def __init__(self):
        super().__init__()
        self.title = 'SOSPEX: SOFIA Spectrum Explorer'
        self.left = 10
        self.top = 10
        self.width = 640
        self.height = 480

        # Get the path of the package
        self.path0, file0 = os.path.split(__file__)
        # Define style
        with open(self.path0+'/yellow-stylesheet.css',"r") as fh:
            self.setStyleSheet(fh.read())
        
        self.initUI()
 
    def initUI(self):
        """ Define the user interface """
        
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # Create main widget
        wid = QWidget()
        self.setCentralWidget(wid)

        # Main layout is horizontal
        mainLayout = QHBoxLayout()

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
        from lines import define_lines
        self.Lines = define_lines()

    def createImagePanel(self):
        """ Panel to display images """

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

        
        # Add widgets to panel
        layout.addWidget(self.itabs)
        

    def addSpectrum(self,b):
        from graphics import SpectrumCanvas
        ''' Add a tab with an image '''
        t = QWidget()
        t.layout = QVBoxLayout(t)
        t.setSizePolicy(QSizePolicy.Ignored,QSizePolicy.Ignored) # Avoid expansion
        self.stabs.addTab(t, b)
        sc = SpectrumCanvas(t, width=11, height=10.5, dpi=100)
        #ih.setVisible(False)
        # Toolbar
        toolbar = QToolBar()
        # Add actions to toolbar
        toolbar.addAction(self.sliceAction)
        toolbar.addAction(self.cutAction)
        toolbar.addAction(self.maskAction)
        toolbar.addSeparator()

        
        # Navigation toolbar
        sc.toolbar = NavigationToolbar(sc, self)

        foot = QWidget()
        foot.layout = QHBoxLayout(foot)
        foot.layout.addWidget(toolbar)
        foot.layout.addWidget(sc.toolbar)

        t.layout.addWidget(sc)
        t.layout.addWidget(foot)
        self.stabs.resize(self.stabs.minimumSizeHint())  # Avoid expansion
        # connect image and histogram to  events
        sid1=sc.mpl_connect('button_release_event', self.onDraw2)
        sid2=sc.mpl_connect('scroll_event',self.onWheel2)
        return t,sc,sid1,sid2

    def addImage(self,b):
        from graphics import ImageCanvas, ImageHistoCanvas
        ''' Add a tab with an image '''
        t = QWidget()
        t.layout = QVBoxLayout(t)
        t.setSizePolicy(QSizePolicy.Ignored,QSizePolicy.Ignored) # Avoid expansion
        self.itabs.addTab(t, b)
        ic = ImageCanvas(t, width=11, height=10.5, dpi=100)
        ih = ImageHistoCanvas(t, width=11, height=0.5, dpi=100)
        ih.setVisible(False)
        ic.toolbar = NavigationToolbar(ic, self)

        # Toolbar
        toolbar = QToolBar()
        toolbar.addAction(self.levelsAction)
        toolbar.addAction(self.blinkAction)
        toolbar.addAction(self.contoursAction)
        toolbar.addAction(self.momentAction)
        toolbar.addAction(self.cropAction)
        toolbar.addAction(self.cloudAction)
        toolbar.addSeparator()
        #toolbar.addWidget(self.apertureAction)        

        # Foot
        foot = QWidget()
        foot.layout = QHBoxLayout(foot)
        foot.layout.addWidget(toolbar)
        foot.layout.addWidget(ic.toolbar)
        
        #ic.toolbar.pan('on')
        t.layout.addWidget(ic)
        t.layout.addWidget(ih)
        t.layout.addWidget(foot)
        self.itabs.resize(self.itabs.minimumSizeHint())  # Avoid expansion
        # connect image and histogram to  events
        cidh=ih.mySignal.connect(self.onChangeIntensity)
        cid1=ic.mpl_connect('button_release_event', self.onDraw)
        cid2=ic.mpl_connect('scroll_event',self.onWheel)
        cid3=ic.mpl_connect('motion_notify_event', self.onMotion)
        return t,ic,ih,cidh,cid1,cid2,cid3

    def removeTab(self, itab):
        print('removing tab no ',itab)
        widget = self.itabs.widget(itab)
        if widget is not None:
            widget.deleteLater()
        self.itabs.removeTab(itab)
        # Disconnect and remove canvases
        ima = self.ici[itab]
        his = self.ihi[itab]
        hcid = self.ihcid[itab]
        c1 = self.icid1[itab]
        c2 = self.icid2[itab]
        c3 = self.icid3[itab]
        his.mpl_disconnect(hcid)
        ima.mpl_disconnect(c1)
        ima.mpl_disconnect(c2)
        ima.mpl_disconnect(c3)
        self.ici.remove(ima)
        self.ihi.remove(his)
        self.ihcid.remove(hcid)
        self.icid1.remove(c1)
        self.icid2.remove(c2)
        self.icid3.remove(c3)
        ima = None
        his = None

    def removeSpecTab(self, stab):
        print('removing tab no ',stab)
        widget = self.stabs.widget(stab)
        if widget is not None:
            widget.deleteLater()
        self.stabs.removeTab(stab)
        # Disconnect and remove canvases
        spec = self.sci[stab]
        c1 = self.scid1[stab]
        c2 = self.scid2[stab]
        spec.mpl_disconnect(c1)
        spec.mpl_disconnect(c2)
        self.sci.remove(spec)
        self.scid1.remove(c1)
        self.scid2.remove(c2)
        spec = None

        
    def onITabChange(self, itab):
        ''' When tab changes check if latest update of ellipse are implemented '''
        #print("current index is ", itab)
        if itab < len(self.ici):
            ima = self.ici[itab]
            if len(self.stabs) > 1:
                # Check if vertices are activated correctly
                istab = self.stabs.currentIndex()
                nap = len(self.stabs)-1
                n = istab-1  # aperture number
                # Activate interactor (toogle on) and disactivate
                for iap in range(nap):
                    ap = ima.photApertures[iap]
                    if iap == n:
                        ap.showverts = True
                    else:
                        ap.showverts = False
                    ap.updateMarkers()
                    ap.line.set_visible(ap.showverts)
                ima.changed = True
            if ima.changed:
                ima.fig.canvas.draw_idle()
                ima.changed = False
            if self.blink == 'select':
                # Select 2nd tab and start blinking until blink status changes ...
                self.btab[1] = itab
                self.blink = 'on'
                self.timer.start(1000)
            if self.contours == 'select':
                self.ctab[1] = itab
                self.contours = 'on'
                self.drawContours()
        
    def onChangeIntensity(self, event):
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        ih = self.ihi[itab]
        # apply intensity limits to the relative figure
        ic.image.set_clim(ih.limits)
        ic.fig.canvas.draw_idle()

    def onHelp(self, event):
        import webbrowser
        #, markdown2,io
        # create new html help files with
        # jupyter nbconvert Help.ipynb --to html --template basic
        # otherwise download as HTML from notebook
        
        #html = markdown2.markdown_path(self.path0+'/../README.md')
        #path = os.path.abspath('temp.html')
        #url = 'file://' + path
        #with open(path, 'w') as f:
        #    f.write(html)
        #webbrowser.open(url)

        #output = io.StringIO(html)
        #webbrowser.open('https://github.com/darioflute/sospex/blob/python3/README.md')
        print(self.path0+'/readme.html')
        webbrowser.open('file://'+os.path.abspath(self.path0+'/help/Help.html'))

    def onMotion(self, event):
        """ Update spectrum when moving an aperture on the image """


        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        if ic.toolbar._active == 'PAN':
            return

        # Grab aperture in the flux image to compute the new fluxes
        istab = self.stabs.currentIndex()
        #print('motion ',itab,istab)
        if istab > 0:
            sc = self.sci[istab]
            s = self.specCube
            # I should find a way to know if the aperture has changed
            if itab != 0:
                self.updateAperture()
            aperture = self.ici[0].photApertures[istab-1].aperture
            path = aperture.get_path()
            transform = aperture.get_patch_transform()
            npath = transform.transform_path(path)
            inpoints = s.points[npath.contains_points(s.points)]
            xx,yy = inpoints.T
        
            fluxAll = np.nansum(s.flux[:,yy,xx], axis=1)
            if s.instrument == 'GREAT':
                sc.updateSpectrum(fluxAll)
            elif s.instrument == 'FIFI-LS':
                ufluxAll = np.nansum(s.uflux[:,yy,xx], axis=1)
                expAll = np.nansum(s.exposure[:,yy,xx], axis=1)
                sc.updateSpectrum(fluxAll,uf=ufluxAll,exp=expAll)

    def updateAperture(self):

        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        ic0 = self.ici[0]
        nap = self.stabs.currentIndex()-1

        aper = ic.photApertures[nap]
        aper0 = ic0.photApertures[nap]

        if aper.type == 'Polygon':
            verts = aper.poly.get_xy()
            adverts = np.array([(ic.wcs.all_pix2world(x,y,1)) for (x,y) in verts])                
            verts = [(ic0.wcs.all_world2pix(ra,dec,1)) for (ra,dec) in adverts]
            aper0.poly.set_xy(verts)
        elif aper.type == 'Ellipse' or aper.type == 'Circle':
            x0,y0 = aper.ellipse.center
            w0    = aper.ellipse.width
            h0    = aper.ellipse.height
            angle = aper.ellipse.angle
            ra0,dec0 = ic.wcs.all_pix2world(x0,y0,1)
            ws = w0 * ic.pixscale; hs = h0 * ic.pixscale
            x0,y0 = ic0.wcs.all_world2pix(ra0,dec0,1)
            w0 = ws/ic0.pixscale; h0 = hs/ic0.pixscale
            aper0.ellipse.center = x0,y0
            aper0.ellipse.width = w0
            aper0.ellipse.height = h0
            aper0.ellipse.angle = angle
        if aper.type == 'Rectangle' or aper.type == 'Square':
            x0,y0 = aper.rect.get_xy()
            w0    = aper.rect.get_width()
            h0    = aper.rect.get_height()
            angle = aper.rect.angle
            ra0,dec0 = ic.wcs.all_pix2world(x0,y0,1)
            ws = w0 * ic.pixscale; hs = h0 * ic.pixscale
            x0,y0 = ic0.wcs.all_world2pix(ra0,dec0,1)
            w0 = ws/ic0.pixscale; h0 = hs/ic0.pixscale
            aper0.rect.set_xy((x0,y0))
            aper0.rect.set_width(w0)
            aper0.rect.set_height(h0)
            aper0.rect.angle = angle
        
            
    def onDraw(self,event):

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
        if istab != 0 and len(ici) > 0:
            aper = ic.photApertures[istab-1]
            #apertype = aper.__class__.__name__
            if aper.type == 'Ellipse' or aper.type == 'Circle':
                x0,y0 = aper.ellipse.center
                w0    = aper.ellipse.width
                h0    = aper.ellipse.height
                angle = aper.ellipse.angle
                ra0,dec0 = ic.wcs.all_pix2world(x0,y0,1)
                ws = w0 * ic.pixscale; hs = h0 * ic.pixscale
                for ima in ici:
                    x0,y0 = ima.wcs.all_world2pix(ra0,dec0,1)
                    w0 = ws/ima.pixscale; h0 = hs/ima.pixscale
                    ap = ima.photApertures[istab-1]
                    ap.ellipse.center = x0,y0
                    ap.ellipse.width = w0
                    ap.ellipse.height = h0
                    ap.ellipse.angle = angle
                    ap.updateMarkers()
                    ima.changed = True
            if aper.type == 'Rectangle' or aper.type == 'Square':
                x0,y0 = aper.rect.get_xy()
                w0    = aper.rect.get_width()
                h0    = aper.rect.get_height()
                angle = aper.rect.angle
                ra0,dec0 = ic.wcs.all_pix2world(x0,y0,1)
                ws = w0 * ic.pixscale; hs = h0 * ic.pixscale
                for ima in ici:
                    x0,y0 = ima.wcs.all_world2pix(ra0,dec0,1)
                    w0 = ws/ima.pixscale; h0 = hs/ima.pixscale
                    ap = ima.photApertures[istab-1]
                    ap.rect.set_xy((x0,y0))
                    ap.rect.set_width(w0)
                    ap.rect.set_height(h0)
                    ap.rect.angle = angle
                    ap.updateMarkers()
                    ima.changed = True
            elif aper.type == 'Polygon':
                verts = aper.poly.get_xy()
                adverts = np.array([(ic.wcs.all_pix2world(x,y,1)) for (x,y) in verts])                
                for ima in ici:
                    verts = [(ima.wcs.all_world2pix(ra,dec,1)) for (ra,dec) in adverts]
                    ap = ima.photApertures[istab-1]
                    ap.poly.set_xy(verts)
                    ap.updateMarkers()
                    ima.changed = True

            

    def onWheel(self,event):
        ''' enable zoom with mouse wheel and propagate changes to other tabs '''
        eb = event.button
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
        print('zooming by a factor ',factor)
        new_width = (curr_xlim[1]-curr_xlim[0])*factor*0.5
        new_height= (curr_ylim[1]-curr_ylim[0])*factor*0.5
        x = [curr_x0-new_width,curr_x0+new_width]
        y = [curr_y0-new_height,curr_y0+new_height]
        ic.axes.set_xlim(x)
        ic.axes.set_ylim(y)
        ic.fig.canvas.draw_idle()
        ici = self.ici.copy()
        ici.remove(ic)
        ra,dec = ic.wcs.all_pix2world(x,y,1)
        for ima in ici:
            x,y = ima.wcs.all_world2pix(ra,dec,1)
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
            elif QMessageBox.No:
                self.sb.showMessage("Redshift value unchanged ", 2000)
                sc.spectrum.redshift = self.specCube.redshift
                for annotation in sc.annotations:
                    annotation.remove()
                sc.zannotation.remove()
                sc.drawSpectrum()
                sc.fig.canvas.draw_idle()
            else:
                pass           

        # Deselect pan & zoom options on mouse release
        if sc.toolbar._active == "PAN":
            sc.toolbar.pan()
        if sc.toolbar._active == "ZOOM":
            sc.toolbar.zoom()

            
            
    def onWheel2(self,event):
        """ Wheel moves right/left the slice defined on spectrum """

        sc = self.sci[self.spectra.index('All')]
        #print(event.button)
        if sc.regionlimits is not None:
            eb = event.button
            xmin,xmax = sc.regionlimits
            dx = (xmax-xmin) * 0.5
            
            #w = self.specCube.wave
            #dw = np.mean(w[1:]-w[:-1])        

            # Increment region limits
            if eb == 'up':
                xmin += dx
                xmax += dx
            elif eb == 'down':
                xmin -= dx
                xmax -= dx
            else:
                pass        
            # redraw images
            self.slice = 'on'
            if sc.xunit == 'THz':
                c = 299792458.0  # speed of light in m/s
                xmin, xmax = c/xmax*1.e-6, c/xmin*1.e-6  # Transform in THz as expected by onSelect
            self.onSelect(xmin,xmax)


    def createSpectralPanel(self):
        """ Panel to plot spectra """

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
        
        # Status bar
        self.sb = QStatusBar()
        self.sb.showMessage("Welcome to SOSPEX !", 10000)
        
        # Add widgets to panel
        banner = QWidget()
        banner.layout = QHBoxLayout(banner)
        banner.layout.addWidget(self.tb)
        banner.layout.addWidget(self.sb)

        layout.addWidget(self.stabs)
        layout.addWidget(banner)

        
    def createToolbar(self):
        """ Toolbar with main commands """

        # Toolbar definition
        self.tb = QToolBar()
        self.tb.setMovable(True)
        self.tb.setObjectName('toolbar')

        # Actions
        self.helpAction = self.createAction(self.path0+'/icons/help.png','Help','Ctrl+q',self.onHelp)
        self.quitAction = self.createAction(self.path0+'/icons/exit.png','Quit program','Ctrl+q',self.fileQuit)
        self.startAction = self.createAction(self.path0+'/icons/new.png','Load new observation','Ctrl+s',self.newFile)
        self.levelsAction = self.createAction(self.path0+'/icons/levels.png','Adjust image levels','Ctrl+L',self.changeVisibility)
        self.blink = 'off'
        self.blinkAction = self.createAction(self.path0+'/icons/blink.png','Blink between 2 images','Ctrl+B',self.blinkImages)
        self.momentAction = self.createAction(self.path0+'/icons/map.png','Compute moment 0','Ctrl+m',self.zeroMoment)
        self.contours = 'off'
        self.contoursAction = self.createAction(self.path0+'/icons/contours.png','Overlap contours','Ctrl+c',self.overlapContours)
        self.apertureAction = self.createApertureAction()
        self.cutAction = self.createAction(self.path0+'/icons/cut.png','Cut part of the cube','Ctrl+k',self.cutCube)
        self.cropAction = self.createAction(self.path0+'/icons/crop.png','Crop the cube','Ctrl+K',self.cropCube)
        self.sliceAction = self.createAction(self.path0+'/icons/slice.png','Select a slice of the cube','Ctrl+K',self.sliceCube)
        self.maskAction =  self.createAction(self.path0+'/icons/mask.png','Mask a slice of the cube','Ctrl+m',self.maskCube)
        self.cloudAction = self.createAction(self.path0+'/icons/cloud.png','Download image from cloud','Ctrl+D',self.selectDownloadImage)


        # Add buttons to the toolbar

        self.spacer = QWidget()
        self.spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        
        ##self.tb.addWidget(self.spacer)
        self.tb.addAction(self.startAction)
        self.tb.addAction(self.quitAction)
        self.tb.addAction(self.helpAction)
        self.tb.addWidget(self.apertureAction)        

        #self.tb.addAction(self.levelsAction)
        #self.tb.addAction(self.blinkAction)
        #self.tb.addAction(self.contoursAction)
        #self.tb.addAction(self.momentAction)
        #self.tb.addAction(self.cropAction)
        #self.tb.addAction(self.cloudAction)
        ##self.tb.addSeparator()
        #self.tb.addAction(self.sliceAction)
        #self.tb.addAction(self.cutAction)
        #self.tb.addAction(self.maskAction)



        
    
    def createAction(self,icon,text,shortcut,action):
        act = QAction(QIcon(icon),text, self)
        act.setShortcut(shortcut)
        act.triggered.connect(action)
        return act


    
    def createApertureAction(self):
        """ Create combo box for choosing an aperture """

        self.apertures = [['apertures','square','rectangle'],
                     ['circle','ellipse','polygon']]

        self.model = QStandardItemModel()
        for d in self.apertures:                
            row = []
            for text in d:
                item = QStandardItem(QIcon(self.path0+'/icons/'+text+'.png'),"")
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
        self.apView.setIconSize(QSize(24,24))

        
        apertureAction = QComboBox()
        apertureAction.setToolTip("Choose an aperture\n")
        #apertureAction.SizeAdjustPolicy(QComboBox.AdjustToContentsOnFirstShow)
        apertureAction.SizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
        apertureAction.setIconSize(QSize(24,24))
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


    def onPolySelect(self, verts):
        from apertures import PolygonInteractor
        
        self.PS.set_active(False)
        # 1 vertices in RA,Dec coords
        itab = self.itabs.currentIndex()
        ic0 = self.ici[itab]

        adverts = np.array([(ic0.wcs.all_pix2world(x,y,1)) for (x,y) in verts])
        # Save aperture with vertices in ra,dec coordinates
        n = len(self.photoApertures)
        self.photoApertures.append(photoAperture(n,'polygon',adverts))

        for ic in self.ici:
            # First adjust vertices to astrometry (they are in xy coords)
            verts = [(ic.wcs.all_world2pix(ra,dec,1)) for (ra,dec) in adverts]
            poly  = PolygonInteractor(ic.axes, verts)
            ic.photApertures.append(poly)
        self.PS = None

        self.drawNewSpectrum(n)

    def drawNewSpectrum(self, n):        
        """ Add tab with the flux inside the aperture """

        from specobj import Spectrum

        apname = "Ap{:d}".format(n)
        self.spectra.append(apname)
        t,sc,scid1,scid2 = self.addSpectrum(apname)
        self.stabi.append(t)
        self.sci.append(sc)
        self.scid1.append(scid1)
        self.scid2.append(scid2)

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
            spec = Spectrum(s.wave, fluxAll, instrument=s.instrument, redshift=s.redshift, l0=s.l0 )
        elif s.instrument == 'FIFI-LS':
            ufluxAll = np.nansum(s.uflux[:,yy,xx], axis=1)
            expAll = np.nansum(s.exposure[:,yy,xx], axis=1)
            spec = Spectrum(s.wave, fluxAll, uflux= ufluxAll,
                            exposure=expAll, atran = s.atran, instrument=s.instrument,
                            redshift=s.redshift, baryshift = s.baryshift, l0=s.l0)
        sc.compute_initial_spectrum(spectrum=spec)
        self.specZoomlimits = [sc.xlimits,sc.ylimits]
        sc.cid = sc.axes.callbacks.connect('xlim_changed' and 'ylim_changed', self.doZoomSpec)
        # Start the span selector to show only part of the cube
        sc.span = SpanSelector(sc.axes, self.onSelect, 'horizontal', useblit=True,
                               rectprops=dict(alpha=0.5, facecolor='LightSalmon'), button=1)
        sc.span.set_visible(False)

        # Select new tab
        self.stabs.setCurrentIndex(n+1)

        # Deselect all other apertures (to do)
        for i in range(n):
            pass
            
        
    def onRectSelect(self, eclick, erelease):
        'eclick and erelease are the press and release events'

        from apertures import EllipseInteractor, RectangleInteractor
        
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
        #print("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))
        #print(" The button you used were: %s %s" % (eclick.button, erelease.button))

        if self.selAp == 'square' or self.selAp == 'rectangle':
            #x0 = (x1+x2)*0.5
            #y0 = (y1+y2)*0.5
            x0=x1;y0=y1
            w  = np.abs(x2-x1)
            h  = np.abs(y2-y1)
        else:
            x0 = (x1+x2)*0.5
            y0 = (y1+y2)*0.5
            w  = np.abs(x2-x1)
            h  = np.abs(y2-y1)

        
        itab = self.itabs.currentIndex()
        ic0 = self.ici[itab]
        r0,d0 = ic0.wcs.all_pix2world(x0,y0,1)
        ws = w*ic0.pixscale; hs = h*ic0.pixscale

        
        n = len(self.photoApertures)
        if self.selAp == 'square':
            self.RS.set_active(False)
            self.RS = None
            # Define square
            data = [r0,d0,ws]
            self.photoApertures.append(photoAperture(n,'square',data))
            for ic in self.ici:
                x0,y0 = ic.wcs.all_world2pix(r0,d0,1)
                w = ws/ic.pixscale; h = hs/ic.pixscale
                square = RectangleInteractor(ic.axes, (x0,y0), w)
                ic.photApertures.append(square)
        elif self.selAp == 'rectangle':
            self.RS.set_active(False)
            self.RS = None
            # Define rectangle
            data = [r0,d0,ws,hs]
            self.photoApertures.append(photoAperture(n,'rectangle',data))
            for ic in self.ici:
                x0,y0 = ic.wcs.all_world2pix(r0,d0,1)
                w = ws/ic.pixscale; h = hs/ic.pixscale
                rectangle = RectangleInteractor(ic.axes, (x0,y0), w, h)
                ic.photApertures.append(rectangle)
        elif self.selAp == 'circle':
            self.ES.set_active(False)
            self.ES = None
            pass
            # Define circle
            data = [r0,d0,ws]
            self.photoApertures.append(photoAperture(n,'circle',data))
            for ic in self.ici:
                x0,y0 = ic.wcs.all_world2pix(r0,d0,1)
                w = ws/ic.pixscale; h = hs/ic.pixscale
                circle = EllipseInteractor(ic.axes, (x0,y0), w)
                ic.photApertures.append(circle)
        elif self.selAp == 'ellipse':
            self.ES.set_active(False)
            self.ES = None
            # Define ellipse
            data = [r0,d0,ws,hs]
            self.photoApertures.append(photoAperture(n,'ellipse',data))
            for ic in self.ici:
                x0,y0 = ic.wcs.all_world2pix(r0,d0,1)
                w = ws/ic.pixscale; h = hs/ic.pixscale
                ellipse = EllipseInteractor(ic.axes, (x0,y0), w, h)
                ic.photApertures.append(ellipse)

        self.drawNewSpectrum(n)

    
    def chooseAperture(self, i):
        """ Choosing an aperture """
        index  = self.apView.selectionModel().currentIndex()
        i = index.row()
        j = index.column()
        #print ('new selection: ', i,j, self.apertures[i][j])
        self.selAp = self.apertures[i][j]
        if self.selAp == 'ellipse':
            self.sb.showMessage("You chose an "+self.selAp, 1000)
        elif self.selAp == 'apertures':
            self.sb.showMessage("Choose an aperture shape ", 1000)
        else:
            self.sb.showMessage("You chose a "+self.selAp, 1000)

        if len(self.ici) == 0:
            self.sb.showMessage("Start by opening a new image ", 1000)
            self.apertureAction.setCurrentIndex(0)
            return
            
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        if self.selAp == 'polygon':
            self.PS = PolygonSelector(ic.axes, self.onPolySelect, lineprops=dict(linestyle='-',color='g'),
                                      useblit=True,markerprops=dict(marker='o',mec='g'),vertex_select_radius=15)
        elif self.selAp == 'rectangle':
            self.RS = RectangleSelector(ic.axes, self.onRectSelect,
                                        drawtype='box', useblit=True,
                                        button=[1, 3],  # don't use middle button
                                        minspanx=5, minspany=5,
                                        spancoords='pixels',
                                        rectprops = dict(facecolor='g', edgecolor = 'g',alpha=0.8, fill=False),
                                        lineprops = dict(color='g', linestyle='-',linewidth = 2, alpha=0.8),
                                        interactive=False)
            self.RS.state.add('center')
        elif self.selAp == 'square':
            self.RS = RectangleSelector(ic.axes, self.onRectSelect,
                                        drawtype='box', useblit=True,
                                        button=[1, 3],  # don't use middle button
                                        minspanx=5, minspany=5,
                                        spancoords='pixels',
                                        rectprops = dict(facecolor='g', edgecolor = 'g',alpha=0.8, fill=False),
                                        lineprops = dict(color='g', linestyle='-',linewidth = 2, alpha=0.8),
                                        interactive=False)
            self.RS.state.add('square')
            self.RS.state.add('center')
        elif self.selAp == 'ellipse':
            self.ES = EllipseSelector(ic.axes, self.onRectSelect,
                                      drawtype='line', useblit=True,
                                      button=[1, 3],  # don't use middle button
                                      minspanx=5, minspany=5,
                                      spancoords='pixels',
                                      rectprops = dict(facecolor='g', edgecolor = 'g',alpha=0.8, fill=False),
                                      lineprops = dict(color='g', linestyle='-',linewidth = 2, alpha=0.8),
                                      interactive=False)
            self.ES.state.add('center')
        elif self.selAp == 'circle':
            self.ES = EllipseSelector(ic.axes, self.onRectSelect,
                                      drawtype='line', useblit=True,
                                      button=[1, 3],  # don't use middle button
                                      minspanx=5, minspany=5,
                                      spancoords='pixels',
                                      rectprops = dict(facecolor='g', edgecolor = 'g',alpha=0.8, fill=False),
                                      lineprops = dict(color='g', linestyle='-',linewidth = 2, alpha=0.8),
                                      interactive=False)
            self.ES.state.add('square')
            self.ES.state.add('center')
        
        if self.selAp != 'apertures':
            ic.fig.canvas.draw_idle()
        #put back to the 0-th item
        self.apertureAction.setCurrentIndex(0)

    def selectDownloadImage(self):
        """ Select the image to download """

        selectDI = QInputDialog()
        selectDI.setStyleSheet("* { font-size: 14pt; }")
        selectDI.setOption(QInputDialog.UseListViewForComboBoxItems)
        selectDI.setWindowTitle("Select image to download")
        selectDI.setLabelText("Selection")
        imagelist = ['sdss-u','sdss-g','sdss-r','sdss-i','sdss-z',
                     'panstarrs-g','panstarrs-r','panstarrs-i','panstarrs-z','panstarrs-y',
                     'wise1','wise2','wise3','wise4',
                     'first','sumss']
        selectDI.setComboBoxItems(imagelist)
        select = selectDI.exec_()
        if select == QDialog.Accepted:
            self.downloadImage(selectDI.textValue())

    def downloadImage(self, band):
        """ Download an image covering the cube """
        from cloud import cloudImage

        # Compute center and size of image (in arcmin)
        nz,ny,nx = np.shape(self.specCube.flux)
        lon,lat = self.specCube.wcs.celestial.all_pix2world(ny//2,nx//2, 0)
        xsize = nx * self.specCube.pixscale /60. #size in arcmin
        ysize = ny * self.specCube.pixscale /60. #size in arcmin
        print('center: ',lon,lat,' and size: ',xsize,ysize)

        # Compute center and size (arcmin) of the displayed image 
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        x = ic.axes.get_xlim()
        y = ic.axes.get_ylim()
        ra,dec = ic.wcs.all_pix2world(x,y,1)
        lon = np.mean(ra)
        lat = np.mean(dec)
        xsize = np.abs(ra[0]-ra[1])*np.cos(lat*np.pi/180.)*60.
        ysize = np.abs(dec[0]-dec[1])*60.
        print('center: ',lon,lat,' and size: ',xsize,ysize)
        
        print('Band selected is: ',band)
        self.downloadedImage = cloudImage(lon,lat,xsize,ysize,band)
        print('image downloaded')
        
        # Open tab and display the image
        if self.downloadedImage.data is not None:
            self.bands.append(band)
            t,ic,ih,h,c1,c2,c3 = self.addImage(band)
            self.tabi.append(t)
            self.ici.append(ic)
            self.ihi.append(ih)
            self.ihcid.append(h)
            self.icid1.append(c1)
            self.icid2.append(c2)
            self.icid3.append(c3)
            
            ic.compute_initial_figure(image=self.downloadedImage.data,wcs=self.downloadedImage.wcs,title=band)
            # Callback to propagate axes limit changes among images
            ic.cid = ic.axes.callbacks.connect('xlim_changed' and 'ylim_changed', self.doZoomAll)
            ih = self.ihi[self.bands.index(band)]
            clim = ic.image.get_clim()
            ih.compute_initial_figure(image=self.downloadedImage.data,xmin=clim[0],xmax=clim[1])
            
            # Add existing apertures
            self.addApertures(ic)
            
            # Align with spectral cube
            ic0 = self.ici[0]
            x = ic0.axes.get_xlim()
            y = ic0.axes.get_ylim()
            ra,dec = ic0.wcs.all_pix2world(x,y,1)
            x,y = ic.wcs.all_world2pix(ra,dec,1)            
            ic.axes.set_xlim(x)
            ic.axes.set_ylim(y)
            ic.changed = True
        else:
            self.sb.showMessage("The selected survey does not cover the displayed image", 2000)


    def addApertures(self, ic):
        """ Add apertures already defined on new image """
        from apertures import PolygonInteractor, EllipseInteractor, RectangleInteractor

        ic0 = self.ici[0]
        for aper in ic0.photApertures:
            #apertype = aper.__class__.__name__
            if aper.type == 'Ellipse' or aper.type == 'Circle':
                x0,y0 = aper.ellipse.center
                w0    = aper.ellipse.width
                h0    = aper.ellipse.height
                angle = aper.ellipse.angle
                ra0,dec0 = ic0.wcs.all_pix2world(x0,y0,1)
                ws = w0 * ic0.pixscale; hs = h0 * ic0.pixscale
                # Add ellipse
                x0,y0 = ic.wcs.all_world2pix(ra0,dec0,1)
                w0 = ws/ic.pixscale; h0 = hs/ic.pixscale
                ellipse = EllipseInteractor(ic.axes, (x0,y0),w0,h0,angle)
                ellipse.type = aper.type
                ellipse.showverts = aper.showverts
                ic.photApertures.append(ellipse)
            elif aper.type == 'Rectangle' or aper.type == 'Square':
                x0,y0 = aper.rect.get_xy()
                w0    = aper.rect.get_width()
                h0    = aper.rect.get_height()
                print(type(h0))
                angle = aper.rect.angle
                ra0,dec0 = ic0.wcs.all_pix2world(x0,y0,1)
                ws = w0 * ic0.pixscale; hs = h0 * ic0.pixscale
                # Add rectangle
                x0,y0 = ic.wcs.all_world2pix(ra0,dec0,1)
                w0 = ws/ic.pixscale; h0 = hs/ic.pixscale
                rectangle = RectangleInteractor(ic.axes, (x0,y0),w0,h0,angle)
                rectangle.type = aper.type
                rectangle.showverts = aper.showverts
                rectangle.rect.set_xy((x0,y0))
                rectangle.updateMarkers()
                ic.photApertures.append(rectangle)                
            elif aper.type == 'Polygon':
                verts = aper.poly.get_xy()
                adverts = np.array([(ic0.wcs.all_pix2world(x,y,1)) for (x,y) in verts])
                verts = [(ic.wcs.all_world2pix(ra,dec,1)) for (ra,dec) in adverts]
                # Add polygon
                poly = PolygonInteractor(ic.axes,verts)                
                poly.showverts = aper.showverts
                ic.photApertures.append(poly)
        
        

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

        

    def cutCube(self):
        """ Cut part of the cube """
        self.sb.showMessage("Drag the mouse over the slice of the cube to save ", 2000)

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
        print('Size is ',size,nx,ny)
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
            elif QMessageBox.No:
                self.sb.showMessage("Cropping aborted ", 2000)
            else:
                pass

    def cropCube2D(self,center,size):
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
        # Create a grid of points
        nz,ny,nx = np.shape(self.specCube.flux)
        xi = np.arange(nx); yi = np.arange(ny)
        xi,yi = np.meshgrid(xi,yi)
        self.points = np.array([np.ravel(xi),np.ravel(yi)]).transpose()
        
        

        
    def saveCube(self):
        """ Save a cut/cropped cube """ # TODO
        from astropy.io import fits
        
        # Dialog to save file
        fd = QFileDialog()
        fd.setNameFilters(["Fits Files (*.fits)","All Files (*)"])
        fd.setOptions(QFileDialog.DontUseNativeDialog)
        fd.setViewMode(QFileDialog.List)
        #fd.setFileMode(QFileDialog.ExistingFile)

        if (fd.exec()):
            fileName = fd.selectedFiles()
            print(fileName[0])
            outfile = fileName[0]
        
            # Reusable header
            header = self.specCube.wcs.to_header()
            header.remove('WCSAXES')
            header['CRPIX3'] = (self.specCube.crpix3,'Reference pixel')
            header['CRVAL3'] = (self.specCube.crval3,'Reference pixel value')
            header['CDELT3'] = (self.specCube.cdelt3,'Increment')
            header['NAXIS3'] = (self.specCube.n,'3rd dimension')
            header['INSTRUME'] = (self.specCube.instrument, 'Instrument')
            
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
                header['NAXIS'] = (3,'Number of axis')
                
                # Primary header
                hdu = fits.PrimaryHDU()
                hdu.header.extend(header)
                
                # Extensions
                hdu1 = self.addExtension(self.specCube.flux,'FLUX','Jy',header)
                hdu2 = self.addExtension(self.specCube.eflux,'ERROR','Jy',header)
                hdu3 = self.addExtension(self.specCube.uflux,'UNCORRECTED_FLUX','Jy',header)
                hdu4 = self.addExtension(self.specCube.euflux,'UNCORRECTED_ERROR','Jy',header)
                hdu5 = self.addExtension(self.specCube.wave,'WAVELENGTH','um',None)
                hdu6 = self.addExtension(self.specCube.x,'X',None,None)
                hdu7 = self.addExtension(self.specCube.y,'Y',None,None)
                hdu8 = self.addExtension(self.specCube.atran,'TRANSMISSION',None,None)
                hdu9 = self.addExtension(self.specCube.response,'RESPONSE',None,None)
                hdu10 = self.addExtension(self.specCube.exposure,'EXPOSURE_MAP',None,header)
                hdul = fits.HDUList([hdu, hdu1, hdu2, hdu3, hdu4, hdu5, hdu6, hdu7, hdu8, hdu9, hdu10])            
                #hdul.info()    
                hdul.writeto(outfile,overwrite=True) # clobber true  allows rewriting
                hdul.close()
            elif self.specCube.instrument == 'GREAT':
                header['OBJECT'] = (self.specCube.objname, 'Object Name')
                c = 299792458.0  # speed of light in m/s 
                header['VELO-LSR'] = self.specCube.redshift * c
                header['RESTFREQ'] = self.specCube.header['RESTFREQ']
                header['CUNIT3'] = ('m/s','Velocity unit')
                eta_fss=0.97
                eta_mb =0.67
                calib = 971.
                factor = calib*eta_fss*eta_mb
                temperature = self.specCube.flux / factor  # Transform flux into temperature
                # Primary header
                hdu = fits.PrimaryHDU(temperature)
                header['NAXIS'] = (3,'Number of axis')
                hdu.header.extend(header)
                hdul = fits.HDUList([hdu])
                hdul.writeto(outfile,overwrite=True) # clobber true  allows rewriting
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
        """ Cut part of the cube """
        self.sb.showMessage("Define slice of the cube ", 1000)
        self.slice = 'on'
        sc = self.sci[self.spectra.index('All')]
        sc.span.set_visible(True)


    def maskCube(self):
        """ Mask a slice of the cube """
        self.sb.showMessage("Drag your mouse over the spectrum to mask part of the cube or click over to unmask", 2000)
        
    def zeroMoment(self):
        """ Compute and display zero moment of flux """

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

        band = 'M0'
        # Open tab and display the image
        self.bands.append(band)
        t,ic,ih,h,c1,c2,c3 = self.addImage(band)
        self.tabi.append(t)
        self.ici.append(ic)
        self.ihi.append(ih)
        self.ihcid.append(h)
        self.icid1.append(c1)
        self.icid2.append(c2)
        self.icid3.append(c3)
        
        ic.compute_initial_figure(image=self.M0,wcs=self.specCube.wcs,title=band)
        # Callback to propagate axes limit changes among images
        ic.cid = ic.axes.callbacks.connect('xlim_changed' and 'ylim_changed', self.doZoomAll)
        ih = self.ihi[self.bands.index(band)]
        clim = ic.image.get_clim()
        ih.compute_initial_figure(image=self.M0,xmin=clim[0],xmax=clim[1])

        # Add apertures
        self.addApertures(ic)
        
        # Align with spectral cube
        ic0 = self.ici[0]
        x = ic0.axes.get_xlim()
        y = ic0.axes.get_ylim()
        ra,dec = ic0.wcs.all_pix2world(x,y,1)
        x,y = ic.wcs.all_world2pix(ra,dec,1)            
        ic.axes.set_xlim(x)
        ic.axes.set_ylim(y)
        ic.changed = True

        
    def overlapContours(self):
        """ Compute contours and overlap them on images """

        if self.contours == 'off':
            self.ctab = [self.itabs.currentIndex(),0]
            self.sb.showMessage("Select another tab to overlap the current image contours", 2000)
            self.contours = 'select'
        else:
            self.sb.showMessage('Contours deleted ', 1000)
            self.contours = 'off'
            # Delete contours
            try:
                for coll in self.imageContours.collections:
                    coll.remove()
                ic = self.ici[self.ctab[1]]
                ic.fig.canvas.draw_idle()
            except:
                print('No contours to delete')

    def drawContours(self):
        """ Draw contours of image in ctab[0] over image in ctab[1] """

        ic0 = self.ici[self.ctab[0]]
        ic1 = self.ici[self.ctab[1]]
        image = ic0.oimage # is this needed (we have to save the original image in the canvas)
        minimum = np.nanmin(image)
        imedian = np.nanmedian(image)
        maximum = np.nanmax(image)
        if self.bands[self.ctab[0]] == 'Cov':
            levels = np.arange(minimum,maximum,(maximum-minimum)/8)
        else:
            sdev = np.nanstd(image-imedian)
            levels = imedian + np.array([1,2,3,5,10,15,20]) * sdev
            #levels  = np.arange(imedian+sdev,maximum,sdev)
        self.imageContours = ic1.axes.contour(image,levels, colors='cyan',transform=ic1.axes.get_transform(ic0.wcs))
        ic1.fig.canvas.draw_idle() # this shouldn't be necessary
        print('Overlap contours of ',self.bands[self.ctab[0]],' over image ',self.bands[self.ctab[1]])
        
        
    def newFile(self):
        """ Display a new image """

        from specobj import specCube, Spectrum

        fd = QFileDialog()
        fd.setNameFilters(["Fits Files (*.fits)","All Files (*)"])
        fd.setOptions(QFileDialog.DontUseNativeDialog)
        fd.setViewMode(QFileDialog.List)
        fd.setFileMode(QFileDialog.ExistingFile)

        if (fd.exec()):
            fileName= fd.selectedFiles()
            print(fileName[0])
            # Ask to save analysis (cubes, spectra) TODO
            # Read the spectral cube
            self.specCube = specCube(fileName[0])
            # Delete pre-existing image tabs
            try:
                # Remove tabs, image and histo canvases and disconnect them
                # The removal is done in reversed order to get all the tabs
                for itab in reversed(range(len(self.ici))):
                    self.removeTab(itab)
            except:
                pass
            # Delete spectral tabs
            try:
                for stab in reversed(range(len(self.sci))):
                    self.removeSpecTab(stab)
            except:
                pass

            # Initialize
            self.tabi = []
            self.ici  = []
            self.ihi  = []
            self.ihcid = []
            self.icid1 = []
            self.icid2 = []
            self.icid3 = []
            self.stabi = []
            self.sci  = []
            self.scid1 = []
            self.scid2 = []

            # Open new tabs and display it
            if self.specCube.instrument == 'FIFI-LS':
                self.bands = ['Flux','uFlux','Exp']
                self.spectra = ['All']
            elif self.specCube.instrument == 'GREAT':
                self.bands = ['Flux']
                self.spectra = ['All']
            else:
                self.spectra = []
                self.bands = []
            print ("bands are ", self.bands)
            self.photoApertures = []
            for b in self.bands:
                t,ic,ih,h,c1,c2,c3 = self.addImage(b)
                self.tabi.append(t)
                self.ici.append(ic)
                self.ihi.append(ih)
                self.ihcid.append(h)
                self.icid1.append(c1)
                self.icid2.append(c2)
                self.icid3.append(c3)
            for s in self.spectra:
                t,sc,scid1,scid2 = self.addSpectrum(s)
                self.stabi.append(t)
                self.sci.append(sc)
                self.scid1.append(scid1)
                self.scid2.append(scid2)
            # Compute initial images
            for ima in self.bands:
                ic = self.ici[self.bands.index(ima)]
                if ima == 'Flux':
                    image = np.nanmedian(self.specCube.flux, axis=0)
                elif ima == 'uFlux':
                    image = np.nanmedian(self.specCube.uflux, axis=0)
                elif ima == 'Exp':
                    image = np.nansum(self.specCube.exposure, axis=0)
                else:
                    pass
                print('size of image is ',np.shape(image))
                ic.compute_initial_figure(image=image,wcs=self.specCube.wcs,title=ima)
                # Callback to propagate axes limit changes among images
                ic.cid = ic.axes.callbacks.connect('xlim_changed' and 'ylim_changed', self.doZoomAll)
                ih = self.ihi[self.bands.index(ima)]
                clim = ic.image.get_clim()
                ih.compute_initial_figure(image=image,xmin=clim[0],xmax=clim[1])
                x = ic.axes.get_xlim()
                y = ic.axes.get_ylim()
                self.zoomlimits = [x,y]
            # Compute initial spectra
            spectrum = self.spectra[0]
            sc = self.sci[self.spectra.index(spectrum)]
            fluxAll = np.nansum(self.specCube.flux, axis=(1,2))
            s = self.specCube
            if s.instrument == 'GREAT':
                spec = Spectrum(s.wave, fluxAll, instrument=s.instrument, redshift=s.redshift, l0=s.l0 )
            elif self.specCube.instrument == 'FIFI-LS':
                ufluxAll = np.nansum(s.uflux, axis=(1,2))
                expAll = np.nansum(s.exposure, axis=(1,2))
                spec = Spectrum(s.wave, fluxAll, uflux= ufluxAll,
                                exposure=expAll, atran = s.atran, instrument=s.instrument,
                                redshift=s.redshift, baryshift = s.baryshift, l0=s.l0)
            print("Compute initial spectrum")
            sc.compute_initial_spectrum(spectrum=spec)
            self.specZoomlimits = [sc.xlimits,sc.ylimits]
            sc.cid = sc.axes.callbacks.connect('xlim_changed' and 'ylim_changed', self.doZoomSpec)
            # Start the span selector to show only part of the cube
            sc.span = SpanSelector(sc.axes, self.onSelect, 'horizontal', useblit=True,
                                   rectprops=dict(alpha=0.5, facecolor='LightSalmon'), button=1)
            sc.span.set_visible(False)
                
            # Re-initiate variables
            self.contours = 'off'
            self.blink = 'off'
            self.slice = 'off'

    def onSelect(self, xmin, xmax):
        """ Consider only a slice of the cube when computing the image """

        if self.slice == 'on':
            # Find indices of the shaded region
            print('xmin, xmax ',xmin,xmax)
            sc = self.sci[self.spectra.index('All')]
            if sc.xunit == 'THz':
                c = 299792458.0  # speed of light in m/s
                xmin, xmax = c/xmax*1.e-6, c/xmin*1.e-6

            print('xmin, xmax ',xmin,xmax)
            indmin, indmax = np.searchsorted(self.specCube.wave, (xmin, xmax))
            indmax = min(len(self.specCube.wave) - 1, indmax)
            print('indmin, indmax', indmin,indmax)
            sc.regionlimits = [xmin,xmax]


            # Draw region on spectrum (All) and hide span selector
            sc.shadeSpectrum()
            sc.fig.canvas.draw_idle()
            sc.span.set_visible(False)

            # Update images (flux, uflux, coverage)
            if self.specCube.instrument == 'GREAT':
                imas = ['Flux']
            elif self.specCube.instrument == 'FIFI-LS':
                imas = ['Flux','uFlux','Exp']
            
            x,y = self.zoomlimits
            for ima in imas:
                ic = self.ici[self.bands.index(ima)]
                ih = self.ihi[self.bands.index(ima)]
                if ima == 'Flux':
                    image = np.nanmedian(self.specCube.flux[indmin:indmax,:,:], axis=0)
                elif ima == 'uFlux':
                    image = np.nanmedian(self.specCube.uflux[indmin:indmax,:,:], axis=0)
                elif ima == 'Exp':
                    image = np.nansum(self.specCube.exposure[indmin:indmax,:,:], axis=0)
                else:
                    pass
                ic.showImage(image)
                # Set image limits to pre-existing values
                ic.axes.set_xlim(x)
                ic.axes.set_ylim(y)
                ic.changed = True
                # Update histogram
                clim = ic.image.get_clim()
                ih.axes.clear()
                ih.compute_initial_figure(image=image,xmin=clim[0],xmax=clim[1])
                ih.fig.canvas.draw_idle()
            self.slice = 'off'
                        
            
    def doZoomAll(self, event):
        ''' propagate limit changes to all images '''
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        if ic.axes == event: # only consider axes on screen (not other tabs)
            self.zoomAll(itab)

    def zoomAll(self, itab):


        # Update total spectrum
        s = self.specCube
        spectrum = self.spectra.index('All')
        sc = self.sci[spectrum]
        
        ic = self.ici[itab]
        if ic.toolbar._active == 'ZOOM':
            ic.toolbar.zoom()  # turn off zoom
        x = ic.axes.get_xlim()
        y = ic.axes.get_ylim()
        ra,dec = ic.wcs.all_pix2world(x,y,1)

        # If not in the flux image, compute values for the flux image
        band = self.bands.index('Flux')
        if itab != band:
            ic = self.ici[band]
            x,y = ic.wcs.all_world2pix(ra,dec,1)            
        self.zoomlimits = [x,y]

        # Compute limits for new total spectrum
        x0 = int(np.min(x)); x1 = int(np.max(x))
        y0 = int(np.min(y)); y1 = int(np.max(y))
        if x0 < 0: x0=0
        if y0 < 0: y0=0
        if x1 >= s.nx: x1 = s.nx-1
        if y1 >= s.ny: y1 = s.ny-1

        # Set new values in other image tabs
        ici = self.ici.copy()
        ici.remove(ic)
        for ima in ici:
            x,y = ima.wcs.all_world2pix(ra,dec,1)
            ima.axes.callbacks.disconnect(ima.cid)
            ima.axes.set_xlim(x)
            ima.axes.set_ylim(y)
            ima.changed = True
            ima.cid = ima.axes.callbacks.connect('xlim_changed' and 'ylim_changed', self.doZoomAll)

        print('itab',itab,'band',band,'limits ',x,y)

        
        fluxAll = np.nansum(self.specCube.flux[:,y0:y1,x0:x1], axis=(1,2))
        if s.instrument == 'GREAT':
            sc.updateSpectrum(fluxAll)
        elif self.specCube.instrument == 'FIFI-LS':
            ufluxAll = np.nansum(s.uflux[:,y0:y1,x0:x1], axis=(1,2))
            expAll = np.nansum(s.exposure[:,y0:y1,x0:x1], axis=(1,2))
            sc.updateSpectrum(fluxAll,uf=ufluxAll,exp=expAll)
            
            
    def doZoomSpec(self,event):
        """ In the future impose the same limits to all the spectral tabs """
        stab = self.itabs.currentIndex()
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

        
    def changeVisibility(self):
        """ Hide/show the histogram of image intensities """
        try:
            itab = self.itabs.currentIndex()
            ih = self.ihi[itab]
            state = ih.isVisible()
            for ih in self.ihi:
                ih.setVisible(not state)
        except:
            self.sb.showMessage("First choose a cube (press arrow) ", 1000)
            

    def onSTabChange(self, stab):
        print('Spectral tab changed to ', stab)

        # This should activate an aperture (put dots on it and/or change color)
        if len(self.stabs) > 1:
            itab  = self.itabs.currentIndex()
            ic = self.ici[itab]
            nap = len(self.stabs)-1
            istab = self.stabs.currentIndex()
            n = istab-1  # aperture number
            # Activate interactor (toogle on) and disactivate
            for iap in range(nap):
                ap = ic.photApertures[iap]
                if iap == n:
                    ap.showverts = True
                else:
                    ap.showverts = False
                ap.line.set_visible(ap.showverts)
            ic.fig.canvas.draw_idle()
            
        

        
if __name__ == '__main__':
    app = QApplication(sys.argv)
    gui = GUI()
    # Adjust geometry to size of the screen
    screen_resolution = app.desktop().screenGeometry()
    width = screen_resolution.width()
    gui.setGeometry(width*0.025, 0, width*0.95, width*0.5)
    gui.hsplitter.setSizes ([width*0.38,width*0.5])
    sys.exit(app.exec_())
