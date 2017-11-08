#!/usr/bin/env python
import sys,os
import numpy as np
from PyQt5.QtWidgets import (QApplication, QWidget, QMainWindow, QTabWidget, QHBoxLayout,
                             QGroupBox, QVBoxLayout, QSizePolicy, QStatusBar, QSplitter,
                             QToolBar, QAction, QFileDialog, QTableView, QComboBox, QAbstractItemView,
                             QToolButton)
from PyQt5.QtGui import QIcon, QStandardItem, QStandardItemModel
from PyQt5.QtCore import Qt, QSize, pyqtSignal

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

 
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
        hsplitter = QSplitter(Qt.Horizontal)
        
        # Create main panels
        self.createImagePanel()
        self.createSpectralPanel()

        # Add panels to splitter
        hsplitter.addWidget(self.imagePanel)
        hsplitter.addWidget(self.spectralPanel)

        # Add panels to main layout
        mainLayout.addWidget(hsplitter)
        wid.setLayout(mainLayout)
        self.show()


    def createImagePanel(self):
        """ Panel to display images """

        #self.imagePanel = QGroupBox("")
        self.imagePanel = QWidget()
        layout = QVBoxLayout(self.imagePanel)
        
        # Tabs with images        
        self.itabs = QTabWidget()
        self.itabs.setTabsClosable(True)
        self.itabs.tabCloseRequested.connect(self.removeTab)
        self.itabs.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding)
        self.itabs.currentChanged.connect(self.onITabChange)  # things to do when changing tab
        self.tabi = []
        self.ici  = []
        self.ihi  = []
        self.ihcid = []
        self.icid1 = []
        self.icid2 = []
        
        # Add widgets to panel
        layout.addWidget(self.itabs)
        #self.imagePanel.setLayout(layout)
        

    def addImage(self,b):
        from graphics import ImageCanvas, ImageHistoCanvas
        ''' Add a tab with an image '''
        t = QWidget()
        t.layout = QVBoxLayout(t)
        self.itabs.addTab(t, b)
        ic = ImageCanvas(t, width=11, height=10.5, dpi=100)
        ih = ImageHistoCanvas(t, width=11, height=0.5, dpi=100)
        cidh=ih.mySignal.connect(self.onChangeIntensity)
        #ih.setVisible(False)
        ic.toolbar = NavigationToolbar(ic, self)
        #ic.toolbar.pan('on')
        t.layout.addWidget(ic)
        t.layout.addWidget(ih)
        t.layout.addWidget(ic.toolbar)
        #t.setLayout(t.layout)
        # connect image to draw events
        cid1=ic.mpl_connect('button_release_event', self.onDraw)
        #ih.mpl_connect('button_release_event', self.onChangeIntensity)
        cid2=ic.mpl_connect('scroll_event',self.onWheel)
        return t,ic,ih,cidh,cid1,cid2

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
        his.mpl_disconnect(hcid)
        ima.mpl_disconnect(c1)
        ima.mpl_disconnect(c2)
        self.ici.remove(ima)
        self.ihi.remove(his)
        self.ihcid.remove(hcid)
        self.icid1.remove(c1)
        self.icid2.remove(c2)
        ima = None
        his = None

        
    def onITabChange(self, itab):
        ''' When tab changes check if latest update of ellipse are implemented '''
        #print("current index is ", itab)
        if itab < len(self.ici):
            ima = self.ici[itab]
            if ima.changed:
                #canvas = ima.arcell[0].figure.canvas
                #canvas.draw_idle()
                ima.fig.canvas.draw_idle()
                ima.changed = False
            #if self.blink == 'select':
            #    # Select 2nd tab and start blinking until blink status changes ...
            #    self.btab[1] = itab
            #    self.blink = 'on'
            #    self.timer.start(1000)

        
    def onChangeIntensity(self, event):
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        ih = self.ihi[itab]
        # apply intensity limits to the relative figure
        ic.image.set_clim(ih.limits)
        ic.fig.canvas.draw_idle()

    def onDraw(self,event):
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        print('Index is ', ic)

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
        new_width = (curr_xlim[1]-curr_xlim[0])*factor*0.5
        new_height= (curr_ylim[1]-curr_ylim[0])*factor*0.5
        x = [curr_x0-new_width,curr_x0+new_width]
        y = [curr_y0-new_height,curr_y0+new_height]
        ic.axes.set_xlim(x)
        ic.axes.set_ylim(y)
        ic.fig.canvas.draw_idle()
        ici = self.ici.copy()
        ici.remove(ic)
        ra,dec = ic.wcsn.all_pix2world(x,y,1)
        for ima in ici:
            x,y = ima.wcs.all_world2pix(ra,dec,1)
            ima.axes.set_xlim(x)
            ima.axes.set_ylim(y)
            ima.changed = True


    def createSpectralPanel(self):
        """ Panel to plot spectra """

        #self.spectralPanel = QGroupBox("")
        self.spectralPanel = QWidget()
        layout = QVBoxLayout(self.spectralPanel)

        # Toolbar
        self.createToolbar()

        # Tabs with plots        
        self.stabs = QTabWidget()
        self.stabs.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding)
        self.stabs.currentChanged.connect(self.onSTabChange)  # things to do when changing tab

        # Status bar
        self.sb = QStatusBar()
        self.sb.showMessage("Welcome to SOSPEX !", 10000)
        
        # Add widgets to panel
        layout.addWidget(self.tb)
        layout.addWidget(self.stabs)
        layout.addWidget(self.sb)

        
    def createToolbar(self):
        """ Toolbar with main commands """

        # Toolbar definition
        self.tb = QToolBar()
        self.tb.setMovable(True)
        self.tb.setObjectName('toolbar')

        # Actions
        startAction = self.createAction(self.path0+'/icons/new.png','Load new observation','Ctrl+s',self.newFile)
        levelsAction = self.createAction(self.path0+'/icons/levels.png','Adjust image levels','Ctrl+L',self.changeVisibility)
        self.apertureAction = self.createApertureAction()
        cutAction = self.createAction(self.path0+'/icons/cut.png','Cut part of the cube','Ctrl+k',self.cutCube)
        cropAction = self.createAction(self.path0+'/icons/crop.png','Crop the cube','Ctrl+K',self.cropCube)
        sliceAction = self.createAction(self.path0+'/icons/slice.png','Select a slice of the cube','Ctrl+K',self.sliceCube)
        quitAction = self.createAction(self.path0+'/icons/exit.png','Quit program','Ctrl+q',self.fileQuit)

        # Add buttons to the toolbar

        self.spacer = QWidget()
        self.spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        
        self.tb.addWidget(self.spacer)
        self.tb.addAction(startAction)
        self.tb.addAction(levelsAction)
        #self.tb.addSeparator()
        self.tb.addWidget(self.apertureAction)        
        self.tb.addAction(sliceAction)
        self.tb.addAction(cutAction)
        self.tb.addAction(cropAction)
        self.tb.addAction(quitAction)



        
    
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
                item = QStandardItem(QIcon('icons/'+text+'.png'),"")
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


    def chooseAperture(self, i):
        """ Choosing an aperture """
        index  = self.apView.selectionModel().currentIndex()
        i = index.row()
        j = index.column()
        #print ('new selection: ', i,j, self.apertures[i][j])
        if self.apertures[i][j] == 'ellipse':
            self.sb.showMessage("You chose an "+self.apertures[i][j], 1000)
        elif self.apertures[i][j] == 'apertures':
            self.sb.showMessage("Choose an aperture shape ", 1000)
        else:
            self.sb.showMessage("You chose a "+self.apertures[i][j], 1000)
        #put back to the 0-th item
        self.apertureAction.setCurrentIndex(0)


        
    def fileQuit(self):
        """ Quitting the program """
        self.close()


    def cutCube(self):
        """ Cut part of the cube """
        self.sb.showMessage("Cut part of the cube ", 1000)

    def cropCube(self):
        """ Cut part of the cube """
        self.sb.showMessage("Crop the cube ", 1000)

    def sliceCube(self):
        """ Cut part of the cube """
        self.sb.showMessage("Define slice of the cube ", 1000)

        
    def newFile(self):
        """ Display a new image """

        from specobj import specCube 

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
            # Delete pre-existing tabs
            try:
                # Remove tabs, image and histo canvases and disconnect them
                # The removal is done in reversed order to get all the tabs
                for itab in reversed(range(len(self.ici))):
                    self.removeTab(itab)
            except:
                pass    
            # Open new tabs and display it
            if self.specCube.instrument == 'FIFI-LS':
                self.bands = ['Flux','uFlux','Exp']
            elif self.specCube.instrument == 'GREAT':
                self.bands = ['Flux']
            else:
                self.bands = []
            print ("bands are ", self.bands)
            for b in self.bands:
                t,ic,ih,h,c1,c2 = self.addImage(b)
                self.tabi.append(t)
                self.ici.append(ic)
                self.ihi.append(ih)
                self.ihcid.append(h)
                self.icid1.append(c1)
                self.icid2.append(c2)
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
                ic.cid = ic.axes.callbacks.connect('xlim_changed' and 'ylim_changed', self.zoomAll)
                ih = self.ihi[self.bands.index(ima)]
                clim = ic.image.get_clim()
                ih.compute_initial_figure(image=image,xmin=clim[0],xmax=clim[1])
            self.imagePanel.update()

    def zoomAll(self, event):
        ''' propagate limit changes to all images '''
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        if ic.axes == event: # only consider axes on screen (not other tabs)
            if ic.toolbar._active == 'ZOOM':
                ic.toolbar.zoom()  # turn off zoom
            x = ic.axes.get_xlim()
            y = ic.axes.get_ylim()
            ra,dec = ic.wcs.all_pix2world(x,y,1)
            ici = self.ici.copy()
            ici.remove(ic)
            for ima in ici:
                x,y = ima.wcs.all_world2pix(ra,dec,1)
                ima.axes.set_xlim(x)
                ima.axes.set_ylim(y)
                ima.changed = True

    def changeVisibility(self):
        """ Hide/show the histogram of image intensities """
        try:
            itab = self.itabs.currentIndex()
            ih = self.ihi[itab]
            state = ih.isVisible()
            for ih in self.ihi:
                ih.setVisible(not state)
            #ic = self.ici[itab]
            #ici = self.ici.copy()
            #ici.remove(ic)
            #for ima in ici:
            #    ima.changed = True
        except:
            self.sb.showMessage("First choose a cube (press arrow) ", 1000)
            

    def onSTabChange(self):
        print('Spectral tab changed')
        
        
if __name__ == '__main__':
    app = QApplication(sys.argv)
    gui = GUI()
    # Adjust geometry to size of the screen
    screen_resolution = app.desktop().screenGeometry()
    width = screen_resolution.width()
    gui.setGeometry(100, 100, width*0.9, width*0.4)
    sys.exit(app.exec_())
