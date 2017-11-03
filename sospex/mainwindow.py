#!/usr/bin/env python
import sys,os
from PyQt5.QtWidgets import (QApplication, QWidget, QMainWindow, QTabWidget, QHBoxLayout,
                             QGroupBox, QVBoxLayout, QSizePolicy, QStatusBar, QSplitter,
                             QToolBar, QAction, QFileDialog)
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import Qt
 
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
        with open(self.path0+'/yellow.stylesheet',"r") as fh:
            self.setStyleSheet(fh.read())
        
        self.initUI()
 
    def initUI(self):
        """ Define the user interface """
        
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # Create main widget
        wid = QWidget(self)
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

        self.imagePanel = QGroupBox("")
        layout = QVBoxLayout()
        
        # Tabs with images        
        self.itabs = QTabWidget()
        self.itabs.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding)
        self.itabs.currentChanged.connect(self.onITabChange)  # things to do when changing tab

        # Add widgets to panel
        layout.addWidget(self.itabs)
        self.imagePanel.setLayout(layout)
        

    def onITabChange(self):
        print('Image tab changed')
        
        
    

    def createSpectralPanel(self):
        """ Panel to plot spectra """

        self.spectralPanel = QGroupBox("")
        layout = QVBoxLayout()

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
        self.spectralPanel.setLayout(layout)

    def createToolbar(self):
        """ Toolbar with main commands """

        # Toolbar definition
        self.tb = QToolBar()
        self.tb.setMovable(True)
        self.tb.setObjectName('toolbar')

        # Actions
        startAction = self.createAction(self.path0+'/icons/start.png','Open new file','Ctrl+s',self.newFile)
        quitAction = self.createAction(self.path0+'/icons/exit.png','Quit program','Ctrl+q',self.fileQuit)

        # Add buttons to the toolbar
        self.tb.addAction(startAction)
        self.tb.addAction(quitAction)

    def createAction(self,icon,text,shortcut,action):
        act = QAction(QIcon(icon),text, self)
        act.setShortcut(shortcut)
        act.triggered.connect(action)
        return act
        
    def fileQuit(self):
        """ Quitting the program """
        self.close()

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
            self.specCube = specCube(fileName[0])
        
    def onSTabChange(self):
        print('Spectral tab changed')
        
        
if __name__ == '__main__':
    app = QApplication(sys.argv)
    gui = GUI()
    # Adjust geometry to size of the screen
    screen_resolution = app.desktop().screenGeometry()
    width = screen_resolution.width()
    gui.setGeometry(100, 100, width*0.9, width*0.5)
    sys.exit(app.exec_())
