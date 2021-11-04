from PyQt5.QtWidgets import (QDialog, QPushButton, QGroupBox, QHBoxLayout,                              QVBoxLayout, QGridLayout, QRadioButton, QLabel,                              QButtonGroup, QCheckBox, QListWidget, QSizePolicy,                             QListWidgetItem, QLineEdit, QDialogButtonBox, QFormLayout)from PyQt5.QtCore import QSize, pyqtSignalfrom PyQt5.QtGui import QIcon# Dialogsclass ApertureParams(QDialog):    """Simple dialog to enter center and radius of a circular aperture"""    def __init__(self, parent=None):        super().__init__(parent)        self.setupUI()            def setupUI(self):        self.createFormGroupBox()        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)        buttonBox.accepted.connect(self.OK)        buttonBox.rejected.connect(self.Cancel)        mainLayout = QVBoxLayout()        mainLayout.addWidget(self.formGroupBox)        mainLayout.addWidget(buttonBox)        self.setLayout(mainLayout)        self.setWindowTitle('Aperture parameters')        self.resize(300,100)            def createFormGroupBox(self):        # Two editable fields for center and radius        self.formGroupBox = QGroupBox("Aperture parameters")        layout = QFormLayout()        self.qcenter = QLineEdit()        self.qcenter.setFixedWidth(200)        self.qcenter.setPlaceholderText("hh:mm:ss.ss +dd:mm:ss.s")        layout.addRow(QLabel("Center:"), self.qcenter)        self.qradius = QLineEdit()        self.qradius.setFixedWidth(200)        self.qradius.setPlaceholderText("Aperture radius in arcsec")        layout.addRow(QLabel("Radius:"), self.qradius)        self.formGroupBox.setLayout(layout)            def OK(self):        self.done(1)            def Cancel(self):        self.done(0)            def save(self):        try:            radius = float(self.qradius.text())        except:            radius = -1  # impossible        position = self.qcenter.text()        # Transform sky coords into ra, dec        from astropy.coordinates import SkyCoord        import astropy.units as u        try:            c = SkyCoord(position, unit=(u.hourangle, u.deg), frame='icrs')            ra = c.ra            dec = c.dec        except:            ra = 0            dec = 0        return ra, dec, radiusclass ContParams(QDialog):    """ Dialog window to define parameters of the continuum fit """        def __init__(self, k, parent=None):        super().__init__(parent)        if k == 1:            self.k = 0        elif k == 5:            self.k = 1        elif k == 9:            self.k = 2        else:            self.k = 0                    self.setupUI()    def setupUI(self):                self.function = self.createGroup('Continuum function', ['Constant', 'Slope'])        self.boundary = self.createGroup('Continuum boundary', ['None', 'Non negative'])        self.kernel   = self.createGroup('Kernel pixels', ['1', '5', '9'], default=self.k)        self.regions = self.createGroup('No of regions', ['16', '32' ,'64', '128', '256', '512'])        self.emlines = self.createGroup('No of emission lines', ['0', '1', '2'])        self.ablines = self.createGroup('No of absorption lines', ['0'])        self.models = self.createGroup('Line model',['Gauss', 'Voigt'])        # OK/Cancel line        hgroup = QGroupBox()        hbox = QHBoxLayout()        self.button1 = QPushButton("OK")        self.button1.clicked.connect(self.OK)        self.button2 = QPushButton("Cancel")        self.button2.clicked.connect(self.Cancel)        hbox.addWidget(self.button1)         hbox.addWidget(self.button2)        hgroup.setLayout(hbox)                # Help        label = QLabel("After defining the parameters, click and drag twice\n " +                       "on the spectrum to define the two continuum regions")        #         grid = QVBoxLayout()        line1 = QHBoxLayout()        line1.addWidget(self.function)        line1.addWidget(self.boundary)        line2 = QHBoxLayout()        line2.addWidget(self.kernel)        line2.addWidget(self.regions)        line3 = QHBoxLayout()        line3.addWidget(self.emlines)        line3.addWidget(self.ablines)        line3.addWidget(self.models)        grid.addLayout(line1)        grid.addLayout(line2)        grid.addLayout(line3)        grid.addWidget(label)        grid.addWidget(hgroup)        self.setLayout(grid)                self.setWindowTitle('Fitting parameters')        self.resize(400,300)            def createGroup(self, title, items, default=0):            """ creates a group of radio buttons  """        group = QGroupBox(title)        group.buttons = QButtonGroup()        vbox = QHBoxLayout()        buttons = []        i = 0        for item in items:            buttons.append(QRadioButton(item))            group.buttons.addButton(buttons[-1], i)            vbox.addWidget(buttons[-1])            i += 1        vbox.addStretch(1)        # Set 1st option as default        buttons[default].setChecked(True)        group.setLayout(vbox)        return group    def OK(self):        self.done(1)    def save(self):        function  = self.function.buttons.checkedButton().text()        boundary  = self.boundary.buttons.checkedButton().text()        kernel    = self.kernel.buttons.checkedButton().text()        regions   = self.regions.buttons.checkedButton().text()        emlines   = self.emlines.buttons.checkedButton().text()        ablines   = self.ablines.buttons.checkedButton().text()        model    = self.models.buttons.checkedButton().text()        return function, boundary, kernel, regions, emlines, ablines, model                def Cancel(self):        self.done(0)        class ContFitParams(QDialog):    """ Dialog window to define type of the continuum fit """    def __init__(self, options, parent=None):        super().__init__(parent)        self.options = options        self.setupUI()    def setupUI(self):                hgroup = QGroupBox()        hbox = QHBoxLayout()        self.button1 = QPushButton("OK")        self.button1.clicked.connect(self.OK)        self.button2 = QPushButton("Cancel")        self.button2.clicked.connect(self.Cancel)        hbox.addWidget(self.button1)         hbox.addWidget(self.button2)        hgroup.setLayout(hbox)                   self.group = QGroupBox('Options')        self.group.buttons = QButtonGroup()        vbox = QVBoxLayout()        buttons = []        i = 0        for option in self.options:            buttons.append(QRadioButton(option))            self.group.buttons.addButton(buttons[-1], i)            vbox.addWidget(buttons[-1])            i += 1        vbox.addStretch(1)        # Set 1st option as default        buttons[0].setChecked(True)        self.group.setLayout(vbox)        grid = QGridLayout()        grid.addWidget(self.group,0,0)        grid.addWidget(hgroup, 1, 0)        self.setLayout(grid)        self.setWindowTitle('Fitting the continuum')        self.resize(400,300)            def OK(self):        self.done(1)    def save(self):        option  = self.group.buttons.checkedButton().text()        return option                def Cancel(self):        self.done(0)        class FitCubeDialog(QDialog):    """Dialog to fit the cube."""        def __init__(self, options, moments=False, lines=False, parent=None):        super().__init__(parent)        self.options = options        self.moments = moments        self.lines = lines        self.setupUI()    def setupUI(self):                hgroup = QGroupBox()        grid = QGridLayout()        # OK/Cancel box        hbox = QHBoxLayout()        self.button1 = QPushButton("OK")        self.button1.clicked.connect(self.OK)        self.button2 = QPushButton("Cancel")        self.button2.clicked.connect(self.Cancel)        hbox.addWidget(self.button1)         hbox.addWidget(self.button2)        hgroup.setLayout(hbox)          ibox = 0        vsize = 200        # Continuum        if self.moments | self.lines:            self.continuum = self.createGroup('', self.options)            # Check buttons            checkbuttons = QButtonGroup(self)            self.checkcontinuum = QCheckBox('Continuum')            self.continuum.setEnabled(False)            checkbuttons.addButton(self.checkcontinuum)            self.checkcontinuum.stateChanged.connect(self.toggleCGroupBox)            grid.addWidget(self.checkcontinuum, ibox, 0)            ibox += 1        else:            self.continuum = self.createGroup('Continuum', self.options)        grid.addWidget(self.continuum, ibox, 0)        # Moments        if self.moments:            if 'Fit region' in set(self.options):                self.momentsbox = self.createGroup('', ['Region','All'])            else:                self.momentsbox = self.createGroup('', ['All'])               self.cbmoments = QCheckBox("Moments")            checkbuttons.addButton(self.cbmoments)            self.momentsbox.setEnabled(False)            self.cbmoments.stateChanged.connect(self.toggleMGroupBox)            ibox += 1            grid.addWidget(self.cbmoments, ibox, 0)            ibox += 1            grid.addWidget(self.momentsbox, ibox, 0)            vsize += 50        # Lines        if self.lines:            if 'Fit region' in set(self.options):                self.linesbox = self.createGroup('Lines', ['Region','All'])            else:                 self.linesbox = self.createGroup('Lines', ['All'])                           self.linesbox.setEnabled(False)            self.cblines = QCheckBox("Lines")            checkbuttons.addButton(self.cblines)            ibox += 1            grid.addWidget(self.cblines, ibox,0)            ibox += 1            grid.addWidget(self.linesbox, ibox, 0)            vsize += 50            self.cblines.stateChanged.connect(self.toggleGroupBox)        ibox += 1         grid.addWidget(hgroup, ibox, 0)        self.setLayout(grid)        self.setWindowTitle('Fitting actions')        self.resize(400,vsize)            def toggleGroupBox(self, state):        if state > 0:            self.linesbox.setEnabled(True)        else:            self.linesbox.setEnabled(False)                def toggleMGroupBox(self, state):        if state > 0:            self.momentsbox.setEnabled(True)        else:            self.momentsbox.setEnabled(False)                def toggleCGroupBox(self, state):        if state > 0:            self.continuum.setEnabled(True)        else:            self.continuum.setEnabled(False)            def createGroup(self, title, items, default=0):            """ creates a group of radio buttons  """        group = QGroupBox(title)        group.buttons = QButtonGroup()        vbox = QVBoxLayout()        buttons = []        i = 0        for item in items:            buttons.append(QRadioButton(item))            group.buttons.addButton(buttons[-1], i)            vbox.addWidget(buttons[-1])            i += 1        vbox.addStretch(1)        # Set 1st option as default        buttons[default].setChecked(True)        group.setLayout(vbox)        return group    def OK(self):        self.done(1)    def save(self):        if self.continuum.isEnabled():            continuum  = self.continuum.buttons.checkedButton().text()        else:            continuum = None        if self.moments:            if self.momentsbox.isEnabled():                moments  = self.momentsbox.buttons.checkedButton().text()            else:                moments = None        else:            moments = None        if self.lines:            if self.linesbox.isEnabled():                lines = self.linesbox.buttons.checkedButton().text()            else:                lines = None        else:            lines = None        return continuum, moments, lines                def Cancel(self):        self.done(0)                class guessParams(QDialog):    """ Dialog window to define guess parameters of continuum and lines fit """    def __init__(self, parent=None):        super().__init__()        self.setupUI()    def setupUI(self):        self.continuum = self.createGroup('Continuum', ['Constant', 'Slope', 'Fixed'], default=0)        self.emission = self.createGroup('Emission lines', ['0', '1', '2', '3'], default=1)        self.absorption = self.createGroup('Absorption lines', ['0', '1', '2', '3'], default=0)        self.function = self.createGroup('Function', ['Gaussian','Voigt'], default=0)        hgroup = QGroupBox()        hbox = QHBoxLayout()        self.button1 = QPushButton("OK")        self.button1.clicked.connect(self.OK)        self.button2 = QPushButton("Cancel")        self.button2.clicked.connect(self.Cancel)        hbox.addWidget(self.button1)        hbox.addWidget(self.button2)        hgroup.setLayout(hbox)                        grid = QVBoxLayout()        line1 = QHBoxLayout()        line1.addWidget(self.continuum)        line1.addWidget(self.function)        line2 = QHBoxLayout()        line2.addWidget(self.emission)        line2.addWidget(self.absorption)        grid.addLayout(line1)        grid.addLayout(line2)        grid.addWidget(hgroup)        self.setLayout(grid)                        self.setWindowTitle('Guess parameters')        self.resize(400, 300)    def createGroup(self, title, items, default=0):        """Creates a group of radio buttons."""        group = QGroupBox(title)        group.buttons = QButtonGroup()        hbox = QHBoxLayout()        buttons = []        i = 0        for item in items:            buttons.append(QRadioButton(item))            group.buttons.addButton(buttons[-1], i)            hbox.addWidget(buttons[-1])            i += 1        hbox.addStretch(1)        # Set 1st option as default        buttons[default].setChecked(True)        group.setLayout(hbox)        return group    def OK(self):        self.done(1)    def save(self):        continuum = self.continuum.buttons.checkedButton().text()        emission = self.emission.buttons.checkedButton().text()        absorption = self.absorption.buttons.checkedButton().text()        function = self.function.buttons.checkedButton().text()        return continuum, emission, absorption, function    def Cancel(self):        self.done(0)class SlicerDialog(QDialog):    """Dialog window to define type of slicer."""    def __init__(self, parent=None):        super().__init__(parent)        self.setupUI()    def setupUI(self):                hgroup = QGroupBox()        hbox = QHBoxLayout()        self.button1 = QPushButton("OK")        self.button1.clicked.connect(self.OK)        self.button2 = QPushButton("Cancel")        self.button2.clicked.connect(self.Cancel)        hbox.addWidget(self.button1)         hbox.addWidget(self.button2)        hgroup.setLayout(hbox)           # Group defining type of slice        self.group = QGroupBox('Show')        self.group.buttons = QButtonGroup()        vbox = QVBoxLayout()        buttons = []        buttons.append(QRadioButton('Channel'))        self.group.buttons.addButton(buttons[-1], 0)        vbox.addWidget(buttons[-1])        buttons.append(QRadioButton('Cube slice'))        self.group.buttons.addButton(buttons[-1], 0)        vbox.addWidget(buttons[-1])        buttons.append(QRadioButton('None'))        self.group.buttons.addButton(buttons[-1], 0)        vbox.addWidget(buttons[-1])        vbox.addStretch(1)        # Set 1st option as default        buttons[0].setChecked(True)        self.group.setLayout(vbox)        # Define the grid        grid = QGridLayout()        grid.addWidget(self.group,0,0)        grid.addWidget(hgroup, 1, 0)        self.setLayout(grid)        self.setWindowTitle('Slicer selection')        self.resize(300,200)            def OK(self):        self.done(1)    def save(self):        option  = self.group.buttons.checkedButton().text()        return option                def Cancel(self):        self.done(0)class cmDialog(QDialog):    dirSignal = pyqtSignal(str)        def __init__(self, cmlist, stlist, clist, currentCM, currentST, currentCC, parent=None):        super().__init__()        import os        path0, file0 = os.path.split(__file__)        self.setWindowTitle('Colors & Stretch')        layout = QVBoxLayout()        label1 = QLabel("Color maps")        self.list = QListWidget(self)        iconSize = QSize(144,10)        self.list.setIconSize(iconSize)        self.list.setSizePolicy(QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding))        self.list.setMaximumSize(QSize(160,150))        self.cmlist = cmlist        for cm in cmlist:            #item = QListWidgetItem(self.list)            #item.setText(cm)            #item.setIcon(QIcon(path0+"/icons/"+cm+".png"))            #item = QListWidgetItem(QIcon(path0+"/icons/"+cm+".png"),'',self.list)            #item.setSizeHint(iconSize)            QListWidgetItem(QIcon(os.path.join(path0,"icons",cm+".png")),'',self.list)        # For some reason this does not work when in the stylesheet        stylesheet = "QListWidget::item {"\                     +"border-style: solid;"\                     +"border-width:1px;" \                     +"border-color:transparent;"\                     +"background-color: transparent;"\                     +"color: white;"\                     +"}"\                     +"QListWidget::item:selected {"\                     +"border-style: solid;" \                     +"border-width:1px;" \                     +"border-color:black;" \                     +"background-color: transparent;"\                     +"color: white;"\                     +"}"                    self.list.setStyleSheet(stylesheet)        n = cmlist.index(currentCM)        self.list.setCurrentRow(n)        # Button to reverse color map direction        b1 = QPushButton("Reverse", self)        b1.clicked.connect(self.reverse)                label2 = QLabel("Stretches")                self.slist = QListWidget(self)        self.slist.setSizePolicy(QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding))        self.slist.setMaximumSize(QSize(160,150))        for st in stlist:            QListWidgetItem(QIcon(os.path.join(path0,"icons",st+"_.png")),st,self.slist)        n = stlist.index(currentST)        self.slist.setCurrentRow(n)        label3 = QLabel("Contour color")                self.clist = QListWidget(self)        iconSize = QSize(144,10)        self.clist.setIconSize(iconSize)        self.clist.setSizePolicy(QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding))        self.clist.setMaximumSize(QSize(160,100))        for col in clist:            QListWidgetItem(QIcon(os.path.join(path0,"icons",col+".png")),'',self.clist)        n = clist.index(currentCC[0])        self.clist.setCurrentRow(n)        self.clist.setStyleSheet(stylesheet)                # Color 2        label4 = QLabel("Contour 2 color")                self.clist2 = QListWidget(self)        self.clist2.setIconSize(iconSize)        self.clist2.setSizePolicy(QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding))        self.clist2.setMaximumSize(QSize(160,100))        for col in clist:            QListWidgetItem(QIcon(os.path.join(path0,"icons",col+".png")),'',self.clist2)        n = clist.index(currentCC[1])        self.clist2.setCurrentRow(n)        self.clist2.setStyleSheet(stylesheet)                # Button with OK to close dialog        b2 = QPushButton("OK",self)        b2.clicked.connect(self.end)        # Layout        layout.addWidget(label1)        layout.addWidget(self.list)        layout.addWidget(b1)        layout.addWidget(label2)        layout.addWidget(self.slist)        layout.addWidget(label3)        layout.addWidget(self.clist)        layout.addWidget(label4)        layout.addWidget(self.clist2)        layout.addWidget(b2)        self.setLayout(layout)    def end(self):        self.close()            def reverse(self):        self.dirSignal.emit('color map reversed')