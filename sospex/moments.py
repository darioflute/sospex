import numpy as np
from PyQt5.QtCore import pyqtSignal,QObject
from PyQt5.QtWidgets import (QDialog, QPushButton, QGroupBox, QHBoxLayout, QVBoxLayout,
                             QGridLayout, QRadioButton, QButtonGroup, QLabel, QCheckBox)
from matplotlib.lines import Line2D
from matplotlib.artist import Artist
import matplotlib.transforms as transforms
import multiprocessing as mp
from lmfit import Parameters, minimize

#class Guess(object):
#    """ class to define a spectrum """
#    def __init__(self, pts, xpts, yc, offset):
#
#        self.min = pts[0]
#        self.max = pts[3]
#        delta = xpts[2]-xpts[1]
#        self.mean  = xpts[1]+delta*0.5
#        self.sigma = delta/5.
#        self.amp  = yc
#        self.offset = offset
#
#        # initialize boundaries
#        if self.amp > 0:
#            self.amplims = [0, self.amp*1.2]
#        else:
#            self.amplims = [self.amp*1.2, 0]
#
#        self.meanlims = [self.mean-self.sigma, self.mean+self.sigma]
#        self.sigmalims = [0,self.sigma*1.5]
#        if self.offset > 0:1
#            self.offlims = [0,2*self.offset]
#        else:
#            self.offlims = [None,None]

class SegmentsSelector(QObject):

    def __init__(self, ax, fig, callback, color='#7ec0ee', zD = True):
        super().__init__()

        self.x = []
        self.y = []
        self.line1 = None
        self.line2 = None
        self.color = color
        self.fig = fig
        self.ax = ax
        self.callback = callback
        self.zeroDeg = zD

        self.__ID1 = self.fig.canvas.mpl_connect('motion_notify_event', self.__motion_notify_callback)
        self.__ID2 = self.fig.canvas.mpl_connect('button_press_event', self.__button_press_callback)
        self.__ID3 = self.fig.canvas.mpl_connect('button_release_event', self.__button_release_callback)

    def __motion_notify_callback(self, event):
        if event.inaxes:
            #ax = event.inaxes
            x, y = event.xdata, event.ydata
            if (event.button == None or event.button == 1):
                if self.line1 != None: # Move line around
                    if self.line2 == None:
                        if self.zeroDeg: self.y[0]=y
                        self.line1.set_data([self.x[0], x],
                                            [self.y[0], y])
                    else:
                        if self.zeroDeg:
                            self.y=[y,y,y,y]
                            self.line1.set_data([self.x[0], self.x[1]],
                                                [self.y[0], self.y[1]])
                        self.line2.set_data([self.x[2], x],
                                            [self.y[2], y])
                    self.fig.canvas.draw_idle()


    def __button_release_callback(self, event):
        if event.inaxes:
            x, y = event.xdata, event.ydata
            #ax = event.inaxes
            if event.button == 1:
                if self.line2 == None:  # Segment 1 completed
                    self.x.append(x)
                    self.y.append(y)
                    if self.zeroDeg:
                        self.y[-2]=self.y[-1]
                    self.line1.set_data([self.x[0], self.x[1]],
                                        [self.y[0], self.y[1]])
                    self.fig.canvas.draw_idle()
                    self.fig.canvas.mpl_disconnect(self.__ID1)
                    self.line2 = 'start'
                else:
                    self.x.append(x)
                    self.y.append(y)
                    if self.zeroDeg:
                        self.y[-1]=self.y[-2]
                        m = 0.
                    else:
                        # Adjust to the same slope between first and last point
                        m = (self.y[3]-self.y[0])/(self.x[3]-self.x[0])
                    for i in range(4):
                        self.y[i] = self.y[0]+m*(self.x[i]-self.x[0])
                    self.line1.set_data([self.x[0], self.x[1]],
                                        [self.y[0], self.y[1]])
                    self.line2.set_data([self.x[2], self.x[3]],
                                        [self.y[2], self.y[3]])
                    self.fig.canvas.draw_idle()
                    self.xy = [(i,j) for (i,j) in zip(self.x,self.y)]
                    # Disconnect
                    self.fig.canvas.mpl_disconnect(self.__ID1) 
                    self.fig.canvas.mpl_disconnect(self.__ID2) 
                    self.fig.canvas.mpl_disconnect(self.__ID3) 
                    # Callback function, pass the vertices
                    self.callback(self.xy)
                    # Remove lines
                    self.remove()

    def __button_press_callback(self, event):
        if event.inaxes:
            x, y = event.xdata, event.ydata
            ax = event.inaxes
            if event.button == 1:
                if self.line2 == None:
                    if self.line1 == None:  # If you press the left button, single click
                        self.line1 = Line2D([x, x],
                                            [y, y],
                                            marker='o',
                                            color=self.color)
                        self.start_point = [x,y]
                        self.previous_point =  self.start_point
                        self.x=[x]
                        self.y=[y]
                        self.__ID1 = self.fig.canvas.mpl_connect('motion_notify_event', self.__motion_notify_callback)
                        ax.add_line(self.line1)
                        # add a segment
                        self.fig.canvas.draw_idle()
                else:
                    if self.zeroDeg:
                        self.y = [y,y]
                        self.line1.set_data([self.x[0], self.x[1]],
                                            [self.y[0], self.y[1]])
                    self.line2 = Line2D([x, x],
                                        [y, y],
                                        marker='o',
                                        color=self.color)
                    self.start_point = [x,y]
                    self.previous_point =  self.start_point
                    self.x.append(x)
                    self.y.append(y)
                    self.__ID1 = self.fig.canvas.mpl_connect('motion_notify_event', self.__motion_notify_callback)
                        
                    ax.add_line(self.line2)
                    self.fig.canvas.draw()

    def remove(self):
        """ Remove lines from plot """
        try:
            self.line1.remove()
            self.line2.remove()
        except:
            print('no lines to remove')


class SegmentsInteractor(QObject):
    """
    A continuum interactor.
    """
    
    showverts = True
    epsilon = 5  # max pixel distance to count as a vertex hit
    mySignal = pyqtSignal(str)
    modSignal = pyqtSignal(str)

    def __init__(self, ax, verts, zeroDeg=True, color='#7ec0ee'):
        super().__init__()

        self.ax = ax
        self.type = 'Continuum'
        # color = 'skyblue'
        # color = '#7ec0ee'
        self.zeroDeg = zeroDeg
        x, y = zip(*verts)
        self.xy = [(i,j) for (i,j) in zip(x,y)]
        self.computeSlope()        
        self.line1 = Line2D(x[:2],y[:2],color=color,linewidth=2, animated = True)
        self.line2 = Line2D(x[2:],y[2:],color=color,linewidth=2, animated = True)
        trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
        self.xl1a = Line2D([x[0], x[0]], [0, 1], color=color, linewidth=1,
                           animated=True, transform=trans)
        self.xl1b = Line2D([x[1], x[1]], [0, 1], color=color, linewidth=1,
                           animated=True, transform=trans)
        self.xl2a = Line2D([x[2], x[2]], [0, 1], color=color, linewidth=1,
                           animated=True, transform=trans)
        self.xl2b = Line2D([x[3], x[3]], [0, 1], color=color, linewidth=1,
                           animated=True, transform=trans)
        self.canvas = ax.figure.canvas
        self.line = Line2D(x, y, marker='o', linestyle=None, linewidth=0., markerfacecolor=color, animated=True)                
        self.artists = [self.line1, self.line2, self.xl1a, self.xl1b, self.xl2a, self.xl2b, self.line]
        for artist in self.artists:
            self.ax.add_line(artist)
        self.cid = self.line1.add_callback(self.si_changed)
        self._ind = None  # the active vert
        self.connect()

    def computeSlope(self):
        xg,yg = zip(*self.xy)
        xg = np.array(xg); yg = np.array(yg)
        if self.zeroDeg:
            self.slope = 0
        else:
            self.slope = (yg[3]-yg[0])/(xg[3]-xg[0])
        self.intcpt = yg[0]-self.slope*xg[0]

    def connect(self):
        self.cid_draw = self.canvas.mpl_connect('draw_event',
                                                self.draw_callback)
        self.cid_press = self.canvas.mpl_connect('button_press_event',
                                                 self.button_press_callback)
        self.cid_release = self.canvas.mpl_connect('button_release_event',
                                                   self.button_release_callback)
        self.cid_motion = self.canvas.mpl_connect('motion_notify_event',
                                                  self.motion_notify_callback)
        self.cid_key = self.canvas.mpl_connect('key_press_event',
                                               self.key_press_callback)
        self.canvas.draw_idle()

    def disconnect(self):
        self.canvas.mpl_disconnect(self.cid_draw)
        self.canvas.mpl_disconnect(self.cid_press)
        self.canvas.mpl_disconnect(self.cid_release)
        self.canvas.mpl_disconnect(self.cid_motion)
        self.canvas.mpl_disconnect(self.cid_key)
        try:
            self.line1.remove()
            self.xl1a.remove()
            self.xl1b.remove()
        except:
            print('no line 1')
        try:
            self.line2.remove()
            self.xl2a.remove()
            self.xl2b.remove()
        except:
            print('no line 2')
        try:
            self.line.remove()
        except:
            print('no markers')
        self.canvas.draw_idle()
        self.aperture = None
        
    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        for artist in self.artists:
            self.ax.draw_artist(artist)

    def si_changed(self, line1):
        'this method is called whenever the line1 object is called'
        # only copy the artist props to the line (except visibility)
        vis = self.line.get_visible()
        Artist.update_from(self.line, line1)
        self.line.set_visible(vis)  

    def get_ind_under_point(self, event):
        'get the index of the point if within epsilon tolerance'
        # Distance is computed in pixels on the screen
        xy = self.ax.transData.transform(self.xy)
        x, y = zip(*xy)
        x = np.array(x); y = np.array(y)
        d = np.hypot(x - event.x, y - event.y)
        indseq, = np.nonzero(d == d.min())
        ind = indseq[0]
        # print('distance is ',d[ind])
        if d[ind] >= self.epsilon:
            ind = None
        return ind

    def key_press_callback(self, event):
        'whenever a key is pressed'
        if not event.inaxes:
            return
        if event.key == 't':
            self.showverts = not self.showverts
            self.line.set_visible(self.showverts)
            if not self.showverts:
                self._ind = None
        elif event.key == 'd':
            ind = self.get_ind_under_point(event)
            if ind is not None:
                self.mySignal.emit('segments deleted')
        self.canvas.draw_idle()

    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        if not self.showverts:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        self._ind = self.get_ind_under_point(event)

    def button_release_callback(self, event):
        'whenever a mouse button is released'
        if not self.showverts:
            return
        if event.button != 1:
            return
        self._ind = None

    def motion_notify_callback(self, event):
        'on mouse movement'
        if not self.showverts:
            return
        if self._ind is None:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        # Rebuild line collection
        x,y = zip(*self.xy)
        x = np.asarray(x)
        y = np.asarray(y)
        # Check order of x (increasing if wavelength, decreasing if frequency)
        if x[3] > x[0]:
            case = 'wave'
        else:
            case = 'freq'
        # Update point
        x_, y_ = event.xdata, event.ydata
        y[self._ind] = y_
        # How avoiding overstepping ...
        if case == 'wave':
            if self._ind == 0:
                if x_ < x[1]:
                    x[0] = x_
            elif self._ind == 1:
                if (x_ > x[0]) & (x_ < x[2]):
                    x[1] = x_
            elif self._ind == 2:
                if (x_ > x[1]) & (x_ < x[3]):
                    x[2] = x_
            elif self._ind == 3:
                if (x_ > x[2]):
                    x[3] = x_
        else:
            if self._ind == 0:
                if x_ > x[1]:
                    x[0] = x_
            elif self._ind == 1:
                if (x_ < x[0]) & (x_ > x[2]):
                    x[1] = x_
            elif self._ind == 2:
                if (x_ < x[1]) & (x_ > x[3]):
                    x[2] = x_
            elif self._ind == 3:
                if (x_ < x[2]):
                    x[3] = x_
        if self.zeroDeg:
            m = 0
        else:
            if self._ind < 2:
                m = (y[3]-y[self._ind])/(x[3]-x[self._ind])
            else:
                m = (y[self._ind]-y[0])/(x[self._ind]-x[0])    
        for i in range(4):
            y[i] = y[self._ind]+m*(x[i]-x[self._ind])
            self.xy[i] = (x[i],y[i])
        # Update segments and markers
        self.updateLinesMarkers()
        self.redraw()
        # Notify callback
        self.modSignal.emit('continuum guess modified')

    def updateLinesMarkers(self):
        x, y = zip(*self.xy)
        self.xl1a.set_xdata([x[0], x[0]])
        self.xl1b.set_xdata([x[1], x[1]])
        self.xl2a.set_xdata([x[2], x[2]])
        self.xl2b.set_xdata([x[3], x[3]])        
        self.line1.set_data(zip(*self.xy[:2]))
        self.line2.set_data(zip(*self.xy[2:]))
        self.line.set_data(zip(*self.xy))

    def redraw(self):
        self.canvas.restore_region(self.background)
        for artist in self.artists:
            self.ax.draw_artist(artist)
        self.canvas.update()
        self.canvas.flush_events()

    def switchUnits(self):
        """ Redraw segments in new units """
        # Rebuild line collection
        x, y = zip(*self.xy)
        x = np.asarray(x)
        y = np.asarray(y)
        c = 299792458.0  # speed of light in m/s
        x = c/x * 1.e-6  # um to THz or viceversa
        for i in range(4):
            self.xy[i] = (x[i],y[i])
        # Update segments and markers
        self.updateLinesMarkers()
        # self.redraw()
        
# Dialogs

class ContParams(QDialog):
    """ Dialog window to define parameters of the continuum fit """
    
    def __init__(self, k, parent=None):
        super().__init__(parent)
        if k == 1:
            self.k = 0
        elif k == 5:
            self.k = 1
        elif k == 9:
            self.k = 2
        else:
            self.k = 0            
        self.setupUI()

    def setupUI(self):        
        self.function = self.createGroup('Continuum function', ['Constant', 'Slope'])
        self.boundary = self.createGroup('Continuum boundary', ['None', 'Non negative'])
        self.kernel   = self.createGroup('Kernel pixels', ['1', '5', '9'], default=self.k)
        self.regions = self.createGroup('No of regions', ['4', '16' ,'64', '128', '256'])
        #self.emlines = self.createGroup('No of emission lines', ['0', '1', '2', '3'])
        #self.ablines = self.createGroup('No of absorption lines', ['0', '1', '2'])
        self.emlines = self.createGroup('No of emission lines', ['0', '1', '2'])
        self.ablines = self.createGroup('No of absorption lines', ['0'])
        # OK/Cancel line
        hgroup = QGroupBox()
        hbox = QHBoxLayout()
        self.button1 = QPushButton("OK")
        self.button1.clicked.connect(self.OK)
        self.button2 = QPushButton("Cancel")
        self.button2.clicked.connect(self.Cancel)
        hbox.addWidget(self.button1) 
        hbox.addWidget(self.button2)
        hgroup.setLayout(hbox)        
        # Help
        label = QLabel("After defining the parameters, click and drag twice\n" +
                       "on the spectrum to define the two continuum regions")
        # Grid
        grid = QGridLayout()
        grid.addWidget(self.function, 0, 0)
        grid.addWidget(self.boundary, 0, 1)
        grid.addWidget(self.kernel, 1, 0)
        grid.addWidget(self.regions, 1, 1)
        grid.addWidget(self.emlines, 2, 0)
        grid.addWidget(self.ablines, 2, 1)
        grid.addWidget(label, 3, 0 , 1, 2)
        grid.addWidget(hgroup, 4, 0, 1, 2)
        self.setLayout(grid)
        self.setWindowTitle('Fitting parameters')
        self.resize(400,300)

    def createGroup(self, title, items, default=0):    
        """ creates a group of radio buttons  """
        group = QGroupBox(title)
        group.buttons = QButtonGroup()
        vbox = QHBoxLayout()
        buttons = []
        i = 0
        for item in items:
            buttons.append(QRadioButton(item))
            group.buttons.addButton(buttons[-1], i)
            vbox.addWidget(buttons[-1])
            i += 1
        vbox.addStretch(1)
        # Set 1st option as default
        buttons[default].setChecked(True)
        group.setLayout(vbox)
        return group

    def OK(self):
        self.done(1)

    def save(self):
        function  = self.function.buttons.checkedButton().text()
        boundary  = self.boundary.buttons.checkedButton().text()
        kernel    = self.kernel.buttons.checkedButton().text()
        regions   = self.regions.buttons.checkedButton().text()
        emlines   = self.emlines.buttons.checkedButton().text()
        ablines   = self.ablines.buttons.checkedButton().text()
        return function, boundary, kernel, regions, emlines, ablines
            
    def Cancel(self):
        self.done(0)
        
class ContFitParams(QDialog):
    """ Dialog window to define type of the continuum fit """
    def __init__(self, options, parent=None):
        super().__init__(parent)
        self.options = options
        self.setupUI()

    def setupUI(self):        
        hgroup = QGroupBox()
        hbox = QHBoxLayout()
        self.button1 = QPushButton("OK")
        self.button1.clicked.connect(self.OK)
        self.button2 = QPushButton("Cancel")
        self.button2.clicked.connect(self.Cancel)
        hbox.addWidget(self.button1) 
        hbox.addWidget(self.button2)
        hgroup.setLayout(hbox)   
        
        self.group = QGroupBox('Options')
        self.group.buttons = QButtonGroup()
        vbox = QVBoxLayout()
        buttons = []
        i = 0
        for option in self.options:
            buttons.append(QRadioButton(option))
            self.group.buttons.addButton(buttons[-1], i)
            vbox.addWidget(buttons[-1])
            i += 1
        vbox.addStretch(1)
        # Set 1st option as default
        buttons[0].setChecked(True)
        self.group.setLayout(vbox)

        grid = QGridLayout()
        grid.addWidget(self.group,0,0)
        grid.addWidget(hgroup, 1, 0)
        self.setLayout(grid)
        self.setWindowTitle('Fitting the continuum')
        self.resize(400,300)
        
    def OK(self):
        self.done(1)

    def save(self):
        option  = self.group.buttons.checkedButton().text()
        return option
            
    def Cancel(self):
        self.done(0)
        
class FitCubeDialog(QDialog):
    """Dialog to fit the cube."""
    
    def __init__(self, options, moments=False, lines=False, parent=None):
        super().__init__(parent)
        self.options = options
        self.moments = moments
        self.lines = lines
        self.setupUI()

    def setupUI(self):        
        hgroup = QGroupBox()
        grid = QGridLayout()
        # OK/Cancel box
        hbox = QHBoxLayout()
        self.button1 = QPushButton("OK")
        self.button1.clicked.connect(self.OK)
        self.button2 = QPushButton("Cancel")
        self.button2.clicked.connect(self.Cancel)
        hbox.addWidget(self.button1) 
        hbox.addWidget(self.button2)
        hgroup.setLayout(hbox)  
        ibox = 0
        vsize = 200
        # Continuum
        if self.moments | self.lines:
            self.continuum = self.createGroup('', self.options)
            # Check buttons
            checkbuttons = QButtonGroup(self)
            self.checkcontinuum = QCheckBox('Continuum')
            self.continuum.setEnabled(False)
            checkbuttons.addButton(self.checkcontinuum)
            self.checkcontinuum.stateChanged.connect(self.toggleCGroupBox)
            grid.addWidget(self.checkcontinuum, ibox, 0)
            ibox += 1
        else:
            self.continuum = self.createGroup('Continuum', self.options)
        grid.addWidget(self.continuum, ibox, 0)
        # Moments
        if self.moments:
            if 'Fit region' in set(self.options):
                self.momentsbox = self.createGroup('', ['Region','All'])
            else:
                self.momentsbox = self.createGroup('', ['All'])   
            self.cbmoments = QCheckBox("Moments")
            checkbuttons.addButton(self.cbmoments)
            self.momentsbox.setEnabled(False)
            self.cbmoments.stateChanged.connect(self.toggleMGroupBox)
            ibox += 1
            grid.addWidget(self.cbmoments, ibox, 0)
            ibox += 1
            grid.addWidget(self.momentsbox, ibox, 0)
            vsize += 50
        # Lines
        if self.lines:
            if 'Fit region' in set(self.options):
                self.linesbox = self.createGroup('Lines', ['Region','All'])
            else:
                 self.linesbox = self.createGroup('Lines', ['All'])               
            self.linesbox.setEnabled(False)
            self.cblines = QCheckBox("Lines")
            checkbuttons.addButton(self.cblines)
            ibox += 1
            grid.addWidget(self.cblines, ibox,0)
            ibox += 1
            grid.addWidget(self.linesbox, ibox, 0)
            vsize += 50
            self.cblines.stateChanged.connect(self.toggleGroupBox)
        ibox += 1 
        grid.addWidget(hgroup, ibox, 0)
        self.setLayout(grid)
        self.setWindowTitle('Fitting actions')
        self.resize(400,vsize)
        

    def toggleGroupBox(self, state):
        if state > 0:
            self.linesbox.setEnabled(True)
        else:
            self.linesbox.setEnabled(False)
            
    def toggleMGroupBox(self, state):
        if state > 0:
            self.momentsbox.setEnabled(True)
        else:
            self.momentsbox.setEnabled(False)
            
    def toggleCGroupBox(self, state):
        if state > 0:
            self.continuum.setEnabled(True)
        else:
            self.continuum.setEnabled(False)
        

    def createGroup(self, title, items, default=0):    
        """ creates a group of radio buttons  """
        group = QGroupBox(title)
        group.buttons = QButtonGroup()
        vbox = QVBoxLayout()
        buttons = []
        i = 0
        for item in items:
            buttons.append(QRadioButton(item))
            group.buttons.addButton(buttons[-1], i)
            vbox.addWidget(buttons[-1])
            i += 1
        vbox.addStretch(1)
        # Set 1st option as default
        buttons[default].setChecked(True)
        group.setLayout(vbox)
        return group

    def OK(self):
        self.done(1)

    def save(self):
        if self.continuum.isEnabled():
            continuum  = self.continuum.buttons.checkedButton().text()
        else:
            continuum = None
        if self.moments:
            if self.momentsbox.isEnabled():
                moments  = self.momentsbox.buttons.checkedButton().text()
            else:
                moments = None
        else:
            moments = None
        if self.lines:
            if self.linesbox.isEnabled():
                lines = self.linesbox.buttons.checkedButton().text()
            else:
                lines = None
        else:
            lines = None
        return continuum, moments, lines
            
    def Cancel(self):
        self.done(0)

class SlicerDialog(QDialog):
    """Dialog window to define type of slicer."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUI()

    def setupUI(self):        
        hgroup = QGroupBox()
        hbox = QHBoxLayout()
        self.button1 = QPushButton("OK")
        self.button1.clicked.connect(self.OK)
        self.button2 = QPushButton("Cancel")
        self.button2.clicked.connect(self.Cancel)
        hbox.addWidget(self.button1) 
        hbox.addWidget(self.button2)
        hgroup.setLayout(hbox)   
        # Group defining type of slice
        self.group = QGroupBox('Show')
        self.group.buttons = QButtonGroup()
        vbox = QVBoxLayout()
        buttons = []
        buttons.append(QRadioButton('Channel'))
        self.group.buttons.addButton(buttons[-1], 0)
        vbox.addWidget(buttons[-1])
        buttons.append(QRadioButton('Cube slice'))
        self.group.buttons.addButton(buttons[-1], 0)
        vbox.addWidget(buttons[-1])
        buttons.append(QRadioButton('None'))
        self.group.buttons.addButton(buttons[-1], 0)
        vbox.addWidget(buttons[-1])
        vbox.addStretch(1)
        # Set 1st option as default
        buttons[0].setChecked(True)
        self.group.setLayout(vbox)
        # Define the grid
        grid = QGridLayout()
        grid.addWidget(self.group,0,0)
        grid.addWidget(hgroup, 1, 0)
        self.setLayout(grid)
        self.setWindowTitle('Slicer selection')
        self.resize(300,200)
        
    def OK(self):
        self.done(1)

    def save(self):
        option  = self.group.buttons.checkedButton().text()
        return option
            
    def Cancel(self):
        self.done(0)


# Functions for multiprocessing continuum fit and moment computation

def residuals(p,x,data=None,eps=None):
    #unpack parameters
    v = p.valuesdict()
    q = v['q']
    # define model
    try:
        m = v['m']
        model = m*x+q
    except:
        model = q
    if data is None:
        return model
    else:
        if eps is None:
            return model-data 
        else:
            return (model-data)/eps

def fitContinuum(p,slope,intcp,posCont,m,w,ff):
    # Take the mean for each wavelength
    f = np.nanmean(ff,axis=1)
    mf = np.isnan(f)
    m[mf] = 0
    if np.sum(m) > 5:
        # Define parameters
        fit_params = Parameters()
        if slope != 0:
            fit_params.add('m',value=slope)
        if posCont:
            fit_params.add('q',value=intcp, min=0.)
        else:
            fit_params.add('q',value=intcp)
        out = minimize(residuals,fit_params,args=(w[m],),kws={'data':f[m]},method='Nelder')
        pars = out.params
    else:
        pars = None
        pass
    return p, pars

def fiteContinuum(p,slope,intcp,posCont,m,w,ff,ee):
    # Skip if exposure is zero
    se = np.nansum(ee)
    if se == 0:
        pars = None
        return p, pars
    # Take the mean for each wavelength
    f = np.nanmean(ff,axis=1)
    mf = np.isnan(f)
    m[mf] = 0
    if np.sum(m) > 5:
        # Define parameters
        fit_params = Parameters()
        if slope != 0:
            fit_params.add('m',value=slope)
        if posCont:
            fit_params.add('q',value=intcp, min=0.)
        else:
            fit_params.add('q',value=intcp)
        e = 1./np.sqrt(np.nanmean(ee,axis=1))  # The error scale with the inverse of the sqrt of the exposure
        out = minimize(residuals,fit_params,args=(w[m],),kws={'data':f[m],'eps':e[m]},method='Nelder')
        pars = out.params
    else:
        pars = None
        pass
    return p, pars


def computeMoments(p,m,w,dw,f):
    """ compute moments on a spatial pixel """
    mf = np.isnan(f)
    m[mf] = 0
    # compute the error on f
    if np.sum(m) > 5:
        f = f[m]
        w = w[m]
        dw = dw[m]
        
        #sf = medfilt(f,5)
        #df = f-sf
        # The dispersion is computed using the difference
        # between two consecutive values. Assuming they have the
        # same dispersion, the dispersion of the difference is
        # sqrt(2) higher.

        pos = f[1:] > 0 # Consider only positive values
        df = f[1:]-f[:-1]
        df = df[pos]
        med = np.nanmedian(df)
        mad = np.nanmedian(np.abs(med))
        sigma = 3 *  mad/np.sqrt(2.) # 3 sigma value
        # Consider only values greater than continuum
        #if np.isfinite(sigma):
        #    ms = f > (-sigma)
        #else:
        #    ms = np.isfinite(f)
        ms = f > 0
        # Compute also negative intensity
        if np.sum(ms) > 5:
            c = 299792458. # m/s
            w_ = w[ms]
            dw_ = dw[ms]
            Snu = f[ms]
            pos = Snu > 0
            #pos = Snu > (-sigma)
            Slambda = c*Snu[pos]/(w_[pos]*w_[pos])*1.e6   # [Jy * Hz / um]
            w_  = w_[pos]
            dw_ = dw_[pos]
            M0 = np.nansum(Slambda*dw_) # [Jy Hz] 
            M1 = np.nansum(w_*Slambda*dw_)/M0 # [um]
            M2 = np.nansum(np.power(w_-M1,2)*Slambda*dw_)/M0 # [um*um]
            SD = np.sqrt(M2)
            M3 = np.nansum(np.power(w_-M1,3)*Slambda*dw_)/M0/np.power(SD,3)
            M4 = np.nansum(np.power(w_-M1,4)*Slambda*dw_)/M0/np.power(SD,4)-3. # Relative to Gaussian which is 3
            M0 *= 1.e-26 # [W/m2]  (W/m2 = Jy*Hz*1.e-26)
        else:
            M0 = np.nan
            M1 = M2 = M3 = M4 = np.nan
    else:
        M0 = np.nan
        M1 = M2 = M3 = M4 = np.nan
        sigma = np.nan
            
    return p, M0, M1, M2, M3, M4, sigma

def multiComputeMoments(m,w,f,c,moments,points):

    # To avoid forking error in MAC OS-X
    try:
        mp.set_start_method('spawn')
    except RuntimeError:
        pass

    # Define noise
    n3,n2,n1 = np.shape(moments)
    noise = np.zeros((n2,n1))

    
    # Compute dw
    dw = [] 
    dw.append([w[1]-w[0]])
    dw.append(list((w[2:]-w[:-2])*0.5))
    dw.append([w[-1]-w[-2]])
    dw = np.concatenate(dw)

    with mp.Pool(processes=mp.cpu_count()) as pool:
        res = [pool.apply_async(computeMoments, (p,m[:,p[1],p[0]],w,dw,f[:,p[1],p[0]]-c[:,p[1],p[0]])) for p in points]
        results = [r.get() for r in res]

        
    for p, M0, M1, M2, M3, M4, sigma in results:
        i,j = p
        moments[0][j,i] = M0
        moments[1][j,i] = M1
        moments[2][j,i] = M2
        moments[3][j,i] = M3
        moments[4][j,i] = M4
        noise[j,i] = sigma
            
    return moments, noise
    

def multiFitContinuum(m, w, f, c, c0, w0, points, slope, intcp, posCont, kernel, exp=None):
    # To avoid forking error in MAC OS-X
    try:
        mp.set_start_method('spawn')
    except RuntimeError:
        pass
    if kernel == 1:
        ik = np.array([0])
        jk = np.array([0])
    elif kernel == 5:
        ik = np.array([-1,0,0,0,1])
        jk = np.array([0,-1,0,1,0])
    elif kernel == 9:
        ik = np.array([-1,-1,-1,0,0,0,1,1,1])
        jk = np.array([-1,0,1,-1,0,1,-1,0,1])
    else:
        print ('unsupported kernel, use one pixel only')
        ik = np.array([0])
        jk = np.array([0])
    if exp is None:
        with mp.Pool(processes=mp.cpu_count()) as pool:
            res = [pool.apply_async(fitContinuum, (p,slope,intcp,posCont,m[:,p[1],p[0]],w,
                                                   f[:,p[1]+ik,p[0]+jk])) for p in points]
            results = [r.get() for r in res]
    else:
        with mp.Pool(processes=mp.cpu_count()) as pool:
            res = [pool.apply_async(fiteContinuum, 
                                    (p, slope, intcp, posCont, m[:,p[1],p[0]], w,
                                     f[:,p[1]+ik,p[0]+jk],exp[:,p[1]+ik,p[0]+jk])) for p in points]
            results = [r.get() for r in res] 
    print('c is writeable ', c.flags)
    c_ = c.copy()
    c0_ = c0.copy()  # continuum at ref wav
    cs_ = c0.copy()  # slope of cont
    for p, pars in results:
        if pars is not None:
            i, j = p
            c_[:, j, i] = residuals(pars, w)
            c0_[j, i] = residuals(pars, w0)
            try:
                cs_[j, i] = pars['m'].value
            except:
                cs_[j, i] = 0.
    return c_, c0_, cs_

# Fit of lines

def fitLines(p, m,w,f,lines):
    """Fit the lines defined in the guess."""
    from lmfit.models import PseudoVoigtModel
    # f is flux without continuum
    # lines is the list of parameters for each line
    mf = np.isnan(f)
    m[mf] = 0
    if np.sum(m) > 5:
        y = f[m]
        x = w[m]
        # Transform into S(lambda)
        c = 299792458. # m/s
        #  W/m2 = Jy*Hz*1.e-26 and c -> c*1.e6 um ...
        y  = c * y / (x * x ) *1.e-20 # Jy * Hz / um --> W/m2 after integration in um
        # Normalization
        norm = np.abs(np.median(y))
        #print('Normalization factor ', norm)
        y /= norm
        # Define lines
        for i, line in enumerate(lines):
            li = 'l' + str(i) + '_'
            vmodel = PseudoVoigtModel(prefix=li)
            if i == 0:
                params = vmodel.make_params()
                model = vmodel
            else:
                params += vmodel.make_params()
                model += vmodel
            x0 = line[0]
            sigma = line[1] / 2.355
            A = line[2] * c / (x0*x0) * 1.e-20 / norm # same units as flux
            params[li+'center'].set(x0, min=(x0 - sigma/2.), max=(x0 + sigma/2.))
            if A > 0:
                params[li+'amplitude'].set(A, min=0.1 * A, max=A * 2)
            else:
                params[li+'amplitude'].set(A, min=2 * A, max=A * 0.1)
            params[li+'sigma'].set(sigma, min=sigma / 2., max=sigma * 2)
            params[li + 'fraction'].set(0.0, vary=False, max = 0.3)  # No Cauchy part (i.e. Gauss)
        # Minimize
        out = model.fit(y, params, x=x, method='nelder')
        #               kws={'data': y, 'eps': e}, method='leastsq')
        # Return lines fitted parameters
        pars = out.params#.valuesdict()
        nlines = len(lines)
        #print("Number of lines fitted: ", nlines)
        linepars = []
        for i in range(nlines):
            li = 'l' + str(i) + '_'
            center = pars[li + 'center'].value  # Observed
            #centerErr = pars[li + 'center'].stderr  # Observed
            sigma = pars[li + 'sigma'] .value   # Observed
            #sigmaErr = pars[li + 'sigma'] .stderr   # Observed
            A =  pars[li+'amplitude'].value
            #Aerr = pars[li+'amplitude'].stderr
            amplitude = A * norm
            #amplitudeErr = amplitude * (Aerr/A + sigmaErr/sigma)
            alpha = pars[li+'fraction'].value
            linepars.append([center, sigma, amplitude, alpha])
        return p, linepars
    else:
        nlines = len(lines)
        return p, np.full((nlines, 4), np.nan)

def multiFitLines(m, w, f, c, lineguesses, linefits, points):

    # To avoid forking error in MAC OS-X
    try:
        mp.set_start_method('spawn')
    except RuntimeError:
        pass

    with mp.Pool(processes=mp.cpu_count()) as pool:
        res = [pool.apply_async(fitLines, 
                                (p, 
                                 m[:, p[1],p[0]], 
                                 w, 
                                 f[:,p[1],p[0]]-c[:,p[1],p[0]], 
                                 lineguesses)
                                ) for p in points]
        results = [r.get() for r in res]

    n = len(lineguesses)
    for p, linepars in results:
        i,j = p
        for k in range(n):
            linefits[k][0][j,i] = linepars[k][0]
            linefits[k][1][j,i] = linepars[k][1]
            linefits[k][2][j,i] = linepars[k][2]
            linefits[k][3][j,i] = linepars[k][3]
            
    return 1

def multiFitLines2(m, w, f, c, lineguesses, linefits, points):

    for p in points:
        res = fitLines(p, 
                       m[:, p[1],p[0]], 
                       w, 
                       f[:,p[1],p[0]]-c[:,p[1],p[0]], 
                       lineguesses)

        n = len(lineguesses)
        pp, linepars = res
        i,j = pp
        print(np.shape(linepars), np.shape(linefits))
        for k in range(n):
            linefits[k][0][j,i] = linepars[k][0]
            linefits[k][1][j,i] = linepars[k][1]
            linefits[k][2][j,i] = linepars[k][2]
            linefits[k][3][j,i] = linepars[k][3]
                
    return 1

def residualsPsf(p, x, y, data=None, err=None):
    '''Residual of a PSF with unknown center'''
    v = p.valuesdict()
    s = v['s']
    A = v['A']
    x0 = v['x0']
    y0 = v['y0']
    dis = np.hypot(x-x0,y-y0)/s
    d2 = np.square(dis.flatten())
    model = A *  np.exp(-0.5 * d2)
    
    if data is None:
        return model
    else:
        if err is None:
            return (model - data.flatten())
        else:
            return (model - data.flatten())/err.flatten()
