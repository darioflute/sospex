import numpy as np
from PyQt5.QtCore import pyqtSignal,QObject
from matplotlib.lines import Line2D

class Guess(object):
    """ class to define a spectrum """
    def __init__(self, pts, xpts, yc, offset):

        self.min = pts[0]
        self.max = pts[3]
        delta = xpts[2]-xpts[1]
        self.mean  = xpts[1]+delta*0.5
        self.sigma = delta/5.
        self.amp  = yc
        self.offset = offset

        # initialize boundaries
        if self.amp > 0:
            self.amplims = [0, self.amp*1.2]
        else:
            self.amplims = [self.amp*1.2, 0]

        self.meanlims = [self.mean-self.sigma, self.mean+self.sigma]
        self.sigmalims = [0,self.sigma*1.5]
        if self.offset > 0:
            self.offlims = [0,2*self.offset]
        else:
            self.offlims = [None,None]


class SegmentsSelector:

    def __init__(self, ax, fig, callback, color='b'):

        self.x = []
        self.y = []
        self.line1 = None
        self.line2 = None
        self.color = color
        self.fig = fig
        self.ax = ax
        self.callback = callback

        self.__ID1 = self.fig.canvas.mpl_connect('motion_notify_event', self.__motion_notify_callback)
        self.__ID2 = self.fig.canvas.mpl_connect('button_press_event', self.__button_press_callback)
        self.__ID3 = self.fig.canvas.mpl_connect('button_release_event', self.__button_release_callback)

    def __motion_notify_callback(self, event):
        if event.inaxes:
            ax = event.inaxes
            x, y = event.xdata, event.ydata
            if (event.button == None or event.button == 1):
                if self.line1 != None: # Move line around
                    if self.line2 == None:
                        self.line1.set_data([self.x[0], x],
                                            [self.y[0], y])
                    else:
                        self.line2.set_data([self.x[2], x],
                                            [self.y[2], y])
                    self.fig.canvas.draw_idle()


    def __button_release_callback(self, event):
        if event.inaxes:
            x, y = event.xdata, event.ydata
            ax = event.inaxes
            if event.button == 1:
                if self.line2 == None:  # Segment 1 completed
                    self.x.append(x)
                    self.y.append(y)
                    self.line1.set_data([self.x[0], self.x[1]],
                                        [self.y[0], self.y[1]])
                    #self.line1 = Line2D([self.x[0], self.x[1]],
                    #                    [self.y[0], self.y[1]],
                    #                    marker='o',
                    #                    color=self.color)
                    self.fig.canvas.draw_idle()
                    self.fig.canvas.mpl_disconnect(self.__ID1)
                    self.line2 = 'start'
                else:
                    self.x.append(x)
                    self.y.append(y)
                    # Adjust to the same slope between first and last point
                    m = (self.y[3]-self.y[0])/(self.x[3]-self.x[0])
                    for i in range(4):
                        self.y[i] = self.y[0]+m*(self.x[i]-self.x[0])
                    self.line1.set_data([self.x[0], self.x[1]],
                                        [self.y[0], self.y[1]])
                    self.line2.set_data([self.x[2], self.x[3]],
                                        [self.y[2], self.y[3]])
                    #self.line2 = Line2D([self.x[2], self.x[3]],
                    #                    [self.y[2], self.y[3]],
                    #                    marker='o',
                    #                    color=self.color)
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
            self.ax.lines.remove()



class SegmentsInteractor(QObject):
    """
    An continuum editor.

    """

    showverts = True
    epsilon = 5  # max pixel distance to count as a vertex hit
    mySignal = pyqtSignal(str)
    modSignal = pyqtSignal(str)

    def __init__(self, ax, verts):
        super().__init__()

        from matplotlib.lines import Line2D
        from matplotlib.artist import Artist
        from matplotlib import collections  as mc

        self.ax = ax
        self.type = 'Continuum'

        x, y = zip(*verts)
        self.xy = [(i,j) for (i,j) in zip(x,y)]
        lines = [[(x[0],y[0]),(x[1],y[1])],[(x[2],y[2]),(x[3],y[3])]]
        self.lc = mc.LineCollection(lines, colors = 'g', linewidths=2)
        self.ax.add_collection(self.lc)

        self.canvas = ax.figure.canvas
        self.line = Line2D(x, y, marker='o', linestyle=None, linewidth=0., markerfacecolor='g', animated=True)                
        self.ax.add_line(self.line)

        self.cid = self.lc.add_callback(self.lc_changed)
        self._ind = None  # the active vert
        self.connect()
        self.aperture = self.continuum


    def connect(self):
        self.cid_draw = self.canvas.mpl_connect('draw_event', self.draw_callback)
        self.cid_press = self.canvas.mpl_connect('button_press_event', self.button_press_callback)
        self.cid_release = self.canvas.mpl_connect('button_release_event', self.button_release_callback)
        self.cid_motion = self.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        self.canvas.draw_idle()

    def disconnect(self):
        self.canvas.mpl_disconnect(self.cid_draw)
        self.canvas.mpl_disconnect(self.cid_press)
        self.canvas.mpl_disconnect(self.cid_release)
        self.canvas.mpl_disconnect(self.cid_motion)
        self.lc.remove()
        self.line.remove()
        self.canvas.draw_idle()
        self.aperture = None
        
    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.lc)
        self.ax.draw_artist(self.line)

    def lc_changed(self, lc):
        'this method is called whenever the polygon object is called'
        # only copy the artist props to the line (except visibility)
        vis = self.line.get_visible()
        Artist.update_from(self.line, self.lc)
        self.line.set_visible(vis)  

    def get_ind_under_point(self, event):
        'get the index of the point if within epsilon tolerance'

        x, y = zip(*self.xy)
        d = np.hypot(x - event.xdata, y - event.ydata)
        indseq, = np.nonzero(d == d.min())
        ind = indseq[0]

        if d[ind] >= self.epsilon:
            ind = None

        return ind


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
        x, y = event.xdata, event.ydata
        self.xy[_ind] = (x,y)

        
        # Rebuild line collection
        x,y = zip(*self.xy)
        if _ind == 0:
            m = (y[3]-y[0])/(x[3]-x[0])
        else:
            m = (y[_ind]-y[0])/(x[_ind]-x[0])

        for i in range(4):
            y[i] = y[0]+m*(x[i]-x[0])
            self.xy[i] = (x[i],y[i])

        
        self.updateMarkers()

        self.canvas.restore_region(self.background)
        self.ax.draw_artist(self.lc)
        self.ax.draw_artist(self.line)
        self.canvas.update()
        self.canvas.flush_events()

        # Notify callback
        self.modSignal.emit('continuum guess modified')

    def updateMarkers(self):
        self.line.set_data(zip(*self.xy))
        
