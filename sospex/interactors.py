import sys
import numpy as np
from PyQt5.QtCore import pyqtSignal, QObject
from matplotlib.lines import Line2D
import matplotlib.transforms as transforms
from matplotlib.patches import Rectangle


class SliderInteractor(QObject):
    """
    A shade with a top slider is shown on the plot.
    Arguments:
        ax      axes associated with the interactor
        x       cursor position of the interactor
        dx      width of the associated range
        epsilon max pixel distance to count as a cursor hit 
    """
    modSignal = pyqtSignal(str)
    
    def __init__(self, ax, x, dx, epsilon=10):
        super().__init__()
        # To avoid crashing with maximum recursion depth exceeded
        sys.setrecursionlimit(10000)  # 10000 is 10x the default value
        self.x = x
        self.dx = dx
        self.epsilon = epsilon
        self.ax = ax
        self.canvas = self.ax.figure.canvas
        # color = '#A9BA9D'  # Laurel green
        # self.region = self.ax.axvspan(x-dx/2., x+dx/2., facecolor='lightgreen', alpha=0.5,
        #                              linewidth=0, zorder=1, animated=True)
        trans = transforms.blended_transform_factory(self.ax.transData, self.ax.transAxes)
        self.region = Rectangle((x - dx * 0.5, 0.), dx, 1, facecolor='lightgreen', alpha=0.3,
                                zorder=1, animated=True, transform=trans)
        # x coords are data, y coords are axes (points stays in the same spot when zooming)
        self.line = Line2D([x], [0.98], marker='^', markersize=10, linestyle=None, linewidth=0,
                           markerfacecolor='lime', animated=True, transform=trans)
        self.ax.add_line(self.line)
        self.ax.add_patch(self.region)
        self.activated = False
        self.connect()
        
    def connect(self):
        self.cid_draw = self.canvas.mpl_connect('draw_event', self.draw_callback)
        self.cid_press = self.canvas.mpl_connect('button_press_event', self.button_press_callback)
        self.cid_release = self.canvas.mpl_connect('button_release_event',
                                                   self.button_release_callback)
        self.cid_motion = self.canvas.mpl_connect('motion_notify_event',
                                                  self.motion_notify_callback)
        # self.cid_key = self.canvas.mpl_connect('key_press_event', self.key_press_callback)
        self.canvas.draw_idle()

    def disconnect(self):
        self.canvas.mpl_disconnect(self.cid_draw)
        self.canvas.mpl_disconnect(self.cid_press)
        self.canvas.mpl_disconnect(self.cid_release)
        self.canvas.mpl_disconnect(self.cid_motion)
        # self.canvas.mpl_disconnect(self.cid_key)
        try:
            self.region.remove()
        except BaseException:
            print('no shaded region')
        try:
            self.line.remove()
        except BaseException:
            print('no markers')
        self.canvas.draw_idle()

    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.line)
        self.ax.draw_artist(self.region)

    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        print('event ', event.inaxes, event.button)
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        self.activate(event)

    def button_release_callback(self, event):
        'whenever a mouse button is released'
        if event.button != 1:
            return
        self.activated = False

    def activate(self, event):
        'activate if within epsilon tolerance from the point'
        # Distance is computed in pixels on the screen
        xy = self.ax.transData.transform((self.x,0))
        x = xy[0]
        x = np.array(x)
        d = np.abs(x - event.x)
        # print('d is ', d)
        if d <= self.epsilon:
            self.activated = True

    def motion_notify_callback(self, event):
        'on mouse movement'
        if self.activated is False:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        self.canvas.restore_region(self.background)
        # Update marker
        x_ = event.xdata
        dx = x_ - self.x
        self.x = x_
        self.redraw(dx)
        # Notify callback to update spatial image
        self.modSignal.emit('slider moved')
        
    def redraw(self, dx):        
        # Update marker
        self.line.set_xdata([self.x])
        # Update shade (set bottom corner x position)
        self.region.set_x(self.x - self.dx/2. + dx)
        self.region.set_width(self.dx)
        # Update segments and markers
        self.ax.draw_artist(self.region)
        self.ax.draw_artist(self.line)
        self.canvas.update()
        self.canvas.flush_events()
        
        
class SliceInteractor(QObject):   
    """
    A shade with two top sliders are shown on the plot.
    Arguments:
        ax      axes associated with the interactor
        xl      left end of the sliced part shaded
        xr      right end of the sliced part shaded
        epsilon max pixel distance to count as a cursor hit 
    """
    modSignal = pyqtSignal(str)
    
    def __init__(self, ax, xl, xr, epsilon=10):
        super().__init__()
        # To avoid crashing with maximum recursion depth exceeded
        sys.setrecursionlimit(10000)  # 10000 is 10x the default value
        self.xl = xl
        self.xr = xr
        self.epsilon = epsilon
        self.ax = ax
        self.canvas = self.ax.figure.canvas
        #self.region = self.ax.axvspan(xl, xr, facecolor='lightgreen', alpha=0.5,
        #                              linewidth=0, zorder=1, animated=True)
        # x coords are data, y coords are axes (points stays in the same spot when zooming)
        trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
        self.region = Rectangle((xl, 0.), xr-xl, 1, facecolor='lightgreen', alpha=0.3,
                                zorder=1, animated=True, transform=trans)
        self.line = Line2D([xl, xr], [0.98, 0.98], marker='^', markersize=10, linestyle=None,
                           linewidth=0,
                           markerfacecolor='lime', animated=True, transform=trans)
        self.ax.add_line(self.line)
        self.ax.add_patch(self.region)
        self.activated = False
        self.ind = None
        self.connect()

    def connect(self):
        self.cid_draw = self.canvas.mpl_connect('draw_event', self.draw_callback)
        self.cid_press = self.canvas.mpl_connect('button_press_event', self.button_press_callback)
        self.cid_release = self.canvas.mpl_connect('button_release_event',
                                                   self.button_release_callback)
        self.cid_motion = self.canvas.mpl_connect('motion_notify_event',
                                                  self.motion_notify_callback)
        self.canvas.draw_idle()

    def disconnect(self):
        self.canvas.mpl_disconnect(self.cid_draw)
        self.canvas.mpl_disconnect(self.cid_press)
        self.canvas.mpl_disconnect(self.cid_release)
        self.canvas.mpl_disconnect(self.cid_motion)
        try:
            self.region.remove()
        except BaseException:
            print('no shaded region')
        try:
            self.line.remove()
        except BaseException:
            print('no markers')
        self.canvas.draw_idle()

    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.line)
        self.ax.draw_artist(self.region)

    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        # print('event ', event.inaxes, event.button)
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        self.activate(event)
        
    def button_release_callback(self, event):
        'whenever a mouse button is released'
        if event.button != 1:
            return
        self.activated = False

    def activate(self, event):
        'identify the endpoint if within epsilon tolerance'
        # Distance is computed in pixels on the screen
        xr = self.ax.transData.transform((self.xr,0))
        xl = self.ax.transData.transform((self.xl,0))
        xr = xr[0]
        xl = xl[0]
        d = np.array([np.abs(xl-event.x), np.abs(xr-event.x)])
        indseq, = np.nonzero(d == d.min())
        ind = indseq[0]
        if d[ind] >= self.epsilon:
            self.ind = None
            self.activated = False
        else:
            self.activated = True
            self.ind = ind

    def motion_notify_callback(self, event):
        'on mouse movement'
        if self.activated is False:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        # Update endpoints
        x_ = event.xdata
        if self.ind == 0:
            if x_ > self.xr:
                return
            else:
                self.xl = x_
        elif self.ind == 1:
            if x_ < self.xl:
                return
            else:
                self.xr = x_
        else:
            return
        self.redraw()
        # Notify callback to update spatial image
        self.modSignal.emit('slicer changed')
        
    def redraw(self):
        self.canvas.restore_region(self.background)
        # Update marker
        self.line.set_xdata([self.xl, self.xr])
        # Update shade (set bottom corner x position)
        self.region.set_x(self.xl)
        self.region.set_width(self.xr - self.xl)
        # Update segments and markers
        self.ax.draw_artist(self.region)
        self.ax.draw_artist(self.line)
        self.canvas.update()
        self.canvas.flush_events()
       
