import sys
import numpy as np
from PyQt5.QtCore import pyqtSignal, QObject
from matplotlib.lines import Line2D
import matplotlib.transforms as transforms


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
        # xlims = self.ax.get_xlim()
        ylims = self.ax.get_ylim()
        dy = np.abs(ylims[1]-ylims[0])
        y = ylims[1]-dy/40.
        print('marker at ', x, y)
        # color = '#A9BA9D'  # Laurel green
        self.region = self.ax.axvspan(x-dx/2., x+dx/2., facecolor='lightgreen', alpha=0.5,
                                      linewidth=0, zorder=1, animated=True)
        # x coords are data, y coords are axes (points stays in the same spot when zooming)
        trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
        self.line = Line2D([x], [0.98], marker='^', markersize=10, linestyle=None, linewidth=0,
                           markerfacecolor='lime', animated=True, transform=trans)
        self.ax.add_line(self.line)
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
        'get the index of the point if within epsilon tolerance'
        # Distance is computed in pixels on the screen
        xy = self.ax.transData.transform((self.x,0))
        x = xy[0]
        x = np.array(x)
        d = np.abs(x - event.x)
        print('d is ', d)
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
        self.line.set_xdata([x_])
        # Update shade (set bottom corner x position)
        xy = self.region.get_xy()
        x, y = zip(*xy)
        xy = [(x_ + dx,y_) for x_,y_ in zip(x,y)]
        self.region.set_xy(xy)
        # Update segments and markers
        self.ax.draw_artist(self.region)
        self.ax.draw_artist(self.line)
        self.canvas.update()
        self.canvas.flush_events()
        # Notify callback to update spatial image
        self.modSignal.emit('slider moved')