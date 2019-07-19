import sys
import numpy as np
from PyQt5.QtCore import pyqtSignal, QObject
from matplotlib.lines import Line2D
import matplotlib.transforms as transforms
from matplotlib.patches import Rectangle, Polygon
from matplotlib.artist import Artist


def dist_point_to_segment(p, s0, s1):
    """
    Get the distance of a point to a segment.

      *p*, *s0*, *s1* are *xy* sequences

    This algorithm from
    http://geomalgorithms.com/a02-_lines.html
    """
    p = np.asarray(p, float)
    s0 = np.asarray(s0, float)
    s1 = np.asarray(s1, float)
    v = s1 - s0
    w = p - s0

    c1 = np.dot(w, v)
    if c1 <= 0:
        return np.hypot(p, s0)

    c2 = np.dot(v, v)
    if c2 <= c1:
        return np.hypot(p, s1)

    b = c1 / c2
    pb = s0 + b * v
    return np.hypot(p, pb)



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
       
        
class DistanceSelector(QObject):

#    def __init__(self, ax, fig, wcs, callback, color='#7ec0ee'):
    def __init__(self, ax, fig, wcs, callback, color='cornflowerblue'):
        super().__init__()

        self.x = []
        self.y = []
        self.color = color
        self.fig = fig
        self.ax = ax
        x1,x2 = self.ax.get_xlim()
        y1,y2 = self.ax.get_ylim()
        self.dx = np.abs(x2 - x1)/30.
        self.dy = np.abs(y2 - y1)/30.
        self.callback = callback
        self.wcs = wcs

        self.__ID2 = self.fig.canvas.mpl_connect('button_press_event', self.__button_press_callback)
        self.__ID3 = self.fig.canvas.mpl_connect('button_release_event', self.__button_release_callback)

    def __motion_notify_callback(self, event):
        if event.inaxes:
            x, y = event.xdata, event.ydata
            if (event.button == None or event.button == 1):
                self.line1.set_data([self.x[0], x], [self.y[0], y])
                self.xline.set_data([self.x[0], x], [self.y[0], self.y[0]])
                self.yline.set_data([x, x], [self.y[0], y])
                # print('x y ', x, y)
                dx = np.abs(x-self.x[0])
                dy = np.abs(y-self.y[0])
                xm = 0.5 * (self.x[0] + x)
                ym = 0.5 * (self.y[0] + y)
                if y > self.y[0]: # northern
                    self.xlab.set_position((xm, self.y[0] - self.dy))
                else:
                    self.xlab.set_position((xm, self.y[0] + self.dy))
                if x < self.x[0]: # western
                    self.ylab.set_position((x - self.dx, ym))
                else:
                    self.ylab.set_position((x + self.dx, ym))
                rad = np.arctan2(y-self.y[0], x-self.x[0])
                angle = np.degrees(rad)
                if angle > 0 and angle < 90:
                    dx = -np.sin(rad)
                    dy = np.cos(rad)
                elif angle >= 90 and angle < 180:
                    dx = np.sin(rad)
                    dy = -np.cos(rad)
                elif angle <=0 and angle > -90:
                    dx = np.sin(rad)
                    dy = -np.cos(rad)
                else:
                    dx = -np.sin(rad)
                    dy = np.cos(rad)
                self.zlab.set_position((xm + self.dx * dx, ym + self.dx * dy))
                pixel = np.array([[x, y]], np.float_)
                world = self.wcs.wcs_pix2world(pixel, 1)                    
                xx = world[0][0]
                yy = world[0][1]
                if angle > 90:
                    angle -= 180
                elif angle < -90:
                    angle += 180
                self.zlab.set_rotation(angle)
                # distances in degs
                dx = np.abs((xx-self.x0)) * np.cos(self.y0 * np.pi/180.)
                dy = np.abs((yy-self.y0))
                dz = np.sqrt(dx*dx+dy*dy)
                # formatting
                def toDMS(x):
                    sx = ""
                    if x > 1:
                        xd_ = int(x)
                        sx += "{0:d}".format(xd_)+":"
                        xm = (x - xd_) * 60.
                    else:
                        xm = x * 60.
                    if xm > 1:
                        xm_ = int(xm)
                        if x > 1:
                            sx += "{0:02d}".format(xm_)+":"
                        else:
                            sx += "{0:d}".format(xm_)+":"
                        xs = (xm - xm_) * 60.
                    else:
                        xs = xm * 60.
                    if x * 60 > 1:
                        sx += "{0:05.2f}".format(xs)
                    else:
                        sx += "{0:.2f}".format(xs)
                    return sx
                sx = toDMS(dx)
                sy = toDMS(dy)
                sz = toDMS(dz)
                if np.abs(angle) < 3:
                    self.xlab.set_text("")                    
                else:
                    self.xlab.set_text(sx)
                if np.abs(np.abs(angle)-90) < 3:
                    self.ylab.set_text("")
                else:
                    self.ylab.set_text(sy)
                self.zlab.set_text(sz)                
                self.fig.canvas.draw_idle()

    def __button_release_callback(self, event):
        if event.inaxes:
            x, y = event.xdata, event.ydata
            if event.button == 1:
                self.x.append(x)
                self.y.append(y)
                self.line1.set_data([self.x[0], self.x[1]], [self.y[0], self.y[1]])
                self.fig.canvas.draw_idle()
                #self.fig.canvas.mpl_disconnect(self.__ID1)
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
            if event.button == 1:
                self.line1 = Line2D([x, x], [y, y], marker='o', color=self.color)
                self.xline = Line2D([x, x], [y, y], linestyle='--', color=self.color)
                self.yline = Line2D([x, x], [y, y], linestyle='--', color=self.color)
                self.xlab = self.ax.annotate("", xy=(x,y), xycoords='data',
                                             ha='center', va='center', color=self.color)
                self.ylab = self.ax.annotate("", xy=(x,y), xycoords='data',
                                             ha='center', va='center',rotation=90, color=self.color)
                self.zlab = self.ax.annotate("", xy=(x,y), xycoords='data',
                                             ha='center', va='center', color=self.color)
                self.x=[x]
                self.y=[y]
                pixel0 = np.array([[self.x[0], self.y[0]]], np.float_)
                world = self.wcs.wcs_pix2world(pixel0, 1)                    
                self.x0 = world[0][0]
                self.y0 = world[0][1]
                self.__ID1 = self.fig.canvas.mpl_connect('motion_notify_event',
                                                         self.__motion_notify_callback)
                self.ax.add_line(self.line1)
                self.ax.add_line(self.xline)
                self.ax.add_line(self.yline)
                # add a segment
                self.fig.canvas.draw_idle()

    def remove(self):
        """ Remove lines from plot """
        try:
            self.line1.remove()
            self.xline.remove()
            self.yline.remove()
            self.xlab.remove()
            self.ylab.remove()
            self.zlab.remove()
        except:
            print('no lines to remove')


class VoronoiInteractor(QObject):
    """
    Draw a Voronoi tessellation with input sites. 
    Update after adding/removing sites.

    Arguments:
        ax       axes used by the interactor
        fig      figure used by the interactor
        sites    points used as input by the Voronoi tessellation

    Key-bindings:

      'v' toggle voronoi markers on and off.  When on,
          site markers can be moved and deleted

      'd' delete the site under point

      'i' insert a site at point.  You must be within epsilon of the
          line connecting two existing vertices
    """

    mySignal = pyqtSignal(str)
    modSignal = pyqtSignal(str)

    def __init__(self, ax, sites, epsilon=10, showsites=False):
        super().__init__()
        from matplotlib.patches import Polygon
        from matplotlib.lines import Line2D
        # from matplotlib.artist import Artist
        self.ax = ax
        self.type = 'Polygon'
        self.epsilon = epsilon
        self.showsites = showsites
        self.poly = Polygon(list(sites), animated=True, fill=False, closed=True, color='red')
        self.ax.add_patch(self.poly)
        self.canvas = self.poly.figure.canvas

        self.sites = sites
        x, y = zip(*self.sites)
        self.line = Line2D(x, y, marker='o', linestyle=None, linewidth=0.,
                           markerfacecolor='red', animated=True)
        self.ax.add_line(self.line)
        #self.canvas = self.line.figure.canvas
        self.line.set_visible(self.showsites)
        self.poly.set_visible(self.showsites)

        # self.sites = self.sites
        self.sites = sites
        self.drawVoronoi()
        #self.cid = self.poly.add_callback(self.poly_changed)
        self._ind = None  # the active vert
        self.connect()
        
    def drawVoronoi(self):
        from scipy.spatial import Voronoi
        self.sites = self.poly.xy
        vor = Voronoi(self.sites)
        # Ridge segments
        self.segments = []
        for simplex in vor.ridge_vertices:
            simplex = np.asarray(simplex)
            if np.all(simplex >= 0):
                x = vor.vertices[simplex, 0]
                y = vor.vertices[simplex, 1]
                segment = Line2D(x, y, color='black', linestyle='-')#, animated=True)
                self.ax.add_line(segment)
                self.segments.append(segment)
        #center = self.sites.mean(axis=0)
        # Ridge rays
        self.rays = []
        for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
            simplex = np.asarray(simplex)
            if np.any(simplex < 0):
                #i = simplex[simplex >= 0][0] # finite end Voronoi vertex
                t = self.sites[pointidx[1]] - self.sites[pointidx[0]]  # tangent
                t = t / np.linalg.norm(t)
                #n = np.array([-t[1], t[0]]) # normal
                #midpoint = self.sites[pointidx].mean(axis=0)
                #far_point = vor.vertices[i] + np.sign(np.dot(midpoint - center, n)) * n * 100
                #x = [vor.vertices[i,0], far_point[0]]
                #y = [vor.vertices[i,1], far_point[1]]
                #ray = Line2D(x, y, color='black', linestyle='--')#, animated=True)
                #self.ax.add_line(ray)
                #self.rays.append(ray)

    def connect(self):
        self.cid_draw = self.canvas.mpl_connect('draw_event', self.draw_callback)
        self.cid_press = self.canvas.mpl_connect('button_press_event', self.button_press_callback)
        self.cid_key = self.canvas.mpl_connect('key_press_event', self.key_press_callback)
        self.cid_release = self.canvas.mpl_connect('button_release_event', self.button_release_callback)
        self.cid_motion = self.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        self.canvas.draw_idle()

    def removeRidges(self):
        for s in self.segments + self.rays:
            s.remove()
        self.poly.remove()
        self.line.remove()
        self.canvas.draw_idle()

    def disconnect(self):
        self.canvas.mpl_disconnect(self.cid_draw)
        self.canvas.mpl_disconnect(self.cid_press)
        self.canvas.mpl_disconnect(self.cid_key)
        self.canvas.mpl_disconnect(self.cid_release)
        self.canvas.mpl_disconnect(self.cid_motion)
        self.aperture = None
        
    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.poly)
        self.ax.draw_artist(self.line)
        #for s in self.segments + self.rays:
        #    s.remove()
        #self.drawVoronoi()
        #self.canvas.blit(self.ax.bbox)
        #self.canvas.udpate()
        #self.canvas.flush_events()

    def poly_changed(self, poly):
        'this method is called whenever the polygon object is called'
        # only copy the artist props to the line (except visibility)
        vis = self.line.get_visible()
        Artist.update_from(self.line, poly)
        #for s in self.segments + self.rays:
        #    s.remove()
        #self.drawVoronoi()
        self.line.set_visible(vis)  # don't use the poly visibility state
        # print('poly changed')

    def get_ind_under_point(self, event):
        'get the index of the vertex under point if within epsilon tolerance'
        # display coords
        xy = np.asarray(self.poly.xy)
        xyt = self.poly.get_transform().transform(xy)
        xt, yt = xyt[:, 0], xyt[:, 1]
        d = np.hypot(xt - event.x, yt - event.y)
        indseq, = np.nonzero(d == d.min())
        ind = indseq[0]
        if d[ind] >= self.epsilon:
            ind = None
        return ind

    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        if not self.showsites:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        self._ind = self.get_ind_under_point(event)

    def button_release_callback(self, event):
        'whenever a mouse button is released'
        if not self.showsites:
            return
        if event.button != 1:
            return
        self._ind = None
        # For some reason, 't' inhibites the changes in Voronoi
        for s in self.segments + self.rays:
            s.remove()
        self.drawVoronoi()
        # Notify callback
        self.modSignal.emit('voronoi modified')

    def key_press_callback(self, event):
        'whenever a key is pressed'
        if not event.inaxes:
            return
        if event.key == 'v':
            self.showsites = not self.showsites
            self.line.set_visible(self.showsites)
            self.poly.set_visible(self.showsites)
            if not self.showsites:
                self._ind = None
            self.canvas.draw_idle()
        elif event.key == 'd':
            ind = self.get_ind_under_point(event)
            if ind is not None:
                if len(self.poly.xy) < 4:  # the minimum polygon has 4 points since the 1st is repeated as final
                    # Delete polygon
                    #self.disconnect()
                    #self.poly = None
                    #self.line = None
                    self.mySignal.emit('voronoi deleted')
                else:
                    #self.poly.xy = np.array([tup
                    #                for i, tup in enumerate(self.poly.xy)
                    #                if i != ind])
                    a = self.poly.get_xy().tolist()
                    del a[ind]
                    # If first point, remove it also at the end and close the loop
                    if ind == 0:
                        a.insert(-1,a[0])
                        del a[-1]
                    self.poly.set_xy(a)
                    print('new lengths ', len(a), len(self.poly.xy))
                    self.line.set_data(zip(*self.poly.xy))
                    self.modSignal.emit(str(ind))  # Pass the index canceled
                for s in self.segments + self.rays:
                    s.remove()
                self.drawVoronoi()
        elif event.key == 'i':
            xys = self.poly.get_transform().transform(self.poly.xy)
            p = event.x, event.y  # display coords
            for i in range(len(xys) - 1):
                s0 = xys[i]
                s1 = xys[i + 1]
                d = dist_point_to_segment(p, s0, s1)
                if d <= self.epsilon:
                    a = list(self.poly.get_xy())
                    a.insert(i, (event.xdata, event.ydata))
                    self.poly.set_xy(a)
                    #self.poly.xy = np.array(
                    #    list(self.poly.xy[:i+1]) +
                    #    [(event.xdata, event.ydata)] +
                    #    list(self.poly.xy[i+1:]))
                    self.line.set_data(zip(*self.poly.xy))
                    break
            self.modSignal.emit('one voronoi site added')
            for s in self.segments + self.rays:
                s.remove()
            self.drawVoronoi()

    def motion_notify_callback(self, event):
        'on mouse movement'
        if not self.showsites:
            return
        if self._ind is None:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        x, y = event.xdata, event.ydata

        self.poly.xy[self._ind] = x, y
        if self._ind == 0:
            self.poly.xy[-1] = x, y
        elif self._ind == len(self.poly.xy) - 1:
            self.poly.xy[0] = x, y
        self.updateMarkers()
        self.canvas.restore_region(self.background)
        self.ax.draw_artist(self.poly)
        self.ax.draw_artist(self.line)
        self.canvas.update()
        self.canvas.flush_events()

    def updateMarkers(self):
        self.line.set_data(zip(*self.poly.xy))
        self.sites = self.poly.xy

class LineInteractor(QObject):
    """
    A Gaussian line interactor.
    Arguments:
        ax axes to which the interactor is associated
        c0 intercept of continuum
        cs slope of continuum
        x0 center of the line
        A  amplitude of the line
        fwhm FWHM of the line
        epsilon max pixel distance to count as a vertex hit
    """
    showverts = True
    mySignal = pyqtSignal(str)
    modSignal = pyqtSignal(str)

    def __init__(self, ax, c0, cs, x0, A, fwhm, n, color='#7ec0ee', epsilon=10):
        super().__init__()
        # To avoid crashing with maximum recursion depth exceeded
        sys.setrecursionlimit(10000)  # 10000 is 10x the default value
        self.epsilon = epsilon
        self.n = n  # ID of the line
        self.ax = ax
        self.fig = ax.figure
        self.canvas = ax.figure.canvas
        self.c0 = c0  # Value of continuum at origin
        self.cs = cs  # Slope of the continuum
        self.type = 'Line'
        self.color = color
        self.x0 = x0
        self.A = A
        self.fwhm = fwhm
        self.computeMarkers()
        self.computeGaussian()
        self.gauss = Polygon(self.verts, animated=True, fill=False, closed=False, color=self.color)
        self.ax.add_patch(self.gauss)
        x, y = zip(*self.xy)
        self.line = Line2D(x, y, marker='o', linestyle=None, linewidth=0.,
                           markerfacecolor=self.color, animated=True)
        self.ax.add_line(self.line)
        self.artists = [self.line, self.gauss]
        self._ind = None  # the active vert
        self.connect()

    def computeMarkers(self):
        'Compute position of markers.'
        x = self.x0 + 0.5 * self.fwhm * np.array([-1, 0, 1])
        y = self.c0 + (x - self.x0) * self.cs + self.A * np.array([0.5, 1., 0.5])
        self.xy = [(i, j) for (i, j) in zip(x, y)]

    def computeGaussian(self):
        'Compute the Gaussian polygon from the position of the markers.'
        self.sigma = self.fwhm / (2 * np.sqrt(2 * np.log(2)))
        # Create an array of x values and compute the value of the Gaussian on it
        x = np.linspace(self.x0 - self.fwhm, self.x0 + self.fwhm, 30)
        dx = (x - self.x0) / self.sigma / np.sqrt(2.)
        y = self.c0 + (x - self.x0) * self.cs + self.A * np.exp(-dx * dx)
        self.verts = [(x_, y_) for x_, y_ in zip(x, y)]
        return

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
            self.line.remove()
        except BaseException:
            print('no markers')
        try:
            self.gauss.remove()
        except BaseException:
            print('no line')
        self.canvas.draw_idle()

    def draw_callback(self, event):
        self.grab_background()
        self.ax.draw_artist(self.gauss)
        self.ax.draw_artist(self.line)

    def safe_draw(self):
        """Temporarily disconnect the draw_event callback to avoid recursion."""
        self.canvas.mpl_disconnect(self.cid_draw)
        self.canvas.draw_idle()
        self.cid_draw = self.canvas.mpl_connect('draw_event', self.draw_callback)

    def grab_background(self):
        self.safe_draw()
        self.background = self.canvas.copy_from_bbox(self.fig.bbox)  # or self.ax.bbox

    def get_ind_under_point(self, event):
        'get the index of the point if within epsilon tolerance'
        # Distance is computed in pixels on the screen
        xy = self.ax.transData.transform(self.xy)
        x, y = zip(*xy)
        x = np.array(x)
        y = np.array(y)
        d = np.hypot(x - event.x, y - event.y)
        indseq, = np.nonzero(d == d.min())
        ind = indseq[0]
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
                self.mySignal.emit('line deleted ' + str(self.n))
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
        # Redrawing
        self.canvas.draw_idle()

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
        # Update markers and Gaussian parameters
        x_, y_ = event.xdata, event.ydata
        x, y = zip(*self.xy)
        if x[-1] > x[0]:  # case wavelength
            if self._ind == 0:
                if x_ < x[1]:
                    self.fwhm = 2 * (x[1] - x_)
            elif self._ind == 1:
                dx = x_ - x[1]
                self.x0 += dx
                dy = y_ - y[1]
                if (self.A > 0) & (dy < -self.A):  # Emission line
                    pass
                elif (self.A < 0) & (dy > -self.A):  # Absorption line
                    pass
                else:
                    self.A += dy
            elif self._ind == 2:
                if x_ > x[1]:
                    self.fwhm = 2 * (x_ - x[1])
        else:
            if self._ind == 0:
                if x_ > x[1]:
                    self.fwhm = 2 * (x_ - x[1])
            elif self._ind == 1:
                dx = x_ - x[1]
                self.x0 += dx
                dy = y_ - y[1]
                if (self.A > 0) & (dy < -self.A):  # Emission line
                    pass
                elif (self.A < 0) & (dy > -self.A):  # Absorption line
                    pass
                else:
                    self.A += dy
            elif self._ind == 2:
                if x_ < x[1]:
                    self.fwhm = 2 * (x[1] - x_)
        self.updateCurves()
        self.redraw()
        # Notify callback
        self.modSignal.emit('line guess modified ' + str(self.n))

    def updateCurves(self):
        self.computeGaussian()
        self.computeMarkers()
        self.line.set_data(zip(*self.xy))
        self.gauss.xy = self.verts
        
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
        self.xy = [(i, j) for (i, j) in zip(x, y)]
        self.line.set_data(zip(*self.xy))
        xy = self.gauss.get_xy()
        x, y = zip(*xy)
        x = np.asarray(x)
        y = np.asarray(y)
        x = c/x * 1.e-6
        xy = [(i, j) for (i, j) in zip(x, y)]
        self.gauss.set_xy(xy)
        # change line parameters
        self.x0 = c / self.x0 * 1.e-6
        self.fwhm *= c / self.x0**2 * 1.e-6

    
class PsfInteractor(QObject):

    epsilon = 5
    showverts = True
    mySignal = pyqtSignal(str)
    modSignal = pyqtSignal(str)
   
    def __init__(self, ax, center, radius, annulus=None):
        super().__init__()
        from matplotlib.patches import Circle
        from matplotlib.lines import Line2D
        # from matplotlib.artist import Artist
        # To avoid crashing with maximum recursion depth exceeded
        import sys
        sys.setrecursionlimit(10000) # 10000 is 10x the default value

        if annulus is None:
            annulus = radius * 0.3
        self.ax = ax
        self.innerCircle = Circle(center, radius, edgecolor='Lime', facecolor='none',
                                   angle=0, fill=False, animated=True)
        self.ax.add_patch(self.innerCircle)
        self.canvas = self.innerCircle.figure.canvas
        self.outerCircle = Circle(center, radius + annulus, edgecolor='Lime',
                                   facecolor='none', angle=0, fill=False,
                                   linestyle='--', animated=True)
        self.ax.add_patch(self.outerCircle)

        # Create a line with center, width, and height points
        self.center = self.innerCircle.center
        self.inRadius = self.innerCircle.radius
        self.outRadius = self.outerCircle.radius

        x0, y0 = self.center
        r0 = self.inRadius
        r1 = self.outRadius

        x = [x0, x0+r0, x0+r1]
        y = [y0, y0, y0]
        self.xy = [(i,j) for (i,j) in zip(x,y)]
        self.line = Line2D(x, y, marker='o', linestyle=None, linewidth=0.,
                           markerfacecolor='g', animated=True)
        self.ax.add_line(self.line)
        self.cid = self.innerCircle.add_callback(self.circles_changed)
        self._ind = None  # the active point
        self.connect()
        self.press = None
        self.lock = None
  
    def connect(self):
        self.cid_draw = self.canvas.mpl_connect('draw_event', self.draw_callback)
        self.cid_press = self.canvas.mpl_connect('button_press_event', self.button_press_callback)
        self.cid_release = self.canvas.mpl_connect('button_release_event', self.button_release_callback)
        self.cid_motion = self.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        self.cid_key = self.canvas.mpl_connect('key_press_event', self.key_press_callback)
        self.canvas.draw_idle()
  
    def disconnect(self):
        self.canvas.mpl_disconnect(self.cid_draw)
        self.canvas.mpl_disconnect(self.cid_press)
        self.canvas.mpl_disconnect(self.cid_release)
        self.canvas.mpl_disconnect(self.cid_motion)
        self.canvas.mpl_disconnect(self.cid_key)
        self.innerCircle.remove()
        self.outerCircle.remove()
        self.line.remove()
        self.canvas.draw_idle()
        self.aperture = None
        
    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.innerCircle)
        self.ax.draw_artist(self.outerCircle)
        self.ax.draw_artist(self.line)

    def circles_changed(self, ellipse):
        'this method is called whenever the polygon object is called'
        # only copy the artist props to the line (except visibility)
        vis = self.line.get_visible()
        Artist.update_from(self.line, ellipse)
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
        x0, y0 = self.innerCircle.center
        r0 = self.inRadius
        r1 = self.outRadius
        # print('center ', x0, y0)
        # print('radii ' ,r0, r1)
        self.press = x0, y0, r0, r1, event.xdata, event.ydata
        self.xy0 = self.xy
        self.lock = "pressed"

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
            #self.disconnect()
            #self.ellipse = None
            #self.line = None
            self.mySignal.emit('psf deleted')
        self.canvas.draw_idle()

    def button_release_callback(self, event):
        'whenever a mouse button is released'
        if not self.showverts:
            return
        if event.button != 1:
            return
        self._ind = None
        self.press = None
        self.lock = "released"
        self.background = None
        # To get other aperture redrawn
        self.canvas.draw_idle()

    def motion_notify_callback(self, event):
        'on mouse movement'
        #if not self.ellipse.contains(event): return
        if not self.showverts:
            return
        if self._ind is None:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        x0, y0, r0, r1, xpress, ypress = self.press
        self.dx = event.xdata - xpress
        self.dy = event.ydata - ypress
        self.update_circles()
        # Redraw ellipse and points
        self.canvas.restore_region(self.background)
        self.ax.draw_artist(self.innerCircle)
        self.ax.draw_artist(self.outerCircle)
        self.ax.draw_artist(self.line)
        self.canvas.update()
        self.canvas.flush_events()
        # Notify callback
        self.modSignal.emit('circles modified')

    def update_circles(self):
        x0, y0, r0, r1, xpress, ypress = self.press
        dx, dy = self.dx, self.dy        
        if self.lock == "pressed":
            if self._ind == 0:
                self.lock = "move"
            else:
                self.lock = "resize"
        elif self.lock == "move":
            if x0+dx < 0:
                xn = x0
                dx = 0
            else:
                xn = x0+dx
            if y0+dy < 0:
                yn = y0
                dy = 0
            else:
                yn = y0+dy
            self.innerCircle.center = xn,yn
            self.outerCircle.center = xn,yn
            # update line
            self.xy = [(i+dx,j+dy) for (i,j) in self.xy0]
            # Redefine line
            self.line.set_data(zip(*self.xy))
        # otherwise rotate and resize
        elif self.lock == 'resize':
            # Avoid to shrink            
            if self._ind == 1:
                if (r0 + dx) > 0:
                    r0 += dx
            elif self._ind == 2:
                if (r1 + dx) > r0:
                    r1 += dx
            if r1 < r0:
                r1 = r0 * 1.2
            # update circles
            self.innerCircle.set_radius(r0)
            self.outerCircle.set_radius(r1)
            self.inRadius = r0
            self.outRadius = r1
            # update points
            self.updateMarkers()
            
    def updateInteractor(self):
        # Redraw ellipse and points
        self.updateMarkers()        
        self.canvas.restore_region(self.background)
        self.ax.draw_artist(self.innerCircle)
        self.ax.draw_artist(self.outerCircle)
        self.ax.draw_artist(self.line)
        self.canvas.update()
        self.canvas.flush_events()

    def updateMarkers(self):
        # update points
        x0,y0  = self.innerCircle.center
        r0 = self.innerCircle.radius
        r1 = self.outerCircle.radius
        x = [x0, x0+r0, x0+r1]
        y = [y0, y0, y0]
        self.xy = [(i,j) for (i,j) in zip(x,y)]
        self.line.set_data(x,y)
