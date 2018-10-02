import numpy as np
from PyQt5.QtCore import pyqtSignal,QObject
from matplotlib.artist import Artist
from matplotlib.mlab import dist_point_to_segment

class photoAperture(QObject):
    """ 
    Photometric aperture 
    
    Data are conserved in RA, Dec units to be passed to
    images with different astrometry
    """

    def __init__(self, n, type, data):
        super().__init__()
        self.n = n
        #self.label = 'Ap'+{:d}.format(n)
        self.type = type

        if type == 'Polygon':
            self.verts = data
        else:
            self.center = (data[0],data[1])
            self.angle = 0.
            if type == 'Square':
                self.width = data[2]
                self.heigth = data[2]
            elif type == 'Rectangle' or type =='Ellipse':
                self.width = data[2]
                self.height = data[3]
            elif type == 'Circle':
                self.radius = data[2]


class PixelInteractor(QObject):

    epsilon = 10
    showverts = True
    mySignal = pyqtSignal(str)
    modSignal = pyqtSignal(str)

    
    def __init__(self,ax,corner,width,angle=0.):
        super().__init__()
        from matplotlib.patches import Rectangle
        from matplotlib.lines import Line2D
        # from matplotlib.artist import Artist
        # To avoid crashing with maximum recursion depth exceeded
        import sys
        sys.setrecursionlimit(10000) # 10000 is 10x the default value

        self.type = 'Pixel'
        height = width
        self.ax = ax
        self.angle  = angle
        self.width  = width
        self.height = width
        self.rect = Rectangle(corner,width,height,edgecolor='Lime',facecolor='none',angle=angle,fill=False,animated=True)
        self.ax.add_patch(self.rect)
        self.canvas = self.rect.figure.canvas

        x,y = self.compute_markers()
        self.line = Line2D(x, y, marker='s', linestyle=None, linewidth=0., markerfacecolor='g', animated=True)
        self.ax.add_line(self.line)

        self.cid = self.rect.add_callback(self.rectangle_changed)
        self._ind = None  # the active point

        self.connect()

        self.aperture = self.rect
        self.press = None
        self.lock = None


    def compute_markers(self):

        # theta0 = self.rect.angle / 180.*np.pi
        w0 = self.rect.get_width()
        # h0 = self.rect.get_height()
        x0,y0 = self.rect.get_xy()
        angle0 = self.rect.angle

        x = [x0+w0/np.sqrt(2.)*np.sin((45.-angle0)*np.pi/180.)]
        y = [y0+w0/np.sqrt(2.)*np.cos((45.-angle0)*np.pi/180.)]

        self.xy = [(x,y)]
        return x, y

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
        self.rect.remove()
        self.line.remove()
        self.canvas.draw_idle()
        self.aperture = None
        
    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.rect)
        self.ax.draw_artist(self.line)


    def rectangle_changed(self, rect):
        'this method is called whenever the polygon object is called'
        # only copy the artist props to the line (except visibility)
        vis = self.line.get_visible()
        Artist.update_from(self.line, rect)
        self.line.set_visible(vis)  

        
    def get_ind_under_point(self, event):
        'get the index of the point if within epsilon tolerance'

        x, y = self.xy[0]
        d = np.hypot(x - event.xdata, y - event.ydata)

        if d >= self.epsilon:
            ind = None
        else:
            ind = 0
            
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
        x0, y0 = self.rect.get_xy()
        w0, h0 = self.rect.get_width(), self.rect.get_height()
        theta0 = self.rect.angle/180*np.pi
        self.press = x0, y0, w0, h0, theta0, event.xdata, event.ydata
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
            self.mySignal.emit('rectangle deleted')

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

        if not self.showverts:
            return
        if self._ind is None:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return

        x0, y0, w0, h0, theta0, xpress, ypress = self.press
        self.dx = event.xdata - xpress
        self.dy = event.ydata - ypress
        self.update_rectangle()

        # Redraw rectangle and points
        self.canvas.restore_region(self.background)
        self.ax.draw_artist(self.rect)
        self.ax.draw_artist(self.line)
        self.canvas.update()
        self.canvas.flush_events()

        # Notify callback
        self.modSignal.emit('rectangle modified')

    def update_rectangle(self):

        x0, y0, w0, h0, theta0, xpress, ypress = self.press
        dx, dy = self.dx, self.dy
        
        if self.lock == "pressed":
            self.lock = "move"
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
            self.rect.set_xy((xn,yn))
            # update line
            self.xy = [(i+dx,j+dy) for (i,j) in self.xy0]
            # Redefine line
            self.line.set_data(zip(*self.xy))
            self.updateMarkers()

    def updateMarkers(self):
        # update points
        x,y = self.compute_markers()
        self.line.set_data(x,y)

    
class EllipseInteractor(QObject):

    epsilon = 5
    showverts = True
    mySignal = pyqtSignal(str)
    modSignal = pyqtSignal(str)
   
    def __init__(self,ax,center,width,height=None,angle=0.):
        super().__init__()
        from matplotlib.patches import Ellipse
        from matplotlib.lines import Line2D
        # from matplotlib.artist import Artist
        # To avoid crashing with maximum recursion depth exceeded
        import sys
        sys.setrecursionlimit(10000) # 10000 is 10x the default value

        if height is None:
            self.type = 'Circle'
            height = width
        else:
            self.type = 'Ellipse'
        self.ax = ax
        self.ellipse = Ellipse(center,width,height,edgecolor='Lime',facecolor='none',angle=angle,fill=False,animated=True)
        self.ax.add_patch(self.ellipse)
        self.canvas = self.ellipse.figure.canvas

        # Create a line with center, width, and height points
        self.center = self.ellipse.center
        self.angle = self.ellipse.angle*180./np.pi
        self.width = self.ellipse.width
        self.height = self.ellipse.height

        ca = np.cos(self.angle); sa = np.sin(self.angle)
        x0, y0 = self.center
        w0 = self.width*0.5
        h0 = self.height*0.5

        x = [x0, x0+w0*ca, x0-h0*sa]
        y = [y0, y0+w0*sa, y0+h0*ca]
        self.xy = [(i,j) for (i,j) in zip(x,y)]
        self.line = Line2D(x, y, marker='o', linestyle=None, linewidth=0., markerfacecolor='g', animated=True)
        self.ax.add_line(self.line)


        self.cid = self.ellipse.add_callback(self.ellipse_changed)
        self._ind = None  # the active point

        self.connect()

        self.aperture = self.ellipse
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
        self.ellipse.remove()
        self.line.remove()
        self.canvas.draw_idle()
        self.aperture = None
        
    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.ellipse)
        self.ax.draw_artist(self.line)

    def ellipse_changed(self, ellipse):
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
        x0, y0 = self.ellipse.center
        w0, h0 = self.ellipse.width, self.ellipse.height
        theta0 = self.ellipse.angle/180*np.pi
        self.press = x0, y0, w0, h0, theta0, event.xdata, event.ydata
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
            self.mySignal.emit('ellipse deleted')
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
        x0, y0, w0, h0, theta0, xpress, ypress = self.press
        self.dx = event.xdata - xpress
        self.dy = event.ydata - ypress
        self.update_ellipse()
        # Redraw ellipse and points
        self.canvas.restore_region(self.background)
        self.ax.draw_artist(self.ellipse)
        self.ax.draw_artist(self.line)
        self.canvas.update()
        self.canvas.flush_events()
        # Notify callback
        self.modSignal.emit('ellipse modified')

    def update_ellipse(self):
        x0, y0, w0, h0, theta0, xpress, ypress = self.press
        dx, dy = self.dx, self.dy        
        if self.lock == "pressed":
            if self._ind == 0:
                self.lock = "move"
            else:
                self.lock = "resizerotate"
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
            self.ellipse.center = xn,yn
            # update line
            self.xy = [(i+dx,j+dy) for (i,j) in self.xy0]
            # Redefine line
            self.line.set_data(zip(*self.xy))
        # otherwise rotate and resize
        elif self.lock == 'resizerotate':
            dtheta = np.arctan2(ypress+dy-y0,xpress+dx-x0)-np.arctan2(ypress-y0,xpress-x0)
            theta_ = (theta0+dtheta) * 180./np.pi
            c, s = np.cos(theta0), np.sin(theta0)
            R = np.matrix('{} {}; {} {}'.format(c, s, -s, c))
            (dx_,dy_), = np.array(np.dot(R,np.array([dx,dy])))
            # Avoid to pass through the center            
            if self._ind == 1:
                w_ = w0+2*dx_  if (w0+2*dx_) > 0 else w0
                if self.type == 'Circle':
                    h_ = w_
                else:
                    h_ = h0
            elif self._ind == 2:
                h_ = h0+2*dy_  if (h0+2*dy_) > 0 else h0
                if self.type == 'Circle':
                    w_ = h_
                else:
                    w_ = w0
            # update ellipse
            self.ellipse.width = w_
            self.ellipse.height = h_
            if self.type == 'Circle':
                theta_ = 0.
            self.ellipse.angle = theta_
            # update points
            self.updateMarkers()

    def updateMarkers(self):
        # update points
        theta_ = self.ellipse.angle*np.pi/180.
        x0,y0  = self.ellipse.center
        w_ = self.ellipse.width
        h_ = self.ellipse.height
        ca = np.cos(theta_); sa = np.sin(theta_)
        x = [x0, x0+w_*0.5*ca, x0-h_*0.5*sa]
        y = [y0, y0+w_*0.5*sa, y0+h_*0.5*ca]
        self.xy = [(i,j) for (i,j) in zip(x,y)]
        self.line.set_data(x,y)


class RectangleInteractor(QObject):

    epsilon = 5
    showverts = True
    mySignal = pyqtSignal(str)
    modSignal = pyqtSignal(str)
    
    def __init__(self,ax,corner,width,height=None,angle=0.):
        super().__init__()
        from matplotlib.patches import Rectangle
        from matplotlib.lines import Line2D
        # from matplotlib.artist import Artist
        # To avoid crashing with maximum recursion depth exceeded
        import sys
        sys.setrecursionlimit(10000) # 10000 is 10x the default value

        if height is None:
            self.type = 'Square'
            height = width
        else:
            self.type = 'Rectangle'
        self.ax = ax
        self.angle  = angle/180.*np.pi
        self.width  = width
        self.height = height
        self.rect = Rectangle(corner,width,height,edgecolor='Lime',facecolor='none',angle=angle,fill=False,animated=True)
        self.ax.add_patch(self.rect)
        self.canvas = self.rect.figure.canvas
        x,y = self.compute_markers()
        self.line = Line2D(x, y, marker='o', linestyle=None, linewidth=0., markerfacecolor='g', animated=True)
        self.ax.add_line(self.line)
        self.cid = self.rect.add_callback(self.rectangle_changed)
        self._ind = None  # the active point
        self.connect()
        self.aperture = self.rect
        self.press = None
        self.lock = None

    def compute_markers(self):

        theta0 = self.rect.angle / 180.*np.pi
        w0 = self.rect.get_width()
        h0 = self.rect.get_height()
        x0,y0 = self.rect.get_xy()
        c, s = np.cos(-theta0), np.sin(-theta0)
        R = np.matrix('{} {}; {} {}'.format(c, s, -s, c))

        x = [0.5*w0, w0, 0.5*w0]
        y = [0.5*h0, 0.5*h0, h0]

        self.xy = []
        x_ = []
        y_ = []
        for dx,dy in zip(x,y):
            (dx_,dy_), = np.array(np.dot(R,np.array([dx,dy])))
            self.xy.append((dx_+x0,dy_+y0))
            x_.append(dx_+x0)
            y_.append(dy_+y0)

        return x_,y_

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
        self.rect.remove()
        self.line.remove()
        self.canvas.draw_idle()
        self.aperture = None
        
    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.rect)
        self.ax.draw_artist(self.line)

    def rectangle_changed(self, rect):
        'this method is called whenever the polygon object is called'
        # only copy the artist props to the line (except visibility)
        vis = self.line.get_visible()
        Artist.update_from(self.line, rect)
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
        x0, y0 = self.rect.get_xy()
        w0, h0 = self.rect.get_width(), self.rect.get_height()
        theta0 = self.rect.angle/180*np.pi
        self.press = x0, y0, w0, h0, theta0, event.xdata, event.ydata
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
            #self.rect = None
            #self.line = None
            self.mySignal.emit('rectangle deleted')
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
        if not self.showverts:
            return
        if self._ind is None:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        x0, y0, w0, h0, theta0, xpress, ypress = self.press
        self.dx = event.xdata - xpress
        self.dy = event.ydata - ypress
        self.update_rectangle()

        # Redraw rectangle and points
        self.canvas.restore_region(self.background)
        self.ax.draw_artist(self.rect)
        self.ax.draw_artist(self.line)
        self.canvas.update()
        self.canvas.flush_events()

        # Notify callback
        self.modSignal.emit('rectangle modified')

    def update_rectangle(self):

        x0, y0, w0, h0, theta0, xpress, ypress = self.press
        dx, dy = self.dx, self.dy
        
        if self.lock == "pressed":
            if self._ind == 0:
                self.lock = "move"
            else:
                self.lock = "resizerotate"
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
            self.rect.set_xy((xn,yn))
            # update line
            self.xy = [(i+dx,j+dy) for (i,j) in self.xy0]
            # Redefine line
            self.line.set_data(zip(*self.xy))
        # otherwise rotate and resize
        elif self.lock == 'resizerotate':
            xc,yc = self.xy0[0] # center is conserved in the markers
            dtheta = np.arctan2(ypress+dy-yc,xpress+dx-xc)-np.arctan2(ypress-yc,xpress-xc)
            theta_ = (theta0+dtheta) * 180./np.pi
            c, s = np.cos(theta0), np.sin(theta0)
            R = np.matrix('{} {}; {} {}'.format(c, s, -s, c))
            (dx_,dy_), = np.array(np.dot(R,np.array([dx,dy])))

            # Avoid to pass through the center            
            if self._ind == 1:
                w_ = w0+2*dx_  if (w0+2*dx_) > 0 else w0
                if self.type == 'Square':
                    h_ = w_
                else:
                    h_ = h0
            elif self._ind == 2:
                h_ = h0+2*dy_  if (h0+2*dy_) > 0 else h0
                if self.type == 'Square':
                    w_ = h_
                else:
                    w_ = w0
            # update rectangle
            self.rect.set_width(w_)
            self.rect.set_height(h_)
            self.rect.angle = theta_
            # update markers
            self.updateMarkers()

    def updateMarkers(self):
        # update points
        x,y = self.compute_markers()
        self.line.set_data(x,y)
        

class PolygonInteractor(QObject):
    """
    An polygon editor.

    Key-bindings

      't' toggle vertex markers on and off.  When vertex markers are on,
          you can move them, delete them

      'd' delete the vertex under point

      'i' insert a vertex at point.  You must be within epsilon of the
          line connecting two existing vertices

    """

    showverts = True
    epsilon = 5  # max pixel distance to count as a vertex hit
    mySignal = pyqtSignal(str)
    modSignal = pyqtSignal(str)

    def __init__(self, ax, verts):
        super().__init__()
        from matplotlib.patches import Polygon
        from matplotlib.lines import Line2D
        # from matplotlib.artist import Artist
        self.ax = ax
        self.type = 'Polygon'
        self.poly = Polygon(list(verts), animated=True, fill=False, closed=True, color='lime')
        self.ax.add_patch(self.poly)
        self.canvas = self.poly.figure.canvas

        x, y = zip(*self.poly.xy)
        self.line = Line2D(x, y, marker='o', linestyle=None, linewidth=0., markerfacecolor='g', animated=True)
        self.ax.add_line(self.line)

        self.cid = self.poly.add_callback(self.poly_changed)
        self._ind = None  # the active vert
        self.connect()
        self.aperture = self.poly

    def connect(self):
        self.cid_draw = self.canvas.mpl_connect('draw_event', self.draw_callback)
        self.cid_press = self.canvas.mpl_connect('button_press_event', self.button_press_callback)
        self.cid_key = self.canvas.mpl_connect('key_press_event', self.key_press_callback)
        self.cid_release = self.canvas.mpl_connect('button_release_event', self.button_release_callback)
        self.cid_motion = self.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        self.canvas.draw_idle()

    def disconnect(self):
        self.canvas.mpl_disconnect(self.cid_draw)
        self.canvas.mpl_disconnect(self.cid_press)
        self.canvas.mpl_disconnect(self.cid_key)
        self.canvas.mpl_disconnect(self.cid_release)
        self.canvas.mpl_disconnect(self.cid_motion)
        self.poly.remove()
        self.line.remove()
        self.canvas.draw_idle()
        self.aperture = None
        
    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.poly)
        self.ax.draw_artist(self.line)
        #self.canvas.blit(self.ax.bbox)
        #self.canvas.udpate()
        #self.canvas.flush_events()

    def poly_changed(self, poly):
        'this method is called whenever the polygon object is called'
        # only copy the artist props to the line (except visibility)
        vis = self.line.get_visible()
        Artist.update_from(self.line, poly)
        self.line.set_visible(vis)  # don't use the poly visibility state

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
                if len(self.poly.xy) < 5:  # the minimum polygon has 4 points since the 1st is repeated as final
                    # Delete polygon
                    #self.disconnect()
                    #self.poly = None
                    #self.line = None
                    self.mySignal.emit('polygon deleted')
                else:
                    self.poly.xy = [tup
                                    for i, tup in enumerate(self.poly.xy)
                                    if i != ind]
                    self.line.set_data(zip(*self.poly.xy))
                    self.mySignal.emit('one vertex of polygon removed')

                
        elif event.key == 'i':
            xys = self.poly.get_transform().transform(self.poly.xy)
            p = event.x, event.y  # display coords
            for i in range(len(xys) - 1):
                s0 = xys[i]
                s1 = xys[i + 1]
                d = dist_point_to_segment(p, s0, s1)
                if d <= self.epsilon:
                    self.poly.xy = np.array(
                        list(self.poly.xy[:i+1]) +
                        [(event.xdata, event.ydata)] +
                        list(self.poly.xy[i+1:]))
                    self.line.set_data(zip(*self.poly.xy))
                    break

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

        # Notify callback
        self.modSignal.emit('polygon modified')

    def updateMarkers(self):
        self.line.set_data(zip(*self.poly.xy))
