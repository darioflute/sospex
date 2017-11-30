import numpy as np

class photoAperture(object):
    """ 
    Photometric aperture 
    
    Data are conserved in RA, Dec units to be passed to
    images with different astrometry
    
    """

    def __init__(self, n, type, data):

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


class DragResizeRotateEllipse:
    #from PyQt5.QtCore import pyqtSignal    
    #updateEllipses = pyqtSignal('QString')

    lock = None
    def __init__(self, ellipse, border_tol=0.6, changed=False):
        self.ellipse = ellipse
        self.ellipse.set_animated(True)
        self.border_tol = border_tol
        self.lock = None
        self.press = None
        self.background = None
        self.changed = changed
        self.axes = ellipse.axes
        canvas = self.axes.figure.canvas        
        # connect to events
        self.ciddraw = canvas.mpl_connect('draw_event', self.on_draw)
        self.cidpress = canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease = canvas.mpl_connect('button_release_event', self.on_release)
        self.cidmotion = canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.canvas = canvas
        self.spot = 0
        
    def on_press(self, event):
        '''on button press it stores some data if mouse is over it'''
        if event.inaxes != self.ellipse.axes: return
        if DragResizeRotateEllipse.lock is not None: return
        contains, attrd = self.ellipse.contains(event)
        if not contains: return

        x0, y0 = self.ellipse.center
        w0, h0 = self.ellipse.width, self.ellipse.height
        theta0 = self.ellipse.angle
        self.lock = "pressed"
        self.press = x0, y0, w0, h0, theta0, event.xdata, event.ydata
        DragResizeRotateEllipse.lock = self
        # draw everything
        self.canvas.draw_idle()
        # register changes
        self.changed = True

    def on_draw(self, event):
        ''' Drawing the ellipse '''
        # first capture the background
        self.background = self.canvas.copy_from_bbox(self.axes.bbox)
        # then draw the ellipse
        self.axes.draw_artist(self.ellipse)

    def on_motion(self, event):
        '''on motion it will act on the ellipse if the mouse is over it'''
        if DragResizeRotateEllipse.lock is not self: return
        if not self.ellipse.contains(event): return

        # store original position, compute coordinates increment, update ellipse
        x0, y0, w0, h0, theta0, xpress, ypress = self.press
        self.dx = event.xdata - xpress
        self.dy = event.ydata - ypress
        self.update_ellipse()

        # restore the background region
        # redraw just the current ellipse
        self.canvas.restore_region(self.background)
        self.axes.draw_artist(self.ellipse)
        # blit just the redrawn area
        #self.canvas.blit(self.axes.bbox)
        self.canvas.update()
        self.canvas.flush_events()

        
    def on_release(self, event):
        '''on release it resets the press data'''
        if DragResizeRotateEllipse.lock is not self:
            return

        self.press = None
        DragResizeRotateEllipse.lock = None
        self.lock = "released"
        self.spot = 0

        # turn off the animation property and reset the background
        self.background = None
        
    def disconnect(self):
        'disconnect all the stored connection ids'
        self.canvas.mpl_disconnect(self.ciddraw)
        self.canvas.mpl_disconnect(self.cidpress)
        self.canvas.mpl_disconnect(self.cidrelease)
        self.canvas.mpl_disconnect(self.cidmotion)

    def update_ellipse(self):
        x0, y0, w0, h0, theta0, xpress, ypress = self.press
        dx, dy = self.dx, self.dy
        bt = self.border_tol
        # Normalized point (to the circle)
        xnorm, ynorm = self.ellipse.get_patch_transform().inverted().transform_point((xpress, ypress))

        # lock into a mode
        if self.lock == "pressed":
            rnorm = np.sqrt(xnorm*xnorm+ynorm*ynorm)
            if rnorm > bt:
                anorm = np.arctan2(ynorm,xnorm)*180./np.pi
                self.lock = "resizerotate"
                th0 = theta0/180.*np.pi
                c, s = np.cos(th0), np.sin(th0)
                R = np.matrix('{} {}; {} {}'.format(c, s, -s, c))
                (dx_,dy_), = np.array(np.dot(R,np.array([dx,dy])))
                if abs(dx_) > 1.2*abs(dy_) and (abs(anorm) < 30.):
                    self.spot = 1
                elif abs(dx_) > 1.2*abs(dy_) and (abs(anorm) > 150.):
                    self.spot = 3
                elif abs(dy_) > 1.2*abs(dx_) and (anorm > 60. and anorm < 120.):
                    self.spot = 2
                elif abs(dy_) > 1.2*abs(dx_) and (anorm < -60. and anorm > -120.):
                    self.spot = 4
                else:
                    self.spot = 0
            else:
                self.lock = "move"
                
        elif self.lock == "move":
            xn = x0+dx; yn =y0+dy
            if xn < 0: xn = x0
            if yn < 0: yn = y0
            self.ellipse.center = (xn,yn)
        elif self.lock == "resizerotate":
            dtheta = np.arctan2(ypress+dy-y0,xpress+dx-x0)-np.arctan2(ypress-y0,xpress-x0)
            dtheta *= 180./np.pi
            theta_ = theta0+dtheta

            anorm = np.arctan2(ynorm,xnorm)*180./np.pi
            th0 = theta0/180.*np.pi
            c, s = np.cos(th0), np.sin(th0)
            R = np.matrix('{} {}; {} {}'.format(c, s, -s, c))
            (dx_,dy_), = np.array(np.dot(R,np.array([dx,dy])))

            if self.spot == 1:
                #print('arc1')
                w_ = w0+2*dx_  if (w0+2*dx_) > 0 else w0 # Avoid flipping
                h_ = h0
            elif self.spot ==3:
                #print('arc3')
                w_ = w0-2*dx_  if (w0-2*dx_) > 0 else w0
                h_ = h0
            elif self.spot == 2:
                #print('arc2')
                h_ = h0+2*dy_  if (h0+2*dx_) > 0 else h0
                w_ = w0
            elif self.spot == 4:
                #print('arc4')
                h_ = h0-2*dy_  if (h0-2*dx_) > 0 else h0
                w_ = w0
            else:
                self.lock = "released"

            if self.lock != "released":
                self.ellipse.width = w_
                self.ellipse.height = h_
                self.ellipse.angle = theta_

from matplotlib.patches import Rectangle
class RectangleSelect:
    ''' Add rectangle to the figure to crop the spectral cube'''
    def __init__(self, parent):
        self.frame = parent
        self.rect = Rectangle((0,0), 0.1, 0.1, facecolor='None', edgecolor='None')
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.rectpatch = self.frame.axes.add_patch(self.rect)
        self.pressed = False
        self.connect()
        
    def connect(self):
        ''' connect to events '''
        self.cidpress = self.rect.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease = self.rect.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cidmotion = self.rect.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def disconnect(self):
        '''disconnect all the stored connection ids'''
        self.rect.figure.canvas.mpl_disconnect(self.cidpress)
        self.rect.figure.canvas.mpl_disconnect(self.cidrelease)
        self.rect.figure.canvas.mpl_disconnect(self.cidmotion)
        self.rectpatch.remove()
        self.frame.canvas.draw_idle()
        
    def on_press(self, event):
        ''' Callback to handle the mouse being clicked and held over the canvas'''
        # Check the mouse press was actually on the canvas 
        if event.xdata is not None and event.ydata is not None:
            # Upon initial press of the mouse record the origin and record the mouse as pressed
            self.pressed = True
            self.rect.set_linestyle('dashed')
            self.rect.set_edgecolor('blue')
            self.x0 = event.xdata
            self.y0 = event.ydata
            print ("pressed ",self.x0,self.y0)
            
    def on_motion(self, event):
        '''Callback to handle the motion event created by the mouse moving over the canvas'''

        # If the mouse has been pressed draw an updated rectangle when the mouse is moved so 
        # the user can see what the current selection is
        if self.pressed:
            # Check the mouse was released on the canvas, and if it wasn't then just leave the width and 
            # height as the last values set by the motion event
            if event.xdata is not None and event.ydata is not None:
                self.x1 = event.xdata
                self.y1 = event.ydata
                
            # Set the width and height and draw the rectangle
            self.rect.set_width(self.x1 - self.x0)
            self.rect.set_height(self.y1 - self.y0)
            self.rect.set_xy((self.x0, self.y0))
            self.frame.canvas.draw()

    def get_xlim(self):
        if self.x1 > self.x0:
            return (self.x0,self.x1)
        else:
            return (self.x1,self.x0)
        
    def get_ylim(self):
        if self.y1 > self.y0:
            return (self.y0,self.y1)
        else:
            return (self.y1,self.y0)

    def get_rect(self):
        return self.rect
 
    def on_release(self, event):
        '''Callback to handle the mouse being released over the canvas'''
        
        # Check that the mouse was actually pressed on the canvas to begin with and this isn't a rouge mouse 
        # release event that started somewhere else
        if self.pressed:
            # Upon release draw the rectangle as a solid rectangle
            self.pressed = False
            self.rect.set_linestyle('solid')

            # Check the mouse was released on the canvas, and if it wasn't then just leave the width and 
            # height as the last values set by the motion event
            if event.xdata is not None and event.ydata is not None:
                self.x1 = event.xdata
                self.y1 = event.ydata

            # Set the width and height and origin of the bounding rectangle
            self.boundingRectWidth =  self.x1 - self.x0
            self.boundingRectHeight =  self.y1 - self.y0
            self.bouningRectOrigin = (self.x0, self.y0)

            # Draw the bounding rectangle
            self.rect.set_width(self.boundingRectWidth)
            self.rect.set_height(self.boundingRectHeight)
            self.rect.set_xy((self.x0, self.y0))
            self.frame.canvas.draw()



class EllipseInteractor(object):

    epsilon = 5
    showverts = True
    
    def __init__(self,ax,center,width,height=None,angle=0.):

        from matplotlib.patches import Ellipse
        from matplotlib.lines import Line2D
        from matplotlib.artist import Artist
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
            self.disconnect()
            self.ellipse = None
            self.line = None

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
            #ca = np.cos(theta_*np.pi/180.); sa = np.sin(theta_*np.pi/180.)
            #x = [x0, x0+w_*0.5*ca, x0-h_*0.5*sa]
            #y = [y0, y0+w_*0.5*sa, y0+h_*0.5*ca]
            #self.xy = [(i,j) for (i,j) in zip(x,y)]
                
        # Redefine line
        #self.line.set_data(zip(*self.xy))


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




class RectangleInteractor(object):

    epsilon = 5
    showverts = True
    
    def __init__(self,ax,corner,width,height=None,angle=0.):

        from matplotlib.patches import Rectangle
        from matplotlib.lines import Line2D
        from matplotlib.artist import Artist
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
            self.disconnect()
            self.rect = None
            self.line = None

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



        
from matplotlib.mlab import dist_point_to_segment
        

class PolygonInteractor(object):
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

    def __init__(self, ax, verts):

        from matplotlib.patches import Polygon
        from matplotlib.lines import Line2D
        from matplotlib.artist import Artist


        self.ax = ax
        self.type = 'Polygon'
        self.poly = Polygon(list(verts), animated=True, fill=False, closed=True, color='g')
        self.ax.add_patch(self.poly)
        self.canvas = self.poly.figure.canvas
        #self.poly = poly

        x, y = zip(*self.poly.xy)
        #self.line = Line2D(x, y, marker='o', markerfacecolor='g', color='g', animated=True)
        self.line = Line2D(x, y, marker='o', linestyle=None, linewidth=0., markerfacecolor='g', animated=True)
        self.ax.add_line(self.line)

        self.cid = self.poly.add_callback(self.poly_changed)
        self._ind = None  # the active vert
        #self.canvas = canvas
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
                    self.disconnect()
                    self.poly = None
                    self.line = None
                else:
                    self.poly.xy = [tup
                                    for i, tup in enumerate(self.poly.xy)
                                    if i != ind]
                    self.line.set_data(zip(*self.poly.xy))

                
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
        self.update_markers()
        #self.line.set_data(zip(*self.poly.xy))

        self.canvas.restore_region(self.background)
        self.ax.draw_artist(self.poly)
        self.ax.draw_artist(self.line)
        self.canvas.update()
        self.canvas.flush_events()

    def update_markers(self):
        self.line.set_data(zip(*self.poly.xy))
        
