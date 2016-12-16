#from matplotlib.patches import Ellipse
#import matplotlib.pyplot as plt
#import matplotlib.transforms as transforms
import numpy as np


class DragResizeRotateEllipse:
    lock = None
    def __init__(self, ellipse, border_tol=0.7, allow_resize=True,
                 fixed_aspect_ratio=True, allow_rotate=True):
        self.ellipse = ellipse
        self.border_tol = border_tol
        self.allow_resize = allow_resize
        self.allow_rotate = allow_rotate
        self.fixed_aspect_ratio = fixed_aspect_ratio
        self.press = None
        self.background = None
    def connect(self):
        'connect to events'
        self.cidpress = self.ellipse.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.ellipse.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.ellipse.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)
    
    def on_press(self, event):
        'on button press it stores some data if mouse is over it'
        if event.inaxes != self.ellipse.axes: return
        if DragResizeRotateEllipse.lock is not None: return
        contains, attrd = self.ellipse.contains(event)
        if not contains: return
        x0, y0 = self.ellipse.center
        w0, h0 = self.ellipse.width, self.ellipse.height
        theta0 = self.ellipse.angle
        self.press = x0, y0, w0, h0, theta0, event.xdata, event.ydata
        DragResizeRotateEllipse.lock = self
        
        # draw everything but the selected ellipse and store the pixel buffer
        canvas = self.ellipse.figure.canvas
        axes = self.ellipse.axes
        self.ellipse.set_animated(True)
        canvas.draw()
        self.background = canvas.copy_from_bbox(self.ellipse.axes.bbox)

        # now redraw just the ellipse
        axes.draw_artist(self.ellipse)

        # and blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_motion(self, event):
        'on motion it will act on the ellipse if the mouse is over it'
        if DragResizeRotateEllipse.lock is not self: return
        if not self.ellipse.contains(event): return
        
        x0, y0, w0, h0, theta0, xpress, ypress = self.press

        self.dx = event.xdata - xpress
        self.dy = event.ydata - ypress
        self.update_ellipse()

        canvas = self.ellipse.figure.canvas
        axes = self.ellipse.axes
        # restore the background region
        canvas.restore_region(self.background)

        # redraw just the current ellipse
        axes.draw_artist(self.ellipse)

        # blit just the redrawn area
        canvas.blit(axes.bbox)
        
    def on_release(self, event):
        'on release it resets the press data'
        if DragResizeRotateEllipse.lock is not self:
            return

        self.press = None
        DragResizeRotateEllipse.lock = None

        # turn off the animation property and reset the background
        self.ellipse.set_animated(False)
        self.background = None

        # redraw the full figure
        self.ellipse.figure.canvas.draw()

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.ellipse.figure.canvas.mpl_disconnect(self.cidpress)
        self.ellipse.figure.canvas.mpl_disconnect(self.cidrelease)
        self.ellipse.figure.canvas.mpl_disconnect(self.cidmotion)

    def update_ellipse(self):
        x0, y0, w0, h0, theta0, xpress, ypress = self.press
        dx, dy = self.dx, self.dy
        bt = self.border_tol
        # Normalized point (to the circle)
        xnorm, ynorm = self.ellipse.get_patch_transform().inverted().transform_point((xpress+dx, ypress+dy))
        rnorm = np.sqrt(xnorm*xnorm+ynorm*ynorm)
        if (rnorm < bt):
            self.ellipse.center = (x0+dx,y0+dy)
        else:
            anorm = np.arctan2(ynorm,xnorm)*180./np.pi
            dtheta = np.arctan2(ypress+dy-y0,xpress+dx-x0)-np.arctan2(ypress-y0,xpress-x0)
            dtheta *= 180./np.pi
            th0 = theta0/180.*np.pi
            c, s = np.cos(th0), np.sin(th0)
            R = np.matrix('{} {}; {} {}'.format(c, s, -s, c))
            (dx_,dy_), = np.array(np.dot(R,np.array([dx,dy])))
            if abs(dx_) > abs(dy_) and (abs(anorm) < 30.):
                self.ellipse.width = w0+dx_*2
            elif abs(dx_) > abs(dy_) and (abs(anorm) > 150.):
                self.ellipse.width = w0-dx_*2
            elif abs(dy_) > abs(dx_) and (anorm > 60. and anorm < 120.):
                self.ellipse.height = h0+dy_*2
            elif abs(dy_) > abs(dx_) and (anorm < -60. and anorm > -120.):
                self.ellipse.height = h0-dy_*2
            else:
                self.ellipse.angle = theta0+dtheta
