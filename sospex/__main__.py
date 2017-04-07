#!/usr/bin/env python
# Using UTF8 encoding
# -*- coding: utf-8 -*-


# In the case of MacOS-X, re-executes with pythonw
import os,sys
if (len(sys.argv) == 1) or (sys.argv[1] != "MAC"):
    if sys.platform == 'darwin' and\
       not sys.executable.endswith('MacOS/Python'):
        print >>sys.stderr,'re-executing using pythonw'
        os.execvp('pythonw',['pythonw',__file__,"MAC"] + sys.argv[1:])

    
# Library imports
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigureCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar
from matplotlib.widgets import Slider
from matplotlib.figure import Figure
from matplotlib.widgets import SpanSelector
from matplotlib.patches import Ellipse
from matplotlib.font_manager import FontProperties

from astropy.io import fits
from astropy.nddata import Cutout2D

import wx, fnmatch, time
import numpy as np
#import wx.lib.platebtn as pbtn
#from wx.lib.buttons import GenBitmapButton
from scipy.spatial import ConvexHull
from wx.lib.dialogs import ScrolledMessageDialog as ScrolledMessage
from wx.lib.agw import pybusyinfo as PBI


# Local imports
from photometry import DragResizeRotateEllipse
from ellipsefit import fit_ellipse
from lines import define_lines
import icons
from fitline import fitContinuum, fitLines, fittedLine
from specobj import specCube, spectrum

        
class MyFrame(wx.Frame):
    def __init__(self, parent, id, title):
        displaySize= wx.DisplaySize()
        self.xysize = (displaySize[0]*.8, displaySize[1]*0.5+80)
        wx.Frame.__init__(self, None, wx.ID_ANY, " SOSPEX  (SOfia SPectrum EXplorer)",size=self.xysize)

        self.topPanel = wx.Panel(self)
        self.Position = (displaySize[0]*0.1,displaySize[1]*0.3)
        self.BackgroundColor = 'white'

        self.panel1 = Panel1(self.topPanel, self.xysize)
        self.toolbar1 = Toolbar1(self.topPanel, self.xysize)
        self.panel2 = Panel2(self.topPanel, self.xysize)
        self.toolbar2 = Toolbar2(self.topPanel, self.xysize)

        # define path as current directory ...
        self.path = os.getcwd()        


        # Menu
        # http://zetcode.com/wxpython/menustoolbars/
        self.menubar = wx.MenuBar()
        fileMenu = wx.Menu()
        fitem = fileMenu.Append(wx.ID_EXIT, 'Quit', 'Quit application')
        self.menubar.Append(fileMenu, '&File')
        self.SetMenuBar(self.menubar)
        self.Bind(wx.EVT_MENU, self.onQuit, fitem)
        
        
        # Select a region in the image with the mouse (and then plot all the spectra  or the average coadded spectrum)
        #http://ebanshi.cc/questions/5266620/matplotlib-draw-a-selection-area-in-the-shape-of-a-rectangle-with-the-mouse
        
        self.sizeFrame()
        
        # States
        self.extSpectrum = None
        self.fitState = 0
        self.displayContours = False
        self.fit  = None 
        self.ufit = None 
        self.shade = False
        self.regionlimits = None
        self.panel1.figure.canvas.mpl_connect('button_release_event', self.onClick)
        self.panel2.figure.canvas.mpl_connect('button_release_event', self.onClick)

        #self.getFile()
        
    def welcomeMessage(self):

        self.panel1.figure.set_canvas(self.panel1.canvas)
        self.panel1.axes = self.panel1.figure.add_subplot(111)
        self.panel1.figure.suptitle('SOFIA Spectrum Explorer', fontsize=20, fontweight='bold')
        self.panel1.axes.text(0.03,0.8, u'To start:', fontweight='bold',style='italic',fontsize=15)
        self.panel1.axes.text(0.03,0.7, u'  Select a file with the double arrow',style='italic',fontsize=15)
        self.panel1.axes.text(0.03,0.6, u'Icons:', fontweight='bold',style='italic',fontsize=15)
        self.panel1.axes.text(0.03,0.5, u'  Hover mouse to know what they do',style='italic',fontsize=15)
        self.panel1.axes.text(0.03,0.4, u'Mouse:', fontweight='bold',style='italic',fontsize=15)
        self.panel1.axes.text(0.03,0.3, u'  Right click to select images and plots',style='italic',fontsize=15)
        self.panel1.axes.text(0.03,0.2, u'Issues/Ideas:', fontweight='bold',style='italic',fontsize=15)
        self.panel1.axes.text(0.03,0.1, u'  Ask Dario Fadda (darioflute@gmail.com)',style='italic',fontsize=15)
        #self.panel1.axes.axis('off')
        self.panel1.axes.xaxis.label.set_visible(False)
        self.panel1.axes.yaxis.label.set_visible(False)
        self.panel1.axes.set_yticklabels([])
        self.panel1.axes.set_xticklabels([])
        self.panel1.axes.xaxis.set_ticks_position('none') 
        self.panel1.axes.yaxis.set_ticks_position('none') 
        
    def getFile(self):
        if wx.Platform != '__WXMAC__':
            defaultFile = "*WXY*.fits"
            wildcard = "Spectrum cube (*WXY*.fits)|*WXY*.fits"
        else:
            defaultFile = "*.fits"
            wildcard = "Spectrum cube (*WXY*.fits)|*.fits"
            
        dlg = wx.FileDialog(
            self, message="Choose a spectrum",
            defaultDir=self.path, 
            defaultFile=defaultFile,
            wildcard=wildcard,
            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST
            )
        if dlg.ShowModal() == wx.ID_OK:
            infile = dlg.GetPath()
            print "You chose the following file:", infile
            # redefine path
            self.path = os.path.dirname(infile)
            try:
                self.spectrum = specCube(infile)
                self.limits = (0, self.spectrum.n)
                self.ijmax = self.panel1.view(self.spectrum, self.limits, 1)
                self.extSpectrum = None
                
                # States
                self.fitState = 0
                self.displayContours = False
                self.fit  = None 
                self.ufit = None 
            
                # Add marker to image and plot spectrum (of flux inside ellipse)
                self.ellipse = Ellipse((self.ijmax[1],self.ijmax[0]), 5, 5, edgecolor='Lime', facecolor='none')
                self.shade = False
                self.drawEllipse()
            
                # Select range in wavelength to average the cube and display in panel 1
                self.regionlimits = None  # initialize regionlimits
                self.span = SpanSelector(self.panel2.ax4, self.onSelect, 'horizontal', useblit=True,
                                         rectprops=dict(alpha=0.5, facecolor='LightSalmon'),button=1)
            
                # Clicks on figure (panel2)
                self.panel1.figure.canvas.mpl_connect('button_release_event', self.onClick)
                self.panel2.figure.canvas.mpl_connect('button_release_event', self.onClick)
            except:
                print "Invalid file"
                
        dlg.Destroy()

        
    def showmsg(self, message, title):
        d = PBI.PyBusyInfo(message, title=title)
        time.sleep(2)
        d = None
        
    def getName(self,path):
        """ list of files matching WXY """    
        wxyfile = fnmatch.filter(os.listdir(path), "*WXY*.fits")
        return path+wxyfile[0]

        
    def onQuit(self, event):
        self.Close()
        
    def sizeFrame(self):
        self.sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.vsizer1 = wx.BoxSizer(wx.VERTICAL)
        self.vsizer1.Add(self.panel1,0,wx.EXPAND|wx.ALL,border=0)
        self.vsizer1.Add(self.toolbar1,0,wx.EXPAND|wx.ALL,border=0)
        self.sizer.Add(self.vsizer1,0,wx.EXPAND|wx.ALL,border=0)
        self.vsizer2 = wx.BoxSizer(wx.VERTICAL)
        self.vsizer2.Add(self.panel2,0,wx.EXPAND|wx.ALL,border=0)
        self.vsizer2.Add(self.toolbar2,0,wx.EXPAND|wx.ALL,border=0)
        self.sizer.Add(self.vsizer2,0,wx.EXPAND|wx.ALL,border=0)
        self.topPanel.SetSizer(self.sizer)
        self.topPanel.Fit()
        self.topPanel.Show()


    def onSelect(self, xmin, xmax):
        indmin, indmax = np.searchsorted(self.spectrum.wave, (xmin, xmax))
        indmax = min(len(self.spectrum.wave) - 1, indmax)
        if xmin < xmax:
            if self.fitState == False:
                print("Data between ", xmin, " and ",xmax)
                print("indmin - indmax ",indmin, indmax)
                self.limits = (indmin, indmax)
                self.shade = True
                try:
                    self.panel2.removeSlider()
                except:
                    pass
                self.panel1.displayMethod = 'Average'
                self.toolbar1.dispBtn.SetBitmap(icons.lambda2.GetBitmap())
                self.panel1.view(self.spectrum,self.limits,0)            
                self.regionlimits = np.array([xmin,xmax])
                self.drawEllipse()
                self.panel1.vbox.Layout()
            elif self.fitState == 1:                
                self.icont1 = (indmin, indmax)
                self.cont1 = np.array([xmin,xmax])
                self.shadeSpectrum()
                self.fitState = 2
            elif self.fitState == 2:
                self.icont2 = (indmin, indmax)
                self.cont2 = np.array([xmin,xmax])
                self.shadeSpectrum()
                x0,y0 = self.ellipse.center
                pixel = np.array([[x0, y0]], np.float_)
                world = self.spectrum.wcs.wcs_pix2world(pixel, 1)
                print "center of ellipse is: ", world
                if self.panel2.displayFlux:
                    self.fit = fittedLine(self.icont1,self.icont2,self.spectrum.wave, \
                                          self.panel2.flux,self.panel2.eflux,\
                                          self.ellipse, world)
                    self.fit = fitContinuum(self.fit)
                    self.panel2.displayLineFit = True
                else:
                    if self.panel2.displayUFlux:
                        self.ufit = fittedLine(self.icont1,self.icont2,self.spectrum.wave, \
                                               self.panel2.uflux,self.panel2.euflux,\
                                               # note: use eflux instead of euflux, which contains also
                                               # uncertainty on response
                                               self.ellipse, world)
                        self.ufit = fitContinuum(self.ufit)
                        self.panel2.displayUncLineFit = True
                    
                # Overplot fit
                if self.fit != None:
                    self.contfitLayer, = self.panel2.ax1.plot(self.fit.wave,self.fit.continuum,color='cyan')
                    self.contfitLayer.set_visible(self.panel2.displayLineFit)
                if self.ufit != None:
                    self.ucontfitLayer, = self.panel2.ax4.plot(self.ufit.wave,self.ufit.continuum,color='palegreen')
                    self.ucontfitLayer.set_visible(self.panel2.displayUncLineFit)
                print "Fit the continuum ..."
                self.fitState = 3
                # Help for fit
                self.showmsg("Click on top of line, then click on icon to fit","Help")

            
    def drawEllipse(self):
        # Add marker to image
        try:
            self.marker.remove()
        except:
            pass
        self.marker = self.panel1.axes.add_patch(self.ellipse)
        self.drr = DragResizeRotateEllipse(self.ellipse)
        self.drr.connect()
        # Plot spectrum
        self.panel2.draw(self.spectrum,self.ellipse)
        self.shadeSpectrum()
        
    def shadeSpectrum(self):
        try:
            self.panel2.region.remove()
        except:
            pass                
        if self.shade == True:
            t = self.regionlimits
            ylim=self.panel2.ax1.get_ylim()
            lowlim = np.zeros(2)+ylim[0]
            upplim = np.zeros(2)+ylim[1]
            self.panel2.region = self.panel2.ax1.fill_between(t,lowlim,upplim,facecolor='Lavender',alpha=0.5,linewidth=0)
        if self.fitState == 1:
            t1 = self.cont1
            ylim=self.panel2.ax1.get_ylim()
            lowlim = np.zeros(2)+ylim[0]
            upplim = np.zeros(2)+ylim[1]
            self.panel2.regionc1 = self.panel2.ax1.fill_between(t1,lowlim,upplim,facecolor='Plum',alpha=0.3,linewidth=0)
        if self.fitState == 2:
            t2 = self.cont2
            ylim=self.panel2.ax1.get_ylim()
            lowlim = np.zeros(2)+ylim[0]
            upplim = np.zeros(2)+ylim[1]
            self.panel2.regionc2 = self.panel2.ax1.fill_between(t2,lowlim,upplim,facecolor='Plum',alpha=0.3,linewidth=0)
        if self.fitState == False:
            try:
                self.panel2.regionc1.remove()
            except:
                pass
            try:
                self.panel2.regionc2.remove()
            except:
                pass
                
            
    def onClick(self, event):
        if event.inaxes == self.panel1.axes:
            print (event.xdata,event.ydata)
            if  event.button == 1:
                self.panel2.draw(self.spectrum,self.ellipse)
                if self.regionlimits is not None: self.shadeSpectrum()
                self.panel2.Layout()
                # if button is 3 --> popupMenu
            elif event.button == 2:
                # Recenter image
                self.panel1.recenterImage(event)
            elif event.button == 3:
                self.panel1.OnRightDown(event)
            else:
                return
                # Error in callback for matplotlib (get the last axes for event.xdata, event.ydata)
                # https://stackoverflow.com/questions/16672530/cursor-tracking-using-matplotlib-and-twinx/16672970#16672970
        elif event.inaxes == self.panel2.ax4:
            if event.button == 1:
                if self.fitState == 3:
                    # Start collecting lines
                    self.xpeak.append(event.xdata) # *(1.+self.spectrum.z))
                    self.ypeak.append(event.ydata)
                    if self.panel2.displayFlux:
                        color='blue'
                    else:
                        color='green'
                    self.panel2.ax1.plot([event.xdata],[event.ydata],'o',color=color)
                else:
                    print "wavelength: ","{:.3f}".format(event.xdata),u"\u03BCm (no baryshift:  ".encode('utf8'), "{:.3f}".format(event.xdata/(1.+self.spectrum.baryshift)),u"\u03BCm )".encode('utf8')
            if event.button == 3:
                self.panel2.OnRightDown(event)
            else:
                return
        else:
            return


    def nextCube(self, event):

        # Choose file
        if wx.Platform != '__WXMAC__':
            defaultFile = "*WXY*.fits"
            wildcard = "Spectrum cube (*WXY*.fits)|*WXY*.fits"
        else:
            defaultFile = "*.fits"
            wildcard = "Spectrum cube (*WXY*.fits)|*.fits"
            
        dlg = wx.FileDialog(
            self, message="Choose a spectrum",
            defaultDir=self.path, 
            defaultFile=defaultFile,
            wildcard=wildcard,
            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST
            )
        if dlg.ShowModal() == wx.ID_OK:
            infile = dlg.GetPath()
            print "You chose the following file:", infile
            self.path = os.path.dirname(infile)  # redefine path name
            self.spectrum = specCube(infile)
            self.limits = (0, self.spectrum.n)
            # remove image and ellipse
            try:
                self.panel1.figure.delaxes(self.panel1.axes)
                self.ellipse.remove()
            except:
                pass
            # draw new image
            self.ijmax = self.panel1.view(self.spectrum, self.limits, 1)
        
            # States
            self.fitState = 0
            self.displayContours = False
            self.fit = None # initialize fit
            self.extSpectrum = None

            # Add marker to image and plot spectrum (of flux inside ellipse)
            self.ellipse = Ellipse((self.ijmax[1],self.ijmax[0]), 5, 5, edgecolor='Lime', facecolor='none')
            self.shade = False
            self.drawEllipse()
            self.panel1.Layout()

            # Select range in wavelength to average the cube and display in panel 1
            self.regionlimits = None  # initialize regionlimits
            self.span = SpanSelector(self.panel2.ax4, self.onSelect, 'horizontal', useblit=True,
                                     rectprops=dict(alpha=0.5, facecolor='LightSalmon'),button=1)
        dlg.Destroy()

        # Forget previous toolbar values
        self.panel1.ntoolbar._views.clear()
        self.panel1.ntoolbar._positions.clear()
        self.panel1.ntoolbar._update_view() 
        self.panel1.ntoolbar.Refresh()
        self.panel2.ntoolbar._views.clear()
        self.panel2.ntoolbar._positions.clear()
        self.panel2.ntoolbar._update_view() 
        self.panel2.ntoolbar.Refresh()

            
    def reloadCube(self, event):
        # Get data        
        filename = self.spectrum.filename
        oldwcs = self.spectrum.wcs
        self.spectrum = specCube(filename)
        self.limits = (0, self.spectrum.n)
        # Compute new ellipse center
        x0,y0 = self.ellipse.center
        pixel = np.array([[x0, y0]], np.float_)
        world = oldwcs.wcs_pix2world(pixel, 1)
        print "world ",world
        pix = self.spectrum.wcs.wcs_world2pix(world,1)
        # Redisplay ...
        # Remove previous axes
        self.panel1.figure.delaxes(self.panel1.axes)
        # Remove and redefine ellipse
        angle = self.ellipse.angle
        width = self.ellipse.width
        height = self.ellipse.height
        self.ellipse.remove()
        self.ellipse = Ellipse((pix[0][0], pix[0][1]), width,height, edgecolor='Lime', facecolor='none',angle=angle)
        # Build new figure
        self.ijmax = self.panel1.view(self.spectrum, self.limits, 1)
        self.panel1.Layout()
        self.drawEllipse()
        # Forget previous toolbar values
        self.panel1.ntoolbar._views.clear()
        self.panel1.ntoolbar._positions.clear()
        self.panel1.ntoolbar._update_view() 
        self.panel1.ntoolbar.Refresh()

    def cutCube(self, event):
        # Dialog window to start cropping
        dlg = wx.MessageDialog(self,"Do you really want to cut it ?", "Question",
                               wx.YES_NO|wx.ICON_QUESTION)
        if dlg.ShowModal() == wx.ID_YES:
            # redefine spectral size (from the span selector)
            limits = self.limits
            l0 = limits[0]
            l1 = limits[1]
            self.spectrum.wave = self.spectrum.wave[l0:l1]
            self.spectrum.crval3 = self.spectrum.wave[0]
            self.spectrum.response = self.spectrum.response[l0:l1]
            self.spectrum.atran = self.spectrum.atran[l0:l1]
            self.spectrum.flux = self.spectrum.flux[l0:l1,:,:]
            self.spectrum.eflux = self.spectrum.eflux[l0:l1,:,:]
            self.spectrum.uflux = self.spectrum.uflux[l0:l1,:,:]
            self.spectrum.euflux = self.spectrum.euflux[l0:l1,:,:]
            self.spectrum.exposure = self.spectrum.exposure[l0:l1,:,:]
            self.spectrum.n = self.spectrum.wave.size
            self.regionlimits = [0,self.spectrum.n]
            self.shade = 'False'
            self.regionlimits = None
            # replot panel 2
            self.panel2.draw(self.spectrum,self.ellipse)
            self.shadeSpectrum()

    def fitGauss(self, event):
        print ("Fitting a line with a Gaussian")
        if self.fitState:
            print ("Do the fit")
            if self.panel2.displayFlux:
                self.fit = fitLines(self.fit, self.xpeak, self.ypeak)
            else:
                if self.panel2.displayUFlux:
                    self.ufit = fitLines(self.ufit, self.xpeak, self.ypeak)                    
            self.fitState = False
            if self.fit != None:
                self.linefitLayer, = self.panel2.ax1.plot(self.fit.wave, self.fit.lines+self.fit.continuum, color='cyan')
                self.linefitLayer.set_visible(self.panel2.displayLineFit)
            if self.ufit != None:
                self.ulinefitLayer, = self.panel2.ax4.plot(self.ufit.wave, self.ufit.lines+self.ufit.continuum, color='palegreen')
                self.ulinefitLayer.set_visible(self.panel2.displayUncLineFit)
            self.panel2.Layout()
            # destroy continuum values
            self.cont1 = None
            self.cont2 = None
            self.xpeak = None
            self.ypeak = None
        else:
            self.fitState = 1
            self.xpeak = []
            self.ypeak = []
            self.showmsg("Select continuum regions dragging the mouse", "Help")
        return

    def saveSpectrum(self, event):
            
        """ 
        Saving the displayed spectrum, ellipse, and fitted lines/continuum
        """
        saveFileDialog = wx.FileDialog(self, "Save As", "", "spectrum.fits", 
                                       "FITS files (*.fits)|*.fits", 
                                       wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        saveFileDialog.ShowModal()
        outfile = saveFileDialog.GetPath()
        saveFileDialog.Destroy()

        # Primary header containing general information
        hdu = fits.PrimaryHDU()
        hdu.header['OBJ_NAME'] = (self.spectrum.objname, 'Object Name')
        hdu.header['FILEGPID'] = (self.spectrum.filegpid, 'File Group ID')
        hdu.header['BARYSHFT'] = (self.spectrum.baryshift, 'Barycentric shift')
        x0,y0 = self.ellipse.center
        pixel = np.array([[x0, y0]], np.float_)
        world = self.spectrum.wcs.wcs_pix2world(pixel, 1)
        hdu.header['RA'] = (world[0][0], 'RA of center of elliptical aperture')
        hdu.header['DEC'] = (world[0][1], 'Dec of center of elliptical aperture')
        hdu.header['ANGLE'] = (self.ellipse.angle, 'Angle of elliptical aperture')
        hdu.header['MAJAX'] = (self.ellipse.width*self.spectrum.pixscale, 'Major axis of elliptical aperture')
        hdu.header['MINAX'] = (self.ellipse.height*self.spectrum.pixscale, 'Minor axis of elliptical aperture')
        hdu.header['RESOLUN'] = (self.spectrum.resolution, 'Spectral resolution')
        hdu.header['ZA_START'] = (self.spectrum.za[0],'Zenith angle [degrees]')
        hdu.header['ZA_END'] = (self.spectrum.za[1],'Zenith angle [degrees]')
        hdu.header['ALTI_STA'] = (self.spectrum.altitude[0],'Equivalent aircraft pressure altitude [feet]')
        hdu.header['ALTI_END'] = (self.spectrum.altitude[1],'Equivalent aircraft pressure altitude [feet]')

        # We should add the parameters of the fits
        
        # Add extensions
        hdu1 = self.addExtension(self.spectrum.wave,'WAVELENGTH','um',None)
        hdu2 = self.addExtension(self.panel2.flux,'FLUX','Jy',None)
        hdu3 = self.addExtension(self.panel2.uflux,'UNCORR_FLUX','Jy',None)
        if self.fit != None:
            hdu4 = self.addExtension(self.fit.continuum,'CONT','Jy',None)
            hdu5 = self.addExtension(self.fit.lines,'LINES','Jy',None)
        if self.ufit != None:
            hdu6 = self.addExtension(self.ufit.continuum,'UNCORR_CONT','Jy',None)
            hdu7 = self.addExtension(self.ufit.lines,'UNCORR_LINES','Jy',None)

        if self.fit == None:
            if self.ufit == None:
                hdul = fits.HDUList([hdu,hdu1,hdu2,hdu3])
            else:
                hdul = fits.HDUList([hdu, hdu1, hdu2, hdu3, hdu6, hdu7])
        else:
            if self.ufit == None:
                hdul = fits.HDUList([hdu, hdu1, hdu2, hdu3, hdu4, hdu5])
            else:
                hdul = fits.HDUList([hdu, hdu1, hdu2, hdu3, hdu4, hdu5, hdu6, hdu7])
                
        # Save file
        hdul.info()    
        hdul.writeto(outfile,clobber=True) # clobber true  allows rewriting
        hdul.close()
        
        return


    def saveCube(self, event):
        """
        Save cube (after cutting or cropping)
        """
        self.currentDirectory = os.getcwd()
        dlg = wx.FileDialog(
            self, message="Save as",
            defaultDir=self.currentDirectory, 
            defaultFile="cube.fits",
            wildcard="FITS files (*.fits)|*.fits",
            style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
            )
        if dlg.ShowModal() == wx.ID_OK:
            outfile = dlg.GetPath()
            print "You chose the following file:", outfile
        dlg.Destroy()
        
        # Reusable header
        header = self.spectrum.wcs.to_header()
        header['CRPIX3'] = (1,'Reference pixel')
        header['CRVAL3'] = (self.spectrum.crval3,'Wavelength at pixel reference')
        header['CDELT3'] = (self.spectrum.cdelt3,'Wavelength increment')
        header['CUNIT3'] = ('um','Wavelength unit')
        header.remove('WCSAXES')

        # Primary header
        hdu = fits.PrimaryHDU()
        hdu.header.extend(header)
        hdu.header['OBJ_NAME'] = (self.spectrum.objname, 'Object Name')
        hdu.header['FILEGPID'] = (self.spectrum.filegpid, 'File Group ID')
        hdu.header['BARYSHFT'] = (self.spectrum.baryshift, 'Barycentric shift')
        hdu.header['RESOLUN'] = (self.spectrum.resolution, 'Spectral resolution')
        hdu.header['ZA_START'] = (self.spectrum.za[0],'Zenith angle [degrees]')
        hdu.header['ZA_END'] = (self.spectrum.za[1],'Zenith angle [degrees]')
        hdu.header['ALTI_STA'] = (self.spectrum.altitude[0],'Equivalent aircraft pressure altitude [feet]')
        hdu.header['ALTI_END'] = (self.spectrum.altitude[1],'Equivalent aircraft pressure altitude [feet]')

        # Extensions
        hdu1 = self.addExtension(self.spectrum.flux,'FLUX','Jy',header)
        hdu2 = self.addExtension(self.spectrum.eflux,'ERROR','Jy',header)
        hdu3 = self.addExtension(self.spectrum.uflux,'UNCORRECTED_FLUX','Jy',header)
        hdu4 = self.addExtension(self.spectrum.euflux,'UNCORRECTED_ERROR','Jy',header)
        hdu5 = self.addExtension(self.spectrum.wave,'WAVELENGTH','um',None)
        hdu6 = self.addExtension(self.spectrum.x,'X',None,None)
        hdu7 = self.addExtension(self.spectrum.y,'Y',None,None)
        hdu8 = self.addExtension(self.spectrum.atran,'TRANSMISSION',None,None)
        hdu9 = self.addExtension(self.spectrum.response,'RESPONSE',None,None)
        hdu10 = self.addExtension(self.spectrum.exposure,'EXPOSURE_MAP',None,header)
        hdul = fits.HDUList([hdu, hdu1, hdu2, hdu3, hdu4, hdu5, hdu6, hdu7, hdu8, hdu9, hdu10])            
        hdul.info()    
        hdul.writeto(outfile,clobber=True) # clobber true  allows rewriting
        hdul.close()


    def addExtension(self,data, extname, unit, hdr):
        hdu = fits.ImageHDU()
        hdu.data = data
        hdu.header['EXTNAME']=(extname)
        if unit !=None: hdu.header['BUNIT']=(unit)
        if hdr != None: hdu.header.extend(hdr)
        return hdu
        
    def uploadImage(self, event):
        """
        Upload existing image on disk
        """
        pass

    def uploadSpectrum(self, event):
        """
        Upload existing spectrum
        """
        self.currentDirectory = os.getcwd()
        dlg = wx.FileDialog(
            self, message="Choose a spectrum",
            defaultDir=self.currentDirectory, 
            defaultFile="*.fits",
            wildcard="Spectrum file (*.fits)|*.fits",
            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST
            )
        if dlg.ShowModal() == wx.ID_OK:
            infile = dlg.GetPath()
            print "You chose the following file:", infile
        dlg.Destroy()
        self.extSpectrum = spectrum(infile)
        # plot over ax1 ?
        self.extspecLayer, = self.panel2.ax1.plot(self.extSpectrum.wave,self.extSpectrum.flux,color='purple')
        self.panel2.displayExtSpec = True
        self.panel2.refreshSpectrum()



        
    def showContours(self, event):
        self.displayContours ^= True
        if self.displayContours:
            self.drawContours()
        else:
            try:
                for coll in self.contour.collections:
                    coll.remove()
            except:
                pass
        self.panel1.Layout()
        return
        
    def drawContours(self):
        limits = self.limits
        if self.panel1.displayMethod == 'Average':
            medianflux = np.nanmean(self.spectrum.flux[limits[0]:limits[1],:,:], axis=0)
        else:
            medianflux = self.spectrum.flux[self.panel2.plane ,:,:]
        # We should check if we are plotting over the wcs,
        # i.e. we are not using an external image as background (add transform ..)
        minima = np.nanmin(medianflux)
        maxima = np.nanmax(medianflux)
        levels = np.arange(minima,maxima,(maxima-minima)/10.)
        self.contour = self.panel1.axes.contour(medianflux, levels, colors = 'cyan')
        
        return

class Toolbar1(wx.Panel):
    def __init__(self, parent, xysize):
        wx.Panel.__init__(self, parent, pos=(0,0), size=wx.Size(xysize[1]*1.2,xysize[1]*0.1))
        self.frame = parent
        self.top = self.GetTopLevelParent()
        self.SetBackgroundColour((240,248,255))

        #self.plotbut = wx.Button(self,-1,"plot", size=(40,30),pos=(10,0))
        #self.plotbut.Bind(wx.EVT_BUTTON,self.plot)
        #self.plotbut.SetToolTipString("Plot")

        self.undoBtn = self.addButton(icons.Undo.GetBitmap(),'Undo',self.Undo,'AliceBlue',(20,0))
        #self.panBtn  = self.addButton(icons.pan2.GetBitmap(),'Pan',self.Pan,'AliceBlue',(60,0))
        #self.zoomBtn  = self.addButton(icons.zoom.GetBitmap(),'Zoom',self.Zoom,'AliceBlue',(100,0))
        self.snapshotBtn  = self.addButton(icons.snapshot.GetBitmap(),'Snapshot',self.Snapshot,'AliceBlue',(60,0))
        self.dispBtn  = self.addButton(icons.lambda2.GetBitmap(),'Display wav plane or range average',self.toggleImage,'AliceBlue',(100,0))
        self.fitEllBtn  = self.addButton(icons.ellipse.GetBitmap(),'Fit ellipse inside selected region',self.top.panel1.fitEllipse,'AliceBlue',(140,0))
        self.cropBtn  = self.addButton(icons.crop.GetBitmap(),'Crop cube on selected image',self.top.panel1.cropCube,'AliceBlue',(180,0))
        #        self.reloadBtn  = self.addButton(icons.reload.GetBitmap(),'Reload original cube',self.top.reloadCube,'AliceBlue')
        self.contoursBtn  = self.addButton(icons.contours.GetBitmap(),'Overplot flux contours',self.top.showContours,'AliceBlue',(220,0))
        self.saveBtn  = self.addButton(icons.save.GetBitmap(),'Save current cube',self.top.saveCube,'AliceBlue',(260,0))
        self.uploadBtn  = self.addButton(icons.upload.GetBitmap(),'Upload image',self.top.uploadImage,'AliceBlue',(300,0))
        self.showHdr = self.addButton(icons.header.GetBitmap(),'Show header',self.showHeader,'AliceBlue',(340,0))
        wx.ToolTip.SetDelay(1000)

    def addButton(self, icon, label, function, color, position):
        button = wx.BitmapButton(self,wx.ID_NEW,bitmap=icon,size=(45,45),pos=position,style=wx.NO_BORDER | wx.BU_EXACTFIT)

        #tip = STT.SuperToolTip(label)
        #tip.SetTarget(button)
        #tip.SetDropShadow(True)
        #tip.ApplyStyle("Office 2007 Blue")
        button.SetToolTipString(label)
        button.Bind(wx.EVT_BUTTON, function)
        button.SetBackgroundColour((240,248,255))


        return button

    def Pan(self,event):
        self.top.panel1.ntoolbar.pan()    

    def Undo(self, event):
        # self.top.panel1.ntoolbar.home()
        self.top.panel1.axes.set_xlim(self.top.panel1.org_xlim)
        self.top.panel1.axes.set_ylim(self.top.panel1.org_ylim)
        self.top.panel1.Layout()


        
    def Zoom(self, event):
        self.top.panel1.ntoolbar.zoom()

    def Snapshot(self, event):
        self.top.panel1.ntoolbar.save_figure()


    def showHeader(self,event):
        #hdrvalues = self.top.spectrum.header.tostring(sep='\n')
        header = self.top.spectrum.header
        hv=header.values()
        hk= header.keys()
        hc=header.comments
        h = []
        cc = False
        for k,v,c in zip(hk,hv,hc):
            if k == 'HISTORY':
                pass
            elif k == 'COMMENT':
                if cc:
                    h.append('        {0:20}'.format(v)+''+c)
                else:
                    h.append('\n        {0:20}'.format(v)+''+c)
                    cc = True
            else:
                if cc:
                    h.append('\n')
                    cc = False
                h.append('{:15s}'.format(k)+ '\t\t=\t {0:15}'.format(v)+' \t\t\t'+c)
    
        s = '\n'.join(h)
        dlg = ScrolledMessage(self,s,'Header',size=(800,400),style=wx.TE_READONLY)
        dlg.ShowModal()
        dlg.Destroy()

    def toggleImage(self, event):
        if self.top.panel1.displayMethod == 'Average':
            self.top.shade = 'False'
            self.top.panel1.displayMethod = 'Plane'
            self.dispBtn.SetBitmap(icons.xbar2.GetBitmap())
            self.top.panel2.addSlider()
        else:
            self.top.shade = 'True'
            self.top.panel1.displayMethod = 'Average'
            self.top.panel2.removeSlider()
            self.dispBtn.SetBitmap(icons.lambda2.GetBitmap())
        print('Method is:', self.top.panel1.displayMethod)
        self.top.shadeSpectrum()
        self.top.panel2.Layout()
        self.top.panel1.refreshImage()


    
class Panel1 (wx.Panel):
    def __init__(self, parent, xysize):
        wx.Panel.__init__(self, parent, pos=(0,0), size=wx.Size(xysize[1]*1.2,xysize[1]*.9))

        # size
        self.xysize = xysize
        self.frame = parent
        self.top = self.GetTopLevelParent()

    
        # Figure
        self.figure=Figure()

        #self.figure.set_tight_layout(True)
        self.canvas=FigureCanvas(self,-1,self.figure)
        self.figure.set_facecolor('AliceBlue')
        self.figure.set_edgecolor('none')


        self.displayImage = 'Flux'
        self.displayMethod = 'Average'
        ## Toolbar
        self.ntoolbar = NavigationToolbar(self.canvas)
        self.ntoolbar.Hide()
        # Mouse wheel zooming (just to prove the concept)
        self.cidPress = self.figure.canvas.mpl_connect('scroll_event', self.onWheel)
        # Sizing
        #self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.vbox=wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas,2,wx.EXPAND,0)
        #self.vbox.Add(self.hbox,0,wx.EXPAND,0)
        self.SetSizer(self.vbox)
        self.frame.Fit()

    def showHeader(self,event):
        #hdrvalues = self.top.spectrum.header.tostring(sep='\n')
        header = self.top.spectrum.header
        hv=header.values()
        hk= header.keys()
        hc=header.comments
        h = []
        cc = False
        for k,v,c in zip(hk,hv,hc):
            if k == 'HISTORY':
                pass
            elif k == 'COMMENT':
                if cc:
                    h.append('        {0:20}'.format(v)+''+c)
                else:
                    h.append('\n        {0:20}'.format(v)+''+c)
                    cc = True
            else:
                if cc:
                    h.append('\n')
                    cc = False
                h.append('{:15s}'.format(k)+ '\t\t=\t {0:15}'.format(v)+' \t\t\t'+c)
    
        s = '\n'.join(h)
        dlg = ScrolledMessage(self,s,'Header',size=(800,400),style=wx.TE_READONLY)
        dlg.ShowModal()
        dlg.Destroy()

        
    def addButton(self, icon, label, function, color):
        button = wx.BitmapButton(self.toolbar,wx.ID_NEW,bitmap=icon,size=(33,33),style=wx.NO_BORDER | wx.BU_EXACTFIT)
        #button = pbtn.PlateButton(self.toolbar, wx.ID_NEW, bmp=icon,style = pbtn.PB_STYLE_DEFAULT | pbtn.PB_STYLE_NOBG)
        ##button = GenBitmapButton(self.toolbar, wx.ID_NEW, bitmap=icon,style = wx.NO_BORDER|wx.BU_EXACTFIT)
        button.SetToolTip(wx.ToolTip(label))
        #button.SetToolTipString(label)
        self.toolbar.AddControl(button, label)
        #button = wx.ToolBar.AddTool(self.toolbar,wx.ID_NEW,icon,shortHelpString=label)
        button.Bind(wx.EVT_BUTTON, function)
        button.SetBackgroundColour(color)
        #self.Bind(wx.EVT_TOOL,function,button)
        return button

    def Pan(self,event):
        self.ntoolbar.pan()    

    def Undo(self, event):
        self.ntoolbar.home()

    def Zoom(self, event):
        self.ntoolbar.zoom()

    def Snapshot(self, event):
        self.ntoolbar.save_figure()
        
    def OnRightDown(self, event):
        self.PopupMenu(PopupMenu1(self), (event.x, self.canvas.Size[1]-event.y))
        
    def __del__(self):
        pass
        
#    def onRadioBox(self,e): 
#        print self.rbox.GetStringSelection(),' is clicked from Radio Box'
#        self.displayImage = self.rbox.GetStringSelection()
#        self.refreshImage()


    def refreshImage(self):    
        self.view(self.top.spectrum,self.top.limits,0)
        self.top.drawEllipse()
        if self.top.displayContours: self.top.drawContours()
        self.Layout()

    def fitEllipse(self, event):
        print ("Fitting the ellipse ....")
        # Find the points inside the existing  ellipse
        ellipse = self.top.ellipse
        spectrum = self.top.spectrum
        x0,y0 = ellipse.center
        path = ellipse.get_path()
        transform = ellipse.get_patch_transform()
        npath = transform.transform_path(path)
        try:
            inpoints = spectrum.points[npath.contains_points(spectrum.points)]
            xx,yy = inpoints.T
            # Select points with flux > 10% current range
            if self.displayMethod == 'Average':
                limits = self.top.limits
                medianflux = np.nanmean(spectrum.flux[limits[0]:limits[1],:,:],axis=0)
            else:
                medianflux = spectrum.flux[self.top.panel2.plane ,:,:]
            flux = medianflux[yy,xx]
            fluxRange = np.nanmax(flux)-np.nanmin(flux)
            threshold = np.nanmin(flux)+0.1*fluxRange
            #print "threshold ",threshold
            index, = np.where(flux >= threshold)
            inpoints = inpoints[index,:]
            # Find the convex hull
            ch = ConvexHull(inpoints)
            hull_indices = ch.vertices
            hull_pts = inpoints[hull_indices, :]
            center, width, height, phi = fit_ellipse(hull_pts.T)
            # Assign new values to the ellipse
            print self.top.ellipse.angle
            self.top.ellipse.center = center
            self.top.ellipse.width = width*2.
            self.top.ellipse.height = height*2.
            self.top.ellipse.angle = phi*180./np.pi-90.
        except:
            self.top.ellipse.width = 5
            self.top.ellipse.height = 5
            self.top.ellipse.angle = 0
        # Update plots
        self.top.drawEllipse()
        self.vbox.Layout()

    def cropCube(self, event):
        print "Cropping the cube using the span tool"
        # Check the current limits
        xlimits = self.axes.get_xlim()
        ylimits = self.axes.get_ylim()
        print "xlimits ", xlimits
        print "ylimits ", ylimits
        center =  ((xlimits[0]+xlimits[1])*0.5,(ylimits[0]+ylimits[1])*0.5)
        #origin = (ylimits[0],xlimits[0])
        size = ((ylimits[1]-ylimits[0]).astype(int),(xlimits[1]-xlimits[0]).astype(int))
        print "center, size: ",center,size
        shape = self.top.spectrum.flux.shape
        print "shape is: ",shape
        if size[1] == shape[1] and size[0] == shape[2]:
            print "No cropping needed"
            return
        else:
            # Dialog window to start cropping
            dlg = wx.MessageDialog(self,"Do you really want to crop it ?", "Question",
                                    wx.YES_NO|wx.ICON_QUESTION)
            if dlg.ShowModal() == wx.ID_YES:
                self.ntoolbar.ToggleTool(self.ntoolbar.wx_ids['Zoom'], False) # turn off zoom tool
                self.ntoolbar.zoom('off')
                ima = self.top.spectrum.flux[0,:,:]
                #print ima.shape
                co = Cutout2D(ima, center, size, wcs=self.top.spectrum.wcs)
                bb = co.bbox_original
                #print "bounding box ", bb
                #print "co data shape ", co.data.shape
                self.top.spectrum.flux = self.top.spectrum.flux[:,bb[0][0]:bb[0][1]+1,bb[1][0]:bb[1][1]+1]
                self.top.spectrum.eflux = self.top.spectrum.eflux[:,bb[0][0]:bb[0][1]+1,bb[1][0]:bb[1][1]+1]
                self.top.spectrum.uflux = self.top.spectrum.uflux[:,bb[0][0]:bb[0][1]+1,bb[1][0]:bb[1][1]+1]
                self.top.spectrum.euflux = self.top.spectrum.euflux[:,bb[0][0]:bb[0][1]+1,bb[1][0]:bb[1][1]+1]
                self.top.spectrum.x = self.top.spectrum.x[bb[1][0]:bb[1][1]+1]
                self.top.spectrum.y = self.top.spectrum.y[bb[0][0]:bb[0][1]+1]
                self.top.spectrum.exposure = self.top.spectrum.exposure[:,bb[0][0]:bb[0][1]+1,bb[1][0]:bb[1][1]+1]
                self.top.spectrum.wcs = co.wcs
                # Recreate grid of points
                xi = np.arange(size[1])
                yi = np.arange(size[0])
                xi,yi = np.meshgrid(xi,yi)
                self.top.spectrum.points = np.array([np.ravel(xi),np.ravel(yi)]).transpose()
                # Modify ellipse center
                x0,y0 = self.top.ellipse.center
                x0 -= bb[1][0] # E-W direction reverses coordinates
                y0 -= bb[0][0]
                dx = self.top.spectrum.flux.shape[2]
                dy = self.top.spectrum.flux.shape[1]
                if x0 < 0 or x0 > dx: x0=dx/2
                if y0 < 0 or y0 > dy: y0=dy/2
                self.top.ellipse.center = (x0, y0)
                # Redisplay ...
                # Remove previous axes
                self.figure.delaxes(self.axes)
                # Remove and redefine ellipse
                angle = self.top.ellipse.angle
                center = (x0,y0)
                width = self.top.ellipse.width
                height = self.top.ellipse.height
                self.top.ellipse.remove()
                self.top.ellipse = Ellipse(center, width,height, edgecolor='Lime', facecolor='none',angle=angle)
                # Build new figure
                self.view(self.top.spectrum,self.top.limits,1)
                self.axes.set_xlim([0., dx])
                self.axes.set_ylim([0., dy])
                #print "ellipse axes" ,self.top.ellipse.axes.get_xlim()
                # update ellipse axes to new axes limits
                #self.top.ellipse.axes.set_xlim([0.,dx])
                #self.top.ellipse.axes.set_ylim([0.,dy])
                #self.top.ellipse.axes.set_title('')  # Remove previous title
                #self.top.ellipse.axes.xaxis.set_ticklabels([])
                #self.top.ellipse.axes.yaxis.set_ticklabels([])
                #print "AXES ",self.axes
                self.top.drawEllipse()
                if self.top.displayContours: self.top.drawContours()
                self.Layout()
                # Forget previous values
                self.ntoolbar._views.clear()
                self.ntoolbar._positions.clear()
                self.ntoolbar._update_view() 
                self.ntoolbar.Refresh()
            dlg.Destroy()        
            
    def onWheel(self,event):
        eb = event.button
        curr_xlim = self.axes.get_xlim()
        curr_ylim = self.axes.get_ylim()
        curr_x0 = (curr_xlim[0]+curr_xlim[1])*0.5
        curr_y0 = (curr_ylim[0]+curr_ylim[1])*0.5
        if eb == 'up':
            factor=0.9
        elif eb == 'down':
            factor=1.1
        new_width = (curr_xlim[1]-curr_xlim[0])*factor*0.5
        new_height= (curr_ylim[1]-curr_ylim[0])*factor*0.5
        self.axes.set_xlim([curr_x0-new_width,curr_x0+new_width])
        self.axes.set_ylim([curr_y0-new_height,curr_y0+new_height])
        self.Layout()

    def recenterImage(self,event):
        curr_xlim = self.axes.get_xlim()
        curr_ylim = self.axes.get_ylim()
        curr_width = (curr_xlim[1]-curr_xlim[0])*0.5
        curr_height= (curr_ylim[1]-curr_ylim[0])*0.5
        self.axes.set_xlim([event.xdata-curr_width,event.xdata+curr_width])
        self.axes.set_ylim([event.ydata-curr_height,event.ydata+curr_height])
        self.Layout()
        
        
    def view(self, spectrum, limits, status):

        self.figure.set_canvas(self.canvas)
        if status == 1:
            self.axes = self.figure.add_subplot(111, projection=spectrum.wcs )
            self.axes.coords[0].set_major_formatter('hh:mm:ss')
        self.axes.clear()
        self.figure.subplots_adjust(bottom=0.14,right=0.85)
        if self.displayImage == 'Flux':
            intensity = spectrum.flux
        elif self.displayImage == 'Unc. Flux':
            intensity = spectrum.uflux
        elif self.displayImage == 'Exposure':
            intensity = spectrum.exposure

        if self.displayMethod == 'Average':
            medianflux = np.nanmedian(intensity[limits[0]:limits[1],:,:],axis=0)
        else:
            medianflux = intensity[self.top.panel2.plane ,:,:]
        if status == 0:
            xlimits = self.axes.get_xlim()
            ylimits = self.axes.get_ylim()
        self.image = self.axes.imshow(medianflux, cmap='gist_heat', origin='lower', interpolation='none')
        self.figure.suptitle(self.displayImage)
        if status == 0:
            self.axes.set_xlim(xlimits)
            self.axes.set_ylim(ylimits)
        # Colorbar
        cbaxes = self.figure.add_axes([0.9,0.1,0.02,0.8])
        self.figure.colorbar(self.image, cax=cbaxes)
        # Sliders to adjust intensity
        self.ax_cmin = self.figure.add_axes([0.1, 0.01, 0.8, 0.01])
        self.ax_cmax = self.figure.add_axes([0.1, 0.04, 0.8, 0.01])
        #self.ax_cmin = self.figure.add_axes([0.85, 0.1, 0.02, 0.8])
        #self.ax_cmax = self.figure.add_axes([0.88, 0.1, 0.02, 0.8])
        self.ax_cmin.clear()
        self.ax_cmax.clear()
        vmin0=np.nanmin(medianflux); vmax0=np.nanmax(medianflux)
        d0 = (vmax0-vmin0)/20.
        self.s_cmin = Slider(self.ax_cmin, 'low', vmin0-d0, vmax0+d0, valinit=vmin0, facecolor='goldenrod')
        self.s_cmax = Slider(self.ax_cmax, 'high', vmin0-d0, vmax0+d0, valinit=vmax0, facecolor='goldenrod')
        self.s_cmin.valtext.set_visible(False)
        self.s_cmax.valtext.set_visible(False)
        self.slider1=self.s_cmin.on_changed(self.updateScale)
        self.slider2=self.s_cmax.on_changed(self.updateScale)

        # find maximum of medianflux
        try:
            ijmax = np.unravel_index(np.nanargmax(medianflux), medianflux.shape)
        except:
            print ("image is all NaNs")
            ijmax = [0,0]
        self.canvas.draw()
        self.canvas.Refresh()

        # Conserve original limits
        self.org_xlim = self.axes.get_xlim()
        self.org_ylim = self.axes.get_ylim()
        
        return ijmax

    def updateScale(self,val):
        _cmin = self.s_cmin.val
        _cmax = self.s_cmax.val
        self.image.set_clim([_cmin, _cmax])
        self.figure.canvas.draw_idle()

        
        

        
    def toggleImage(self, event):
        if self.displayMethod == 'Average':
            self.top.shade = 'False'
            self.displayMethod = 'Plane'
            self.dispBtn.SetBitmap(icons.xbar2.GetBitmap())
            self.top.panel2.addSlider()
        else:
            self.top.shade = 'True'
            self.displayMethod = 'Average'
            self.top.panel2.removeSlider()
            self.dispBtn.SetBitmap(icons.lambda2.GetBitmap())
        print('Method is:', self.displayMethod)
        self.top.shadeSpectrum()
        self.top.panel2.Layout()
        self.refreshImage()



class Toolbar2(wx.Panel):
    def __init__(self, parent, xysize):
        wx.Panel.__init__(self, parent, pos=(0,0), size=wx.Size(xysize[1]*1.2,xysize[1]*0.1))
        self.frame = parent
        self.top = self.GetTopLevelParent()
        self.SetBackgroundColour((240,248,255)) # AliceBlue

        self.nextBtn  = self.addButton(icons.next.GetBitmap(),'Load another cube',self.top.nextCube,'AliceBlue',(0,0))
        self.reloadBtn  = self.addButton(icons.reload.GetBitmap(),'Reload original cube',self.top.reloadCube,'AliceBlue',(40,0))
        self.quitBtn    = self.addButton(icons.quit.GetBitmap(),'Quit',self.top.onQuit,'AliceBlue',(80,0))

        self.undoBtn = self.addButton(icons.Undo.GetBitmap(),'Undo',self.Undo,'AliceBlue',(200,0))
        self.panBtn  = self.addButton(icons.pan2.GetBitmap(),'Pan',self.Pan,'AliceBlue',(240,0))
        self.zoomBtn  = self.addButton(icons.zoom.GetBitmap(),'Zoom',self.Zoom,'AliceBlue',(280,0))
        self.snapshotBtn  = self.addButton(icons.snapshot.GetBitmap(),'Snapshot',self.Snapshot,'AliceBlue',(320,0))
        self.cutBtn    = self.addButton(icons.cut.GetBitmap(),'Cut cube as selected spectrum range',self.top.cutCube,'AliceBlue',(360,0))
        self.fitBtn    = self.addButton(icons.gauss.GetBitmap(),'Fit lines and continuum',self.top.fitGauss,'AliceBlue',(400,0))
        self.saveBtn    = self.addButton(icons.save.GetBitmap(),'Save spectrum/fit',self.top.saveSpectrum,'AliceBlue',(440,0))
        self.uploadBtn    = self.addButton(icons.upload.GetBitmap(),'Upload spectrum',self.top.uploadSpectrum,'AliceBlue',(480,0))

        wx.ToolTip.SetDelay(1000)

    def addButton(self, icon, label, function, color, position):
        button = wx.BitmapButton(self,wx.ID_NEW,bitmap=icon,size=(45,45),pos=position,style=wx.NO_BORDER | wx.BU_EXACTFIT)
        button.SetToolTipString(label)
        button.Bind(wx.EVT_BUTTON, function)
        button.SetBackgroundColour((240,248,255))
        return button

    def Pan(self,event):
        self.top.panel2.ntoolbar.pan()    

    def Undo(self, event):
        self.top.panel2.ntoolbar.home()

    def Zoom(self, event):
        self.top.panel2.ntoolbar.zoom()

    def Snapshot(self, event):
        self.top.panel2.ntoolbar.save_figure()


        
class Panel2 (wx.Panel):
    def __init__(self, parent, xysize):
        wx.Panel.__init__(self, parent, pos=(xysize[1],0), size=wx.Size(xysize[0]-xysize[1]*1.2,xysize[1]*0.9))

        self.frame = parent
        self.top = self.GetTopLevelParent()
        
        # Figure
        self.figure=Figure()
        self.figure.set_facecolor('AliceBlue')
        self.figure.set_edgecolor('none')
        self.figure.set_tight_layout(True)
        self.canvas=FigureCanvas(self,-1,self.figure)
        self.ax1 = self.figure.add_subplot(111)
        self.ax2 = self.ax1.twinx()
        self.ax3 = self.ax1.twinx()
        self.ax4 = self.ax1.twinx()

        # Lines
        self.Lines = define_lines()


        # Toolbar
        self.ntoolbar = NavigationToolbar(self.canvas)
        self.ntoolbar.Hide()
        
        # Size
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas,1,wx.LEFT|wx.TOP|wx.GROW,0)
        self.SetSizer(self.vbox)
        self.frame.Fit()
        
        #self.cidPress = self.figure.canvas.mpl_connect('button_press_event', self.onPress)
        #self.cidrelease = self.figure.canvas.mpl_connect('button_release_event', self.onRelease)
        #self.cidmotion = self.figure.canvas.mpl_connect('motion_notify_event', self.onMotion)

        self.displayFlux = True
        self.displayUFlux = False
        self.displayATran = False
        self.displayExposure = False
        self.displayLines = False
        self.displayLineFit = False
        self.displayUncLineFit = False
        self.displayExtSpec = False
        
    def __del__(self):
        pass


    def addButton(self, icon, label, function, color):
        button = wx.BitmapButton(self.toolbar,wx.ID_NEW,bitmap=icon,size=(33,33),style=wx.NO_BORDER | wx.BU_EXACTFIT)
        #button = pbtn.PlateButton(self.toolbar, wx.ID_NEW, bmp=icon,style = pbtn.PB_STYLE_DEFAULT)
        button.SetToolTip(wx.ToolTip(label))
        self.toolbar.AddControl(button, label)
        button.Bind(wx.EVT_BUTTON, function)
        button.SetBackgroundColour(color)
        return button


    def Pan(self,event):
        self.ntoolbar.pan()    

    def Undo(self, event):
        self.ntoolbar.home()

    def Zoom(self, event):
        self.ntoolbar.zoom()

    def Snapshot(self, event):
        self.ntoolbar.save_figure()
        
        
    def draw(self,spectrum,ellipse):
        x0,y0 = ellipse.center
        path = ellipse.get_path()
        transform = ellipse.get_patch_transform()
        npath = transform.transform_path(path)
        inpoints = spectrum.points[npath.contains_points(spectrum.points)]
        xx,yy = inpoints.T
        self.flux = np.nansum(spectrum.flux[:,yy,xx],axis=1)
        self.eflux = np.sqrt(np.nansum(spectrum.eflux[:,yy,xx]**2,axis=1))
        self.uflux = np.nansum(spectrum.uflux[:,yy,xx],axis=1)
        self.euflux = np.sqrt(np.nansum(spectrum.euflux[:,yy,xx]**2,axis=1))
        expos = np.nansum(spectrum.exposure[:,yy,xx],axis=1).astype(float)
        self.ax1.clear()
        self.ax2.clear()
        self.ax3.clear()
        self.ax4.clear()
        self.fluxLayer, = self.ax1.step(spectrum.wave,self.flux,color='blue')
        self.ax1.grid(True,which='both')
        # Atmosphere transmission shifted by the barycentric shift to be aligned to the corrected flux
        self.atranLayer, = self.ax2.plot(spectrum.wave*(1.+spectrum.baryshift), spectrum.atran,color='red')
        self.ax2.set_ylim([0,1.1])
        self.ax4.tick_params(labelright='off',right='off')
        self.exposureLayer, = self.ax3.step(spectrum.wave, expos, color='orange')
        ymax = np.nanmax(self.flux); ymin = np.nanmin(self.flux)
        yumax = np.nanmax(self.uflux); yumin = np.nanmin(self.uflux)
        if yumax > ymax: ymax=yumax
        if yumin < ymin: ymin=yumin
        self.ax1.set_ylim([ymin, ymax+(ymax-ymin)/8.])
        self.ax3.set_ylim([0.5,np.nanmax(expos)*1.04])
        # Uncorrected flux shifted to be aligned to the corrected flux
        self.ufluxLayer, = self.ax4.step(spectrum.wave*(1.+spectrum.baryshift),self.uflux,color='green')
        # Put same y-axis limits on ax4 as ax1:
        self.ax4.set_ylim(self.ax1.get_ylim())
        self.ax1.set_xlim([np.min(spectrum.wave),np.max(spectrum.wave)])
        self.ax1.set_title(spectrum.objname+" ["+spectrum.filegpid+"] @ "+spectrum.obsdate)
        
        # Line names as annotations
        LinesNames = self.Lines.keys()
        self.annotations = []
        xlim = self.ax1.get_xlim()
        ylim = self.ax1.get_ylim()
        dy = ylim[1]-ylim[0]
        font = FontProperties(family='DejaVu Sans',  size=12)
        for line in LinesNames:
            nline = self.Lines[line][0]
            wline = self.Lines[line][1]*(1.+self.top.spectrum.redshift)
            if (wline < xlim[1] and wline > xlim[0]):
                wdiff = abs(spectrum.wave - wline)
                y = self.flux[(wdiff == wdiff.min())]
                y1 = y  # ylim[0]+dy/2.
                #y2 = np.max([y + 0.1*dy,ylim[1]-0.05*dy])
                if (ylim[1]-(y+0.2*dy)) > ((y-0.2*dy)-ylim[0]):
                    y2 = y+0.2*dy
                else:
                    y2 = y-0.2*dy
                annotation = self.ax1.annotate(nline, xy=(wline, y1), xytext=(wline, y2),color='blue', alpha=0.4,
                                                arrowprops=dict(color='blue',facecolor='y', arrowstyle='-',alpha=0.4),
                                                rotation = 90, fontstyle = 'italic', fontproperties=font, visible=self.displayLines,)
                self.annotations.append(annotation)        

        # Draw fit (if exists)
        if self.top.fit != None:
            fit = self.top.fit
            if self.top.fit != None:
                fit = self.top.fit
                self.top.contfitLayer, = self.ax1.plot(fit.wave,fit.continuum,color='cyan')
                self.top.linefitLayer, = self.ax1.plot(fit.wave, fit.lines+fit.continuum, color='cyan')
            if self.top.ufit != None:
                fit = self.top.ufit    
                self.top.ucontfitLayer, = self.ax4.plot(fit.wave,fit.continuum,color='palegreen')
                self.top.ulinefitLayer, = self.ax4.plot(fit.wave, fit.lines+fit.continuum, color='palegreen')

        # Draw ext spectrum (if exists)
        if self.top.extSpectrum != None:
            s = self.top.extSpectrum
            self.top.extspecLayer, = self.ax1.plot(s.wave,s.flux,color='purple')
            
        self.refreshSpectrum()

    def refreshSpectrum(self):
        self.fluxLayer.set_visible(self.displayFlux)
        self.ufluxLayer.set_visible(self.displayUFlux)
        self.atranLayer.set_visible(self.displayATran)
        if self.displayATran:
            self.ax2.get_yaxis().set_tick_params(labelright='on',right='on')            
            self.ax2.get_yaxis().set_tick_params(which='both', direction='out',colors='red')
        else:
            self.ax2.get_yaxis().set_tick_params(labelright='off',right='off')            
        self.exposureLayer.set_visible(self.displayExposure)
        if self.displayExposure:
            self.ax3.tick_params(labelright='on',right='on',direction='in',pad=-30,colors='orange')
        else:
            self.ax3.get_yaxis().set_tick_params(labelright='off',right='off')
        for annotation in self.annotations:
            annotation.set_visible(self.displayLines)

        if self.top.fit != None:
            self.top.contfitLayer.set_visible(self.displayLineFit)
            self.top.linefitLayer.set_visible(self.displayLineFit)
        if self.top.ufit != None:
            self.top.ucontfitLayer.set_visible(self.displayUncLineFit)
            self.top.ulinefitLayer.set_visible(self.displayUncLineFit)

        if self.top.extSpectrum != None:
            self.top.extspecLayer.set_visible(self.displayExtSpec)
            
        self.Layout()

    def onSlider(self, event):
        self.plane = self.sld.GetValue()
        #        print('Event', self.plane)
        self.top.panel1.view(self.top.spectrum,[0,1],0)
        self.top.drawEllipse()
        if self.top.displayContours: self.top.drawContours()
        self.top.panel1.Layout()

    def addSlider(self):
        self.plane = self.top.spectrum.n/2
        self.sld = wx.Slider(self, -1,  self.plane, 0, self.top.spectrum.n, (-1,-1), (150,-1),wx.SL_LABELS)
        #self.sld.SetThumbLength(5)
        #rect = self.sld.GetRect()
        #self.sld.ClearBackground()
        self.sld.SetBackgroundColour((240,248,255))
        #self.sld.SetBackgroundColour('AliceBlue')
        #self.frame.RefreshRect(rect)
        #w,h = self.sld.GetClientSize()
        #self.sld.SetClientSize((w+1,h))
        #self.sld.SetClientSize((w,h))
        #self.sld.Bind(wx.EVT_SCROLL_THUMBRELEASE,self.onSlider)
        self.sld.Bind(wx.EVT_SLIDER,self.onSlider)
        self.hbox.Add(self.sld,wx.EXPAND,0)
        self.vbox.Layout()
        self.frame.Fit()
        self.top.Fit()

    def removeSlider(self):
        self.hbox.Hide(self.sld)
        self.hbox.Remove(self.sld)
        self.vbox.Layout()
        self.frame.Fit()
        self.top.Fit()

    def OnRightDown(self, event):
#        self.PopupMenu(PopupMenu2(self), (event.x, event.y))
        self.PopupMenu(PopupMenu2(self), (event.x, self.canvas.Size[1]-event.y))
        #        if self.canvas.HasCapture(): self.canvas.ReleaseMouse()



"""
Other classes
"""

# Choice of images on panel 1        
class PopupMenu1(wx.Menu):

    def __init__(self, parent):
        super(PopupMenu1, self).__init__()

        self.parent = parent
        self.addItem('Show flux', self.showFlux, 'Flux')
        self.addItem('Show uncorrected flux', self.showUncFlux, 'Unc. Flux')
        self.addItem('Show exposure', self.showExposure, 'Exposure')
        
    def addItem(self, comment, action, status):
        item = wx.MenuItem(self, wx.NewId(), comment, kind = wx.ITEM_CHECK)
        self.AppendItem(item)
        self.Bind(wx.EVT_MENU, action, item)
        if self.parent.displayImage == status: item.Check(True)

    def showFlux(self, event):
        self.parent.displayImage = 'Flux'
        self.parent.refreshImage()

    def showUncFlux(self, event):
        self.parent.displayImage = 'Unc. Flux'
        self.parent.refreshImage()

    def showExposure(self, event):
        self.parent.displayImage = 'Exposure'
        self.parent.refreshImage()


# Choice of curves on panel 2
class PopupMenu2(wx.Menu):
    def __init__(self, parent):
        super(PopupMenu2, self).__init__()
        self.parent = parent

        self.addItem('Show flux', self.showFlux, self.parent.displayFlux)
        self.addItem('Show uncorrected flux', self.showUncFlux, self.parent.displayUFlux)
        self.addItem('Show atmospheric transmission', self.showAtmTransm, self.parent.displayATran)
        self.addItem('Show exposure', self.showExposure, self.parent.displayExposure)
        self.addItem('Show lines', self.showLines, self.parent.displayLines)
        self.addItem('Show fit', self.showFit, self.parent.displayLineFit)
        self.addItem('Show uncorrected fit', self.showUncFit, self.parent.displayUncLineFit)
        self.addItem('Show external spectrum', self.showExtSpec, self.parent.displayExtSpec)


    def addItem(self, comment, action, status):
        item = wx.MenuItem(self, wx.NewId(), comment, kind = wx.ITEM_CHECK)
        self.AppendItem(item)
        self.Bind(wx.EVT_MENU, action, item)
        if status: item.Check(True)
                            
    def showFlux(self, event):
        self.parent.displayFlux ^= True
        self.parent.refreshSpectrum()

    def showUncFlux(self, event):
        self.parent.displayUFlux ^= True
        self.parent.refreshSpectrum()

    def showAtmTransm(self, event):
        self.parent.displayATran ^= True
        self.parent.refreshSpectrum()

    def showExposure(self, event):
        self.parent.displayExposure ^= True
        self.parent.refreshSpectrum()

    def showLines(self, event):
        self.parent.displayLines ^= True
        self.parent.refreshSpectrum()

    def showFit(self, event):
        self.parent.displayLineFit ^= True
        self.parent.refreshSpectrum()

    def showUncFit(self, event):
        self.parent.displayUncLineFit ^= True
        self.parent.refreshSpectrum()

    def showExtSpec(self, event):
        self.parent.displayExtSpec ^= True
        self.parent.refreshSpectrum()


"""
Main code
"""      
  
class MyApp(wx.App):
    def OnInit(self):
        frame = MyFrame(None, -1, 'frame')
        frame.Show(True)
        # Welcome message on panel 1
        frame.welcomeMessage()
        return True


# In the case of MacOS-X, re-executes with pythonw
#if len(sys.argv) == 1:
#    reexec_with_pythonw()
# Starting the code
app = MyApp(0)
app.MainLoop()
