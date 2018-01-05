import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

class specCube(object):
    """ spectral cube """
    def __init__(self, infile):

        c = 299792458.0  # speed of light in m/s 

        hdl = fits.open(infile)
        header = hdl['PRIMARY'].header
        self.header = header
        self.filename = infile
        try:
            self.instrument = header['INSTRUME']        
        except:
            origin = header['ORIGIN']
            if origin == 'GILDAS Consortium':
                self.instrument = 'GREAT'
        try:
            self.obsdate = header['DATE-OBS']
        except:
            self.obsdate = header['DATE']
        self.wcs = WCS(header).celestial
        self.crpix3 = header['CRPIX3']
        self.crval3 = header['CRVAL3']
        self.cdelt3 = header['CDELT3']


        if self.instrument == 'FIFI-LS':        
            self.objname = header['OBJ_NAME']
            self.filegpid = header['FILEGPID']
            self.baryshift = header['BARYSHFT']
            self.pixscale = header['PIXSCAL']
            self.resolution = header['RESOLUN']
            self.za = (header['ZA_START'], header['ZA_END'])
            self.altitude = (header['ALTI_STA'],header['ALTI_END'])
            try:
                self.redshift = hdl['REDSHIFT'].data
            except:
                self.redshift = 0.0
        elif self.instrument == 'GREAT':
            self.objname = header['OBJECT']
            self.redshift = header['VELO-LSR'] # in m/s
            self.redshift /= c
            self.pixscale,ypixscale = proj_plane_pixel_scales(self.wcs)*3600. # Pixel scale in arcsec
        else:
            print('This is not a standard spectral cube')


        if self.instrument == 'FIFI-LS':
            self.flux = hdl['FLUX'].data
            self.eflux = hdl['ERROR'].data
            self.uflux = hdl['UNCORRECTED_FLUX'].data
            self.euflux = hdl['UNCORRECTED_ERROR'].data
            self.wave = hdl['WAVELENGTH'].data
            self.n = len(self.wave)
            self.l0 = np.nanmedian(self.wave)
            #self.vel = np.zeros(self.n)  # prepare array of velocities
            self.x = hdl['X'].data
            self.y = hdl['Y'].data
            self.atran = hdl['TRANSMISSION'].data
            self.response = hdl['RESPONSE'].data
            self.exposure = hdl['EXPOSURE_MAP'].data
        elif self.instrument == 'GREAT':
            self.n = header['NAXIS3']
            naxes = header['NAXIS']
            if naxes == 4:
                self.flux = (hdl['PRIMARY'].data)[0,:,:,:]
            else:
                self.flux = hdl['PRIMARY'].data
            eta_fss=0.97
            eta_mb =0.67
            calib = 971.
            factor = calib*eta_fss*eta_mb
            self.flux *= factor   # Transformed from temperature to S_nu [Jy]            
            nu0 = header['RESTFREQ']  # MHz
            l0 = c/nu0  # in micron
            vel = -self.cdelt3*self.crpix3+self.cdelt3*np.arange(self.n)+self.crval3
            self.l0 = l0
            #self.vel = vel
            self.wave = l0 + l0*vel/c

            
        hdl.close()
        # Create a grid of points
        
        self.nz,self.ny,self.nx = np.shape(self.flux)
        xi = np.arange(self.nx); yi = np.arange(self.ny)
        xi,yi = np.meshgrid(xi,yi)
        self.points = np.array([np.ravel(xi),np.ravel(yi)]).transpose()


class ExtSpectrum(object):
    """ class for external spectrum """
    def __init__(self, infile):
        hdl = fits.open(infile)
        header = hdl[0].header
        # Assuming flux and wavelength are conserved in the respective extensions
        self.flux = hdl['FLUX'].data # in Jansky
        self.wave = hdl['WAVELENGTH'].data  # in micronmeters
        try:
            self.redshift = header['REDSHIFT']
        except:
            self.redshift = 0.
        hdl.close()
        # If wavelength is in the header, use:
        #self.wave = crval+(np.arange(naxis)-crpix+1)*cdelt

class Spectrum(object):
    """ class to define a spectrum """
    def __init__(self, wave, flux, uflux=None, exposure=None, atran=None, instrument=None, baryshift=None, redshift=None, l0=None, area=None):
        self.wave = wave
        self.flux = flux
        if exposure is not None:
            self.exposure = exposure
        if atran is not None:
            self.atran = atran
        if uflux is not None:
            self.uflux = uflux
        if instrument is not None:
            self.instrument = instrument
        if baryshift is not None:
            self.baryshift=baryshift
        if redshift is not None:
            self.redshift = redshift
        if l0 is not None:
            self.l0 = l0
        if area is not None:
            self.area = area
            
