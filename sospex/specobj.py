import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

class specCube(object):
    """ spectral cube """
    def __init__(self, infile):

        c = 299792458.0  # speed of light in m/s 

        hdl = fits.open(infile)
        header = hdl['PRIMARY'].header
        self.header = header
        self.filename = infile
        self.instrument = header['INSTRUME']        
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
        else:
            print('This is not a standard spectral cube')
            
        self.obsdate = header['DATE-OBS']
        self.wcs = WCS(header).celestial
        self.crpix3 = header['CRPIX3']
        self.crval3 = header['CRVAL3']
        self.cdelt3 = header['CDELT3']

        if self.instrument == 'FIFI-LS':
            self.flux = hdl['FLUX'].data
            self.eflux = hdl['ERROR'].data
            self.uflux = hdl['UNCORRECTED_FLUX'].data
            self.euflux = hdl['UNCORRECTED_ERROR'].data
            self.wave = hdl['WAVELENGTH'].data
            self.n = len(self.wave)
            self.x = hdl['X'].data
            self.y = hdl['Y'].data
            self.atran = hdl['TRANSMISSION'].data
            self.response = hdl['RESPONSE'].data
            self.exposure = hdl['EXPOSURE_MAP'].data
        elif self.instrument == 'GREAT':
            self.n = header['NAXIS3']
            self.flux = hdl['PRIMARY'].data
            eta_fss=0.97
            eta_mb =0.67
            calib = 971.
            factor = calib*eta_fss*eta_mb
            self.flux *= factor   # Transformed from temperature to S_nu [Jy]            
            nu0 = header['RESTFREQ']  # Hz
            l0 = c/nu0*1.e6
            vel = -self.cdelt3*self.crpix3+self.cdelt3*np.arange(self.n)+self.crval3
            self.wave = l0 + l0*vel/c

        hdl.close()
        # Create a grid of points
        nz,ny,nx = np.shape(self.flux)
        xi = np.arange(nx); yi = np.arange(ny)
        xi,yi = np.meshgrid(xi,yi)
        self.points = np.array([np.ravel(xi),np.ravel(yi)]).transpose()


class spectrum(object):
    """ class for external spectrum """
    def __init__(self, infile):
        hdl = fits.open(infile)
        header = hdl[0].header
        self.flux = hdl[0].data
        crval = header['CRVAL']
        crpix = header['CRPIX']
        cdelt = header['CDELT']
        naxis = header['NAXIS1']
        # I should consider the unit here,
        # now I assume everything in um ...
        self.wave = crval+(np.arange(naxis)-crpix+1)*cdelt
