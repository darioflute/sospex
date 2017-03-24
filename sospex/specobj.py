import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

class specCube(object):
    """ spectral cube """
    def __init__(self, infile):
#        wx.BeginBusyCursor()
        hdl = fits.open(infile)
        header = hdl[0].header
        self.filename = infile
        self.wcs = WCS(header).celestial
        self.crval3 = header['CRVAL3']
        self.cdelt3 = header['CDELT3']
        self.objname = header['OBJ_NAME']
        self.filegpid = header['FILEGPID']
        self.baryshift = header['BARYSHFT']
        self.pixscale = header['PIXSCAL']
        self.resolution = header['RESOLUN']
        self.za = (header['ZA_START'], header['ZA_END'])
        self.altitude = (header['ALTI_STA'],header['ALTI_END'])
        self.obsdate = header['DATE-OBS']
        self.flux = hdl['FLUX'].data
        self.eflux = hdl['ERROR'].data
        self.uflux = hdl['UNCORRECTED_FLUX'].data
        self.euflux = hdl['UNCORRECTED_ERROR'].data
        self.wave = hdl['WAVELENGTH'].data
        self.x = hdl['X'].data
        self.y = hdl['Y'].data
        self.atran = hdl['TRANSMISSION'].data
        self.response = hdl['RESPONSE'].data
        self.exposure = hdl['EXPOSURE_MAP'].data
        self.n = self.wave.size
        # Create a grid of points
        xi = np.arange(hdl['FLUX'].header['NAXIS1'])
        yi = np.arange(hdl['FLUX'].header['NAXIS2'])
        xi,yi = np.meshgrid(xi,yi)
        self.points = np.array([np.ravel(xi),np.ravel(yi)]).transpose()
#        wx.EndBusyCursor()

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
