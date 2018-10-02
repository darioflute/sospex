import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

class specCube(object):
    """ spectral cube """
    def __init__(self, infile):
        import time
        t = time.process_time()
        c = 299792458.0  # speed of light in m/s 
        # Option None seems faster than False
        hdl = fits.open(infile,memmap=None)
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

        if self.instrument == 'FIFI-LS':        
            self.wcs = WCS(header).celestial
            self.crpix3 = header['CRPIX3']
            self.crval3 = header['CRVAL3']
            self.cdelt3 = header['CDELT3']
            self.objname = header['OBJ_NAME']
            try:
                self.filegpid = header['FILEGPID']
            except:
                self.filegpid = 'Unknown'
            self.baryshift = header['BARYSHFT']
            self.pixscale = header['PIXSCAL']
            self.resolution = header['RESOLUN']
            self.za = (header['ZA_START'], header['ZA_END'])
            self.altitude = (header['ALTI_STA'],header['ALTI_END'])
            try:
                self.redshift = hdl['REDSHIFT'].data
            except:
                self.redshift = 0.0
            self.flux = hdl['FLUX'].data
            self.eflux = hdl['ERROR'].data
            self.uflux = hdl['UNCORRECTED_FLUX'].data
            self.euflux = hdl['UNCORRECTED_ERROR'].data
            self.wave = hdl['WAVELENGTH'].data
            self.l0 = np.nanmedian(self.wave)
            self.n = len(self.wave)
            #self.vel = np.zeros(self.n)  # prepare array of velocities
            self.x = hdl['X'].data
            self.y = hdl['Y'].data
            #self.atran = hdl['TRANSMISSION'].data
            try:
                utran = hdl['UNSMOOTHED_TRANSMISSION'].data
                w = utran[0,:]
                t = utran[1,:]
                self.atran = np.interp(self.wave,w,t)  # Interpolation at the resolution of the wavelength grid
            except:
                print('The unsmoothed transmission is not available')
                self.atran = hdl['TRANSMISSION'].data   
            self.response = hdl['RESPONSE'].data
            self.exposure = hdl['EXPOSURE_MAP'].data
        elif self.instrument == 'GREAT':
            #self.cmin = header['DATAMIN']
            #self.cmax = header['DATAMAX']
            self.wcs = WCS(header).celestial
            self.crpix3 = header['CRPIX3']
            self.crval3 = header['CRVAL3']
            self.cdelt3 = header['CDELT3']
            self.objname = header['OBJECT']
            self.redshift = header['VELO-LSR'] # in m/s
            self.redshift /= c
            self.pixscale,ypixscale = proj_plane_pixel_scales(self.wcs)*3600. # Pixel scale in arcsec
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
            self.wave = l0 + l0*vel/c
        elif self.instrument == 'PACS':
            """ Case of PACS spectral cubes """
            print('This is a PACS spectral cube')
            self.objname = header['OBJECT']
            try:
                self.redshift = header['REDSHFTV']*1000. # in km/s
                self.redshift /= c
            except:
                self.redshift = 0.
            print('Object is ',self.objname)
            self.flux = hdl['image'].data
            print('Flux read')
            self.exposure = hdl['coverage'].data
            print('Coverage read')
            wave = hdl['wcs-tab'].data
            print('Wvl read')
            nwave = len(np.shape(wave['wavelen']))
            if nwave == 3:
                self.wave = np.concatenate(wave['wavelen'][0])
            else:
                self.wave = np.concatenate(wave['wavelen'])
            self.l0 = np.nanmedian(self.wave)
            self.n = len(self.wave)
            print('Length of wavelength ',self.n)
            header = hdl['IMAGE'].header
            self.header = header
            print('Header ',header)
            hdu = fits.PrimaryHDU(self.flux)
            hdu.header
            hdu.header['CRPIX1']=header['CRPIX1']
            hdu.header['CRPIX2']=header['CRPIX2']
            hdu.header['CDELT1']=header['CDELT1']
            hdu.header['CDELT2']=header['CDELT2']
            hdu.header['CRVAL1']=header['CRVAL1']
            hdu.header['CRVAL2']=header['CRVAL2']
            hdu.header['CTYPE1']=header['CTYPE1']
            hdu.header['CTYPE2']=header['CTYPE2']
            self.wcs = WCS(hdu.header).celestial
            print('astrometry ', self.wcs)
            self.pixscale,ypixscale = proj_plane_pixel_scales(self.wcs)*3600. # Pixel scale in arcsec
            self.crpix3 = 1
            w = self.wave
            self.crval3 = w[0]
            self.cdelt3 = np.median(w[1:]-w[:-1])
        else:
            print('This is not a standard spectral cube')
        hdl.close()

        # index of the ref wavelength
        self.n0 = np.argmin(np.abs(self.wave - self.l0))
        print('ref wavelength at n: ', self.n0)
        # Create a grid of points
        self.nz, self.ny, self.nx = np.shape(self.flux)
        xi = np.arange(self.nx); yi = np.arange(self.ny)
        xi,yi = np.meshgrid(xi, yi)
        self.points = np.array([np.ravel(xi), np.ravel(yi)]).transpose()
           
        elapsed_time = time.process_time() - t
        print('Reading of cube completed in ', elapsed_time,' s')

    def computeExpFromNan(self):
        """Compute an exposure cube from NaN in the flux cube."""
        self.exposure = np.ones(np.shape(self.flux))
        mask = self.flux == np.nan
        self.exposure[mask] = 0
        

class ExtSpectrum(object):
    """ class for external spectrum """
    def __init__(self, infile):
        hdl = fits.open(infile, memmap=False)
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
    def __init__(self, wave, flux, eflux=None, uflux=None, exposure=None, atran=None, instrument=None, baryshift=None, redshift=None, l0=None, area=None):
        self.wave = wave
        self.flux = flux
        if eflux is not None:
            self.eflux = eflux
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
        self.continuum = np.ones(len(wave))
        self.continuum.fill(np.nan)
