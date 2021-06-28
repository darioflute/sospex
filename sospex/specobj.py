import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales


def computeBaryshift(header):
    from astropy.coordinates import FK5, solar_system, UnitSphericalRepresentation, CartesianRepresentation
    from astropy.time import Time
    import astropy.units as u
    import astropy.constants as const

    equinox = header['TELEQUI']
    time = Time(header['DATE-OBS'])
    sc = FK5(header['TELRA'] * u.hourangle, header['TELDEC'] * u.deg, equinox=equinox)
    sc_cartesian = sc.represent_as(UnitSphericalRepresentation).\
    represent_as(CartesianRepresentation)
    _, ev = solar_system.get_body_barycentric_posvel('earth', time)
    helio_vel = sc_cartesian.dot(ev).to(u.km / u.s)
    # Compute solar velocity wrt LSR
    # http://conga.oan.es/~alonso/jparsec/doc/jparsec/ephem/stars/StarEphem.html
    # Location of the LSR in J2000, 18h 03m 50.2s, 30° 00' 16.8", and 19.5 km/s of speed set as radius.
    sunpos = FK5((18+3/60.+50.2/3600.) * u.hourangle, (30+16.8/3600.) * u.deg, equinox='J2000') # Solar apex
    # Precess to current equinox
    #sunpos = sunpos.transform_to(FK5(equinox=equinox))
    #  Component of peculiar motion of the sun relative to the LSR (local star ref system)
    sun_v0 = 13.4 * u.km / u.s   # wikipedia - relative to average velocity of other stars
    sun_cartesian = sunpos.represent_as(UnitSphericalRepresentation).\
    represent_as(CartesianRepresentation)
    sun_vel = sc_cartesian.dot(sun_cartesian) * sun_v0
    print('Heliocentric correction, sun velocity component ', helio_vel, sun_vel)
    vlsr = helio_vel + sun_vel
    speed_of_light = const.c.to(vlsr.unit)
    return vlsr / speed_of_light

def computeBaryshiftAstropy(header):
    # Baryshift from astropy - it's different
    from astropy import units as u
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, EarthLocation
    lat = header['LAT_STA']
    lon = header['LON_STA']
    height = header['ALTI_STA'] * 0.3048 # in meter
    ra = header['TELRA']
    dec = header['TELDEC']
    sofia = EarthLocation.from_geodetic(lat=lat*u.deg, lon=lon*u.deg, height=height*u.m)
    sc = SkyCoord(ra=ra*u.hourangle, dec=dec*u.deg)
    dateobs = header['DATE-OBS']
    barycorr = sc.radial_velocity_correction(kind='barycentric', obstime=Time(dateobs), location=sofia)  
    heliocorr = sc.radial_velocity_correction('heliocentric', obstime=Time(dateobs), location=sofia) 
    print('Barycentric correction ', barycorr.to(u.km/u.s), heliocorr.to(u.km/u.s))
    c = 299792.458
    #self.baryshift = barycorr.to(u.km/u.s).value/c
    return(barycorr.to(u.km/u.s).value/c)


class specCubeAstro(object):
    """ spectral cube - read with AstroPy routines """
    def __init__(self, infile):
        import time
        t = time.process_time()
        
        print('Using Astropy ...')
        # Option None seems faster than False
        # ignore blank to speed up reading of GREAT cubes
        hdl = fits.open(infile, memmap=None, ignore_blank=True)
        header = hdl['PRIMARY'].header
        self.header = header
        self.filename = infile
        try:
            self.instrument = header['INSTRUME']        
        except:
            origin = header['ORIGIN']
            if origin == 'GILDAS Consortium':
                self.instrument = 'GREAT'
            elif origin[0:6] == 'Miriad':
                print('Origin is ', origin)
                self.instrument = 'HI'
        try:
            telescope = self.header['TELESCOP']
            if telescope == 'ALMA':
                self.instrument = 'ALMA'
            elif telescope == 'IRAM30M':
                self.instrument = 'IRAM'
            print('Telescope is ', telescope)
        except:
            print('Unknown telescope')
        try:
            self.obsdate = header['DATE-OBS']
        except:
            self.obsdate = header['DATE']
        try:
            self.observer = header['OBSERVER']
        except:
            self.observer = ''
        # Reading files
        if self.instrument == 'FIFI-LS':     
            self.readFIFI(hdl)
        elif self.instrument == 'GREAT':
            self.readGREAT(hdl)
        elif self.instrument == 'PACS':
            self.readPACS(hdl)
        elif self.instrument == 'SPIRE':
            self.readSPIRE(hdl)
        elif self.instrument == 'FORCAST':
            self.readFORCAST(hdl)
        elif self.instrument == 'HI':
            self.readHI(hdl)
        elif self.instrument == 'MUSE':
            self.readMUSE(hdl)
        elif self.instrument == 'IRAM':
            self.readIRAM(hdl)
        elif self.instrument == 'SITELLE':
            self.readSITELLE(hdl)
        elif self.instrument == 'PCWI':
            self.readPCWI(hdl)
        elif self.instrument == 'HALPHA':
            self.readHalpha(hdl)
        #elif self.instrument == 'OVRO': # Case of Carma data
        #    se
        else:
            print('This is not a supported spectral cube')
        hdl.close()
        # Index of the ref wavelength
        self.n0 = np.argmin(np.abs(self.wave - self.l0))
        print('ref wavelength at n: ', self.n0)
        # Create a grid of points
        self.nz, self.ny, self.nx = np.shape(self.flux)
        if self.n0 <= 0:
            self.n0 = self.nz // 2
        xi = np.arange(self.nx); yi = np.arange(self.ny)
        xi,yi = np.meshgrid(xi, yi)
        # Compute rotation angle
        h1 = self.wcs.to_header()
        try:
            self.crota2 = np.arctan2(-h1["PC2_1"], h1["PC2_2"]) * 180./np.pi
        except:
            self.crota2 = 0.
        print('rotation angle ', self.crota2)
        # Alternative way
        # self.points = np.array([np.ravel(xi), np.ravel(yi)]).transpose()
        self.points = np.c_[np.ravel(xi), np.ravel(yi)]
        # Time used for reading
        elapsed_time = time.process_time() - t
        print('Reading of cube completed in ', elapsed_time,' s')

    def computeExpFromNan(self):
        """Compute an exposure cube from NaN in the flux cube."""
        if self.instrument == 'GREAT':
            # Blank values (highest value is blank)
            try:
                idx = self.flux > self.header['DATAMAX']
                self.flux[idx] = np.nan
                self.exposure = np.asarray(~idx, np.dfloat32)           
            except:
                print('No data max in the header')
                self.exposure = np.asarray(np.isfinite(self.flux), np.float32)
        else:
            self.exposure = np.asarray(np.isfinite(self.flux), np.float32)
       
    def readFIFI(self, hdl):
        print('This is a FIFI-LS spectral cube')
        self.wcs = WCS(self.header).celestial
        self.crpix3 = self.header['CRPIX3']
        self.crval3 = self.header['CRVAL3']
        self.cdelt3 = self.header['CDELT3']
        self.objname = self.header['OBJ_NAME']
        print('Object is ', self.objname)
        try:
            self.filegpid = self.header['FILEGPID']
        except:
            print('File group ID not defined')
            self.filegpid = 'Unknown'
        try:
            self.baryshift = self.header['BARYSHFT']
        except:
            print('No baryshift defined in the cube !')
            self.baryshift = 0.
        self.pixscale = self.header['PIXSCAL']
        self.resolution = self.header['RESOLUN']
        self.za = (self.header['ZA_START'], self.header['ZA_END'])
        self.altitude = (self.header['ALTI_STA'],self.header['ALTI_END'])
        try:
            self.redshift = self.header['REDSHIFT']
        except:
            print('No redshift present')
            self.redshift = 0.0
        self.flux = hdl['FLUX'].data
        self.eflux = hdl['ERROR'].data
        try:
            self.uflux = hdl['UNCORRECTED_FLUX'].data
            self.euflux = hdl['UNCORRECTED_ERROR'].data
        except:
            self.uflux = self.flux.copy()
            self.euflux = self.eflux.copy()
            print('No uncorrected flux defined !')
        self.wave = hdl['WAVELENGTH'].data
        self.n = len(self.wave)
        #print('Size of wave ', self.n)
        #self.vel = np.zeros(self.n)  # prepare array of velocities
        self.x = hdl['X'].data
        #print('X  read')
        self.y = hdl['Y'].data
        #print('Y  read')
        self.channel = self.header['DETCHAN']
        #print('channel ', self.channel)
        if self.channel == 'BLUE':
            self.order = self.header["G_ORD_B"]
        else:
            self.order = '1'
        #print('Channel order ', self.channel, self.order)
        # Read reference wavelength from file group name
        try:
            self.l0 = self.header['RESTWAV']
            print('Restwave is ', self.l0)
        except:
            try:
                if self.channel == 'RED':
                    filegp = self.header['FILEGP_R']
                else:
                    filegp = self.header['FILEGP_B']
                print('file group id ', filegp)
                names = filegp.split('_')
                wr = names[-1]
                print('wr is ',wr)
                self.l0 = float(wr)
            except:
                self.l0 = np.nanmedian(self.wave)
                print('No file group present, assumed central wavelength ', self.l0)
        #print('min wave ', np.nanmin(self.wave))
        try:
            utran = hdl['UNSMOOTHED_TRANSMISSION'].data
            w = utran[0,:]
            t = utran[1,:]
            idx = (w > np.nanmin(self.wave)) & (w < np.nanmax(self.wave))
            print('Atran ', len(idx), len(w))
            self.watran = w[idx]
            self.uatran = t[idx]
            #self.atran = np.interp(self.wave,w,t)  # Interpolation at the resolution of the wavelength grid
        except:
            print('The unsmoothed transmission is not available')
            self.watran = None
            self.uatran = None
        self.atran = hdl['TRANSMISSION'].data
        self.response = hdl['RESPONSE'].data
        nexp = self.header['NEXP']
        exptime = self.header['EXPTIME']
        # Exposure contains number of exposures - if all the exposure last the same time, the 
        # following is correct, otherwise it is an approximation
        self.exposure = hdl['EXPOSURE_MAP'].data.astype(float) * exptime/nexp
          
    def readGREAT(self, hdl):
        #from scipy.special import erf
        print('This is a GREAT spectral cube')
        #self.cmin = header['DATAMIN']
        #self.cmax = header['DATAMAX']
        self.wcs = WCS(self.header).celestial
        self.crpix3 = self.header['CRPIX3']
        self.crval3 = self.header['CRVAL3']
        self.cdelt3 = self.header['CDELT3']
        self.objname = self.header['OBJECT']
        self.redshift = self.header['VELO-LSR'] # in m/s
        c = 299792458.0  # speed of light in m/s 
        self.redshift /= c
        self.pixscale, ypixscale = proj_plane_pixel_scales(self.wcs) * 3600. # Pixel scale in arcsec
        self.n = self.header['NAXIS3']
        naxes = self.header['NAXIS']
        if naxes == 4:
            self.flux = (hdl['PRIMARY'].data)[0,:,:,:]
        else:
            self.flux = hdl['PRIMARY'].data
        # Generate the exposure mask    
        idx = (self.flux > self.header['DATAMAX']) | (self.flux < self.header['DATAMIN'])
        self.flux[idx] = np.nan
        self.exposure = np.asarray(~idx, np.dfloat32)           
        #
        self.bunit = self.header['BUNIT']
        eta_fss=0.97
        eta_mb =0.67
        calib = 971. # approx calibration
        if self.bunit == 'K (Tmb)':
            self.Tb2Jy = calib*eta_fss*eta_mb
        else: # Case of K (Ta*)
            self.Tb2Jy = calib
        nu0 = self.header['RESTFREQ']  # MHz
        l0 = c/nu0  # in micron
        # Compute calibration at reference frequency
        calibnu = 2 * 1.38e-16 * (nu0/(c*100))**2
        print('calibration at nu ', calibnu)
        # Transform in Jy/pixel
        # Compute the beam size at the wavelength
        bmaj = self.header['BMAJ'] * 3600. # Beam major axis in arcsec
        bmin = self.header['BMIN'] * 3600. # Beam minor axis in arcsec
        # Multiply by the flux fraction in the pixel assuming a 2D Gaussian curve                    
        #pixfraction = 0.5 * erf(self.pixscale*0.5/bmaj) * erf(ypixscale*0.5/bmin)
        #print('Beam fraction on pixel ', pixfraction)
        self.npix_per_beam = 1.331 * bmaj * bmin / (self.pixscale * ypixscale)
        # Flux in K /beam. So divide by the number of pixels per beam to have K/pix
        self.flux /= self.npix_per_beam
        print('Tb2Jy ', self.Tb2Jy)
        #self.Tb2Jy *= pixfraction
        # Multiplication is delayed in the code to speed up reading
        # self.flux *= self.Tb2Jy   # Transformed from temperature to S_nu [Jy] per pixel
        pix0 = 0
        vel = self.cdelt3 * (np.arange(self.n) - self.crpix3 + pix0) + self.crval3
        self.l0 = l0
        self.wave = l0 + l0 * vel / c
        
    def readFORCAST(self, hdl): 
        print('This is a FORCAST spectral cube - astropy')
        wcs = WCS(self.header)
        self.wcs = wcs.celestial
        self.crpix3 = self.header['CRPIX3']
        self.crval3 = self.header['CRVAL3']
        self.cdelt3 = self.header['CDELT3']
        self.objname = self.header['OBJECT']
        self.redshift = 0 # in m/s
        self.pixscale, ypixscale = proj_plane_pixel_scales(self.wcs) * 3600. # Pixel scale in arcsec
        self.n = self.header['NAXIS3']
        self.flux = hdl['FLUX'].data
        self.eflux = hdl['ERROR'].data
        exptime = self.header['EXPTIME']
        exp = hdl['EXPOSURE'].data.astype(float) * exptime
        self.exposure = np.broadcast_to(exp, np.shape(self.flux))
        pix0=0
        self.wave = self.cdelt3 * (np.arange(self.n) - self.crpix3 + pix0) + self.crval3
        self.l0 = np.nanmedian(self.wave)
        # self.watran = hdl['WAVEPOS'].data
        # computing watran
        ha = hdl['TRANSMISSION'].header
        self.watran = ha['CDELT1'] * (np.arange(ha['NAXIS1'])-ha['CRPIX1']) + ha['CRVAL1']
        self.uatran = hdl['TRANSMISSION'].data
        #self.baryshift = self.header['WAVSHIFT']*self.header['CDELT3']#/self.l0
        self.baryshift = computeBaryshiftAstropy(self.header)
        c = 299792.458
        print('Barycentric shift [km/s]: ', self.baryshift * c)

    def readPACS(self, hdl):
        """ Case of PACS spectral cubes """
        print('This is a PACS spectral cube')
        self.objname = self.header['OBJECT']
        try:
            self.redshift = self.header['REDSHFTV']*1000. # in km/s
            c = 299792458.0  # speed of light in m/s 
            self.redshift /= c
        except:
            self.redshift = 0.
        print('Object is ',self.objname)
        self.flux = hdl['image'].data
        self.eflux = hdl['error'].data
        try:
            self.exposure = hdl['coverage'].data
            print('Coverage read')
        except:
            print('No coverage available - range observation')
            self.computeExpFromNan()
            print('New exposure computed')
            
        wave = hdl['wcs-tab'].data
        nwave = len(np.shape(wave['wavelen']))
        if nwave == 3:
            self.wave = np.concatenate(wave['wavelen'][0])
        else:
            self.wave = np.concatenate(wave['wavelen'])
        self.l0 = np.nanmedian(self.wave)
        self.n = len(self.wave)
        print('Length of wavelength ',self.n)
        print('Min and max wavelengths: ', np.nanmin(self.wave), np.nanmax(self.wave))
        print('wavelenght shape ', np.shape(self.wave))
        header = hdl['IMAGE'].header
        self.header = header
        # print('Header ',header)
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
        #print('astrometry ', self.wcs)
        self.pixscale, ypixscale = proj_plane_pixel_scales(self.wcs) * 3600. # Pixel scale in arcsec
        self.crpix3 = 0
        w = self.wave
        self.crval3 = w[0]
        self.cdelt3 = np.median(w[1:] - w[:-1])

    def readSPIRE(self, hdl):
        """ Case of SPIRE spectral cubes """
        print('This is a SPIRE spectral cube')
        self.objname = self.header['OBJECT']
        #try:
        #    self.redshift = self.header['REDSHFTV']*1000. # in km/s
        #    c = 299792458.0  # speed of light in m/s 
        #    self.redshift /= c
        #except:
        # Does not exist in original cube, but can be added later
        try:
            self.redshift = self.header['REDSHIFT']
        except:
            self.redshift = 0.
        print('Object is ',self.objname)
        self.flux = hdl['image'].data 
        self.eflux = hdl['error'].data
        self.exposure = hdl['coverage'].data
            
        nz,ny,nx = np.shape(self.flux)
        imah = hdl['image'].header
        crpix3 = imah['CRPIX3']
        crval3 = imah['CRVAL3']
        cdelt3 = imah['CDELT3']
        pix0 = 0
        self.n = nz
        frequency = cdelt3 * (np.arange(self.n) - crpix3 + pix0) + crval3 # Hz
        self.wave = 299792.458 / frequency
        self.l0 = np.nanmedian(self.wave)
        #self.n = len(self.wave)
        #print('Length of wavelength ',self.n)
        print('Min and max wavelengths: ', np.nanmin(self.wave), np.nanmax(self.wave))
        print('wavelenght shape ', np.shape(self.wave))
        #self.header = hdl['image'].header
        self.wcs = WCS(imah).celestial
        self.pixscale, ypixscale = proj_plane_pixel_scales(self.wcs) * 3600. # Pixel scale in arcsec
        pixelstd = (self.pixscale/206264.8)**2 # pixel size in sterad
        self.flux *= 1.e26 * pixelstd  # transformed  in Jy/pixel
        self.eflux *= 1.e26 * pixelstd 
        self.crpix3 = 0
        w = self.wave
        # Flip order to increasing wavelength
        self.wave = self.wave[::-1]
        self.flux = self.flux[::-1,:,:]
        self.eflux = self.eflux[::-1,:,:]
        self.exposure = self.exposure[::-1,:,:]
        self.crval3 = w[0]
        self.cdelt3 = np.median(w[1:] - w[:-1])

        
    def readHI(self, hdl):
        """Case of generic radio cube (in this case HI cubes from Westerbrock)."""
        self.objname = self.header['OBJECT']
        print('Object of HI is ',self.objname)
        self.header = hdl['PRIMARY'].header
        self.redshift = 0.
        data = hdl['PRIMARY'].data
        naxis = self.header['NAXIS']
        if naxis == 3:
            self.flux = data
        elif naxis == 4: # polarization
            self.flux = data[0,:,:,:]
        nz, ny, nx = np.shape(self.flux)
        self.n = nz
        print('nz: ',nz, self.header['NAXIS3'])
        wcs = WCS(self.header)
        self.wcs = wcs.celestial
        self.crpix3 = self.header['CRPIX3']
        self.crval3 = self.header['CRVAL3']
        self.cdelt3 = self.header['CDELT3']
        ctype3 = self.header['CTYPE3']
        if (ctype3 == 'VELO-HEL') or (ctype3 == 'VELO-LSR'):
            pix0 = 0 # it's 1 normally ...
            velocity = self.cdelt3 * (np.arange(self.n) - self.crpix3 + pix0) + self.crval3 # m/s
            # Neutral Hydrogen (HI)
            nu0 = self.header['RESTFREQ']
            print('reference frequency', nu0)
            c = 299792458.0 # m/s
            # self.l0 = 21.1061140542 * 1.e4 #um
            self.l0 = c/nu0 * 1.e6 #um
            self.wave = self.l0 * (1 + velocity/c) #um
        self.pixscale, ypixscale = proj_plane_pixel_scales(self.wcs) * 3600. # Pixel scale in arcsec
        print('scale is ', self.pixscale)
        # Back to wavelength
        w = self.wave
        self.crval3 = w[0]
        self.cdelt3 = np.median(w[1:] - w[:-1])
        
    def readHalpha(self, hdl):
        """Case of Halpha images."""
        self.objname = self.header['OBJECT']
        print('Object of HALPHA is ',self.objname)
        self.header = hdl['PRIMARY'].header
        self.redshift = 0.
        data = hdl['PRIMARY'].data
        naxis = self.header['NAXIS']
        if naxis == 3:
            self.flux = data
        elif naxis == 4: # polarization
            self.flux = data[0,:,:,:]
        nz, ny, nx = np.shape(self.flux)
        self.n = nz
        print('nz: ',nz, self.header['NAXIS3'])
        pix0 = 0
        self.crpix3 = self.header['CRPIX3']
        self.crval3 = self.header['CRVAL3']
        self.cdelt3 = self.header['CDELT3']
        ctype3 = self.header['CTYPE3']
        if ctype3 in ['VELO-HEL','VELO-LSR','VOPT']:
            velocity = self.cdelt3 * (np.arange(self.n) - self.crpix3 + pix0) + self.crval3 # m/s
            nu0 = self.header['RESTFREQ']
            print('reference frequency', nu0)
            c = 299792458.0 # m/s
            self.l0 = c / nu0 * 1.e6 #um
            print('Rest wavelength ', self.l0)
            self.wave = self.l0 * (1 + velocity / c) #um
        self.wcs = WCS(self.header).celestial
        #print('WCS ', self.wcs)
        self.pixscale, ypixscale = proj_plane_pixel_scales(self.wcs) * 3600. # Pixel scale in arcsec
        print('scale is ', self.pixscale)
        # Back to wavelength
        w = self.wave
        self.crval3 = w[0]
        self.cdelt3 = np.median(w[1:] - w[:-1])
        # Rescale flux
        bscale = self.header['BSCALE']
        bzero = self.header['BZERO']
        print('bzero, bscale ', bzero, bscale)
        self.flux = bzero + self.flux * bscale   # Flux in units of erg/cm2/s/A
        # Transform into Jy
        #btype = self.header['BTYPE']
        #if btype == 'erg/cm2/s/A':
        print('Transform into F_lambda [erg/cm2/s/A] into F_nu [Jy] ...')
        c = 299792458.0 # m/s
        wA = w * 1.e4 # Wavelength in A (1 um = 1.e4 A)
        self.flux *=  3.33564095e4 * wA**2
        #self.flux = (w * 1.e4) * self.flux  * (w * 1.e-6)/c * 1.e23 # F_nu

        
    def readIRAM(self, hdl):
        """Case of generic radio cube (in this case HI cubes from Westerbrock)."""
        #self.objname = self.header['OBJECT']
        #print('Object of HI is ',self.objname)
        print('Reading IRAM cube')
        self.objname = ''
        self.header = hdl['PRIMARY'].header
        self.redshift = 0.
        data = hdl['PRIMARY'].data
        naxis = self.header['NAXIS']
        if naxis == 3:
            self.flux = data
        elif naxis == 4: # polarization
            self.flux = data[0,:,:,:]
        nz, ny, nx = np.shape(self.flux)
        # Flux in HERACLES data is expressed in T_MB (Main Beam)
        # From Leroy et al. 2009, and the header values we know that  T_MB = F_eff/B_eff T_A = 0.92/0.58 T_A = 1.59 T_A,
        # In reality, the calibration paper of Kramer
        # (https://safe.nrao.edu/wiki/pub/KPAF/KfpaPipelineReview/kramer_1997_cali_rep.pdf)
        # show in the Table at page 17 that this is a mean value.
        # Once converted the T_MB in T_A, we can convert to Flux in Jy using the Kramer's formula:
        # S_nu/T_A = 3.906 F_eff/eta_A [Jy/K]
        # S_nu = 3.906 B_eff/eta_A * T_MB = 3.906 * .58 / 0.32 * T_MB = 7.08 T_MB, using the header value
        # eta_A from the table at page 17.
        # Using eta_A from https://web-archives.iram.fr/IRAMFR/ARN/aug05/node6.html for 2005, and interpolating at 230Hz
        # we have eta_A = 0.41 --> S_nu = 3.906 * 0.58/ 0.41 * T_MB = 5.52
        # If I use the data at page 17 at 230 GHz, this means:
        # S_nu [Jy] = 3.906 0.86/0.32 T_A[K] = 3.906 0.86/0.32 *T_MB[K]/1.59 = 6.6021 * T_MB[K]
        # The result is in Jy/beam
        self.Tmb2Jy = 5.52
        self.flux *= self.Tmb2Jy
        
        self.n = nz
        print('nz: ',nz, self.header['NAXIS3'])
        wcs = WCS(self.header)
        self.wcs = wcs.celestial
        self.crpix3 = self.header['CRPIX3']
        self.crval3 = self.header['CRVAL3']
        self.cdelt3 = self.header['CDELT3']
        ctype3 = self.header['CTYPE3']
        pix0 = 0
        if ctype3 == 'VELOCITY':
            velocity = self.cdelt3 * (np.arange(self.n) - self.crpix3 + pix0) + self.crval3 # m/s
            # Neutral Hydrogen (HI)
            nu0 = self.header['RESTFREQ']
            print('reference frequency', nu0)
            c = 299792458.0 # m/s
            # self.l0 = 21.1061140542 * 1.e4 #um
            self.l0 = c/nu0 * 1.e6 #um
            #self.wave = self.l0 * (1 + velocity/c) #um
            # Radio convention: https://science.nrao.edu/facilities/vla/docs/manuals/obsguide/modes/line
            self.wave = self.l0 / (1 - velocity/c) #um  
        self.pixscale, ypixscale = proj_plane_pixel_scales(self.wcs) * 3600. # Pixel scale in arcsec
        print('scale is ', self.pixscale)

        # Now transform into Jy/pix
        # Compute the beam size at the wavelength
        bmaj = self.header['BMAJ'] * 3600. # Beam major axis in arcsec
        bmin = self.header['BMIN'] * 3600. # Beam minor axis in arcsec
        # Multiply by the flux fraction in the pixel assuming a 2D Gaussian curve                    
        #pixfraction = 0.5 * erf(self.pixscale*0.5/bmaj) * erf(ypixscale*0.5/bmin)
        #print('Beam fraction on pixel ', pixfraction)
        self.npix_per_beam = 1.331 * bmaj * bmin / (self.pixscale * ypixscale)
        # Flux in Jy /beam. So divide by the number of pixels per beam to have Jy/pix
        self.flux /= self.npix_per_beam

        # Back to wavelength
        w = self.wave
        self.crval3 = w[0]
        self.cdelt3 = np.median(w[1:] - w[:-1])
        
    def readMUSE(self, hdl):
        """MUSE integral field spectrometer at VLT"""
        self.objname = self.header['OBJECT']
        print('Object of MUSE is ', self.objname)
        self.flux = hdl['DATA'].data  # 10**(-20)*erg/s/cm**2/Angstrom
        nz, ny, nx = np.shape(self.flux)
        self.n = nz
        self.header = hdl['DATA'].header
        wcs = WCS(self.header)
        self.wcs = wcs.celestial
        self.crpix3 = self.header['CRPIX3']
        self.crval3 = self.header['CRVAL3']
        self.cdelt3 = self.header['CD3_3']
        #ctype3 = self.header['CTYPE3'].strip()
        pix0=0
        self.wave = self.cdelt3 * (np.arange(self.n) - self.crpix3 + pix0) + self.crval3 # Angstrom
        self.wave *= 1.e-4  # um
        #print('min wave ', np.nanmin(self.wave))
        self.pixscale, ypixscale = proj_plane_pixel_scales(self.wcs) * 3600. # Pixel scale in arcsec
        print('scale is ', self.pixscale)
        self.redshift = 0
        self.l0 = np.nanmedian(self.wave)
        #self.l0 = (self.header['WAVELMIN']+self.header['WAVELMAX']) * 0.5 * 1.e-3 # wav in nm
    
    def readSITELLE(self, hdl):
        """SITELLE integral field spectrometer at CFH"""
        self.objname = self.header['OBJECT']
        print('Object of SITELLE is ', self.objname)
        self.flux = hdl['PRIMARY'].data  # 10**(-17)*erg/s/cm**2/Angstrom
        nz, ny, nx = np.shape(self.flux)
        self.n = nz
        self.header = hdl['PRIMARY'].header
        wcs = WCS(self.header, naxis=2)
        self.wcs = wcs.celestial
        self.crpix3 = self.header['CRPIX3']
        self.crval3 = self.header['CRVAL3']
        self.cdelt3 = self.header['CDELT3']
        self.wave = 1. / (self.crval3 + self.cdelt3 * np.arange(self.n)) * 1.e4 #um
        self.pixscale, ypixscale = proj_plane_pixel_scales(self.wcs) * 3600. # Pixel scale in arcsec
        print('scale is ', self.pixscale)
        self.redshift = 0
        self.l0 = np.nanmedian(self.wave)
    
    def readPCWI(self, hdl):
        """
        Palomar Cosmic Web Imager (http://www.srl.caltech.edu/sal/cosmic-web-imager.html)
        """
        self.objname = self.header['OBJECT']
        print('Object of PCWI is ', self.objname)
        self.flux = hdl[0].data  # 10**(-20)*erg/s/cm**2/Angstrom
        nz, ny, nx = np.shape(self.flux)
        self.n = nz
        self.header = hdl[0].header
        wcs = WCS(self.header)
        self.wcs = wcs.celestial
        self.crpix3 = self.header['CRPIX3']
        self.crval3 = self.header['CRVAL3']
        self.cdelt3 = self.header['CD3_3']
        pix0=0
        self.wave = self.cdelt3 * (np.arange(self.n) + pix0 - self.crpix3) + self.crval3 # Angstrom
        self.wave *= 1.e-4  # um
        self.pixscale, self.ypixscale = proj_plane_pixel_scales(self.wcs) * 3600. # Pixel scale in arcsec
        print('scale is ', self.pixscale, self.ypixscale)
        try:
            self.redshift = self.header['REDSHIFT']
        except:
            self.redshift = 0
        try:
            self.l0 = self.header['RESTWAV']
        except:
            self.l0 = self.header['WAVMID']*1.e-4 
    
    def getResolutionFIFI(self):
        """Compute resolution at reference wavelength for FIFI-LS"""
        if self.instrument != 'FIFI-LS':
            return
        l0 = self.l0
        z = self.redshift
        l = l0 * (1.+z)
        if self.channel == 'RED':
            return 0.062423 * l * l - 6.6595 * l + 647.65
        else:
            if self.order == '1':
                return 0.16864 * l * l - 22.831 * l + 1316.6 
            else:
                return 1.9163 * l * l - 187.35 * l + 5496.9

class specCube(object):
    """ spectral cube - read with fitsio routines"""
    def __init__(self, infile):
        import time
        try:
            import fitsio
        except ImportError:
            print('install fitsio library:  conda install -c conda-forge fitsio')
            exit
        t = time.process_time()
        # Option None seems faster than False
        #hdl = fits.open(infile,memmap=None)
        hdl = fitsio.FITS(infile)
        print('File opened ')
        #header = hdl['PRIMARY'].header
        self.header = hdl[0].read_header()
        self.purifyHeader()
        self.filename = infile
        # Change internal file name to external file name
        self.header['FILENAME']=self.filename
        print('Using FITSIO ...')
        try:
            self.instrument = self.header['INSTRUME'].strip()
            print('Instrument: ', self.instrument)
            if self.instrument == 'OVRO MMA':
                self.telescope = 'OVRO'
                self.instrument = 'MMA'
        except:
            try:
                origin = self.header['ORIGIN'].strip()
                if origin == 'GILDAS Consortium':
                    self.instrument = 'GREAT'
                elif origin[0:6] == 'Miriad':
                    print('Origin is ', origin)
                    self.instrument = 'HI'
            except:
                print('Unknown origin')
            try:
                telescope = self.header['TELESCOP'].strip()
                if telescope == 'ALMA':
                    self.instrument = 'ALMA'
                elif telescope == 'IRAM30M':
                    self.instrument = 'IRAM'
                elif telescope == 'OVRO':
                    self.instrument = 'CARMA'
                print('telescope is ', telescope)
            except:
                print('Unknown telescope')
        try:
            self.obsdate = self.header['DATE-OBS'].strip()
        except:
            self.obsdate = self.header['DATE'].strip()
        try:
            self.observer = self.header['OBSERVER'].strip()
        except:
            self.observer = ''
        # Reading files
        if self.instrument == 'FIFI-LS':     
            self.readFIFI(hdl)
        elif self.instrument == 'GREAT':
            self.readGREAT(hdl)
        elif self.instrument == 'PACS':
            self.readPACS(hdl)
        elif self.instrument == 'FORCAST':
            self.readFORCAST(hdl)
        elif self.instrument == 'HI':
            self.readHI(hdl)
        elif self.instrument == 'IRAM':
            self.readIRAM(hdl)
        elif self.instrument in ['VLA','ALMA','CARMA','MMA']:
            self.readVLA(hdl)
        elif self.instrument == 'MUSE':
            self.readMUSE(hdl)
        elif self.instrument == 'HALPHA':
            self.readHalpha(hdl)
        else:
            print('This is not a supported spectral cube')
        hdl.close()
        # Index of the ref wavelength
        self.n0 = np.argmin(np.abs(self.wave - self.l0))
        print('ref wavelength at n: ', self.n0)
        # Create a grid of points
        self.nz, self.ny, self.nx = np.shape(self.flux)
        if self.n0 <= 0:
            self.n0 = self.nz // 2
        xi = np.arange(self.nx); yi = np.arange(self.ny)
        xi, yi = np.meshgrid(xi, yi)
        # Compute rotation angle
        h1 = self.wcs.to_header()
        try:
            self.crota2 = np.arctan2(-h1["PC2_1"], h1["PC2_2"]) * 180./np.pi
        except:
            self.crota2 = 0.
        print('rotation angle ', self.crota2)
        # Alternative way
        # self.points = np.array([np.ravel(xi), np.ravel(yi)]).transpose()
        self.points = np.c_[np.ravel(xi), np.ravel(yi)]
        # Time used for reading
        elapsed_time = time.process_time() - t
        print('Reading of cube completed in ', elapsed_time,' s')
        
    def purifyHeader(self):
        orig_header = self.header
        header = fits.Header()
        for dict_key in orig_header.keys():
            try:
                header[dict_key] = orig_header[dict_key]
            except:
                pass
        self.header = header

    def computeExpFromNan(self):
        """Compute an exposure cube from NaN in the flux cube."""
        if self.instrument == 'GREAT':
            # Blank values (highest value is blank)
            try:
                idx = self.flux > self.header['DATAMAX']
                #blank = self.header['BZERO']+self.header['BSCALE']*self.header['BLANK']
                #idx = self.flux == blank
                #self.flux[idx] = np.nan
                self.exposure = ~idx            
            except:
                print('No data max in the header')
                self.exposure = np.isfinite(self.flux)
        else:
            self.exposure = np.isfinite(self.flux)
            print('exp map is ', np.shape(self.exposure))
        
    def readFIFI(self, hdl):
        print('This is a FIFI-LS spectral cube (read with FITSIO')
        #self.header.delete('ASSC_AOR')  # this cannot be interpreted by WCS
        self.wcs = WCS(self.header).celestial
        self.crpix3 = self.header['CRPIX3']
        self.crval3 = self.header['CRVAL3']
        self.cdelt3 = self.header['CDELT3']
        self.bunit = self.header['BUNIT']
        self.objname = self.header['OBJ_NAME'].strip()
        print('Name :', self.objname)
        try:
            self.filegpid = self.header['FILEGPID'].strip()
        except:
            print('File group ID not defined')
            self.filegpid = 'Unknown'
        self.baryshift = self.header['BARYSHFT']
        self.pixscale = self.header['PIXSCAL']
        self.resolution = self.header['RESOLUN']
        self.za = (self.header['ZA_START'], self.header['ZA_END'])
        self.altitude = (self.header['ALTI_STA'],self.header['ALTI_END'])
        try:
            self.redshift = self.header['REDSHIFT']
            print('Redshift ', self.redshift)
        except:
            print('No redshift present')
            self.redshift = 0.0
            
        extnames = [f.get_extname() for f in hdl]
        print('ext names ', extnames)
        
        self.flux = hdl[extnames.index('FLUX')].read()
        self.eflux = hdl[extnames.index('ERROR')].read()
        self.uflux = hdl[extnames.index('UNCORRECTED_FLUX')].read()
        #self.uflux = hdl[extnames.index('ERROR')].read()
        self.euflux = hdl[extnames.index('UNCORRECTED_ERROR')].read()
        self.wave = hdl[extnames.index('WAVELENGTH')].read()
        self.n = len(self.wave)
        #self.vel = np.zeros(self.n)  # prepare array of velocities
        self.x = hdl[extnames.index('X')].read()
        self.y = hdl[extnames.index('Y')].read()
        #self.atran = hdl['TRANSMISSION'].data
        self.channel = self.header['DETCHAN'].strip()
        if self.channel == 'BLUE':
            self.order = str(self.header["G_ORD_B"])
        else:
            self.order = '1'
        print('channel ', self.channel, ' order ', self.order)
        # Read reference wavelength from file group name
        try:
            self.l0 = self.header['RESTWAV']
            print('Rest wavelength is ', self.l0)
        except:
            try:
                if self.channel == 'RED':
                    filegp = self.header['FILEGP_R'].strip()
                else:
                    filegp = self.header['FILEGP_B'].strip()
                print('file group id ', filegp)
                names = filegp.split('_')
                wr = names[-1]
                print('wr is ',wr)
                self.l0 = float(wr)
            except:
                self.l0 = np.nanmedian(self.wave)
                print('No file group present, assumed central wavelength ', self.l0)
        print('ref wav ', self.l0)
        
        try:
            utran = hdl['UNSMOOTHED_TRANSMISSION'].read()
            w = utran[0,:]
            t = utran[1,:]
            idx = (w > np.nanmin(self.wave)) & (w < np.nanmax(self.wave))
            print('Atran ', len(idx), len(w))
            self.watran = w[idx]
            self.uatran = t[idx]
            #self.atran = np.interp(self.wave,w,t)  # Interpolation at the resolution of the wavelength grid
        except:
            print('The unsmoothed transmission is not available')
            self.watran = None
            self.uatran = None
        self.atran = hdl['TRANSMISSION'].read()
        #try:
        #    utran = hdl[extnames.index('UNSMOOTHED_TRANSMISSION')].read()
        #    w = utran[0,:]
        #    t = utran[1,:]
        #    self.atran = np.interp(self.wave,w,t)  # Interpolation at the resolution of the wavelength grid
        #except:
        #    print('The unsmoothed transmission is not available')
        #    self.atran = hdl[extnames.index('TRANSMISSION')].read()
        #print('atran read')
        self.response = hdl[extnames.index('RESPONSE')].read()
        nexp = self.header['NEXP']
        exptime = self.header['EXPTIME']
        # Exposure contains number of exposures - if all the exposure last the same time, the 
        # following is correct, otherwise it is an approximation
        self.exposure = hdl[extnames.index('EXPOSURE_MAP')].read().astype(float) * exptime/nexp
        # Baryshift
        #self.baryshift = self.computeBaryshift()
        #c = 299792.458 
        #print('Baryshift ', self.baryshift*c, computeBaryshift(self.header)*c,computeBaryshiftAstropy(self.header)*c)
        try:
            self.baryshift = computeBaryshiftAstropy(self.header)
        except:
            print('Baryshift from header - includes LSR velocity correction')
            pass
          
    def readGREAT(self, hdl):
        from scipy.special import erf
        import time
        print('This is a GREAT spectral cube')
        #self.cmin = header['DATAMIN']
        #self.cmax = header['DATAMAX']
        t = time.process_time()
        self.wcs = WCS(self.header).celestial
        self.crpix3 = self.header['CRPIX3']
        self.crval3 = self.header['CRVAL3']
        self.cdelt3 = self.header['CDELT3']
        self.objname = self.header['OBJECT'].strip()
        self.redshift = self.header['VELO-LSR'] # in m/s
        print('redshift ', self.redshift)
        c = 299792458.0  # speed of light in m/s 
        self.redshift /= c
        self.pixscale, ypixscale = proj_plane_pixel_scales(self.wcs) * 3600. # Pixel scale in arcsec
        self.n = self.header['NAXIS3']
        self.bunit = self.header['BUNIT'].strip()
        eta_fss=0.97
        eta_mb =0.67
        calib = 971.
        nu0 = self.header['RESTFREQ']  # MHz
        # Compute trasformation from K -> Jy using definition of Ta
        # Rayleigh-Jeans approximation is appropriate in GREAT range
        calibnu = 2 * 1.38e-16 * (nu0*1.e6/(c*100))**2
        print('calibration at nu ', calibnu)
        print('rest freq ', nu0)
        if self.bunit == 'K (Tmb)':
            self.Tb2Jy = calib*eta_fss*eta_mb
        else: # Case of K (Ta*)
            self.Tb2Jy = calib
        l0 = c/nu0  # in micron
        # Transform in Jy/pixel
        # Compute the beam size at the wavelength
        bmaj = self.header['BMAJ'] * 3600. # Beam major axis in arcsec
        bmin = self.header['BMIN'] * 3600. # Beam minor axis in arcsec
        xsigma = 2 * bmaj / 2.355 / np.sqrt(2) # sigma from FWHM divided by sqrt(2)
        ysigma = 2 * bmin / 2.355 / np.sqrt(2)
        ## Multiply by the flux fraction in the pixel assuming a 2D Gaussian curve                    
        pixfraction = erf(self.pixscale*0.5/xsigma) * erf(ypixscale*0.5/ysigma)
        print('Beam fraction on pixel ', pixfraction)
        # Transform Jy/beam to Jy/pixel
        self.npix_per_beam = 1.331 * bmaj * bmin / (self.pixscale * ypixscale)
        print('Pixels per beam', self.npix_per_beam)
        print('Tb2Jy ', self.Tb2Jy)
        #self.Tb2Jy *= pixfraction
        naxes = self.header['NAXIS']
        print('no of axes ', naxes)
        if naxes == 4:
            self.flux = (hdl[0].read())[0,:,:,:]
        else:
            self.flux = hdl[0].read()
        # Flux in K, non in K/beam - so  divided by beam and multiply by pix area
        # Generate the exposure mask [remove data exceeding max and min values]  
        #print('Datamax', self.header['DATAMAX'])
        #print('Datamin', self.header['DATAMIN'])
        #idx = (self.flux > self.header['DATAMAX']) | (self.flux < self.header['DATAMIN'])
        #self.flux[idx] = np.nan
        #self.exposure = np.asarray(~idx, np.dfloat32)           
        #
        self.flux /= self.npix_per_beam
        t1 = time.process_time()
        print('Flux reading completed in ', t1-t,' s')
        #self.flux *= self.Tb2Jy   # Transformed from temperature to S_nu [Jy]  Move this to code
        pix0 = 0
        vel = self.cdelt3 * (np.arange(self.n) - self.crpix3 + pix0) + self.crval3
        self.l0 = l0
        print('Reference wavelength ', self.l0)
        self.wave = l0 * (1 + vel / c)
        #t2 = time.process_time()
        #print('GREAT transf completed in ', t2-t1,' s')
        
    def readFORCAST(self, hdl): 
        print('This is a FORCAST spectral cube - fitsio')
        wcs = WCS(self.header)
        self.wcs = wcs.celestial
        self.crpix3 = self.header['CRPIX3']
        self.crval3 = self.header['CRVAL3']
        self.cdelt3 = self.header['CDELT3']
        self.objname = self.header['OBJECT'].strip()
        self.redshift = 0 # in m/s
        self.pixscale, ypixscale = proj_plane_pixel_scales(self.wcs) * 3600. # Pixel scale in arcsec
        self.n = self.header['NAXIS3']
        extnames = [f.get_extname() for f in hdl]
        self.flux = hdl[extnames.index('FLUX')].read()
        try:
            self.eflux = hdl[extnames.index('ERROR')].read()
        except:
            self.eflux = np.sqrt(hdl[extnames.index('VARIANCE')].read())
        exptime = self.header['EXPTIME']
        exp = hdl[extnames.index('EXPOSURE')].read().astype(float) * exptime
        self.exposure = np.broadcast_to(exp, np.shape(self.flux))
        pix0 = 0 
        self.wave = self.cdelt3 * (np.arange(self.n) - self.crpix3 + pix0) + self.crval3
        self.l0 = np.nanmedian(self.wave)
        #self.watran = hdl['WAVEPOS'].read()
        ha = hdl['TRANSMISSION'].read_header()
        self.watran = ha['CDELT1'] * (np.arange(ha['NAXIS1']) - ha['CRPIX1']) + ha['CRVAL1']
        self.uatran = hdl['TRANSMISSION'].read()
        #self.baryshift = self.header['WAVSHIFT']*self.header['CDELT3']#/self.l0
        self.baryshift = computeBaryshiftAstropy(self.header)
        c = 299792.458
        print('Barycentric shift [km/s]: ', self.baryshift * c)
    
    def readPACS(self, hdl):
        """ Case of PACS spectral cubes """
        print('This is a PACS spectral cube')
        extnames = [f.get_extname() for f in hdl]
        print('Extensions ', extnames)
        self.objname = self.header['OBJECT'].strip()
        try:
            self.redshift = self.header['REDSHFTV']*1000. # in km/s
            c = 299792458.0  # speed of light in m/s 
            self.redshift /= c
        except:
            self.redshift = 0.
        print('Object for PACS is ',self.objname)
        self.flux = hdl[extnames.index('image')].read()
        print('Flux read')
        self.eflux = hdl[extnames.index('error')].read()
        print('eflux read')
        try:
            self.exposure = hdl[extnames.index('coverage')].read()
            print('Coverage read')
        except:
            print('No coverage available - range observation')
            self.computeExpFromNan()
            print('New exposure computed')
            print('shape of exp ', np.shape(self.exposure))
            
        wave = hdl[extnames.index('wcs-tab')].read()
        print('Wvl read')
        nwave = len(np.shape(wave['wavelen']))
        if nwave == 3:
            self.wave = np.concatenate(wave['wavelen'][0])
        else:
            self.wave = np.concatenate(wave['wavelen'])
        self.l0 = np.nanmedian(self.wave)
        self.n = len(self.wave)
        print('Length of wavelength ',self.n)
        #print('Min and max wavelengths: ', np.nanmin(self.wave), np.nanmax(self.wave))
        #print(np.shape(self.wave))
        #print('image ext ', extnames.index('image'))
        header = hdl[extnames.index('image')].read_header()
        self.header = header
        #self.purifyHeader()
        #self.wcs = WCS(self.header)
        self.getWCS()
        print('astrometry ', self.wcs)
        self.pixscale, ypixscale = proj_plane_pixel_scales(self.wcs) * 3600. # Pixel scale in arcsec
        self.crpix3 = 1
        w = self.wave
        self.crval3 = w[0]
        self.cdelt3 = np.median(w[1:] - w[:-1])
        
    def readHI(self, hdl):
        """Case of generic radio cube (in this case HI cubes from Westerbrock)."""
        self.objname = self.header['OBJECT'].strip()
        print('Object of HI is ',self.objname)
        #self.header = hdl[0].read_header()
        self.redshift = 0.
        data = hdl[0].read()
        naxis = self.header['NAXIS']
        if naxis == 3:
            self.flux = data
        elif naxis == 4: # polarization
            self.flux = data[0,:,:,:]
        nz, ny, nx = np.shape(self.flux)
        self.n = nz
        print('nz: ',nz, self.header['NAXIS3'])
        wcs = WCS(self.header)
        self.wcs = wcs.celestial
        self.crpix3 = self.header['CRPIX3']
        self.crval3 = self.header['CRVAL3']
        self.cdelt3 = self.header['CDELT3']
        ctype3 = self.header['CTYPE3'].strip()
        pix0 = 0
        if (ctype3 == 'VELO-HEL') or (ctype3 == 'VELO-LSR'):
            velocity = self.cdelt3 * (np.arange(self.n) - self.crpix3 + pix0) + self.crval3 # m/s
            # Neutral Hydrogen (HI)
            nu0 = self.header['RESTFREQ']
            print('reference frequency', nu0)
            c = 299792458.0 # m/s
            # self.l0 = 21.1061140542 * 1.e4 #um
            self.l0 = c/nu0 * 1.e6 #um
            self.wave = self.l0 * (1 + velocity/c) #um
        self.pixscale, ypixscale = proj_plane_pixel_scales(self.wcs) * 3600. # Pixel scale in arcsec
        print('scale is ', self.pixscale)
        # Back to wavelength
        w = self.wave
        self.crval3 = w[0]
        self.cdelt3 = np.median(w[1:] - w[:-1])

    def readHalpha(self, hdl):
        """Case of generic H-alpha cube."""
        self.objname = self.header['OBJECT'].strip()
        print('Object of Halpha is ',self.objname)
        #self.header = hdl[0].read_header()
        self.redshift = 0.
        data = hdl[0].read()
        naxis = self.header['NAXIS']
        if naxis == 3:
            self.flux = data
        elif naxis == 4: # polarization
            self.flux = data[0,:,:,:]
        nz, ny, nx = np.shape(self.flux)
        self.n = nz
        print('nz: ',nz, self.header['NAXIS3'])
        
        self.crpix3 = self.header['CRPIX3']
        self.crval3 = self.header['CRVAL3']
        self.cdelt3 = self.header['CDELT3']
        ctype3 = self.header['CTYPE3'].strip()
        pix0=0
        if ctype3 in ['VELO-HEL', 'VELO-LSR', 'VOPT']:
            velocity = self.cdelt3 * (np.arange(self.n) - self.crpix3 + pix0) + self.crval3 # m/s
            # Neutral Hydrogen (HI)
            nu0 = self.header['RESTFREQ'] # H-alpha
            print('reference frequency', nu0)
            c = 299792458.0 # m/s
            # self.l0 = 21.1061140542 * 1.e4 #um
            self.l0 = c/nu0 * 1.e6 #um
            print('Rest wavelength ', self.l0)
            self.wave = self.l0 * (1 + velocity/c) #um
        self.wcs = WCS(self.header).celestial
        self.pixscale, ypixscale = proj_plane_pixel_scales(self.wcs) * 3600. # Pixel scale in arcsec
        print('scale is ', self.pixscale)
        # Back to wavelength
        w = self.wave
        self.crval3 = w[0]
        self.cdelt3 = np.median(w[1:] - w[:-1])
        # Rescale flux
        bscale = self.header['BSCALE']
        bzero = self.header['BZERO']
        print('bzero, bscale ',bzero, bscale)
        self.flux = bzero + self.flux * bscale   # Flux in units of erg/cm2/s/A
        # Transform into Jy
        #btype = self.header['BTYPE']
        #if btype == 'erg/cm2/s/A':
        print('Transform into Jy ...')
        c = 299792458.0 # m/s
        wA = self.wave * 1.e4
        print(np.nanmean(33356.4095 * wA**2))
        print(np.nanmean(self.flux))
        for f, w in zip(self.flux, wA):
            f *= 33356.4095 * w**2
        #self.flux = (w * 1.e4) * self.flux  * (w * 1.e-6)/c * 1.e-23 # F_nu

    def readIRAM(self, hdl):
        """Case of Heracles observations with IRAM 30M."""
        print('This is an IRAM cube')
        try:
            self.objname = self.header['OBJECT'].strip()
            print('Object of HI is ',self.objname)
        except:
            self.objname = ''
        self.redshift = 0.
        data = hdl[0].read()
        naxis = self.header['NAXIS']
        if naxis == 3:
            self.flux = data
        elif naxis == 4: # polarization
            self.flux = data[0,:,:,:]
        nz, ny, nx = np.shape(self.flux)

        # Flux in HERACLES data is expressed in T_MB (Main Beam)
        # From Leroy et al. 2009, and the header values we know that  T_MB = F_eff/B_eff T_A = 0.92/0.58 T_A = 1.59 T_A,
        # In reality, the calibration paper of Kramer
        # (https://safe.nrao.edu/wiki/pub/KPAF/KfpaPipelineReview/kramer_1997_cali_rep.pdf)
        # show in the Table at page 17 that this is a mean value.
        # Once converted the T_MB in T_A, we can convert to Flux in Jy using the Kramer's formula:
        # S_nu/T_A = 3.906 F_eff/eta_A [Jy/K]
        # S_nu = 3.906 B_eff/eta_A * T_MB = 3.906 * .58 / 0.32 * T_MB = 7.08 T_MB, using the header value
        # eta_A from the table at page 17.
        # Using eta_A from https://web-archives.iram.fr/IRAMFR/ARN/aug05/node6.html for 2005, and interpolating at 230Hz
        # we have eta_A = 0.41 --> S_nu = 3.906 * 0.58/ 0.41 * T_MB = 5.52
        # If I use the data at page 17 at 230 GHz, this means:
        # S_nu [Jy] = 3.906 0.86/0.32 T_A[K] = 3.906 0.86/0.32 *T_MB[K]/1.59 = 6.6021 * T_MB[K]
        # The result is in Jy/beam
        self.Tmb2Jy = 5.52
        self.flux *= self.Tmb2Jy

        self.n = nz
        print('nz: ',nz, self.header['NAXIS3'])
        wcs = WCS(self.header)
        self.wcs = wcs.celestial
        self.crpix3 = self.header['CRPIX3']
        self.crval3 = self.header['CRVAL3']
        self.cdelt3 = self.header['CDELT3']
        ctype3 = self.header['CTYPE3'].strip()
        pix0=0
        if ctype3 == 'VELOCITY':
            velocity = self.cdelt3 * (np.arange(self.n) - self.crpix3 + pix0) + self.crval3 # m/s
            # Neutral Hydrogen (HI)
            nu0 = self.header['RESTFREQ']
            print('reference frequency', nu0)
            c = 299792458.0 # m/s
            # self.l0 = 21.1061140542 * 1.e4 #um
            self.l0 = c/nu0 * 1.e6 #um
            #self.wave = self.l0 * (1 + velocity/c) #um # Optical convention
            # Radio convention: https://science.nrao.edu/facilities/vla/docs/manuals/obsguide/modes/line
            self.wave = self.l0 / (1 - velocity/c) #um  

        self.pixscale, ypixscale = proj_plane_pixel_scales(self.wcs) * 3600. # Pixel scale in arcsec
        print('scale is ', self.pixscale)

        # Now transform into Jy/pix
        # Compute the beam size at the wavelength
        bmaj = self.header['BMAJ'] * 3600. # Beam major axis in arcsec
        bmin = self.header['BMIN'] * 3600. # Beam minor axis in arcsec
        # Multiply by the flux fraction in the pixel assuming a 2D Gaussian curve                    
        #pixfraction = 0.5 * erf(self.pixscale*0.5/bmaj) * erf(ypixscale*0.5/bmin)
        #print('Beam fraction on pixel ', pixfraction)
        self.npix_per_beam = 1.331 * bmaj * bmin / (self.pixscale * ypixscale)
        # Flux in Jy /beam. So divide by the number of pixels per beam to have Jy/pix
        self.flux /= self.npix_per_beam

        # Back to wavelength
        w = self.wave
        self.crval3 = w[0]
        self.cdelt3 = np.median(w[1:] - w[:-1])
        
    def readVLA(self, hdl):
        """Case of VLA cube (from VIVA)."""
        print('Inside read VLA ')
        c = 299792458.0 # m/s
        try:
            self.objname = self.header['OBJECT'].strip()
            print('Object of ', self.instrument, ' is ', self.objname)
        except:
            self.objname = ""
            print('No object name in the header')
        data = hdl[0].read()
        naxis = self.header['NAXIS']
        if naxis == 3:
            self.flux = data
        elif naxis == 4: # polarization
            self.flux = data[0,:,:,:]
        try:
            self.redshift = self.header['REDSHIFT']
            print('Redshift ', self.redshift)
        except:
            self.redshift = 0.
        if self.instrument == 'MMA':
            idx = np.isfinite(self.flux)
            self.flux[~idx] = np.nan            
        nz, ny, nx = np.shape(self.flux)
        self.n = nz
        print('nz: ',nz, self.header['NAXIS3'])
        wcs = WCS(self.header)
        self.wcs = wcs.celestial
        self.crpix3 = self.header['CRPIX3']
        self.crval3 = self.header['CRVAL3']
        self.cdelt3 = self.header['CDELT3']
        ctype3 = self.header['CTYPE3'].strip()
        print('ctype 3 is ', ctype3)
        pix0 =0 
        if ctype3 == 'FREQ':
            if self.instrument in ['VLA','CARMA']:
                nu0 = self.header['RESTFREQ']
            elif self.instrument == 'ALMA':
                nu0 = self.header['RESTFRQ']
            print('reference frequency', nu0, 'Hz')
            # Check if altrval exists and its different from crval3
            try:
                altrval = self.header['ALTRVAL']  # ref velocity in m/s
                vcrval3 = c * (1 - self.crval3/nu0)
                faltrval = nu0 * (1 - altrval/c)  # alt crval3
                if np.abs(vcrval3 - altrval) > 20:
                    self.crval3 = faltrval
                #self.redshift = altrval/c
            except:
                pass
            try:
                altrpix = self.header['ALTRPIX'] # Ref pixel for 3rd dimension
                self.crpix3 = altrpix
            except:
                pass
            freq = self.cdelt3 * (np.arange(self.n) - self.crpix3 + pix0) + self.crval3 
            self.l0 = c/nu0 * 1.e6 #um
            print('Reference wavelength at ', self.l0)
            self.wave = c/freq * 1.e6 #um
            # Reorder wave (and flux)
            idx = np.argsort(self.wave)
            self.wave = self.wave[idx]
            self.flux = self.flux[idx, :, :]
        elif ctype3 in ['VELO-HEL', 'VELO-LSR', 'VRAD','FELO-HEL']:
            velocity = self.cdelt3 * (np.arange(self.n) - self.crpix3 + pix0) + self.crval3 # m/s
            if self.instrument == 'VLA':
                nu0 = self.header['RESTFREQ']
            elif self.instrument in ['ALMA', 'CARMA']:
                nu0 = self.header['RESTFRQ']
            elif self.instrument == 'MMA':
                nu0 = self.header['FREQ0']
            print('reference frequency', nu0)
            c = 299792458.0 # m/s
            # self.l0 = 21.1061140542 * 1.e4 #um
            self.l0 = c/nu0 * 1.e6 #um
            if self.instrument == 'MMA':
                c = 299792.458 # km/s
            try:
                altrval = self.header['ALTRVAL']  # ref frequency
                altrpix = self.header['ALTRPIX']  # ref pixel
                vcrval3 = c * (1. - altrval/nu0)
                #faltrval = nu0 * (1 - altrval/c)  # alt crval3
                #if np.abs(self.crval3 - vcrval3) > 20:
                self.crval3 = vcrval3
                self.crpix3 = altrpix
                velocity = self.cdelt3 * (np.arange(self.n) - self.crpix3 + pix0) + self.crval3 # m/s
            except:
                pass
            self.wave = self.l0 * (1 + velocity/c) #um
        else:
            print('Ctype3 ', ctype3,' is not supported')
        pixscale, ypixscale = proj_plane_pixel_scales(self.wcs) # Pixel scale in arcsec
        # Back to wavelength
        w = self.wave
        self.crval3 = w[0]
        self.cdelt3 = np.median(w[1:] - w[:-1])
        try:
            # Flux conversion from Jy/beam to Jy/pixel
            bmaj = self.header['BMAJ']
            bmin = self.header['BMIN']
            area_pixel = pixscale * ypixscale
            self.npix_per_beam = 1.331 * bmaj * bmin / area_pixel
            self.flux /= self.npix_per_beam
        except:
            self.npix_per_beam = 1
            print('No beam size given in the header')
        self.pixscale = pixscale * 3600.0 # pixel scale in arcsec
        print('scale is ', self.pixscale)
       
    def readMUSE(self, hdl):
        "MUSE integral field spectrometer at VLT"
        self.objname = self.header['OBJECT'].strip()
        print('Object of MUSE is ',self.objname)
        self.flux = hdl[1].read()  # 10**(-20)*erg/s/cm**2/Angstrom
        nz, ny, nx = np.shape(self.flux)
        self.n = nz
        wcs = WCS(self.header)
        self.wcs = wcs.celestial
        self.crpix3 = self.header['CRPIX3']
        self.crval3 = self.header['CRVAL3']
        self.cdelt3 = self.header['CDELT3']
        #ctype3 = self.header['CTYPE3'].strip()
        pix0=0
        self.wave = self.cdelt3 * (np.arange(self.n) - self.crpix3 + pix0) + self.crval3 # Angstrom
        self.wave *= 1.e4  # um
        self.pixscale, ypixscale = proj_plane_pixel_scales(self.wcs) * 3600. # Pixel scale in arcsec
        print('scale is ', self.pixscale)
                
    def getWCS(self):
        #hdu = fits.PrimaryHDU(self.flux)
        hdr = fits.Header()
        #hdu.header
        hdr['CRPIX1']=self.header['CRPIX1']
        hdr['CRPIX2']=self.header['CRPIX2']
        hdr['CDELT1']=self.header['CDELT1']
        hdr['CDELT2']=self.header['CDELT2']
        hdr['CRVAL1']=self.header['CRVAL1']
        hdr['CRVAL2']=self.header['CRVAL2']
        hdr['CTYPE1']=self.header['CTYPE1'].strip()
        hdr['CTYPE2']=self.header['CTYPE2'].strip()
        self.wcs = WCS(hdr).celestial 
    
    def getResolutionFIFI(self):
        """Compute resolution at reference wavelength for FIFI-LS"""
        if self.instrument != 'FIFI-LS':
            return
        l0 = self.l0
        z = self.redshift
        l = l0 * (1.+z)
        if self.channel == 'RED':
            return 0.062423 * l * l - 6.6595 * l + 647.65
        else:
            if self.order == '1':
                return 0.16864 * l * l - 22.831 * l + 1316.6 
            else:
                return 1.9163 * l * l - 187.35 * l + 5496.9

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
    def __init__(self, wave, flux, eflux=None, uflux=None, exposure=None,
                 atran=None, uatran=None, watran=None,
                 instrument=None, baryshift=None, redshift=None, l0=None, area=None, Tb2Jy=None,
                 bunit=None, yunit=None, pixscale=None):
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
        if uatran is not None:
            self.uatran = uatran
        else:
            self.uatran = None
        if watran is not None:
            self.watran = watran
        else:
            self.watran = None
        if instrument is not None:
            self.instrument = instrument
        if baryshift is not None:
            self.baryshift = baryshift
        if redshift is not None:
            self.redshift = redshift
        if l0 is not None:
            self.l0 = l0
        if area is not None:
            self.area = area
        if Tb2Jy is not None:
            self.Tb2Jy = Tb2Jy
        if bunit is not None:
            self.bunit = bunit
        if yunit is None:
            self.yunit = 'Jy/pix'
        else:
            self.yunit = yunit
        if pixscale is None:
            self.pixscale = 1
        else:
            self.pixscale = pixscale
        self.continuum =  np.full(len(wave), np.nan)
