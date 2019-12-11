import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

class specCubeAstro(object):
    """ spectral cube - read with AstroPy routines """
    def __init__(self, infile):
        import time
        t = time.process_time()
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
        elif self.instrument == 'MUSE':
            self.readMUSE(hdl)
        elif self.instrument == 'IRAM':
            self.readIRAM(hdl)
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
        try:
            self.filegpid = self.header['FILEGPID']
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
        except:
            print('No redshift present')
            self.redshift = 0.0
        self.flux = hdl['FLUX'].data
        self.eflux = hdl['ERROR'].data
        self.uflux = hdl['UNCORRECTED_FLUX'].data
        self.euflux = hdl['UNCORRECTED_ERROR'].data
        self.wave = hdl['WAVELENGTH'].data
        self.n = len(self.wave)
        #self.vel = np.zeros(self.n)  # prepare array of velocities
        self.x = hdl['X'].data
        self.y = hdl['Y'].data
        #self.atran = hdl['TRANSMISSION'].data
        self.channel = self.header['DETCHAN']
        if self.channel == 'BLUE':
            self.order = self.header["G_ORD_B"]
        else:
            self.order = '1'
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
        print('min wave ', np.nanmin(self.wave))
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
        from scipy.special import erf
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
        self.bunit = self.header['BUNIT']
        eta_fss=0.97
        eta_mb =0.67
        calib = 971.
        if self.bunit == 'K (Tmb)':
            self.Tb2Jy = calib*eta_fss*eta_mb
        else: # Case of K (Ta*)
            self.Tb2Jy = calib
        nu0 = self.header['RESTFREQ']  # MHz
        l0 = c/nu0  # in micron
        # Transform in Jy/pixel
        # Compute the beam size at the wavelength
        bmaj = self.header['BMAJ'] * 3600. # Beam major axis in arcsec
        bmin = self.header['BMIN'] * 3600. # Beam minor axis in arcsec
        # Multiply by the flux fraction in the pixel assuming a 2D Gaussian curve                    
        pixfraction = 0.5 * erf(self.pixscale*0.5/bmaj) * erf(ypixscale*0.5/bmin)
        print('Beam fraction on pixel ', pixfraction)
        self.Tb2Jy *= pixfraction
        # Multiplication is delayed in the code to speed up reading
        # self.flux *= self.Tb2Jy   # Transformed from temperature to S_nu [Jy] per pixel
        vel = self.cdelt3 * (np.arange(self.n) - self.crpix3 + 1) + self.crval3
        self.l0 = l0
        self.wave = l0 + l0*vel/c
        
    def readFORCAST(self, hdl): 
        print('This is a FORCAST spectral cube')
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
        self.eflux = np.sqrt(hdl['VARIANCE'].data)
        # nz, ny, nx = np.shape(self.flux)
        # print('nz is ',nz)
        exptime = self.header['EXPTIME']
        exp = hdl['EXPOSURE'].data.astype(float) * exptime
        self.exposure = np.broadcast_to(exp, np.shape(self.flux))
        print('shape of exposure is ', np.shape(self.exposure))
        self.wave = self.cdelt3 * (np.arange(self.n) - self.crpix3 + 1) + self.crval3
        self.l0 = np.nanmedian(self.wave)

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
        #print('Flux read')
        try:
            self.exposure = hdl['coverage'].data
            print('Coverage read')
        except:
            print('No coverage available - range observation')
            self.computeExpFromNan()
            print('New exposure computed')
            
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
        print('astrometry ', self.wcs)
        self.pixscale, ypixscale = proj_plane_pixel_scales(self.wcs) * 3600. # Pixel scale in arcsec
        self.crpix3 = 1
        w = self.wave
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
            velocity = self.cdelt3 * (np.arange(self.n) - self.crpix3 + 1) + self.crval3 # m/s
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
        self.n = nz
        print('nz: ',nz, self.header['NAXIS3'])
        wcs = WCS(self.header)
        self.wcs = wcs.celestial
        self.crpix3 = self.header['CRPIX3']
        self.crval3 = self.header['CRVAL3']
        self.cdelt3 = self.header['CDELT3']
        ctype3 = self.header['CTYPE3']
        if ctype3 == 'VELOCITY':
            velocity = self.cdelt3 * (np.arange(self.n) - self.crpix3 + 1) + self.crval3 # m/s
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
        
    def readMUSE(self, hdl):
        "MUSE integral field spectrometer at VLT"
        self.objname = self.header['OBJECT']
        print('Object of MUSE is ',self.objname)
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
        self.wave = self.cdelt3 * (np.arange(self.n) - self.crpix3 + 1) + self.crval3 # Angstrom
        self.wave *= 1.e-4  # um
        #print('min wave ', np.nanmin(self.wave))
        self.pixscale, ypixscale = proj_plane_pixel_scales(self.wcs) * 3600. # Pixel scale in arcsec
        print('scale is ', self.pixscale)
        self.redshift = 0
        self.l0 = np.nanmedian(self.wave)
        #self.l0 = (self.header['WAVELMIN']+self.header['WAVELMAX']) * 0.5 * 1.e-3 # wav in nm
    
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
        try:
            self.instrument = self.header['INSTRUME'].strip()
            print('Instrument: ', self.instrument)
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
        elif self.instrument in ['VLA','ALMA','CARMA']:
            print('Ok on VLA')
            self.readVLA(hdl)
        elif self.instrument == 'MUSE':
            self.readMUSE(hdl)
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
                self.flux[idx] = np.nan
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
        if self.bunit == 'K (Tmb)':
            self.Tb2Jy = calib*eta_fss*eta_mb
        else: # Case of K (Ta*)
            self.Tb2Jy = calib
        nu0 = self.header['RESTFREQ']  # MHz
        print('rest freq ', nu0)
        l0 = c/nu0  # in micron
        # Transform in Jy/pixel
        # Compute the beam size at the wavelength
        bmaj = self.header['BMAJ'] * 3600. # Beam major axis in arcsec
        bmin = self.header['BMIN'] * 3600. # Beam minor axis in arcsec
        xsigma = 2 * bmaj / 2.355 / np.sqrt(2) # sigma from FWHM divided by sqrt(2)
        ysigma = 2 * bmin / 2.355 / np.sqrt(2)
        # Multiply by the flux fraction in the pixel assuming a 2D Gaussian curve                    
        pixfraction = erf(self.pixscale*0.5/xsigma) * erf(ypixscale*0.5/ysigma)
        print('Beam fraction on pixel ', pixfraction)
        self.Tb2Jy *= pixfraction
        naxes = self.header['NAXIS']
        print('no of axes ', naxes)
        if naxes == 4:
            self.flux = (hdl[0].read())[0,:,:,:]
        else:
            self.flux = hdl[0].read()
        t1 = time.process_time()
        print('Flux reading completed in ', t1-t,' s')
        #self.flux *= self.Tb2Jy   # Transformed from temperature to S_nu [Jy]  Move this to code
        vel = self.cdelt3 * (np.arange(self.n) - self.crpix3 + 1) + self.crval3
        self.l0 = l0
        self.wave = l0 + l0 * vel / c
        #t2 = time.process_time()
        #print('GREAT transf completed in ', t2-t1,' s')
        
    def readFORCAST(self, hdl): 
        print('This is a FORCAST spectral cube')
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
        self.eflux = np.sqrt(hdl[extnames.index('VARIANCE')].read())
        # nz, ny, nx = np.shape(self.flux)
        # print('nz is ',nz)
        exptime = self.header['EXPTIME']
        exp = hdl[extnames.index('EXPOSURE')].read().astype(float) * exptime
        self.exposure = np.broadcast_to(exp, np.shape(self.flux))
        print('shape of exposure is ', np.shape(self.exposure))
        self.wave = self.cdelt3 * (np.arange(self.n) - self.crpix3 + 1) + self.crval3
        self.l0 = np.nanmedian(self.wave)

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
        print('Object is ',self.objname)
        self.flux = hdl[extnames.index('image')].read()
        #print('Flux read')
        try:
            self.exposure = hdl[extnames.index('coverage')].read()
            print('Coverage read')
        except:
            print('No coverage available - range observation')
            self.computeExpFromNan()
            print('New exposure computed')
            print('shape of exp ', np.shape(self.exposure))
            
        wave = hdl[extnames.index('wcs-tab')].read()
        #print('Wvl read')
        nwave = len(np.shape(wave['wavelen']))
        if nwave == 3:
            self.wave = np.concatenate(wave['wavelen'][0])
        else:
            self.wave = np.concatenate(wave['wavelen'])
        self.l0 = np.nanmedian(self.wave)
        self.n = len(self.wave)
        #print('Length of wavelength ',self.n)
        #print('Min and max wavelengths: ', np.nanmin(self.wave), np.nanmax(self.wave))
        #print(np.shape(self.wave))
        #print('image ext ', extnames.index('image'))
        header = hdl[extnames.index('image')].read_header()
        self.header = header
        #self.purifyHeader()
        #self.wcs = WCS(self.header)
        self.getWCS()
        #print('astrometry ', self.wcs)
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
        if (ctype3 == 'VELO-HEL') or (ctype3 == 'VELO-LSR'):
            velocity = self.cdelt3 * (np.arange(self.n) - self.crpix3 + 1) + self.crval3 # m/s
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
        self.n = nz
        print('nz: ',nz, self.header['NAXIS3'])
        wcs = WCS(self.header)
        self.wcs = wcs.celestial
        self.crpix3 = self.header['CRPIX3']
        self.crval3 = self.header['CRVAL3']
        self.cdelt3 = self.header['CDELT3']
        ctype3 = self.header['CTYPE3'].strip()
        if ctype3 == 'VELOCITY':
            velocity = self.cdelt3 * (np.arange(self.n) - self.crpix3 + 1) + self.crval3 # m/s
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
        self.redshift = 0.
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
                self.redshift = altrval/c
            except:
                pass
            freq = self.cdelt3 * (np.arange(self.n) - self.crpix3 + 1) + self.crval3 
            self.l0 = c/nu0 * 1.e6 #um
            self.wave = c/freq * 1.e6 #um
            # Reorder wave (and flux)
            idx = np.argsort(self.wave)
            self.wave = self.wave[idx]
            self.flux = self.flux[idx, :, :]
        elif ctype3 in ['VELO-HEL', 'VELO-LSR', 'VRAD']:
            velocity = self.cdelt3 * (np.arange(self.n) - self.crpix3 + 1) + self.crval3 # m/s
            if self.instrument == 'VLA':
                nu0 = self.header['RESTFREQ']
            elif self.instrument in ['ALMA', 'CARMA']:
                nu0 = self.header['RESTFRQ']
            print('reference frequency', nu0)
            c = 299792458.0 # m/s
            # self.l0 = 21.1061140542 * 1.e4 #um
            self.l0 = c/nu0 * 1.e6 #um
            try:
                altrval = self.header['ALTRVAL']  # ref frequency
                altrpix = self.header['ALTRPIX']  # ref pixel
                vcrval3 = c * (1. - altrval/nu0)
                #faltrval = nu0 * (1 - altrval/c)  # alt crval3
                if np.abs(self.crval3 - vcrval3) > 20:
                    self.crval3 = vcrval3
                    self.crpix3 = altrpix
                    self.redshift = vcrval3 / c
                    velocity = self.cdelt3 * (np.arange(self.n) - self.crpix3 + 1) + self.crval3 # m/s
            except:
                pass
            self.wave = self.l0 * (1 + velocity/c) #um
        else:
            print('Ctype3 ', ctype3,' is not supported')
        self.pixscale, ypixscale = proj_plane_pixel_scales(self.wcs) * 3600. # Pixel scale in arcsec
        print('scale is ', self.pixscale)
        # Back to wavelength
        w = self.wave
        self.crval3 = w[0]
        self.cdelt3 = np.median(w[1:] - w[:-1])
       
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
        self.wave = self.cdelt3 * (np.arange(self.n) - self.crpix3 + 1) + self.crval3 # Angstrom
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
                 bunit=None):
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
            self.baryshift=baryshift
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
        self.continuum =  np.full(len(wave), np.nan)
