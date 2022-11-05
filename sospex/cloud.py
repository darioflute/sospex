import urllib, urllib.request
import os
from io import BytesIO
from astropy.io import ascii,fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.utils.data import download_file
from astropy.wcs import WCS
from reproject import reproject_interp
import numpy as np
from html.parser import HTMLParser
from PyQt5.QtWidgets import QFileDialog, QMessageBox
from sospex.specobj import specCube, specCubeAstro

class MyHTMLParser(HTMLParser):

    def __init__(self):
        HTMLParser.__init__(self)
        self.data = []
        self.values = []
        self.recording = 0

    def handle_starttag(self, tag, attrs):
        # Only parse the 'anchor' tag.
        if tag == "a" or tag == "base":
           # Check the list of defined attributes.
           for name, value in attrs:
               # If href is defined, print it.
               if name == "href":
                    #print name, value
                    self.recording = 1
                    self.values.append(value)
                    
    def handle_endtag(self,tag):
        if tag == 'a':
            self.recording -=1 
            #print "Encountered the end of a %s tag" % tag 

    def handle_data(self, data):
        if self.recording:
            self.data.append(data)

class cloudImage(object):
    def __init__(self, lon, lat, xsize, ysize, source, pixscale=None):

        self.lon = lon
        self.lat = lat
        self.xsize = xsize
        self.ysize = ysize
        self.source = source
        self.data = None
        self.wcs = None
        self.crota2 = 0
        if pixscale is not None:
            self.pixscale = pixscale

        if source in ['local image', 'velocity image']:
            self.openLocal()
        elif source == 'local spectral cube':
            self.openLocalCube()
        elif source == 'wise1':
            self.downloadWise(1)
        elif source == 'wise2':
            self.downloadWise(2)
        elif source == 'wise3':
            self.downloadWise(3)
        elif source == 'wise4':
            self.downloadWise(4)
        elif source == '2mass-j':
            self.download2MASS('J')
        elif source == '2mass-h':
            self.download2MASS('H')
        elif source == '2mass-k':
            self.download2MASS('K')
        elif source == 'panstarrs-g':
            self.downloadPanSTARRS('g')
        elif source == 'panstarrs-r':
            self.downloadPanSTARRS('r')
        elif source == 'panstarrs-i':
            self.downloadPanSTARRS('i')
        elif source == 'panstarrs-z':
            self.downloadPanSTARRS('z')
        elif source == 'panstarrs-y':
            self.downloadPanSTARRS('y')
        elif source == 'sdss-u':
            self.downloadSDSS('u')
        elif source == 'sdss-r':
            self.downloadSDSS('g')
        elif source == 'sdss-g':
            self.downloadSDSS('r')
        elif source == 'sdss-i':
            self.downloadSDSS('i')
        elif source == 'sdss-z':
            self.downloadSDSS('z')
        elif source == 'first':
            self.downloadFIRST()
        elif source == 'nvss':
            self.downloadNVSS()
        elif source == 'sumss':
            self.downloadSUMSS()
        else:
            print('Source not supported')

        # possible contours
        self.contours = None

    def download2MASS(self, band):

        c = SkyCoord(ra=self.lon*u.degree, dec=self.lat*u.degree, frame='icrs')
        coords = c.to_string('hmsdms',sep=' ')+' Equ J2000'
        size = int(np.max([self.xsize,self.ysize])*60.) # Size in arcsec
        if size > 1024: size = 1024

        url = 'http://irsa.ipac.caltech.edu/cgi-bin/2MASS/IM/nph-im_pos'
        
        
        post_params = {
            'POS': coords,
            'subsz': size,
            'band': band,
            '':'Submit'
        }
        
        post_args = urllib.parse.urlencode(post_params)#.encode("utf-8") 
        full_url=url + '?' + post_args

        #print(full_url)
        response = urllib.request.urlopen(full_url)
        html = response.read()
        parser = MyHTMLParser()
        parser.feed(html.decode('utf-8'))
        # data= parser.data
        values= parser.values
        parser.close()

        #for d,v in zip(data,values):
        #    print (d,' ',v)
        
        file = None
        base = None
        for v in values:
            if '.fits' in v:
                file = v
            if 'irsa.ipac' in v:
                base = v    

        file = base+file
        #print('file is: ', file)
        
        if '.fits' in file:
            #print('file is: ',file)
            response = urllib.request.urlopen(file)
            output = response.read()
            fitsfile= BytesIO(output)  # Read the downloaded FITS data
            hdulist = fits.open(fitsfile,memmap=False)
            self.header = hdulist['PRIMARY'].header
            self.data = hdulist['PRIMARY'].data
            hdulist.close()
            self.wcs = WCS(self.header).celestial
        else:
            self.data = None
            self.wcs = None
            print('Coordinates out of the 2MASS survey')


    def openLocalCube(self):
        """ Open local spectral cube."""
        
        # Open a dialog
        fd = QFileDialog()
        fd.setLabelText(QFileDialog.Accept, "Import")
        fd.setNameFilters(["Fits Files (*.fit*)","WXY fits files (*WXY*.fits*)", "All Files (*)"])
        fd.setOptions(QFileDialog.DontUseNativeDialog)
        fd.setViewMode(QFileDialog.List)
        fd.setFileMode(QFileDialog.ExistingFile)

        if (fd.exec()):
            filenames= fd.selectedFiles()
            cube_file = filenames[0]
            print("File selected is: ", filenames[0])
            try:           
                print('opening ', cube_file)
                try:
                    self.data = specCube(cube_file)
                except:
                    self.data = specCubeAstro(cube_file)
                self.wcs = self.data.wcs
                newpixel = self.data.pixscale
                pixel = self.pixscale
                ratio = (pixel/newpixel)**2
                print('area ratio: ', ratio)
                self.data.flux *= ratio                    
                # Check if coordinates are inside the image
                # This should be modified to check if the cube is inside the axes
                #x, y = self.wcs.wcs_world2pix(self.lon, self.lat, 0)
                #print('x y ', x, y)
                #if x >= 0 and x <= self.data.nx and y >= 0 and y  <= self.data.ny:
                #    print('Source inside the FITS cube')
                #    # At this point the fluxes are rescaled to the pixel size
                #    # of the original cube
                #    newpixel = self.data.pixscale
                #    pixel = self.specCube.pixscale
                #    ratio = (newpixel/pixel)**2
                #    self.data.flux *= ratio                    
                #else:
                #    self.data = None
                #    self.wcs = None
                #    print('The selected  FITS is not a valid file')
            except:
                self.data = None
                self.wcs = None
                print('The selected  FITS is not a valid file')
                

    def openLocal(self):
        """ Open a local FITS / check if coordinates fall into the FITS """
        
        # Open a dialog
        fd = QFileDialog()
        fd.setLabelText(QFileDialog.Accept, "Import")
        fd.setNameFilters(["Fits Files (*.fit*)","All Files (*)"])
        fd.setOptions(QFileDialog.DontUseNativeDialog)
        fd.setViewMode(QFileDialog.List)
        fd.setFileMode(QFileDialog.ExistingFile)

        if (fd.exec()):
            filenames= fd.selectedFiles()
            image_file = filenames[0]
            print("File selected is: ", filenames[0])
            try:           
                print('opening ',image_file)
                hdulist = fits.open(image_file,memmap=False)
                hdu = hdulist['PRIMARY']
                hdulist.info()
                header = hdulist['PRIMARY'].header
                self.header = header
                try:
                    instrument = header['INSTRUME']
                    self.source = instrument
                    try:
                        if instrument == 'PACS':
                            w = header['WAVELNTH']
                            self.source += str(int(w // 1))
                    except:
                        pass
                    try:
                        if instrument == 'IRAC':
                            fc = header['FLUXCONV']
                            #print('fc ', fc)
                            if fc == 0.1069:
                                self.source += '1'
                            elif fc == 0.1382:
                                self.source += '2'
                            elif fc == 0.5858:
                                self.source += '3'
                            elif fc == 0.2026:
                                self.source += '4'
                    except:
                        pass   
                    if instrument == 'MIPS':
                        channel = header['CHNLNUM']
                        if channel == 1:
                            self.source +='24'
                        elif channel == 2:
                            self.source +='70'
                        elif channel == 3:
                            self.source +='160'
                except:
                    self.source = 'unknown'
                print('Instrument ',self.source)
                naxis = header['NAXIS']
                if naxis == 2:
                    self.data = hdulist['PRIMARY'].data
                elif naxis == 3:
                    self.data = hdulist['PRIMARY'].data[0,:,:]
                    try:
                        del header['CTYPE4']
                    except:
                        pass
                elif naxis > 3:
                    self.data = hdulist['PRIMARY'].data[0,0,:,:]
                elif naxis == 0:
                    # Switch to image extensions if exists.
                    try:
                        hdu = hdulist['IMAGE']
                        header = hdu.header
                        self.data = hdu.data
                    except:
                        print('This is not an image')
                hdulist.close()
                # NaN for infinite values
                idx = np.isfinite(self.data)
                if np.sum(~idx) > 0:
                    self.data[~idx] = np.nan
                try:
                    # Transform FK4 into ICRS
                    if header['EPOCH'] == 1950:
                        print('This image is in epoch 1950')
                        from astropy import units as u
                        from astropy.coordinates import SkyCoord
                        cold = SkyCoord(ra=header['crval1']*u.degree, dec=header['crval2']*u.degree, frame='fk4')
                        cnew = cold.transform_to('icrs')
                        header['CRVAL1'] = cnew.ra.value
                        header['CRVAL2'] = cnew.dec.value
                        header['EPOCH'] = 2000
                except:
                    pass
                self.wcs = WCS(header).celestial
                # Check if coordinates are inside the image
                x, y = self.wcs.wcs_world2pix(self.lon,self.lat,0)
                ny,nx = np.shape(self.data)
                print('WCS - prerot ', self.wcs)
                print('Data - prerot ', np.shape(self.data))
                if x >= 0 and x <= nx and y >= 0 and y  <= ny:
                    print('Source inside the FITS image')
                    self.rotateImage(image_file, header)
                    print('WCS ', self.wcs)
                    print('Data ', np.shape(self.data))
                else:
                    self.data = None
                    self.wcs = None
                    print('Coordinates out of the selected FITS image')
            except:
                self.data = None
                self.wcs = None
                print('The selected  FITS is not a valid file')

    def purifyHeader(self, orig_header):
        header = fits.Header()
        for dict_key in orig_header.keys():
            try:
                header[dict_key] = orig_header[dict_key]
            except:
                pass
        return header

    def openLocalFitsio(self):
        """ Open a local FITS / check if coordinates fall into the FITS """
        import fitsio        
        # Open a dialog
        fd = QFileDialog()
        fd.setLabelText(QFileDialog.Accept, "Import")
        fd.setNameFilters(["Fits Files (*.fits)","All Files (*)"])
        fd.setOptions(QFileDialog.DontUseNativeDialog)
        fd.setViewMode(QFileDialog.List)
        fd.setFileMode(QFileDialog.ExistingFile)

        if (fd.exec()):
            filenames= fd.selectedFiles()
            image_file = filenames[0]
            print("File selected is: ", filenames[0])
            try:           
                print('opening ',image_file)
                hdulist = fitsio.FITS(image_file)
                hdu = hdulist[0]
                #hdulist.info()
                h = hdu.read_header()
                header = self.purifyHeader(h)
                try:
                    instrument = header['INSTRUME']
                    self.source = instrument
                    try:
                        if instrument == 'PACS':
                            w = header['WAVELNTH']
                            self.source += str(int(w // 1))
                    except:
                        pass
                except:
                    self.source = 'unknown'
                print('Instrument ',self.source)
                naxis = header['NAXIS']
                print('Number of axes ', naxis)
                if naxis == 2:
                    self.data = hdulist[0].read()
                elif naxis > 2:
                    self.data = hdulist[0].read()[0,0,:,:]
                elif naxis == 0:
                    # Switch to image extensions if exists.
                    extnames = [f.get_extname() for f in hdulist]
                    print('ext names ', extnames)
                    if 'IMAGE' in extnames:
                        hdu = hdulist[extnames.index('IMAGE')]
                        h = hdu.read_header()
                        header = self.purifyHeader(h)
                        self.data = hdu.data.read()
                    elif 'image' in extnames:
                        hdu = hdulist[extnames.index('image')]
                        h = hdu.read_header()
                        header = self.purifyHeader(h)
                        self.data = hdu.data.read()
                    else:
                        print('This is not an image')
                        return
                hdulist.close()
                # NaN for infinite values
                idx = np.isfinite(self.data)
                self.data[~idx] = np.nan
                #print('header: ',header)
                self.wcs = WCS(header).celestial
                # Check if coordinates are inside the image
                x,y = self.wcs.wcs_world2pix(self.lon,self.lat,0)
                print('x y ',x,y)
                ny,nx = np.shape(self.data)
                if x >= 0 and x <= nx and y >= 0 and y  <= ny:
                    print('Source inside the FITS image')
                    self.rotateImage(image_file, header)
                    print('WCS ', self.wcs)
                    print('Data ', np.shape(self.data))
                else:
                    self.data = None
                    self.wcs = None
                    print('Coordinates out of the selected FITS image')
            except:
                self.data = None
                self.wcs = None
                print('The selected  FITS is not a valid file')

    def rotateImage(self, image_file, header):
        """ Rotate the FITS image """
        
        h1 = self.wcs.to_header()
        try:
            self.crota2 = np.arctan2(-h1["PC2_1"], h1["PC2_2"]) * 180./np.pi
        except:
            self.crota2 = 0.
        print('Rotation ', self.crota2)
        try:
            print('h1 ',h1)
            pc11=h1["PC1_1"]
            pc12=h1["PC1_2"]
            pc21=h1["PC2_1"]
            pc22=h1["PC2_2"]
            pc1 = -np.hypot(pc11,pc12)
            pc2 = np.hypot(pc21,pc22)
            rota = np.arctan(pc21/pc22)    
            crota = rota * 180. / np.pi 
        except:
            crota = 0
        print('crota ', crota)
        if (np.abs(crota) > 5.):
            flags = QMessageBox.Yes | QMessageBox.No
            question = "The image is rotated by {:.2f} degs. \n".format(crota)
            question += "Do you want to rotate the image to N top, E left ?"
            msgBox = QMessageBox()
            msgBox.setText(question)
            msgBox.setStandardButtons(flags)
            response = msgBox.exec()
            if response == QMessageBox.Yes:
                if h1["CRVAL2"] < 0:
                    pc2 = -pc2
                h1.update(pc1_1=pc1, pc1_2=0.0, pc2_1=0.0, pc2_2=pc2, orientat=0.)
                h1['NAXIS'] = 2
                # Compute the size of the rectangle containing the rotated rectangle
                rota = np.abs(rota)
                ny, nx = np.shape(self.data)
                n1 = int(np.rint(nx*np.sin(rota)+ny*np.cos(rota)))
                n2 = int(np.rint(nx*np.cos(rota)+ny*np.sin(rota)))
                h1['NAXIS1'] = n1
                h1['NAXIS2'] = n2
                crpix1 = header['CRPIX1']
                crpix2 = header['CRPIX2']
                h1['crpix1'] = crpix1+(n1-nx)/2.
                h1['crpix2'] = crpix2+(n2-ny)/2.
                print("Rotating the image ....")
                # Generate hdu
                hdu = fits.PrimaryHDU(self.data)
                hdu.header.extend(header)
                # Reproject
                array, footprint = reproject_interp(hdu, h1)
                print('Shape of rotated image: ', np.shape(array))  
                # Rotated data and wcs
                self.data = array
                self.wcs = WCS(h1)
                self.header = h1
                self.crota2 = 0. # there is no rotation left
                # Save the rotated image
                self.saveRotatedFits(image_file)
            else:
                message = "The image will appear rotated by: {:.2f} ".format(self.crota2)
                print(message)
                print('wcs ', self.wcs)
                print('data ', np.shape(self.data))
                #self.sb.showMessage(message)
        else:
            message = 'No rotation needed. Rotation angle: {:.2f}'.format(crota)
            #self.sb.showMessage(message)
            print(message)
        #except
        #    message = "Rotation failed"
        #    self.sb.showMessage(message)
        #    print(message)
        print('End of rotation module')

    def saveRotatedFits(self, name_orig):
        """ Save the downloaded FITS image """
        
        filename, file_extension = os.path.splitext(name_orig)
        #fileroot = os.path.basename(filename)
        #print('file root is ', fileroot)
        print('Saving ',filename+'_NE.fits')
        # Dialog to save file
        fd = QFileDialog()
        fd.setLabelText(QFileDialog.Accept, "Save as")
        fd.setNameFilters(["Fits Files (*.fits)","All Files (*)"])
        fd.setOptions(QFileDialog.DontUseNativeDialog)
        fd.setViewMode(QFileDialog.List)
        fd.selectFile(filename+'_NE.fits')        
        if (fd.exec()):
            fileName = fd.selectedFiles()
            outfile = fileName[0]
            filename, file_extension = os.path.splitext(outfile)
            # basename = os.path.basename(filename)            
            # Primary header
            image = self.data
            wcs   = self.wcs
            header = wcs.to_header()
            header.remove('WCSAXES')
            header['INSTRUME'] = (self.source, 'Instrument')
            hdu = fits.PrimaryHDU(image)
            hdu.header.extend(header)
            hdul = fits.HDUList([hdu])
            hdul.writeto(outfile,overwrite=True) # clobber true  allows rewriting
            hdul.close()

                
    def downloadWise(self,band):
        """ Download a Wise image """

        url = 'http://irsa.ipac.caltech.edu/ibe/search/wise/allsky/4band_p3am_cdd?POS='+\
                            "{0:.3f}".format(float(self.lon))+",{0:.3f}".format(float(self.lat))
        print(url)
        request = urllib.request.Request(url)
        response = urllib.request.urlopen(request)
        table = response.read()
        t = BytesIO(table)
        data = ascii.read(t,format='ipac')
        bands = data['band']
        coadds = data['coadd_id']
        bands = np.array(bands)
        coadds = np.array(coadds)
        coadd = coadds[bands == band]
        print(coadd)

        params = { 'coadd_id': coadd[0],'band': band,}
        params['coaddgrp'] = params['coadd_id'][:2]
        params['coadd_ra'] = params['coadd_id'][:4]
        path = str.format(
            '{coaddgrp:s}/{coadd_ra:s}/{coadd_id:s}/{coadd_id:s}-w{band:1d}-int-3.fits',
            **params)

        # We again specify center and size in WISE pixels (what is the pixel of WISE ? 1arcmin=45pix)
        #size=30,45arcsec
        #s1 = "{0:.0f}".format(pixscale[0]*60*ny*45)
        #s2 = "{0:.0f}".format(pixscale[1]*60*nx*45)
        s1 = "{0:.0f}".format(self.xsize*45)
        s2 = "{0:.0f}".format(self.ysize*45)

        url = 'http://irsa.ipac.caltech.edu/ibe/data/wise/allsky/4band_p3am_cdd/' + path+\
              "?center={0:.3f}".format(float(self.lon))+",{0:.3f}".format(float(self.lat))+"&size="+s1+','+s2+"pix&gzip=false"
        image_file = download_file(url,cache=True)

        hdulist = fits.open(image_file,memmap=False)
        #hdulist.info()
        self.header = hdulist['PRIMARY'].header
        self.data = hdulist['PRIMARY'].data
        hdulist.close()
        self.wcs = WCS(self.header)


    def downloadFIRST(self):
        """ Download data from the VLA FIRST survey """

        c = SkyCoord(ra=self.lon*u.degree, dec=self.lat*u.degree, frame='icrs')
        coords = c.to_string('hmsdms',sep=' ')
        #print ('coords of the center are: ', coords)
        Equinox = 'J2000'
        isize = np.max([self.xsize,self.ysize]) # Size in arcmin
        itype ='FITS Image'

        url='https://third.ucllnl.org/cgi-bin/firstcutout'
        post_params = { 
            'RA'  : coords,
            'Equinox': Equinox,
            'ImageSize': isize,
            'ImageType': itype,
            '.submit': " Extract the Cutout "
        }
        post_args = urllib.parse.urlencode(post_params).encode("utf-8")        
        request = urllib.request.Request(url, post_args)
        response = urllib.request.urlopen(request)
        output = response.read()
        fitsfile= BytesIO(output)  # Read the downloaded FITS data

        # If a FITS file was downloaded pass image and header, otherwise failes
        try:
            hdulist = fits.open(fitsfile,memmap=False)
            self.header = hdulist['PRIMARY'].header
            self.data = hdulist['PRIMARY'].data
            hdulist.close()
            self.wcs = WCS(self.header).celestial
        except:
            self.data = None
            self.wcs = None
            print('Coordinates out of the FIRST VLA survey')


    def downloadNVSS(self):
        """ Download data from the Sydney University Molonglo Sky Survey """

        c = SkyCoord(ra=self.lon*u.degree, dec=self.lat*u.degree, frame='icrs')
        cra = c.ra.hms
        cdec = c.dec.dms
        ra = "{:.0f} {:.0f} {:.1f}".format(cra[0],cra[1],cra[2])
        if self.lat > 0:
            dec = "+{:.0f} {:.0f} {:.1f}".format(cdec[0],cdec[1],cdec[2])
        else:
            dec = "-{:.0f} {:.0f} {:.1f}".format(-cdec[0],-cdec[1],-cdec[2])
        Equinox = 'J2000'
        xsize = self.xsize/60. # Size in degs
        ysize = self.ysize/60.
        fieldsize = "{:.1f} {:.1f}".format(xsize,ysize)
        itype ='application/octet-stream'
        
        #print('ra ',ra)
        #print('dec ',dec)
        url='http://www.cv.nrao.edu/cgi-bin/postage.pl'
        #url='http://www.cv.nrao.edu/nvss/postage.shtml'
        post_params = {
            'Equinox': Equinox,
            'PolType': 'I',
            'ObjName':'',
            'RA'  : ra,
            'Dec' : dec,
            'Size': fieldsize,
            'Cells':'15.0 15.0',
            'MAPROJ':'SIN',
            'rotate': '0',
            'Type': itype,
            '': "Submit!"
        }

        post_args = urllib.parse.urlencode(post_params)#.encode("utf-8") 
        full_url=url + '?' + post_args

        try:
            response = urllib.request.urlopen(full_url)
            output = response.read()
            fitsfile= BytesIO(output)  # Read the downloaded FITS data
            hdulist = fits.open(fitsfile,memmap=False)
            hdulist.info()
            self.header = hdulist['PRIMARY'].header
            self.data = hdulist['PRIMARY'].data[0,0,:,:]
            hdulist.close()
            self.wcs = WCS(self.header).celestial
        except:
            self.data = None
            self.wcs = None
            print('Coordinates out of the NVSS survey')

            
    def downloadSUMSS(self):
        """ Download data from the Sydney University Molonglo Sky Survey """

        c = SkyCoord(ra=self.lon*u.degree, dec=self.lat*u.degree, frame='icrs')
        cra = c.ra.hms
        cdec = c.dec.dms
        ra = "{:.0f} {:.0f} {:.1f}".format(cra[0],cra[1],cra[2])
        if self.lat > 0:
            dec = "+{:.0f} {:.0f} {:.1f}".format(cdec[0],cdec[1],cdec[2])
        else:
            dec = "-{:.0f} {:.0f} {:.1f}".format(-cdec[0],-cdec[1],-cdec[2])
        Equinox = 'J2000'
        xsize = self.xsize/60. # Size in degs
        ysize = self.ysize/60.
        fieldsize = "{:.1f} {:.1f}".format(xsize,ysize)
        itype ='FITS'
        
        
        url='http://www.astrop.physics.usyd.edu.au/cgi-bin/postage.pl'
        post_params = { 
            'RA'  : ra,
            'DEC' : dec,
            'fieldsize': fieldsize,
            'scale':'11 11',
            'equinox': Equinox,
            'projection':'SIN',
            'rotation': '0',
            'imagetypes': itype,
            'action': "submit"
        }
        post_args = urllib.parse.urlencode(post_params).encode("utf-8")        
        request = urllib.request.Request(url, post_args)
        response = urllib.request.urlopen(request)
        html = response.read()

        parser = MyHTMLParser()
        parser.feed(html.decode('utf-8'))
        # data= parser.data
        values= parser.values
        parser.close()
        # Retrieve file
        file = values[0]
        # Check if fits file:
        if '.fits' in file:
            #print('file is: ',file)
            response = urllib.request.urlopen(file)
            output = response.read()
            fitsfile= BytesIO(output)  # Read the downloaded FITS data
            hdulist = fits.open(fitsfile,memmap=False)
            self.header = hdulist['PRIMARY'].header
            self.data = hdulist['PRIMARY'].data
            hdulist.close()
            self.wcs = WCS(self.header).celestial
        else:
            self.data = None
            self.wcs = None
            print('Coordinates out of the SUMSS survey')

    def downloadPanSTARRS(self, band):
        """ Download Pan STARRS images """
        url='http://ps1images.stsci.edu/cgi-bin/ps1cutouts'
        
        # scale: 60 arcsec = 240 pix, 1" = 4 pixels
        c = SkyCoord(ra=self.lon*u.degree, dec=self.lat*u.degree, frame='icrs')
        pos = c.to_string('hmsdms',sep=' ')

        fieldsize = int(np.max([self.xsize,self.ysize])*240.) # in pixels
        if fieldsize > 6000: fieldsize = 6000  # current limit
        
        gchecked=''
        rchecked=''
        ichecked=''
        zchecked=''
        ychecked=''
        
        if band == 'y':
            ychecked = 'checked'
        elif band == 'g':
            gchecked = 'checked'
        elif band == 'r':
            rchecked = 'checked'
        elif band == 'i':
            ichecked = 'checked'
        elif band == 'z':
            zchecked = 'checked'
            
        
        post_params = { 
            'pos'  : pos,
            'color':'checked',
            'g': gchecked,
            'r': rchecked,
            'i': ichecked,
            'z': zchecked,
            'y': ychecked,
            'data':'checked',
            'mask':'',
            'wt':'',
            'exp':'',
            'expwt':'',
            'num':'',
            'size': fieldsize,
            'submit': "Submit"
        }

        post_args = urllib.parse.urlencode(post_params).encode("utf-8")        
        request = urllib.request.Request(url, post_args)
        response = urllib.request.urlopen(request)
        html = response.read()

        parser = MyHTMLParser()
        parser.feed(html.decode('utf-8'))
        # data= parser.data
        values= parser.values
        parser.close()

        file = None
        for v in values:
            if '.'+band+'.unconv.fits' in v:
                file = 'http:'+v
                print('File: ', file)
                
        if file is not None:
            #idx = data.index('FITS-cutout')
            #file = 'http:'+values[idx]
            try:
                response = urllib.request.urlopen(file)
                content= response.read()
                fitsfile = BytesIO(content)
                hdulist = fits.open(fitsfile,memmap=False)
                self.header = hdulist['PRIMARY'].header
                self.data = hdulist['PRIMARY'].data
                hdulist.close()
                self.wcs = WCS(self.header).celestial
            except:
                print('Problems with connection')
        else:
            self.data = None
            self.wcs = None
            print('Coordinates out of the PanSTARRS survey')


    def downloadSDSS(self,band):
        """ Download a script to build an image using SWarp - also needs wget and bzip2 installed"""
        #url='https://dr12.sdss.org/mosaics/'

        import bz2

        # Other way:
        # Build a search with:
        # https://dr12.sdss.org/fields/raDec?ra=00+41+37.8+&dec=-09+20+33
        #
        # Select a FITS image
        # Perhaps cut it if too big
        c = SkyCoord(ra=self.lon*u.degree, dec=self.lat*u.degree, frame='icrs')
        cra = c.ra.hms
        cdec = c.dec.dms
        url0="https://dr12.sdss.org/fields/raDec?ra="
        ras = "{:02.0f}+{:02.0f}+{:.1f}".format(cra[0],cra[1],cra[2])
        if self.lat > 0:
            url = url0+ras+"+&dec={:02.0f}+{:02.0f}+{:02.0f}".format(cdec[0],cdec[1],cdec[2])
        else:
            url = url0+ras+"+&dec=-{:02.0f}+{:02.0f}+{:02.0f}".format(-cdec[0],-cdec[1],-cdec[2])

        print('request ', url)
        request = urllib.request.Request(url)
        response = urllib.request.urlopen(request)
        html = response.read()
        
        parser = MyHTMLParser()
        parser.feed(html.decode('utf-8'))
        # data= parser.data
        values= parser.values
        parser.close()

        file = None
        for v in values:
            if 'frame-'+band+'-' in v:
                file = 'https://dr14.sdss.org/'+v
                print('file ', file)
                
        if file is not None:
            #idx = data.index(band+'-band FITS')
            #print("Value is:", values[idx-1])
            #file = 'https://dr12.sdss.org/'+values[idx-1]
            try:
                response = urllib.request.urlopen(file)
                content= response.read()        
                fitsfile = BytesIO(bz2.decompress(content))
                hdulist=fits.open(fitsfile,memmap=False)
                self.header = hdulist['PRIMARY'].header
                self.data = hdulist['PRIMARY'].data
                print('shape ', np.shape(self.data))
                hdulist.close()
                self.wcs = WCS(self.header)
                h1 = self.wcs.to_header()
                self.crota2 = np.arctan2(-h1["PC2_1"], h1["PC2_2"]) * 180./np.pi
                print('rotation angle ', self.crota2)
            except:
                print('Problems with retrieving file.')
        else:
            self.data = None
            self.wcs = None
            print('Coordinates out of the SDSS survey')
  

        
            
