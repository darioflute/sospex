import urllib, urllib.request
from io import StringIO,BytesIO
from astropy.io import ascii,fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.utils.data import download_file
from astropy.wcs import WCS
import numpy as np
from html.parser import HTMLParser
from PyQt5.QtWidgets import QFileDialog

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
    def __init__(self, lon, lat, xsize, ysize, source):

        self.lon = lon
        self.lat = lat
        self.xsize = xsize
        self.ysize = ysize
        self.source = source
        self.data = None
        self.wcs = None

        if source == 'local':
            image_file = self.openLocal()
        elif source == 'wise1':
            image_file = self.downloadWise(1)
        elif source == 'wise2':
            image_file = self.downloadWise(2)
        elif source == 'wise3':
            image_file = self.downloadWise(3)
        elif source == 'wise4':
            image_file = self.downloadWise(4)
        elif source == '2mass-j':
            image_file = self.download2MASS('J')
        elif source == '2mass-h':
            image_file = self.download2MASS('H')
        elif source == '2mass-k':
            image_file = self.download2MASS('K')
        elif source == 'panstarrs-g':
            image_file = self.downloadPanSTARRS('g')
        elif source == 'panstarrs-r':
            image_file = self.downloadPanSTARRS('r')
        elif source == 'panstarrs-i':
            image_file = self.downloadPanSTARRS('i')
        elif source == 'panstarrs-z':
            image_file = self.downloadPanSTARRS('z')
        elif source == 'panstarrs-y':
            image_file = self.downloadPanSTARRS('y')
        elif source == 'sdss-u':
            image_file = self.downloadSDSS('u')
        elif source == 'sdss-r':
            image_file = self.downloadSDSS('g')
        elif source == 'sdss-g':
            image_file = self.downloadSDSS('r')
        elif source == 'sdss-i':
            image_file = self.downloadSDSS('i')
        elif source == 'sdss-z':
            image_file = self.downloadSDSS('z')
        elif source == 'first':
            image_file = self.downloadFIRST()
        elif source == 'nvss':
            image_file = self.downloadNVSS()
        elif source == 'sumss':
            image_file = self.downloadSUMSS()
        else:
            print('Source not supported')

        # possible contours
        self.contours = None

    def download2MASS(self,band):

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
        data= parser.data
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
            hdulist = fits.open(fitsfile)
            header = hdulist['PRIMARY'].header
            self.data = hdulist['PRIMARY'].data
            hdulist.close()
            self.wcs = WCS(header).celestial
        else:
            self.data = None
            self.wcs = None
            print('Coordinates out of the 2MASS survey')


    def openLocal(self):
        """ Open a local FITS / check if coordinates fall into the FITS """

        # Open a dialog
        fd = QFileDialog()
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
                hdulist = fits.open(image_file)
                hdulist.info()
                header = hdulist['PRIMARY'].header
                self.data = hdulist['PRIMARY'].data
                hdulist.close()
                print(header)
                self.wcs = WCS(header)#.celestial
                print('wcs is ',self.wcs)
                # Check if coordinates are inside the image
                print('data shape is',np.shape(self.data))
                x,y = self.wcs.all_world2pix(self.lon,self.lat,1)
                print('x y ',x,y)
                ny,nx = np.shape(self.data)
                if x >= 0 and x< nx and y >= 0 and y  <= ny:
                    print('Source inside the FITS image')
                else:
                    self.data = None
                    self.wcs = None
                    print('Coordinates out of the selected FITS image')
            except:
                self.data = None
                self.wcs = None
                print('The selected  FITS is not a valid file')


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

        hdulist = fits.open(image_file)
        #hdulist.info()
        header = hdulist['PRIMARY'].header
        self.data = hdulist['PRIMARY'].data
        hdulist.close()
        self.wcs = WCS(header)


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
            hdulist = fits.open(fitsfile)
            header = hdulist['PRIMARY'].header
            self.data = hdulist['PRIMARY'].data
            hdulist.close()
            self.wcs = WCS(header).celestial
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
            hdulist = fits.open(fitsfile)
            hdulist.info()
            header = hdulist['PRIMARY'].header
            self.data = hdulist['PRIMARY'].data[0,0,:,:]
            hdulist.close()
            self.wcs = WCS(header).celestial
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
        data= parser.data
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
            hdulist = fits.open(fitsfile)
            header = hdulist['PRIMARY'].header
            self.data = hdulist['PRIMARY'].data
            hdulist.close()
            self.wcs = WCS(header).celestial
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
        data= parser.data
        values= parser.values
        parser.close()

        file = None
        for v in values:
            if '.'+band+'.unconv.fits' in v:
                file = 'http:'+v
                
        if file is not None:
            #idx = data.index('FITS-cutout')
            #file = 'http:'+values[idx]
            response = urllib.request.urlopen(file)
            content= response.read()
            fitsfile = BytesIO(content)
            hdulist = fits.open(fitsfile)
            header = hdulist['PRIMARY'].header
            self.data = hdulist['PRIMARY'].data
            hdulist.close()
            self.wcs = WCS(header).celestial
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
        ras = "{:02.0f}+{:.0f}+{:.1f}".format(cra[0],cra[1],cra[2])
        if self.lat > 0:
            url = url0+ras+"+&dec=+{:02.0f}+{:.0f}+{:.0f}".format(cdec[0],cdec[1],cdec[2])
        else:
            url = url0+ras+"+&dec=-{:02.0f}+{:.0f}+{:.0f}".format(-cdec[0],-cdec[1],-cdec[2])

        request = urllib.request.Request(url)
        response = urllib.request.urlopen(request)
        html = response.read()
        
        parser = MyHTMLParser()
        parser.feed(html.decode('utf-8'))
        data= parser.data
        values= parser.values
        parser.close()

        file = None
        for v in values:
            if 'frame-'+band+'-' in v:
                file = 'https://dr12.sdss.org/'+v
                
        if file is not None:
            #idx = data.index(band+'-band FITS')
            #print("Value is:", values[idx-1])
            #file = 'https://dr12.sdss.org/'+values[idx-1]
            response = urllib.request.urlopen(file)
            content= response.read()        
            fitsfile = BytesIO(bz2.decompress(content))
            hdulist=fits.open(fitsfile)
            header = hdulist['PRIMARY'].header
            self.data = hdulist['PRIMARY'].data
            hdulist.close()
            self.wcs = WCS(header)
        else:
            self.data = None
            self.wcs = None
            print('Coordinates out of the SDSS survey')
  
    
        
            
