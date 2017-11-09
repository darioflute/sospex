import urllib, urllib.request
from io import StringIO,BytesIO
from io import StringIO,BytesIO
from astropy.io import ascii,fits
from astropy.utils.data import download_file
from astropy.wcs import WCS
import numpy as np


class cloudImage(object):
    def __init__(self, lon, lat, xsize, ysize, source):

        self.lon = lon
        self.lat = lat
        self.xsize = xsize
        self.ysize = ysize
        self.source = source

        if source == 'wise1':
            image_file = self.downloadWise(1)
        elif source == 'wise2':
            image_file = self.downloadWise(2)
        elif source == 'wise3':
            image_file = self.downloadWise(3)
        elif source == 'wise4':
            image_file = self.downloadWise(4)
        else:
            print('Source not supported')

        # possible contours
        self.contours = None
        
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

        # Then, using the columns of the table, we can get an image
        params = { 'coadd_id': coadd[0],'band': 1,}
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
