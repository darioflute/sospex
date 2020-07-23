import json
import io
import numpy as np
from collections import OrderedDict
from PyQt5.QtWidgets import QFileDialog

class MyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(MyEncoder, self).default(obj)
        
def exportContours(self):
    """Export defined contours."""
    itab = self.itabs.currentIndex()
    ic0 = self.ici[itab]
    ih0 = self.ihi[itab]
    if ic0.contour is not None:
        if ic0.title in ['Flux','Coverage map','Flux [no atm. corr.]']:
            source = self.specCube.instrument
        else:
            source = ic0.title
        info = [
                ('source', source),
                ('levels', ih0.levels)
                ]
        data = OrderedDict(info)
        # Open a dialog
        fd = QFileDialog()
        fd.setLabelText(QFileDialog.Accept, "Export as")
        fd.setNameFilters(["Json Files (*.json)","All Files (*)"])
        fd.setOptions(QFileDialog.DontUseNativeDialog)
        fd.setViewMode(QFileDialog.List)
        if (fd.exec()):
            filenames= fd.selectedFiles()
            filename = filenames[0]
            if filename[-5:] != '.json':
                filename += '.json'               
            print("Exporting contour levels to file: ", filename)
            with io.open(filename, mode='w') as f:
                str_= json.dumps(data, indent=2, separators=(',', ': '),
                                 ensure_ascii=False, cls=MyEncoder)
                f.write(str_)
            self.sb.showMessage("Aperture exported in file "+filename, 3000)
    else:
        return

    
def importContours(self):
    """Import defined contours.""" 
    # Open a dialog
    fd = QFileDialog()
    fd.setLabelText(QFileDialog.Accept, "Import")
    fd.setNameFilters(["Json Files (*.json)","All Files (*)"])
    fd.setOptions(QFileDialog.DontUseNativeDialog)
    fd.setViewMode(QFileDialog.List)
    fd.setFileMode(QFileDialog.ExistingFile)
    if (fd.exec()):
        filenames= fd.selectedFiles()
        print("Loading contour levels from file: ", filenames[0])
        with open(filenames[0],'r') as f:
            data = json.load(f)
    if self.contours == 'on':
        self.removeContours()
    self.drawContours(data['levels']) 
    self.contours = 'on'
    

def computeAreaPolygon(verts):
    """ 
    Compute area of polygon as sum of trapezoids.
    Coordinates of vertices are passed in degrees.
    The output area is expressed in square arcsec.
    """
    x = verts[:,0] * np.cos(verts[:,1] * np.pi / 180.) * 3600
    y = verts[:,1] * 3600
    x -= np.nanmin(x)
    y -= np.nanmin(y)
    area = 0
    x_ = x[0]
    y_ = y[0]
    for xi, yi in zip(x, y):
        area += (xi - x_) * (yi + y_) * 0.5
        x_ = xi
        y_ = yi
    return area


def exportAperture(self):
    """Export an aperture in Json format."""
    # Check if tab with aperture is open
    istab = self.stabs.currentIndex()
    instrument = self.specCube.instrument
    if istab > 1:
        nap = istab-1
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        istab = self.stabs.currentIndex()
        sc = self.sci[istab]        
        aperture = ic.photApertures[nap]
        type = aperture.type
        info = [
                ('waveUnit', 'micrometers'),
                ('amplitudeUnit', 'Jy'),
                ('fluxUnit', 'W/m2'),
                ('raUnit', 'deg'),
                ('decUnit', 'deg'),
                ('areaUnit', 'sq arcsec'),
                ('type', aperture.type),
                ('redshift', self.specCube.redshift)
                ]
        if type == 'Polygon':
            verts = aperture.poly.get_xy()
            adverts = np.array([(ic.wcs.wcs_pix2world(x,y,0)) for (x,y) in verts])
            area = computeAreaPolygon(adverts)
            info.extend([
                ('area', area),
                ('verts', adverts.tolist())
                ])
        elif type in ['Square', 'Rectangle']:
            x0,y0 = aperture.rect.get_xy()
            r0,d0 = ic.wcs.wcs_pix2world(x0,y0,0)
            width = aperture.rect.get_width() * ic.pixscale
            height = aperture.rect.get_height() * ic.pixscale
            area = width * height  # area in sq arcsec
            info.extend([
                ('area', area),
                ('width', width),
                ('height', height),
                ('angle', aperture.rect.angle - ic.crota2),
                ('ra0', r0.tolist()),
                ('dec0', d0.tolist())
                ])
        elif type in['Ellipse', 'Circle']:
            x0,y0 = aperture.ellipse.center
            r0,d0 = ic.wcs.wcs_pix2world(x0,y0,0)
            ax1 = aperture.ellipse.width * ic.pixscale * 0.5
            ax2 = aperture.ellipse.height * ic.pixscale * 0.5
            area = np.pi * ax1 * ax2  # area in sq arcsec
            info.extend([
                    ('area', area),
                    ('width',  ax1 * 2),
                    ('height', ax2 * 2),
                    ('angle',  aperture.ellipse.angle - ic.crota2),
                    ('ra0', r0.tolist()),
                    ('dec0', d0.tolist())
                    ])
        else:
            data = {}
        # Add information about line guesses and fitted parameters 
        if type in ['Polygon', 'Square', 'Rectangle', 'Ellipse', 'Circle']:
            # Add guess of continuum
            if sc.guess is not None:            
                g = sc.guess
                xg, yg = zip(*g.xy)
                info.append(('xg', xg))
                info.append(('yg', yg))
            # Add guesses of lines
            if sc.lguess is not None:
                lineguesses = []
                for line in sc.lines:
                    lineguesses.append([line.x0, line.fwhm, line.A])
                info.append(('lines', lineguesses))
            try:
                nlines = len(sc.aplines)
                info.append(('nlines',nlines))
                info.append(('model', sc.function))
                data = OrderedDict(info)
                # Add fitted parameters for each line
                for i, apline in enumerate(sc.aplines):
                    if sc.function == 'Gaussian':
                        c0, ec0, slope, x, ex, A, eA, sigma, esigma = apline
                        # Compute FWHM
                        FWHM = 2 * np.sqrt(2*np.log(2)) * sigma
                        eFWHM = 2 * np.sqrt(2*np.log(2)) * esigma
                        c = 299792458. # m/s
                        FWHMv = c * FWHM / x / 1000.
                        eFWHMv = FWHMv / sigma * esigma
                        # Compute intensity of line in W/m2
                        # The flux is integrated in um, so has to be divided
                        # by um, then c/x is expressed in m. 1.e-26 factor to
                        # transform Jy in W/m2/Hz
                        #jy2wm2 = c / (x * x) * 1.e-20 
                        flux = A #* jy2wm2
                        eflux = eA #* jy2wm2
                        A = A * 1.e20  * x * x / c
                        eA = eA * 1.e20  * x * x / c
                        line = 'line '+str(i+1)
                        data[line] = {
                                'instrument': instrument,
                                'continuum [Jy]': c0,
                                'errContinuum [Jy]': ec0,                                
                                'slope of continuum': slope,
                                'center [um]': x,
                                'errCenter [um]': ex,
                                'amplitude [Jy]': A,
                                'errAmplitude [Jy]': eA,
                                'sigma [um]': sigma,
                                'errSigma [um]': esigma,
                                'FWHM [um]': FWHM,
                                'errFWHM [um]': eFWHM,
                                'FWHM [km/s]': FWHMv,
                                'errFWHM [km/s]': eFWHMv,
                                'Flux [W/m2]': flux,
                                'errFlux [W/m2]': eflux
                                }
                    else:
                        c0, ec0, slope, x, ex, A, eA, sigma, esigma, alpha = apline
                        # Compute FWHM
                        FWHM = 2 * np.sqrt(2*np.log(2)) * sigma
                        eFWHM = 2 * np.sqrt(2*np.log(2)) * esigma
                        c = 299792458. # m/s
                        FWHMv = c * FWHM / x / 1000.
                        eFWHMv = FWHMv / sigma * esigma
                        # Compute intensity of line in W/m2
                        #jy2wm2 = c / (x * x) * 1.e-20 
                        flux = A #* jy2wm2
                        eflux = eA #* jy2wm2   
                        A = A * 1.e20  * x * x / c
                        eA = eA * 1.e20  * x * x / c
                        line = 'line '+str(i+1)
                        data[line] = {
                                'instrument': instrument,
                                'continuum [Jy]': c0,
                                'errContinuum [Jy]': ec0,
                                'slope of continuum': slope,
                                'center [um]': x,
                                'errCenter [um]': ex,
                                'amplitude [Jy]': A,
                                'errAmplitude [Jy]': eA,
                                'sigma [um]': sigma,
                                'errSigma [um]': esigma,
                                'alpha': alpha,
                                'FWHM [um]': FWHM,
                                'errFWHM [um]': eFWHM,
                                'FWHM [km/s]': FWHMv,
                                'errFWHM [km/s]': eFWHMv,
                                'Flux [W/m2]': flux,
                                'errFlux [W/m2]': eflux
                                }                       
            except:
                info.append(('nlines', 0))
                data = OrderedDict(info)
        # Open a dialog
        fd = QFileDialog()
        fd.setWindowTitle('Export aperture/fit')
        fd.setLabelText(QFileDialog.Accept, "Export as")
        fd.setNameFilters(["Json Files (*.json)","All Files (*)"])
        fd.setOptions(QFileDialog.DontUseNativeDialog)
        fd.setViewMode(QFileDialog.List)
        if (fd.exec()):
            filenames= fd.selectedFiles()
            filename = filenames[0]
            if filename[-5:] != '.json':
                filename += '.json'               
            print("Exporting aperture to file: ", filename)
            with io.open(filename, mode='w') as f:
                str_= json.dumps(data, indent=2, separators=(',', ': '),
                                 ensure_ascii=False, cls=MyEncoder)
                f.write(str_)
            self.sb.showMessage("Aperture exported in file "+filename, 3000)
    else:
        self.sb.showMessage("To export an aperture, select the tab with the desired aperture ",
                            3000)

def importAperture(self):
    """Import an aperture defined in a Json file."""
    from sospex.interactors import SegmentsInteractor, InteractorManager
    # Open a dialog
    fd = QFileDialog()
    fd.setWindowTitle('Import aperture/fit')
    fd.setLabelText(QFileDialog.Accept, "Import")
    fd.setNameFilters(["Json Files (*.json)","All Files (*)"])
    fd.setOptions(QFileDialog.DontUseNativeDialog)
    fd.setViewMode(QFileDialog.List)
    fd.setFileMode(QFileDialog.ExistingFile)
    if (fd.exec()):
        filenames= fd.selectedFiles()
        print("Loading aperture from file: ", filenames[0])
        with open(filenames[0],'r') as f:
            data = json.load(f)
        # Decoding info and opening new tab
        try:
            type = data['type']
            if type == 'Polygon':
                adverts = data['verts']
                self.drawNewPolygonAperture(adverts)
            else:
                itab = self.itabs.currentIndex()
                ic = self.ici[itab]
                r0 = data['ra0']
                d0 = data['dec0']
                w = data['width']
                h = data['height']
                angle = data['angle'] + ic.crota2
                self.drawNewAperture(type,r0,d0,w,h,angle)
            # Import lines
            if data['nlines'] > 0:
                istab = self.stabs.currentIndex()
                sc = self.sci[istab]
                sc.function = data['model']
                # Update redshift
                sc.spectrum.redshift = data['redshift']
                self.specCube.redshift = data['redshift']
                # Draw segments
                x = data['xg']
                y = data['yg']
                lines = data['lines']
                sc.xguess = [x]
                # Add lines
                if len(lines) > 0:
                    sc.lguess = lines
                # Plot continuum guess
                if y[3] != y[0]:
                    self.zeroDeg = False
                else:
                    self.zeroDeg = True
                xy = [(i,j) for (i,j) in zip(x,y)]
                sc.guess = SegmentsInteractor(sc.axes, xy, self.zeroDeg)
                sc.guess.modSignal.connect(self.onModifiedGuess)
                sc.guess.mySignal.connect(self.onRemoveContinuum)
                interactors = [sc.guess]
                print('Plotted segment interactor')
                # Plot lines
                if sc.lguess is not None:
                    sc.emslines = 0
                    sc.abslines = 0
                    ax0s = []
                    ex0s = []
                    afwhms = []
                    efwhms = []
                    aAs = []
                    eAs = []
                    for line in lines:
                        if line[2] >= 0:
                            sc.emslines += 1
                            ex0s.append(line[0]/( 1 + sc.spectrum.redshift))
                            efwhms.append(line[1])
                            eAs.append(line[2])
                        else:
                            sc.abslines += 1
                            ax0s.append(line[0]/( 1 + sc.spectrum.redshift))
                            afwhms.append(line[1])
                            aAs.append(line[2])
                    if sc.emslines > 0:
                        print('There are ', sc.emslines, ' emission lines')
                        ex0s = np.array(ex0s)
                        efwhms = np.array(efwhms)
                        eAs = np.array(eAs)
                        sc.lines = self.addLines(sc.emslines, x, 'emission', 
                                                  x0s=ex0s, fwhms=efwhms, As=eAs)
                    else:
                        sc.lines = []
                    if sc.abslines > 0:
                        print('There are ', sc.abslines, ' absorption lines')
                        ax0s = np.array(ax0s)
                        afwhms = np.array(afwhms)
                        aAs = np.array(aAs)
                        sc.lines.extend(self.addLines(sc.abslines, x, 'absorption',
                                                  nstart=sc.emslines, 
                                                  x0s=ax0s, fwhms=afwhms, As=aAs))
                else:
                    sc.lines = []
                # Plot guesses
                interactors.extend(sc.lines)
                sc.interactorManager = InteractorManager(sc.axes, interactors)
                # 
                sc.drawSpectrum()
                # Do the fit
                self.fitApLines()
        except:
            self.sb.showMessage("The file is not a valide aperture file.", 3000)
    else:
        self.sb.showMessage("To import an aperture, select a tab ", 3000)
        
def exportGuesses(self):
    "Export tessellation and guesses of continuum and lines."
    istab = self.spectra.index('Pix')
    sc = self.sci[istab]
    if sc.guess is None:
        return
    try:
        model = sc.model
    except:
        model = ''
    info = [('ncells', self.ncells),
            ('waveUnit', 'micrometers'),
            ('fluxUnit', 'Jy/pixel'),
            ('raUnit', 'deg'),
            ('decUnit', 'deg'),
            ('redshift', self.specCube.redshift),
            ('wavref', self.specCube.l0),
            ('kernel', self.kernel),
            ('model', model)
            ]
    if self.ncells == 1:
        xg,yg = zip(*sc.guess.xy)
        info.append(('x', xg))
        info.append(('y', yg))
        if sc.lguess is not None:
            info.append(('lines', sc.lguess))                
        data = OrderedDict(info)  
    elif self.ncells > 1:
        # Record sites
        imtab = self.bands.index('Flux')
        ic = self.ici[imtab]
        x, y = zip(*self.sites)
        #print("sites ", self.sites)
        ra, dec = ic.wcs.wcs_pix2world(x, y, 0)
        data = OrderedDict(info)
        # Check deleted lines
        idx = [i for i, x in enumerate(sc.lines) if x == None]
        idx.reverse()
        g = sc.guess
        xg, yg = zip(*g.xy)
        print('no of cells ', self.ncells)
        for i in range(self.ncells):
            xg = sc.xguess[i]
            lines = []
            if sc.lguess is not None:
                for line in sc.lguess:
                    lines.append(line[i])
                # Remove deleted line guesses
                for j in idx:
                    del lines[j]
            data[i] = {
                    'ra': ra[i],
                    'dec': dec[i],
                    'x': xg,
                    'y': yg,
                    'lines': lines
                    }
            data.move_to_end(i, last=True)  # Move element to the end
    # Open a dialog
    fd = QFileDialog()
    fd.setWindowTitle('Export guesses')
    fd.setLabelText(QFileDialog.Accept, "Export as")
    fd.setNameFilters(["Json Files (*.json)","All Files (*)"])
    fd.setOptions(QFileDialog.DontUseNativeDialog)
    fd.setViewMode(QFileDialog.List)
    if (fd.exec()):
        filenames= fd.selectedFiles()
        filename = filenames[0]
        if filename[-5:] != '.json':
            filename += '.json'              
        print("Exporting guesses to file: ", filename)
        with io.open(filename, mode='w') as f:
            str_= json.dumps(data, indent=2, separators=(',', ': '),
                             ensure_ascii=False, cls=MyEncoder)
            f.write(str_)
        self.sb.showMessage("Guesses exported to file "+filename, 3000)
    
def importGuesses(self):
    "Import previously defined tessellation and guesses of continuum and lines."
    from sospex.interactors import SegmentsInteractor, InteractorManager
    # Open a dialog
    fd = QFileDialog()
    fd.setWindowTitle('Import guesses')
    fd.setLabelText(QFileDialog.Accept, "Import as")
    fd.setNameFilters(["Json Files (*.json)","All Files (*)"])
    fd.setOptions(QFileDialog.DontUseNativeDialog)
    fd.setViewMode(QFileDialog.List)
    if (fd.exec()):
        filenames= fd.selectedFiles()
        filename = filenames[0]
        with open(filename) as f:
            data = json.load(f, object_pairs_hook=OrderedDict)
        self.ncells = data['ncells']
        #stab = self.stabs.currentIndex()        
        istab = self.spectra.index('Pix')
        sc = self.sci[istab]
        self.stabs.setCurrentIndex(istab)
        sc.spectrum.redshift = data['redshift']
        sc.spectrum.l0 = data['wavref']
        model = data['model']
        if model != '':
            sc.model = model
        # Update 
        self.onDraw2(1)
        imtab = self.bands.index('Flux')
        ic = self.ici[imtab]
        self.itabs.setCurrentIndex(imtab)
        # Delete previous guess and select new one
        if sc.guess is not None:
            self.onRemoveContinuum('segments deleted')
        if self.ncells == 1:
            x = data['x']
            y = data['y']
            lines = data['lines'][0]
            sc.xguess = [x]
            # Add lines
            if len(lines) > 0:
                sc.lguess = lines
        else:
            # Case with tessellation
            ra = []
            dec = []
            sc.lguess = []
            sc.xguess = []
            for i in range(self.ncells):
                g = data[str(i)]
                sc.xguess.append(g['x'])
                ra.append(g['ra'])
                dec.append(g['dec'])
                sc.lguess.append(g['lines'])
            sc.lguess = [list(i) for i in zip(*sc.lguess)]  # Transposing
            x, y = ic.wcs.wcs_world2pix(ra, dec, 0)
            # Round to 2 decimal figures to avoid crazy Voronoi vertices values
            self.sites = [(round(i,2),round(j,2)) for (i,j) in zip(x,y)]
            #print('sites in ', self.sites)
            x = data['0']['x']
            y = data['0']['y']
            lines = data['0']['lines']
        # Plot continuum guess
        if y[3] != y[0]:
            self.zeroDeg = False
        else:
            self.zeroDeg = True
        xy = [(i,j) for (i,j) in zip(x,y)]
        sc.guess = SegmentsInteractor(sc.axes, xy, self.zeroDeg)
        sc.guess.modSignal.connect(self.onModifiedGuess)
        sc.guess.mySignal.connect(self.onRemoveContinuum)
        interactors = [sc.guess]
        # Initialize continuum
        self.openContinuumTab()
        self.positiveContinuum = False  # Default value
        # Plot lines
        if sc.lguess is not None:
            sc.emslines = 0
            sc.abslines = 0
            for line in lines:
                if line[2] >= 0:
                    sc.emslines += 1
                else:
                    sc.abslines += 1
            if sc.emslines > 0:
                sc.lines = self.addLines(sc.emslines, x, 'emission', x0s=line[0],
                                         fwhms=line[1], As=line[2])
            else:
                sc.lines = []
            if sc.abslines > 0:
                sc.lines.append(self.addLines(sc.abslines, x, 'absorption', nstart=sc.emslines, 
                                              x0s=line[0], fwhms=line[1], As=line[2]))
        # Then plot guess and activate tessellation
        interactors.extend(sc.lines)
        sc.interactorManager = InteractorManager(sc.axes, interactors)
        #sc.fig.canvas.draw_idle()
        # Create Voronoi sites, KDTree, plot Voronoi ridges on image
        self.removeVI()
        if self.ncells > 1:
            self.createVI()
        else:
            ic.fig.canvas.draw_idle()
        # Set kernel
        self.setKernel(data['kernel'])
    else:
        print('Try again after exporting a set of guesses.')