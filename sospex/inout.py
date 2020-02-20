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


def exportAperture(self):
    """Export an aperture in Json format."""
    # Check if tab with aperture is open
    istab = self.stabs.currentIndex()
    if istab > 1:
        nap = istab-1
        itab = self.itabs.currentIndex()
        ic = self.ici[itab]
        aperture = ic.photApertures[nap]
        type = aperture.type
        print('type is ', type)
        print('rotation angle ', ic.crota2)
        if type == 'Polygon':
            verts = aperture.poly.get_xy()
            adverts = np.array([(ic.wcs.wcs_pix2world(x,y,0)) for (x,y) in verts])                
            data = {
                'type': aperture.type,
                'verts': adverts.tolist()
            }
        elif (type == 'Square') | (type == 'Rectangle'):
            x0,y0 = aperture.rect.get_xy()
            r0,d0 = ic.wcs.wcs_pix2world(x0,y0,0)
            data = {
                'type': aperture.type,
                'width': aperture.rect.get_width()*ic.pixscale,
                'height': aperture.rect.get_height()*ic.pixscale,
                'angle': aperture.rect.angle - ic.crota2,
                'ra0': r0.tolist(),
                'dec0': d0.tolist()
                }
        elif  (type == 'Ellipse') | (type == 'Circle'):
            x0,y0 = aperture.ellipse.center
            r0,d0 = ic.wcs.wcs_pix2world(x0,y0,0)
            data = {
                    'type': aperture.type,
                    'width':  aperture.ellipse.width*ic.pixscale,
                    'height': aperture.ellipse.height*ic.pixscale,
                    'angle':  aperture.ellipse.angle - ic.crota2,
                    'ra0': r0.tolist(),
                    'dec0': d0.tolist()
                    }
        else:
            data = {}
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
            print("Exporting aperture to file: ", filename)
            with io.open(filename, mode='w') as f:
                str_= json.dumps(data,indent=2,sort_keys=True,separators=(',',': '),
                                 ensure_ascii=False)
                # print(str_)
                f.write(str_)
            self.sb.showMessage("Aperture exported in file "+filename, 3000)
    else:
        self.sb.showMessage("To export an aperture, select the tab with the desired aperture ",
                            3000)

def importAperture(self):
    """Import an aperture defined in a Json file."""
    # Open a dialog
    fd = QFileDialog()
    fd.setLabelText(QFileDialog.Accept, "Import")
    fd.setNameFilters(["Json Files (*.json)","All Files (*)"])
    fd.setOptions(QFileDialog.DontUseNativeDialog)
    fd.setViewMode(QFileDialog.List)
    fd.setFileMode(QFileDialog.ExistingFile)
    if (fd.exec()):
        filenames= fd.selectedFiles()
        print("Loading aperture from file: ", filenames[0])
        with open(filenames[0],'r') as f:
            data = json.load(f)#, object_pairs_hook=OrderedDict)
        # Decoding info and opening new tab
        try:
            type = data['type']
            if type == 'Polygon':
                adverts = data['verts']
                #verts = np.array([(ic.wcs.all_world2pix(r,d,1)) for (r,d) in adverts])       
                self.drawNewPolygonAperture(adverts)
            else:
                itab = self.itabs.currentIndex()
                ic = self.ici[itab]
                print('rotation angle ', ic.crota2)
                r0 = data['ra0']
                d0 = data['dec0']
                w = data['width']
                h = data['height']
                angle = data['angle'] + ic.crota2
                self.drawNewAperture(type,r0,d0,w,h,angle)
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
    info = [('ncells', self.ncells),
            ('waveUnit', 'micrometers'),
            ('fluxUnit', 'Jy/pixel'),
            ('raUnit', 'deg'),
            ('decUnit', 'deg'),
            ('redshift', self.specCube.redshift),
            ('wavref', self.specCube.l0),
            ('kernel', self.kernel)
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
            self.emslines = 0
            self.abslines = 0
            for line in lines:
                if line[2] >= 0:
                    self.emslines += 1
                else:
                    self.abslines += 1
            if self.emslines > 0:
                sc.lines = self.addLines(self.emslines, x, 'emission', x0s=line[0], fwhms=line[1], As=line[2])
            else:
                sc.lines = []
            if self.abslines > 0:
                sc.lines.append(self.addLines(self.abslines, x, 'absorption', nstart=self.emslines, 
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