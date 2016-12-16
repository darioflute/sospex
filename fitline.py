import numpy as np
from lmfit import Parameters, minimize
from lmfit.models import QuadraticModel


class fittedLine (object):
    """ fitted lines + continuum """
    def __init__(self,c1,c2,w,f,ef,ellipse,center):
        # limits of the continuum
        a = c1[0]
        b = c2[1]
        # Raw data 
        self.wave = w[a:b]
        self.flux = f[a:b]
        self.eflux = ef[a:b]
        self.sig = (w[c2[0]]-w[c1[1]])/4.
        self.wavec = np.append( w[c1[0]:c1[1]], w[c2[0]:c2[1]])
        self.fluxc = np.append(f[c1[0]:c1[1]], f[c2[0]:c2[1]])
        self.efluxc = np.append(ef[c1[0]:c1[1]], ef[c2[0]:c2[1]])

        # Fit results
        self.continuum=None
        self.continuumPars=None
        self.lines=None
        self.linesPars=None
        self.ellipse = ellipse # ellipse used to compute flux
        self.center = center # center of ellipse in world coordinates
    # Methods to compute line intensity, FWHM, EW, ...
        


def fitContinuum(fit):
    """
    Fit a 2nd deg polynomium to the continuum values in the two shaded regions around the line
    """

    # Model
    cont_mod = QuadraticModel()
    # Fit
    pars = cont_mod.guess(fit.fluxc,x=fit.wavec)
    out = cont_mod.fit(fit.fluxc,pars,x=fit.wavec,weights=fit.efluxc)
    fitc = cont_mod.eval(out.params, x=fit.wave)
    fit.continuum = fitc
    fit.continuumPars = out.params
    # Output
    return fit
    
def residuals(p,x,data=None,eps=None):
    # unpack parameters:
    #  extract .value attribute for each parameter
    v = p.valuesdict()
    model = x*0.
    n = len(v)/3
    for i in range(n):
        li = 'l'+str(i+1)+'_'
        lc = v[li+'center']
        la = v[li+'amplitude']
        ls = v[li+'sigma']
        model += la/(np.sqrt(2*np.pi)*ls)*np.exp(-(x-lc)*(x-lc)*0.5/(ls*ls))

    if data is None:
        return model
    else:
        if eps is None:
            return (model - data)
        else:
            return (model - data)/eps



def fitLines(fit, xpeak, ypeak):
    """
    Fit the lines after continuum subtraction
    """
    print "Fitting the lines ... "
    for i in range(len(xpeak)):
        print "Lines: ", i, xpeak[i], ypeak[i]

            
    # Define parameters
    fit_params = Parameters()
    # guess of sigma based on wavelength
    sig = fit.sig
    
    for i in range(len(xpeak)):
        li = 'l'+str(i+1)+'_'
        xp=xpeak[i]
        idx = np.argmin(np.abs(fit.wave-xp))
        fit_params.add(li+'center', value=xp,min=(xp-0.02), max=(xp+0.02))
        dy = (fit.flux[idx]-fit.continuum[idx])*np.sqrt(2.*np.pi)*sig
        print "dy ",dy
        if dy > 0:
            fit_params.add(li+'amplitude', value=dy, min = (dy*0.5),max=(dy*2))
        else:
            fit_params.add(li+'amplitude', value=dy, max = (dy*0.5),min=(dy*2))
        fit_params.add(li+'sigma', value=sig,min=sig*0.25,max=sig*1.5)

        x = fit.wave #/ (1.+self.galaxies[self.ngal].z)    
        y = fit.flux-fit.continuum
        e = fit.eflux  

        # output these data in a file
        #values = np.vstack((x, y)).T
        #np.savetxt('spectrum.txt',values,delimiter = "\t")

        out = minimize(residuals, fit_params, args=(x,), kws={'data':y,'eps':e},method='Nelder')
        yfit = residuals(out.params, x)
        fit.lines = yfit
        fit.linesPars = out.params
        
        # Now identify lines and store them in dictionary
        par = out.params.valuesdict()
        print "parameters: ", par
        nlines = len(par)/3
        print "Number of lines fitted: ", nlines
#        self.panel2.displayUncLineFit = True
    return fit

    
#def saveFit (self):
#    """
#    Save the fit (continuum, line, ellipse values in the header)
#    Continuum saved as an extension, line as another extension, wavelenght as 3rd extension
#    Also, wavelenght should be saved in the header (normal way to guarantee other displays will work)
#    """
#    print "Saving the fit in a FITS file ..."
