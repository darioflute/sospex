import numpy as np
import multiprocessing as mp
from lmfit import Parameters, minimize


def biweight(data, axis=None):
    """Biweight estimation of location and scale according to
    Understanding Robust and Exploratory Data Analysis, Hoaghley, Mosteller, & Tukey,
    formula 4. Tuning constant as at page 417. Location on page 421.
    See also Beers, Flynn, & Gebhardt, 1990, AJ 100:32.
    """
    import numpy as np
    np.seterr(divide='ignore', invalid='ignore')
    
    c1 = 6.0
    c2 = 9.0
    data = np.asarray(data)
    M = np.nanmedian(data, axis=axis, keepdims=True)
    S = np.nanmedian(np.abs(data-M), axis=axis, keepdims=True)
        
    if np.nanmean(S) < 0.0001:
        return M.squeeze(), S.squeeze()
    
    u = (data - M) / (c1 * S)
    mask = np.abs(u) > 1.0
    u = (1 - u * u) ** 2
    u[mask] = 0
    s3 = np.nansum((data - M) * u, axis=axis)
    s4 = np.nansum(u, axis=axis)
    if axis is not None:
        s34 = np.expand_dims(s3 / s4, axis=axis)
    else:
        s34 = s3 / s4
    M += s34
    
    u = (data - M ) / (c2 * S)
    mask = np.abs(u) >= 1.
    u = u * u
    u1 = 1 - u
    u1[mask] = 0 
    s1 = np.nansum((data - M)**2 * u1**4, axis=axis)
    s2 = np.nansum((1 - 5 * u) * u1, axis=axis)
    ndata = np.ma.size(data, axis=axis)
    S = np.sqrt( s1 / (ndata-1) ) / ( np.abs(s2) / ndata)
        
    return M.squeeze(), S

def weightedMedian(data, weights, axis=None):
    """
        Median of a weighted array.
        The function works also on a 2D or 3D array one an axis for median computation is defined.
    """
    import numpy as np
    np.seterr(divide='ignore', invalid='ignore')
    # Transform into numpy arrays and sort
    a = np.asarray(data)
    w = np.asarray(weights)
    ids = np.argsort(a, axis=axis)
    a = np.take_along_axis(a, ids, axis=axis)
    w = np.take_along_axis(w, ids, axis=axis)
    
    # Compute midpoint and find the median
    wcum = np.cumsum(w, axis=axis)
    wmid = np.sum(w, axis=axis, keepdims=True)*0.5
    mask = wcum <= wmid
    if axis is None:
        idx = np.sum(mask) - 1
        if np.max(w) > wmid:
            wmed = (a[w == np.max(w)])[0]
        else:
            if wcum[idx]  == wmid:
                wmed = (a[idx] + a[idx+1]) * 0.5
            else:
                wmed = a[idx+1]
    else:
        wmid = wmid.squeeze()
        idx = np.sum(mask, axis=axis, keepdims=True) - 1
        icum = np.take_along_axis(wcum, idx, axis=axis).squeeze()
        wmed = np.take_along_axis(a, idx+1, axis=axis).squeeze()
        w1 = np.take_along_axis(a, idx, axis=axis).squeeze()
        id = icum == wmid
        if np.sum(id) > 0:
            wmed[id] = (wmed[id]+w1[id])*0.5
        # Case of weight > midpoint
        id = np.nanmax(w, axis=axis) > wmid
        if np.sum(id) > 0:
            imax = np.expand_dims(np.nanargmax(w, axis=axis), axis=axis)
            amax = np.take_along_axis(a, imax, axis=axis)
            wmed[id] = amax[id]
        
    return wmed
    
        

# Functions for multiprocessing continuum fit and moment computation

def residuals(p,x,data=None,eps=None):
    #unpack parameters
    v = p.valuesdict()
    q = v['q']
    # define model
    try:
        m = v['m']
        model = m*x+q
    except:
        model = q
    if data is None:
        return model
    else:
        if eps is None:
            return model-data 
        else:
            return (model-data)/eps

def fitContinuum(p,slope,intcp,posCont,m,w,ff):
    # Take the mean for each wavelength
    f = np.nanmean(ff,axis=1)
    mf = np.isnan(f)
    m[mf] = 0
    if np.sum(m) > 5:
        # Define parameters
        fit_params = Parameters()
        if slope != 0:
            fit_params.add('m',value=slope)
        if posCont:
            fit_params.add('q',value=intcp, min=0.)
        else:
            fit_params.add('q',value=intcp)
        out = minimize(residuals,fit_params,args=(w[m],),kws={'data':f[m]},method='Nelder')
        pars = out.params
    else:
        pars = None
        pass
    return p, pars

def fiteContinuum(p,slope,intcp,posCont,m,w,ff,ee):
    # Skip if exposure is zero
    se = np.nansum(ee)
    if se == 0:
        pars = None
        return p, pars
    # Take the mean for each wavelength
    f = np.nanmean(ff,axis=1)
    mf = np.isnan(f)
    m[mf] = 0
    if np.sum(m) > 5:
        # Define parameters
        fit_params = Parameters()
        if slope != 0:
            fit_params.add('m',value=slope)
        if posCont:
            fit_params.add('q',value=intcp, min=0.)
        else:
            fit_params.add('q',value=intcp)
        e = 1./np.sqrt(np.nanmean(ee,axis=1))  # The error scale with the inverse of the sqrt of the exposure
        out = minimize(residuals,fit_params,args=(w[m],),kws={'data':f[m],'eps':e[m]},method='Nelder')
        pars = out.params
    else:
        pars = None
        pass
    return p, pars


def computeMoments(p,m,w,dw,f):
    """ compute moments on a spatial pixel """
    mf = np.isnan(f)
    m[mf] = 0
    # compute the error on f
    if np.sum(m) > 5:
        f = f[m]
        w = w[m]
        dw = dw[m]
        
        #sf = medfilt(f,5)
        #df = f-sf
        # The dispersion is computed using the difference
        # between two consecutive values. Assuming they have the
        # same dispersion, the dispersion of the difference is
        # sqrt(2) higher.
        pos = f[1:] > 0 # Consider only positive values
        df = f[1:]-f[:-1]
        df = df[pos]
        med = np.nanmedian(df)
        # Divide by sqrt(2.) since I did a subtraction 
        mad = np.nanmedian(np.abs(med))/np.sqrt(2.) * 1.48 # MAD for Gaussian distributions
        sigma = -5.0 *  mad# n sigma value
        # Consider only values greater than continuum
        #if np.isfinite(sigma):
        #    ms = f > (-sigma)
        #else:
        #    ms = np.isfinite(f)
        # ms = f > sigma
        ms = f > 0
        # Compute also negative intensity
        if np.sum(ms) > 5:
            c = 299792458. # m/s
            w_ = w[ms]
            dw_ = dw[ms]
            Snu = f[ms]
            pos = Snu > sigma
            Slambda = c*Snu[pos]/(w_[pos]*w_[pos])*1.e6   # [Jy * Hz / um]
            w_  = w_[pos]
            dw_ = dw_[pos]
            # Compute only for values > sigma
            M0 = np.nansum(Slambda*dw_) # [Jy Hz] 
            M1 = np.nansum(w_*Slambda*dw_)/M0 # [um]
            M2 = np.nansum(np.power(w_-M1,2)*Slambda*dw_)/M0 # [um*um]
            SD = np.sqrt(M2)
            M3 = np.nansum(np.power(w_-M1,3)*Slambda*dw_)/M0/np.power(SD,3)
            M4 = np.nansum(np.power(w_-M1,4)*Slambda*dw_)/M0/np.power(SD,4)-3. # Relative to Gaussian which is 3
            #w_ = w[ms]
            #dw_ = dw[ms]
            #Snu = f[ms]
            # Compute M0 only for positive values (or mildly negative)
            #pos = Snu > sigma
            #Slambda = c*Snu[pos]/(w_[pos]*w_[pos])*1.e6   # [Jy * Hz / um]
            #w_  = w_[pos]
            #dw_ = dw_[pos]
            #M0 = np.nansum(Slambda*dw_) # [Jy Hz]             
            M0 *= 1.e-26 # [W/m2]  (W/m2 = Jy*Hz*1.e-26)
        else:
            M0 = np.nan
            M1 = M2 = M3 = M4 = np.nan
    else:
        M0 = np.nan
        M1 = M2 = M3 = M4 = np.nan
        sigma = np.nan
            
    return p, M0, M1, M2, M3, M4, sigma

def multiComputeMoments(m,w,f,c,moments,points):

    # To avoid forking error in MAC OS-X
    try:
        mp.set_start_method('spawn')
        print('started spawing')
    except RuntimeError:
        pass

    # Define noise
    n3,n2,n1 = np.shape(moments)
    noise = np.zeros((n2,n1))

    
    # Compute dw
    dw = [] 
    dw.append([w[1]-w[0]])
    dw.append(list((w[2:]-w[:-2])*0.5))
    dw.append([w[-1]-w[-2]])
    dw = np.concatenate(dw)

    with mp.Pool(processes=mp.cpu_count()) as pool:
        res = [pool.apply_async(computeMoments, (p,m[:,p[1],p[0]],w,dw,f[:,p[1],p[0]]-c[:,p[1],p[0]])) for p in points]
        results = [r.get() for r in res]

        
    for p, M0, M1, M2, M3, M4, sigma in results:
        i,j = p
        moments[0][j,i] = M0
        moments[1][j,i] = M1
        moments[2][j,i] = M2
        moments[3][j,i] = M3
        moments[4][j,i] = M4
        noise[j,i] = sigma
            
    return moments, noise
    

def multiFitContinuum(m, w, f, c, c0, w0, points, slope, intcp, posCont, kernel, exp=None):
    # To avoid forking error in MAC OS-X
    try:
        mp.set_start_method('spawn')
        print('started spawing')
    except RuntimeError:
        pass
    if kernel == 1:
        ik = np.array([0])
        jk = np.array([0])
    elif kernel == 5:
        ik = np.array([-1,0,0,0,1])
        jk = np.array([0,-1,0,1,0])
    elif kernel == 9:
        ik = np.array([-1,-1,-1,0,0,0,1,1,1])
        jk = np.array([-1,0,1,-1,0,1,-1,0,1])
    else:
        print ('unsupported kernel, use one pixel only')
        ik = np.array([0])
        jk = np.array([0])
    if exp is None:
        with mp.Pool(processes=mp.cpu_count()) as pool:
            res = [pool.apply_async(fitContinuum, (p,slope,intcp,posCont,m[:,p[1],p[0]],w,
                                                   f[:,p[1]+ik,p[0]+jk])) for p in points]
            results = [r.get() for r in res]
    else:
        with mp.Pool(processes=mp.cpu_count()) as pool:
            res = [pool.apply_async(fiteContinuum, 
                                    (p, slope, intcp, posCont, m[:,p[1],p[0]], w,
                                     f[:,p[1]+ik,p[0]+jk],exp[:,p[1]+ik,p[0]+jk])) for p in points]
            results = [r.get() for r in res] 
    print('c is writeable ', c.flags)
    c_ = c.copy()
    c0_ = c0.copy()  # continuum at ref wav
    cs_ = c0.copy()  # slope of cont
    for p, pars in results:
        if pars is not None:
            i, j = p
            c_[:, j, i] = residuals(pars, w)
            c0_[j, i] = residuals(pars, w0)
            try:
                cs_[j, i] = pars['m'].value
            except:
                cs_[j, i] = 0.
    return c_, c0_, cs_

# Fit of lines

def fitLines(p, m, w, f, lines, model):
    """Fit the lines defined in the guess."""
    from lmfit.models import PseudoVoigtModel
    # f is flux without continuum
    # lines is the list of parameters for each line
    mf = np.isnan(f)
    m[mf] = 0
    if np.sum(m) > 5:
        y = f[m]
        x = w[m]
        # Transform into S(lambda)
        c = 299792458. # m/s
        #  W/m2 = Jy*Hz*1.e-26 and c -> c*1.e6 um ...
        y  = c * y / (x * x ) *1.e-20 # Jy * Hz / um --> W/m2 after integration in um
        
        # Find normalization
        amplitudes = []
        for line in lines:
            sigma = line[1] / 2.355
            x0 = line[0]
            A = line[2] * c / x0**2 * 1.e-20  # same units as flux
        amplitudes.append(A)
        amplitudes = np.array(amplitudes)
        norm = np.nanmax(np.abs(amplitudes))
        y /= norm  
        # Define lines
        for i, line in enumerate(lines):
            li = 'l' + str(i) + '_'
            vmodel = PseudoVoigtModel(prefix=li)
            if i == 0:
                params = vmodel.make_params()
                model = vmodel
            else:
                params += vmodel.make_params()
                model += vmodel
            x0 = line[0]
            sigma = line[1] / 2.355
            A = line[2] * c / (x0*x0) * 1.e-20  # same units as flux
            A /= norm
            params[li+'center'].set(x0, min=(x0 - sigma/2.), max=(x0 + sigma/2.))
            #params[li+'amplitude'].set(A)
            if A > 0:
                params[li+'amplitude'].set(A, min=0, max=1.1)
            elif A == 0:
                params[li+'amplitude'].set(A, vary=False)
            else:
                params[li+'amplitude'].set(A, max=-1.1, min=0)# min=2 * A)
                
            params[li+'sigma'].set(sigma, min=sigma / 2., max=sigma * 2)
            if model == 'Gauss':
                    params[li + 'fraction'].set(0.0, vary=False)
            else:
                params[li + 'fraction'].set(0.4, vary=True, min=0, max = 0.45)  # No Cauchy part (i.e. Gauss)
        # Minimize
        out = model.fit(y, params, x=x, method='leastsq')
        # Return lines fitted parameters
        pars = out.params#.valuesdict()
        nlines = len(lines)
        linepars = []
        for i in range(nlines):
            li = 'l' + str(i) + '_'
            center = pars[li + 'center'].value  # Observed
            centerErr = pars[li + 'center'].stderr  # Observed
            if centerErr is None:
                centerErr = np.nan
            sigma = pars[li + 'sigma'] .value   # Observed
            sigmaErr = pars[li + 'sigma'].stderr   # Observed
            if sigmaErr is None:
                sigmaErr = np.nan
            A =  pars[li+'amplitude'].value * norm
            Aerr = pars[li+'amplitude'].stderr
            if Aerr is None:
                Aerr = np.nan
            else:
                Aerr *= norm
            alpha = pars[li+'fraction'].value
            linepars.append([center, sigma, A, alpha, centerErr, sigmaErr, Aerr])
        return p, linepars
    else:
        nlines = len(lines)
        return p, np.full((nlines, 7), np.nan)

def multiFitLines(m, w, f, c, lineguesses, model, linefits, points):

    # To avoid forking error in MAC OS-X
    try:
        mp.set_start_method('spawn')
        print('started spawing')
    except RuntimeError:
        pass

    with mp.Pool(processes=mp.cpu_count()) as pool:
        res = [pool.apply_async(fitLines, 
                                (p, 
                                 m[:, p[1],p[0]], 
                                 w, 
                                 f[:,p[1],p[0]]-c[:,p[1],p[0]], 
                                 lineguesses, model)
                                ) for p in points]
        results = [r.get() for r in res]

    n = len(lineguesses)
    for p, linepars in results:
        i,j = p
        for k in range(n):
            for l in range(7):
                linefits[k][l][j,i] = linepars[k][l]
            
    return 1

def multiFitLinesSingle(m, w, f, c, lineguesses, model, linefits, points):

    for p in points:
        res = fitLines(p, 
                       m[:, p[1],p[0]], 
                       w, 
                       f[:,p[1],p[0]]-c[:,p[1],p[0]], 
                       lineguesses, model)

        n = len(lineguesses)
        pp, linepars = res
        i,j = pp
        #print(np.shape(linepars), np.shape(linefits))
        for k in range(n):
            for l in range(7):
                linefits[k][l][j,i] = linepars[k][l]
                
    return 1

def residualsPsf(p, x, y, data=None, err=None):
    '''Residual of a PSF with unknown center'''
    v = p.valuesdict()
    s = v['s']
    A = v['A']
    x0 = v['x0']
    y0 = v['y0']
    dis = np.hypot(x-x0,y-y0)/s
    d2 = np.square(dis.flatten())
    model = A *  np.exp(-0.5 * d2)
    
    if data is None:
        return model
    else:
        if err is None:
            return (model - data.flatten())
        else:
            return (model - data.flatten())/err.flatten()

def histoImage(image, percent, xmin=None, xmax=None):
    #ima = image.ravel()
    #mask = np.isfinite(ima)
    #ima = ima[mask]
    ima = image[np.isfinite(image)]
    nh = len(ima)
    if nh > 0:
        ima = np.sort(ima)
        # Take uniq to avoid repeated values (which could be border values)
        u, indices = np.unique(ima, return_index=True)
        ima = ima[indices]
        nh = len(ima)
        s = np.size(ima)
        if percent is None:
            percent = 99.0
        p1 = (100.-percent)/200.
        p2 = 1. - p1
        smin = int(s*p1)
        smax = min(int(s*p2)-1,s-1)
        nbins=256
        imedian = np.median(ima)
        sdev = np.nanmedian(np.abs(ima - np.nanmedian(ima))) * 1.4826 # Gauss distr
        imin = np.min(ima)
        imax = np.max(ima)
        epsilon = sdev/3.            
        # Define the interval containing 99% of the values
        if xmin == None:
            xmin = ima[smin]
        if xmax == None:
            xmax = ima[smax]
        smin = int(s*0.005)
        smax = min(int(s*0.995), s-1)
        hmin = ima[smin]
        hmax = ima[smax]
        # Avoid excessively lower flux
        if hmin > (imedian - 1.5 * sdev):
            hmin = imedian - 1.5 * sdev
        if hmax < (imedian + 6 * sdev):
            hmax = imedian + 6 * sdev        
    else:
        nbins = 0
        xmin = 0.
        xmax = 0.
        imedian = 0.
        sdev = 0.
        imin = 0.
        imax = 0.
        epsilon = 0.
        hmin = 0.
        hmax = 0.
        
    return ima, nbins, xmin, xmax, hmin, hmax, imedian, imin, imax, sdev, epsilon, nh
    
# Functions to fit lines inside an aperture
def contResiduals(p, x, data=None, eps=None):
    # unpack parameters:
    #  extract .value attribute for each parameter
    v = p.valuesdict()
    intcp = v['intercept']
    try:
        slope = v['slope']
        model = intcp + slope * x
    except BaseException:
        model = intcp
    if data is None:
        return model
    else:
        if eps is None:
            return (model - data)
        else:
            return (model - data) / eps


def linesGaussResiduals(p, x, data=None, eps=None):
    # unpack parameters:
    # extract .value attribute for each parameter
    v = p.valuesdict()
    model = 0
    n = len(v) // 3
    for i in range(n):
        li = 'l' + str(i) + '_'
        lc = v[li + 'center']
        la = v[li + 'amplitude']
        ls = v[li + 'sigma']
        model += la / (np.sqrt(2 * np.pi) * ls) * np.exp(-(x - lc) * (x - lc) * 0.5 / (ls * ls))
    if data is None:
        return model
    else:
        if eps is None:
            return (model - data)
        else:
            return (model - data) / eps
        
def linesVoigtResiduals(p, x, data=None, eps=None):
    """Pseudo-Voigt function."""
    v = p.valuesdict()
    model = 0
    n = len(v) // 4
    for i in range(n):
        li = 'l' + str(i) + '_'
        xc = v[li + 'center']
        A = v[li + 'amplitude']
        sigma = v[li + 'sigma']
        alpha = v[li + 'alpha']
        sigmag = sigma/np.sqrt(2*np.log(2.))
        xx2 = (x - xc)**2
        gauss = np.exp(-xx2 / (2 * sigmag**2)) / (np.sqrt(2*np.pi) * sigmag)
        cauchy = sigma / np.pi / (xx2 + sigma**2)
        model += A * ((1-alpha)* gauss + alpha*cauchy)
    if data is None:
        return model
    else:
        if eps is None:
            return (model - data)
        else:
            return (model - data) / eps
   
def fitApertureContinuum(sc):
    """Fit the continuum defined in the guess."""
    from lmfit import Parameters, minimize
    
    slope = sc.guess.slope
    intcpt = sc.guess.intcpt
    print('Guess values ', intcpt, slope)
    xg, yg = zip(*sc.guess.xy)
    xg = np.array(xg)
    wc = sc.spectrum.wave
    fc = sc.spectrum.flux
    try:
        ec = sc.spectrum.eflux
    except:
        ec = None
    #c = sp.gal.c
    idx = ((wc > xg[0]) & (wc < xg[1]) ) | ((wc > xg[2]) & (wc < xg[3])) & np.isfinite(fc)
    x = wc[idx]
    y = fc[idx]
    if ec is not None:
        e = ec[idx]
    # Definition of the model
    fit_params = Parameters()
    fit_params.add('intercept', value=intcpt)
    if slope != 0:
        fit_params.add('slope', value=slope)
    print('Fit continuum ... ')
    if ec is None:
        out = minimize(contResiduals, fit_params, args=(x,), kws={'data':y}, method='leastsq')
    else:
        out = minimize(contResiduals, fit_params, args=(x,), kws={'data':y,'eps':e}, method='leastsq')
    # out = minimize(contResiduals, fit_params, args=(x,), kws={'data': y, 'eps': e}, method='Nelder')
    par = out.params#.valuesdict()
    p1 = par['intercept']
    ic = p1.value
    eic = p1.stderr
    print('error on cont intercept ', eic)
    try:
        p2 = par['slope']
        s = p2.value
        es = p2.stderr
        print('error on slope ', es)
    except:
        s = 0.0
        es = 0.0
    return ic, eic, s, es


def fitApertureLines(sc, intercept, slope):
    """Fit the lines defined in the guess."""
    from lmfit import Parameters, minimize
    # Select the input values for the fit
    #z = sc.spectrum.redshift
    wc = sc.spectrum.wave
    fc = sc.spectrum.flux
    try:
        ec = sc.spectrum.eflux
    except:
        print('No errors given')
        ec = None
#    cc = sp.gal.c
    xg, yg = zip(*sc.guess.xy)
    xg = np.array(xg)
#    xg *= (1. + z)  # Back to observed
    idx = (wc > xg[0]) & (wc < xg[3]) & np.isfinite(fc)
    x = wc[idx]
    y = fc[idx]
    if ec is not None:
        e = ec[idx]
    intc, eintc = intercept
    slop, eslop = slope
    continuum = intc + slop * x

    print('variables defined')
    # Transform F_nu into F_lambda to do integration        
    # Normalization (biggest amplitude)
    amplitudes = []
    c = 299792458. # m/s
    for line in sc.lines:
        sigma = line.fwhm / 2.355
        x0 = line.x0
        A = line.A * (np.sqrt(2*np.pi) * sigma) * c / x0**2 * 1.e-20  # same units as flux
        amplitudes.append(A)
        print('x0, sigma, A, flux', x0, sigma, line.A, A)
    amplitudes = np.array(amplitudes)
    norm = np.nanmax(np.abs(amplitudes))
    print('Fit normalization ', norm)
    y *= c / x**2 * 1.e-20 / norm
    if ec is not None:
        e *= c / x**2 * 1.e-20 / norm
    continuum *= c / x**2 * 1.e-20 / norm
    # Define the model
    fit_params = Parameters()
    # Define lines
    for i, line in enumerate(sc.lines):
        li = 'l' + str(i) + '_'
        x0 = line.x0 #* (1. + z)
        sigma = line.fwhm / 2.355 #* (1. + z)
        fit_params.add(li + 'center', value=x0, min=(x0 -  sigma * 0.5), max=(x0 +  sigma * 0.5))
        A = line.A * (np.sqrt(2*np.pi) * sigma) * c / x0**2 * 1.e-20 / norm
        #print('A ', A)
        if A > 0:
            fit_params.add(li + 'amplitude', value=A, min=0., max=A * 2)
        else:
            fit_params.add(li + 'amplitude', value=A, max=0., min=A * 2)
        fit_params.add(li + 'sigma', value=sigma, min=sigma * 0.5, max=sigma * 2)
        if sc.function == 'Voigt':
            fit_params.add(li + 'alpha', value=0.2, max=0.6)
            
    # Minimize
    if ec is None:
        kws = {'data': y - continuum}
    else:
        kws = {'data': y - continuum, 'eps': e}
    if sc.function == 'Voigt':
        npars = 4
        out = minimize(linesVoigtResiduals, fit_params, args=(x,), kws=kws, method='leastsq')
    else:
        npars = 3
        out = minimize(linesGaussResiduals, fit_params, args=(x,), kws=kws, method='leastsq')
    # Return lines fitted parameters
    pars = out.params#.valuesdict()
    #print('out pars ', pars)
    nlines = len(pars) // npars
    
    linepars = []
    for i in range(nlines):
        li = 'l' + str(i) + '_'
        center = pars[li + 'center'].value  # Observed
        centerErr = pars[li + 'center'].stderr  # Observed
        sigma = pars[li + 'sigma'].value   # Observed
        sigmaErr = pars[li + 'sigma'].stderr   # Observed
        A =  pars[li+'amplitude'].value * norm
        Aerr = pars[li+'amplitude'].stderr * norm
        c0 = intc + slop * center # Continuum at line center
        if slop == 0:
            ec0 = eintc
        else:
            ec0 = c0 * (eintc/intc + eslop/slop)
        if sc.function == 'Voigt':
            alpha = pars[li + 'alpha'].value
            #factor = (1-alpha)/np.sqrt(np.pi/np.log(2)) + alpha/np.pi
            #amplitude = A / sigma *  factor
            #amplitudeErr = amplitude * (Aerr / A + sigmaErr / sigma)
            linepars.append([c0, ec0, slop, center, centerErr, A, Aerr, sigma, sigmaErr, alpha])
        else:
            #amplitude = A  / (np.sqrt(2 * np.pi) * sigma)
            #amplitudeErr = amplitude * (Aerr / A + sigmaErr / sigma)
            linepars.append([c0, ec0, slop, center, centerErr, A, Aerr, sigma, sigmaErr])
    return linepars

