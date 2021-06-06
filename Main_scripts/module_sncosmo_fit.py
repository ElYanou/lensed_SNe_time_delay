import sncosmo
import pickle as pkl
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.ticker import AutoMinorLocator
from scipy.optimize import minimize
from scipy.optimize import differential_evolution
import copy
import pycs3.gen.lc_func
import pycs3.gen.lc
import logging
loggerformat='%(message)s'
logging.basicConfig(format=loggerformat,level=logging.INFO)


################# PyCS to sncosmo ###############################

def PyCS_to_sncosmo(lc_, Filter, zpsys = 'ab', zp = 25, f_func = None, ferr_func = None):
    """
    Function that can be used to go from PyCS3 light curve object to sncosmo light curve object.
    
    : param Filter : name of the filter used for observation as string
    
    : param zpsys : name of the magnitude system convention as sting
    
    : param zp : zero point to scale flux
    
    These 3 parameters are only needed by sncosmo and are irrelevant for PyCS3 analysis. Notice that it will always be better to directly use 'load_1' function with the direct observations since computational approximation can arise.
    
    : param f_func : function to transform magnitude to flux
    
    : param ferr_func : function to transform flux, magnitude and error on magnitude to error on flux
    
    """
    lc = lc_.copy()

    lc.resetshifts()   # It makes no sense to work in sncosmo with an optimized liht curve

    time = lc.getjds()
    mag = lc.getmags()
    magerr = lc.getmagerrs()
    
    if f_func is None:
        flux = 10**(-0.4 * mag)
    else:
        flux = f_func(mag)
        
    if ferr_func is None:
        fluxerr = np.abs(10**(-0.4 * (mag - magerr)) - flux)
    else:
        fluxerr = ferr_func(flux,mag,magerr)
    
    filt = [Filter for i in range(len(flux))]
    magsys = [zpsys for i in range(len(flux))]
    
    
    obs = Table([time, flux, fluxerr, filt, magsys, zp*np.ones(len(flux))], names=('time', 'flux', 'fluxerr', 'band', 'zpsys', 'zp'))
    
    return obs
    

def sncosmo_to_PyCS(obs, obj = "Unknown", m_func = None, merr_func = None):
    """
    Function that can be used to go from sncosmo light curve object to PyCS3 light curve object.
    
    """
    
    
    time = obs['time']
    flux = obs['flux']
    fluxerr = obs['fluxerr']
    
    if m_func is None:
        mag = -2.5*np.log10(flux)
    else:
        mag = m_func(flux)
    
    if merr_func is None:
        mag_err = np.abs(-2.5 * np.log10(flux + flux_err) - mag)
    else:
        mag_err = merr_func(mag,flux,fluxerr)
        
    mag_err[mag_err<10e-4] = 10e-4  # Avoid too close to 0 errors
    
    l = pycs3.gen.lc_func.factory(np.array(time), np.array(mag), np.array(mag_err), object=obj)

    return l
   
    
################# Loaders ###############################
def merger(parameters, data_array):
    """
    Function that return a single merged light curve from the array given
    
    : param parameters : time shifts and flux scale factor in order {delta_t_i, delta_m_i}
    
    : param data_array : array containing the light curve to merge
    
    """
    
    if len(parameters) != 2*(len(data_array)-1):
        raise RuntimeError("I can't merge if the parameters' array doesn't contain time shifts and flux scale factor for every curve except the first one !")
    if len(data_array) == 1:
        raise RuntimeError("I can't merge if you give a single light curve !")
    
    lenght = len(data_array) - 1
    
    obs = copy.deepcopy(data_array[0])
    
    time = obs['time']
    flux = obs['flux']
    fluxerr = obs['fluxerr']
    filt = obs['band']
    magsys = obs['zpsys']
    zps = obs['zp']
    
    
    for i in range(lenght):
        obs2 = copy.deepcopy(data_array[i+1])

        time = np.concatenate([time, obs2['time'] - parameters[2*i]])
        flux = np.concatenate([flux, obs2['flux'] / parameters[2*i + 1]])
        fluxerr = np.concatenate([fluxerr,obs2['fluxerr'] / parameters[2*i + 1]])
        filt = np.concatenate([filt, obs2['band']])
        magsys = np.concatenate([magsys, obs2['zpsys']])
        zps = np.concatenate([zps, obs2['zp']])
    
    obs1 = Table([time, flux, fluxerr, filt, magsys, zps], names=('time', 'flux', 'fluxerr', 'band', 'zpsys', 'zp'))   
        
    return obs1
    
    
def f_to_m(f, t, ferr = None):
    """
    Flux to magnitude with flux under 0 filtered out
    """
    
    
    ind = np.where(f <= 0)[0]
    f2 = np.delete(f,ind)
    m = -2.5*np.log10(f2)
    time = np.delete(t,ind)
    
    if ferr is not None:
        merr = np.abs(-2.5*np.log10(f+ferr)-m)
        merr[merr<10e-4] = 10e-4
        return m, time, merr
        
    else:
        return m, time

def clean_flux(obs):
    """
    Returns a copy of the light curve obs with flux under 0 filtered away 
    """
    
    ind = np.where(obs['flux']<0)[0]
    
    t = np.delete(obs['time'],ind)
    f = np.delete(obs['flux'],ind)
    fe = np.delete(obs['fluxerr'],ind)
    b = np.delete(obs['band'],ind)
    zy = np.delete(obs['zpsys'],ind)
    zp = np.delete(obs['zp'],ind)
    
    obs1 = Table([t, f, fe, b, zy, zp], names=('time', 'flux', 'fluxerr', 'band', 'zpsys', 'zp'))
    
    return obs1
  
  
def load_1(file_name, Filter, zpsys = 'ab', zp = 25):
    """
    Function that returns a light curve with the format used by sncosmo.
    
    : param file_name : name of the file containing the light curve time, flux and error on the flux. The 2 first rows are ignored and the file is in .dat format since our data set was built like this.
    
    : param Filter : name of the filter used for observation as string
    
    : param zpsys : name of the magnitude system convention as sting
    
    : param zp : zero point to scale flux 
    
    """
    
    
    data = np.loadtxt(file_name, skiprows = 2)
    l = len(data[:,1])
    
    filt = [Filter for i in range(l)]
    magsys = [zpsys for i in range(l)]
    
    obs = Table([data[:,0], data[:,1], data[:,2], filt, magsys, zp*np.ones(l)], names=('time', 'flux', 'fluxerr', 'band', 'zpsys', 'zp'))
    obs = clean_flux(obs)
    
    return obs
      
    
    
def load_multiple(file_names, Filters, zpsys, zp):
    """
    Function that returns multiple light curves as a single object with the format used by sncosmo. The parameters are the same as for load_1 but must be contained in arrays with the same length as file_names.
    
    """
    
    if (len(file_names) != len(Filters)) or (len(file_names) != len(zpsys)) or (len(file_names) != len(zp)):
        raise RuntimeError("The array arguments must be of the same size !")
    arr = []
    p = []
    
    for i in range(len(file_names)):
        arr.append(load_1(file_names[i], z, Filters[i], zpsys[i], zp[i]))
        p.append(0) 
        p.append(1)     
    
    obs1 = merger(p,arr)
    
    return obs1  
    

################# Objective functions ###############################

# take care this function works only for templates with time shift and magnitude shift as 2nd and 3rd parameter of the model !!!
def objective(parameters, models, data_array, z):
    """
    Function that computes a global chi square. It is defined as the sum of each chi square for each light curve given the same template for each.
    
    : param parameters : should contain in order {delta_t_i, delta_m_i} and model.parameters 
    
    : param models : the model to use for the fit 
    
    : param data_array : list or tuple of light curves
    
    : param z : known redshift
    
    """
    
    if len(parameters) != 2*(len(data_array)-1) + len(models.parameters) - 1:
        raise RuntimeError("Parameters' array has not the right size, check the documentation of this function to see what is expected !")
    
    chi2 = 0  
    lenght = len(data_array) - 1
    model = copy.deepcopy(models)
    
    model.parameters = np.concatenate(([z],parameters[2 * lenght:]))
    t0 = parameters[2 * lenght]
    c = parameters[2 * lenght + 1]
    

    for i in range(lenght+1):
        
        if i == 0:
            model.parameters[1] = t0
            model.parameters[2] = c
        else:
            model.parameters[1] = t0 + parameters[2*(i-1)]
            model.parameters[2] = c * parameters[2*(i-1) + 1]
        
               
        chi2 += (1/len(data_array[i])) * sncosmo.chisq(data_array[i], model)
        
    
    return chi2
        
    
    
  
def objective_merge(parameters, models, data_array, z):
    """
    Function that computes a chi square on the merged light curves. 
    
    : param parameters : should contain in order {delta_t_i, delta_m_i} and model.parameters 
    
    : param models : the model to use for the fit 
    
    : param data_array : list or tuple of light curves
    
    : param z : known redshift
    
    """
    
    if len(parameters) != 2*(len(data_array)-1) + len(models.parameters) - 1:
        raise RuntimeError("Parameters' array has not the right size, check the documentation of this function to see what is expected !")
        
    lenght = len(data_array) - 1    
    model = copy.deepcopy(models)
    
    model.parameters = np.concatenate(([z],parameters[2 * lenght:]))
    obs1 = merger(parameters[:2*lenght], data_array)
    
    chi2 = sncosmo.chisq(obs1, model)
    
    return chi2
    
# the 2 objective_merge are the same except for the size of the array containing the parameters and the fact that in merge2 the model is excepted to have its parameters already attached   
def objective_merge2(parameters, models, data_array):
    """
    Same as objective_merge but doesn't optimize the parameters of the model
    
    : param parameters : should contain in order {delta_t_i, delta_m_i} 
    
    : param models : the model to use for the fit with its parameters attached
    
    : param data_array : list or tuple of light curves
    
    """

    if len(parameters) != 2*(len(data_array)-1):
        raise RuntimeError("Parameters' array has not the right size, check the documentation of this function to see what is expected !")
        
    model = copy.deepcopy(models)   
    obs1 = merger(parameters, data_array)
    
    chi2 = sncosmo.chisq(obs1, model)
    
    return chi2
    
# the model is supposed to already have its parameters attached   
def objective_time(param, obs, model):
    """
    Compute the chi square of several given light curves with time axis rescaling, should be used with merge method. 
    
    : param param : float, time axis rescaling
    
    : param model : the model to use for the fit with its parameters attached
    
    : param obs : list or tuple of light curves
    
    """
    
    chi2 = 0
    
    for i in range(len(obs)):
        lc = copy.deepcopy(obs[i])
    
        lc['time'] = lc['time']*param
        chi2 += sncosmo.chisq(lc,model)
    
    return chi2
    
def objective_time2(param, obs, model, g):
    """
    Compute the chi square of several given light curves with time axis rescaling and time and flux scale factor. It is linked to the no merge method. Should only be used with models that have time and flux scale factor as 2nd and 3rd param. 
    
    : param param : float, time axis rescaling
    
    : param model : the model to use for the fit with its parameters attached
    
    : param obs : list or tuple of light curves
    
    : param g : should contain in order {delta_t_i, delta_m_i} 
    
    """
    
    chi2 = 0
    mod = copy.deepcopy(model)
    
    for i in range(len(obs)):
        lc = copy.deepcopy(obs[i])
        if i != 0:
            model.parameters[1] += g[2*(i-1)] * param
            model.parameters[2] = g[2*(i-1) + 1] * param
        
        lc['time'] = lc['time'] * param
        chi2 += sncosmo.chisq(lc,mod)
    
    return chi2

################# Plotter ###############################

def plot_as_PyCS(model, t, array, delta_t, delta_m, Filter, at=1):  
    """
    Returns a plot of the light curves and the spline as PyCS3 does.
    
    : param model : model to use to draw the spline with parameters attached 
    
    : param t : time range of the light curves, should be a linspace of at least 1000 points
    
    : param array : list or tuple of light curves
    
    : param delta_t : time shifts
    
    : param delta_m : flux scale factors
    
    : param Filter : which filter do you want to display, we can only display one filter at the time despite the light curves might be composed of several 
    
    : param at : time axis scale factor
    
    """
    
    if len(array)-1 != len(delta_t) or len(array)-1 != len(delta_m):
        raise RuntimeError("You should give me time shift and flux scale factors arrays of size len(lc_list)-1 !")
    
    c = ['darkorange', 'royalblue', 'seagreen', 'purple', 'brown', 'magenta', 'orange']
    mod, time = f_to_m(model.bandflux(Filter,t*at,25,'ab'),t)
    
    
    fig, ax = plt.subplots()
    ax.yaxis.set_minor_locator(AutoMinorLocator(10))    
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    
    ax.plot(time,mod,color='k',zorder=15)
    
    for i in range(len(array)):
        
        ind = np.where(array[i]['band']!=Filter)[0]
        obs = copy.deepcopy(array[i])
        f = np.delete(obs['flux'],ind)
        ferr = np.delete(obs['fluxerr'],ind)
        ti = np.delete(obs['time'],ind)
        
        s1, t1, err1 = f_to_m(f,ti,ferr)
        
        if i == 0:
            ax.errorbar(t1, s1, err1, ecolor="#BBBBBB", elinewidth=0.5, fmt=".", markersize=6, color=c[i],capsize = 8)
        else:
            ax.errorbar(t1 - delta_t[i-1], s1 + 2.5*np.log10(delta_m[i-1]), err1, ecolor="#BBBBBB", elinewidth=0.5, fmt=".", markersize=6, color=c[i],capsize = 8)


    ax.invert_yaxis()
    ax.set_xlabel('HJD - 2400000.5 (day)',fontsize=12)
    ax.set_ylabel('Magnitude (relative)',fontsize=12)


def template_visual_research(template_list, obs1, z, at = 1):
    """
    Returns plots of light curves as fitted by each templates given.
    
    : param template_list : list or tuple of templates' name
    
    : param obs1 : a light curve, should be your reference one
    
    : param z : redshift
    
    : param at : time axis scale factor
    
    """
    
    obs = copy.deepcopy(obs1)
    obs['time'] = obs['time'] * at   
    
    for i in template_list:
        source = sncosmo.get_source(i)
        model = sncosmo.Model(source=source)
        model.set(z=z)
        param_name = ['t0'] + model.source.param_names
        
        _ , fitted_model = sncosmo.fit_lc(obs, model, param_name)
        
        sncosmo.plot_lc(obs,fitted_model)
    
        print('For model', i, 'chi square is', sncosmo.chisq(obs,fitted_model))

    if len(template_list) >= 19:
        raise RuntimeError("I won't let you continue, you have too much graph to show (more than 19), your pc might explode !")
        
def time_scale_visual_research(template, obs1, z, at_list):  
    """
    Returns plots of light curves as fitted with each time axis scale factor and a given template, should be used to look for right bounds on at.
    
    : param template : template to use
    
    : param obs1 : a light curve, should be your reference one
    
    : param z : redshift
    
    : param at_list : time axis scale factor list or tuple
    
    """
    
    source = sncosmo.get_source(template)
    model = sncosmo.Model(source=source)
    model.set(z=z)
    param_name = ['t0'] + model.source.param_names
    
    for i in at_list:
        obs = copy.deepcopy(obs1)
        obs['time'] = obs['time'] * i
        
        _ , fitted_model = sncosmo.fit_lc(obs, model, param_name)
        
        sncosmo.plot_lc(obs,fitted_model)
    
        print('For time scaling', i, 'chi square is', sncosmo.chisq(obs,fitted_model))

    if len(at_list) >= 19:
        raise RuntimeError("I won't let you continue, you have too much graph to show (more than 19), your pc might explode !")
        
def time_scale_research(template, obs1, z, at_list):  
    """
    Returns plots of the chi square as a function of the scale factor, you should have already a good guess on at to be efficient
    
    : param template : template to use
    
    : param obs1 : a light curve, should be your reference one
    
    : param z : redshift
    
    : param at_list : time axis scale factor list or tuple
    
    """

    source = sncosmo.get_source(template)
    model = sncosmo.Model(source=source)
    model.set(z=z)
    param_name = ['t0'] + model.source.param_names
    chi = []
    
    for i in at_list:
        obs = copy.deepcopy(obs1)
        obs['time'] = obs['time'] * i
        
        _ , fitted_model = sncosmo.fit_lc(obs, model, param_name)
        
        chi.append(sncosmo.chisq(obs,fitted_model))

    plt.plot(at_list,chi,color='orangered')   
    plt.xlabel('Value of the time axis rescaling', fontsize=14)
    plt.ylabel(r'$\chi^2$ with optimal param', fontsize=14)
    plt.grid(True)

    


################# Optimizers ###############################


def optimize_simple(data_array, model, z):
    """
    optimize time shifts and flux scale factors for 2 parameters models, notice you have no guess for this optimizer
    
    : param data_array : list or tuple of light curves
    
    : param model : model to use
    
    : param z : redshift
    
    """
    
    if len(model.parameters) != 3:
        raise RuntimeError("You're trying to use a method made for 2 parameters models !")
    
    model.set(z = z)
    param_name = ['t0'] + model.source.param_names
    
    result, fitted_model = sncosmo.fit_lc(data_array[0], model, param_name)
    delta_t = []
    delta_m = []
    t0, m0 = fitted_model.parameters[1], fitted_model.parameters[2]
    
    for i in range(len(data_array) - 1):
        result2, fitted_model2 = sncosmo.fit_lc(data_array[i + 1], model, param_name)
        t, m = fitted_model2.parameters[1], fitted_model2.parameters[2]
        
        delta_t.append(t - t0)
        delta_m.append(m / m0)
        
    model.parameters = [z,t0,m0]    
    return delta_t, delta_m
      
               


# faster than opt_time_merge since we use minimize, a bit different since you optimize every parameter at the same time. Take care this method isn't made to be ran with time scaled axis but you can manually do it and rescale the time shift at the end.
def optimize_merge(data_array, model, z, guess):
    sol = minimize(objective_merge, guess, args=(model, data_array, z)) 
    model.parameters = np.concatenate(([z], sol.x[(len(data_array)-1)*2:]))
    return sol.x
    
 
    
def opt_time_nomerge(guess, obs, model, z, t_in, nit=1, th1 = 30, th2 = 1, solu = None, bounds_t = None):
    """
    optimize time shifts and flux scale factors for 2 parameters models with no merge method
    
    : param guess : should contain in order a guess of {delta_t_i, delta_m_i} and model.parameters 
    
    : param obs : tuple or list of light curves
    
    : param z : redshift
    
    : param model : model to use
    
    : param t_in : guess of time axis scale factor
    
    : param nit : number of run
    
    : param t_in : guess of time axis scale factor
    
    : param th1 : threshold for time shifts, bounds are set to (guess-th1,guess+th1)
    
    : param th2 : threshold for flux scale factor, bounds are set to (guess-th2,guess+th2), you should avoid negative values in the interval
    
    : param solu : if given doesn't optimize the time axis scale factor and use the given value
    
    : param bounds_t : bounds of time axis scale factor, should be really tight
    
    """

    if len(guess) != 2*(len(obs)-1) + len(model.parameters) - 1:
        raise RuntimeError("Guesses' array has not the right size, check the documentation of this function to see what is expected !")
        
    g = guess
    lenght = len(obs) - 1
    
    mod = copy.deepcopy(model)
    
    mod.parameters = np.concatenate(([z],guess[2 * lenght:]))
    if solu is None:
        ts = t_in
    else:
        ts = solu
    
    for i in range(nit):
        lc = copy.deepcopy(obs)
        sc = copy.deepcopy(ts)
        k = copy.deepcopy(g)
        
        for i in range(len(obs)):
            if i != 0:
                k[2*(i-1)] = ts * k[2*(i-1)]
            lc[i]['time'] = lc[i]['time'] * ts
        
        bounds = [(k[0]-ts*th1,k[0]+ts*th1),(k[1]-th2,k[1]+th2),(k[2]-ts*th1,k[2]+ts*th1),(k[3]-th2,k[3]+th2),(k[4]-ts*th1,k[4]+ts*th1),(k[5]-th2,k[5]+th2),(k[6]-ts*th1,k[6]+ts*th1),(k[7]-th2,k[7]+th2),(k[8]-abs(k[8]),k[8]+abs(k[8])),(k[9]-abs(k[9]),k[9]+abs(k[9])),(k[10]-abs(k[10]),k[10]+abs(k[10])),(k[11]-abs(k[11]),k[11]+abs(k[11]))]
        
        g = minimize(objective, k, args=(model, lc, z), bounds = bounds).x
        mod.parameters = np.concatenate(([z],g[2 * lenght:]))
        
        for i in range(lenght):              
            g[2*i] = g[2*i] / ts
            
        if solu is None:      
            if bounds_t is None:
                bounds_t = [(1/(z+1)**2-0.01, 1/(z+1)**2+0.01)]  
            sol = differential_evolution(objective_time2, bounds_t, args=(obs,mod,g))
            ts = sol.x[0]
        else:
            ts = solu        
    
    model.parameters = np.concatenate(([z],g[2 * lenght:]))
    
    return g, sc
       
    
def opt_time_merge(guesst, guessp, obs, model, z, t_in, nit=1, th1 = 30, th2 = 1, solu = None, bounds_t = None):
    """
    optimize time shifts and flux scale factors for 2 parameters models with the merge method
    
    : param guesst : should contain a guess of {delta_t_i, delta_m_i} 
    
    : param guessp : should contain a guess of model.parameters 
    
    : param obs : tuple or list of light curves
    
    : param z : redshift
    
    : param model : model to use
    
    : param t_in : guess of time axis scale factor
    
    : param nit : number of run
    
    : param t_in : guess of time axis scale factor
    
    : param th1 : threshold for time shifts, bounds are set to (guess-th1,guess+th1)
    
    : param th2 : threshold for flux scale factor, bounds are set to (guess-th2,guess+th2), you should avoid negative values in the interval
    
    : param solu : if given doesn't optimize the time axis scale factor and use the given value
    
    : param bounds_t : bounds of time axis scale factor, should be really tight
    
    """

    if len(guesst) != 2*(len(obs)-1):
        raise RuntimeError("Guesses' array has not the right size, check the documentation of this function to see what is expected !")
    
    gt = guesst
    gp = guessp
    lenght = len(obs) - 1 
    bounds = []
    mod = copy.deepcopy(model)
     
    
    mod.parameters = np.concatenate(([z],gp))
    if solu is None:
        ts = t_in
    else:
        ts = solu
        
    for i in range(nit):
        sc = copy.deepcopy(ts)
        lc = copy.deepcopy(obs)

        k = copy.deepcopy(gt)
        
        for i in range(lenght+1):
            if i != 0:
                k[2*(i-1)] = ts * k[2*(i-1)]
                bounds.append(((gt[2*(i-1)]-th1)*ts, (gt[2*(i-1)]+th1)*ts))
                bounds.append((gt[2*(i-1)+1]-th2, gt[2*(i-1)+1]+th2))
            
            lc[i]['time'] = lc[i]['time'] * ts
        try:
            gt = differential_evolution(objective_merge2, bounds, args=(mod, lc)).x
            lcs = merger(k,lc)
            
            param_name = ['t0'] + model.source.param_names
            _, mod = sncosmo.fit_lc(lcs, mod, param_name)
        
        except Exception:
            print("I couldn't find an optimal solution !")
        
        gp = mod.parameters[1:]
        for i in range(lenght):              
            gt[2*i] = gt[2*i] / ts
        
        if solu is None:      
            if bounds_t is None:
                bounds_t = [(1/(z+1)**2-0.01, 1/(z+1)**2+0.01)]  
            sol = differential_evolution(objective_time, bounds_t, args=(obs,mod))
            ts = sol.x[0]
        else:
            ts = solu
    
    model.parameters = np.concatenate(([z],gp))
       
    return gt, gp, sc
    
    
    

      
    

