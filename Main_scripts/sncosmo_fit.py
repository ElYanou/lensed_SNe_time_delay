#!/usr/bin/env python3
import os
import shutil
import re
from multiprocessing import Pool
from module_sncosmo_fit import *
from param_sncosmo_fit import *


################################## Multirun process ######################################

def fit_delay(data_array, model, startpts, k, optimizer, t_in, nit = 1, th1 = 30, th2 = 1, solu = None, bounds_t = None):
    """
    Compute the time delays
    
    : param k : guess on the model's paramters
    
    """
    
    z = model.get('z')
    if optimizer == optimize_simple:
        if solu is None:
            raise RuntimeError("You must specify a time axis scale solution to use this optimizer ! If you don't want to rescale it just set to 1.")
            
        lc = copy.deepcopy(data_array)
        for i in lc:
            i['time'] = i['time'] * solu
                
        fitted_delay, flux_scale = optimize_simple(lc, model, z)
        fitted_delay = np.array(fitted_delay) / solu
        at = solu
        
    if optimizer == optimize_merge:
        if solu is None:
            raise RuntimeError("You must specify a time axis scale solution to use this optimizer ! If you don't want to rescale it just set to 1.")
            
        guess = np.concatenate((startpts, k))
        guess[0] = guess[0] * solu
        guess[2] = guess[2] * solu
        guess[4] = guess[4] * solu
        guess[6] = guess[6] * solu
        lc = copy.deepcopy(data_array)
        for i in lc:
            i['time'] = i['time'] * solu
            
        final = optimize_merge(lc, model, z, guess)    
        fitted_delay = [final[0] / solu,final[2] / solu,final[4] / solu,final[6] / solu]  
        flux_scale = [final[1],final[3],final[5],final[7]]  
        at = solu 
        
    if optimizer == opt_time_nomerge:
        guess = np.concatenate((startpts, k))
        final, at = opt_time_nomerge(guess, data_array, model, z, t_in, nit, th1, th2, solu, bounds_t)
        fitted_delay = [final[0],final[2],final[4],final[6]]
        flux_scale = [final[1],final[3],final[5],final[7]]
        
    if optimizer == opt_time_merge: 
        final, _ , at = opt_time_merge(startpts, k, data_array, model, z, t_in, nit, th1, th2, solu, bounds_t)
        fitted_delay = [final[0],final[2],final[4],final[6]]
        flux_scale = [final[1],final[3],final[5],final[7]]
            
    return fitted_delay, flux_scale, at


def run(data_array, model, startpts, optimizer, file_name, folder_name, write_name, t_in, nit = 1, th1 = 30, th2 = 1, solu = None, bounds_t = None):    
    """
    Compute time-delay difference
    
    : param data_array : tuple or list of light curves
    
    : param startpts : should contain a guess of {delta_t_i, delta_m_i} 
    
    : param optimizer : optimizer to use
    
    : param model : model to ue with redshift attached
    
    : param file_name : name of the file containing the true delays
    
    : param file_name : path of the folder containing the true delays files
    
    : param write_name : path to a folder to store the time delay difference
    
    # the following parameters are not useful for each optimizer
    
    : param t_in : guess of time axis scale factor
    
    : param nit : number of run
    
    : param t_in : guess of time axis scale factor
    
    : param th1 : threshold for time shifts, bounds are set to (guess-th1,guess+th1)
    
    : param th2 : threshold for flux scale factor, bounds are set to (guess-th2,guess+th2), you should avoid negative values in the interval
    
    : param solu : if given doesn't optimize the time axis scale factor and use the given value
    
    : param bounds_t : bounds of time axis scale factor, should be really tight
    
    """

    
    z = model.get('z')
    
    if z == 0:
        raise RuntimeError("You must set z in your model given in argument !")
 
    print('I will do ', file_name)
    
    if solu is not None:
        t_in = solu
        
    obs1_bis = copy.deepcopy(data_array[0])
    obs1_bis['time'] = obs1_bis['time'] * t_in
      
    param_name = ['t0'] + model.source.param_names
    
    _ , fitted_model = sncosmo.fit_lc(obs1_bis, model, param_name)
    k = fitted_model.parameters[1:]
    
                   
    try:
        fitted_delay, _, _ = fit_delay(data_array, model, startpts, k, optimizer, t_in, nit, th1, th2, solu, bounds_t)     
        
        file = open(folder_name + '/' + file_name + '/phot.ratios','r')
        file_delays = file.readlines()
        stockage = np.array([])
        
        
    
        for line in file_delays:
            line = re.findall(r"[-+]?\d*\.\d+|\d+", line)
            stockage = np.concatenate((stockage,np.array([float(line[-1])])))
    
        true_delays = np.array([])
        true_mu = np.array([])
      
        for j in range(2,len(stockage)):
            if j % 2 == 1:
                true_mu = np.concatenate((true_mu,np.array([stockage[j]])))
            else: 
                true_delays = np.concatenate((true_delays,np.array([stockage[j]])))
    
        
        seter = true_delays - fitted_delay
        
        writer = write_name + file_name
    
        np.savetxt(writer, seter)
        
    except:
        print('light curve ', file_name, ' won t be used !')
        
    
    
        
        
   
            
def run_mul(files):
    
    source = sncosmo.get_source(template_p)
    
    model = sncosmo.Model(source=source)

    model.set(z = z_p)
   
    obs1 = load_1(folder_p + '/' + files + '/phot_s1_WFC3_IR_' + filt_dat + '.dat',Filter_p,'ab',25)
    obs2 = load_1(folder_p + '/' + files + '/phot_s2_WFC3_IR_' + filt_dat + '.dat',Filter_p,'ab',25)
    obs3 = load_1(folder_p + '/' + files + '/phot_s3_WFC3_IR_' + filt_dat + '.dat',Filter_p,'ab',25)
    obs4 = load_1(folder_p + '/' + files + '/phot_s4_WFC3_IR_' + filt_dat + '.dat',Filter_p,'ab',25)
    obsx = load_1(folder_p + '/' + files + '/phot_sx_WFC3_IR_' + filt_dat + '.dat',Filter_p,'ab',25)
    
    np.random.seed()
    un = np.random.uniform(low=-10.0, high=10.0, size = 4)
    
    pts = startpts_p
    
    pts[0] += un[0]
    pts[2] += un[1]
    pts[4] += un[2]
    pts[6] += un[3]
   
    run(np.array([obs1, obs2, obs3, obs4, obsx],dtype=object), model, pts, optimizer_p, files, folder_p, path_stockage + '/', t_in_p, nit_p, th1_p, th2_p, solu_p, bounds_t_p)
 
    
def multiproc(files):
    
    pool = Pool(os.cpu_count())                         
    pool.map(run_mul, files)
    
    pool.close()
    pool.join()
    
    delay_s2 = np.array([])
    delay_s3 = np.array([])
    delay_s4 = np.array([])
    delay_sx = np.array([])


    for i in os.listdir(path_stockage):
        seted = np.array([])
        file = open(path_stockage + '/' + i,'r')
        file_delays = file.readlines()
    
        for line in file_delays:
            seted = np.concatenate((seted,np.array([float(line)])))
    
        delay_s2 = np.append(delay_s2, seted[0])
        delay_s3 = np.append(delay_s3, seted[1])
        delay_s4 = np.append(delay_s4, seted[2])
        delay_sx = np.append(delay_sx, seted[3])
    
    final = np.vstack((delay_s2,delay_s3,delay_s4,delay_sx))

    writer = open(final_path,"wb")
    pkl.dump(final,writer)
    writer.close()
    
    

    print("We're done !")


def runner(file, folder = path_stockage):
    os.mkdir(folder)

    multiproc(file)

    shutil.rmtree(folder)


################################### Computation ##############################################

# multi process
if mul_processing is True:
    filess = os.listdir(folder_p)

    filess.remove(name_0) # we have no comparison with true time delay for it

    runner(filess)
    print("Multi-processing is over, i'll do next task")

# run optimizer over 1 set of lc and plot it 
if single_run is True:
    if solu_p is not None:
        t_in_p = solu_p
        
    obs = []
    source = sncosmo.get_source(template_p)    
    model = sncosmo.Model(source=source)
    model.set(z = z_p)

    for i in file_list_p:
        obs.append(load_1(i,Filter_p))
    
    obs1_bis = copy.deepcopy(obs[0])
    obs1_bis['time'] = obs1_bis['time'] * t_in_p
      
    param_name = ['t0'] + model.source.param_names
    
    _ , fitted_model = sncosmo.fit_lc(obs1_bis, model, param_name)
    k = fitted_model.parameters[1:]
    
    x, y, at = fit_delay(obs, model, startpts_p, k, optimizer_p, t_in_p, nit_p, th1_p, th2_p, solu_p, bounds_t_p)
    print("Time delays are :", x)
    plot_as_PyCS(model, t_p, obs, x, y, Filter_p, at)
    print("Close figure to start next task !")
    plt.show()

# template research
if template_search is True:
    obs1 = load_1(file_list_p[0],Filter_p)
    
    template_visual_research(template_list_p, obs1, z_p, solu_p)
    print("Close figures to start next task !")
    plt.show()
    

# time axis scale factor visual check
if at_visual_search is True:
    obs1 = load_1(file_list_p[0],Filter_p)
    
    time_scale_visual_research(template_p, obs1, z_p, at_list_p)
    print("Close figures to start next task !")
    plt.show()
    
# time axis scale factor research
if at_search is True:
    obs1 = load_1(file_list_p[0],Filter_p)
    
    time_scale_research(template_p, obs1, z_p, at_biglist_p)
    print("Close figures to start next task !")
    plt.show()

print("No further tasks ! Feel free to modify the code to implement new ones !")


























