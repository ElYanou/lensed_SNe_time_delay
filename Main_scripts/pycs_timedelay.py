#!/usr/bin/env python3
import os
import re
import shutil
from multiprocessing import Pool
import pickle as pkl
import numpy as np 
import pycs3.gen.lc_func
import pycs3.gen.splml
import pycs3.spl.topopt
import pycs3.gen.util

import logging
loggerformat='%(message)s'
logging.basicConfig(format=loggerformat,level=logging.INFO)

from param_pycs_timedelay import *

filess = os.listdir(sim_path)

filess.remove(file0_name) # we have no comparison with true time delay for it

def true_delay_get(num):
    """
    Loads the true time delay of light curve num, this is hardcoded and might be changed depending on the 
    dataset and the way your light curve are called.
    """


    file = open(data_folder + '/phot_dolphot_2_1.0000_True_1_' + num + '/phot.ratios','r')
    file_delays = file.readlines()
    stockage = np.array([])
    
    
    for line in file_delays:
        line = re.findall(r"[-+]?\d*\.\d+|\d+", line)
        stockage = np.concatenate((stockage,np.array([float(line[-1])])))
    
    return stockage


def run(files, write_name, kn, nb_cycle, degree, numbml, sttpos, magstt, do_you_want_ml, do_you_display = False):
    """
    Computes the time delay differences given a non optimized set of light curves as pkl file. 
    """
    
    lcs1 = pycs3.gen.util.readpickle(sim_path + '/' + files)
    
    np.random.seed()
    startpts = sttpos + np.random.uniform(low=-10.0, high=10.0, size=len(lcs1))
    
    pycs3.gen.lc_func.applyshifts(lcs1, startpts, magstt)
    
    if do_you_want_ml:
        if numbml == 1:
            for l in lcs1[1:]:
                pycs3.gen.polyml.addtolc(l, nparams=degree+1, autoseasonsgap=600.0)  # polynomial spline
        else:
            for l in lcs1[1:]:  
                curve_length = l.jds[-1] - l.jds[0]
                mlbokeps = np.floor(curve_length / numbml)
            
                pycs3.gen.splml.addtolc(l, n = numbml, bokeps = mlbokeps)   # full spline
    
    try:
        
        spline = pycs3.spl.topopt.opt_fine(lcs1, nit = nb_cycle, knotstep=kn,bokeps=kn / 3.0, stabext=100, verbose = False)
        all_delays = pycs3.gen.lc_func.getdelays(lcs1, to_be_sorted=False)

        
        am = files[-3:]
        if am[1] == "_":
            am = am[2]          
        elif am[0] == "_":
            am = am[1:]
    
        stockage = true_delay_get(am)
    
        true_delays = np.array([])
        true_mu = np.array([])
        
        for j in range(2,len(stockage)):
            if j % 2 == 1:
                true_mu = np.concatenate((true_mu,np.array([stockage[j]])))
            else: 
                true_delays = np.concatenate((true_delays,np.array([stockage[j]])))
    

        seter = true_delays + all_delays[:4]
        
        writers = write_name + files
    
        np.savetxt(writers, seter)
        
        
        if do_you_display:
            pycs3.gen.lc_func.display(lcs1, [spline], figsize = (12,9), showlegend =False, title=r"$\mathrm{SN\,Refsdal}$", nicefont=True, capsize = 15, showdelays =False, hidecolourbar = True)
        
    
    except Exception:
        print('light curve ', files, ' won t be used !')
        

def multi(files):
    
    run(files, writer + '/', kn_p, nb_cycle_p, degree_p, numbml_p, sttpos_p, magstt_p, do_you_want_ml_p, do_you_display_p)


def multiproc(files):
    if do_you_display_p:
        n = 1
    else :
        n = os.cpu_count()
    
    pool = Pool(n)                         
    pool.map(multi, files)
    
    pool.close()
    pool.join()
    
    delay_s2 = np.array([])
    delay_s3 = np.array([])
    delay_s4 = np.array([])
    delay_sx = np.array([])


    for i in os.listdir(writer):
        seted = np.array([])
        file = open(writer + '/' + i,'r')
        file_delays = file.readlines()
    
        for line in file_delays:
            seted = np.concatenate((seted,np.array([float(line)])))
    
        delay_s2 = np.append(delay_s2, seted[0])
        delay_s3 = np.append(delay_s3, seted[1])
        delay_s4 = np.append(delay_s4, seted[2])
        delay_sx = np.append(delay_sx, seted[3])
    
    final = np.vstack((delay_s2,delay_s3,delay_s4,delay_sx))

    writers = open(final_name,"wb")
    pkl.dump(final,writers)
    writers.close()   

def runner(file, folder = writer):
    os.mkdir(folder)

    multiproc(file)

    shutil.rmtree(folder)


runner(filess)











