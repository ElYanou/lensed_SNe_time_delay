#!/usr/bin/env python3
import os
import numpy as np  
import pycs3.gen.lc_func
import pycs3.gen.mrg
import pycs3.gen.util

import logging
loggerformat='%(message)s'
logging.basicConfig(format=loggerformat,level=logging.INFO)

from param_lc_maker import *


abs_path = os.path.abspath(path_dir+Filter)
print(path_dir)
if os.path.exists(abs_path) == False:
    os.mkdir(path_dir+Filter)

# --------------------------------------------------------------------------------------------- #
    files = os.listdir(data_dir)
    supernovae = ['s1','s2','s3','s4','sx']  # it might be changed depending on which data you use 
   
    for j in files: 
        lcs1 = []  # contains the 5 light curves 
        for i in supernovae :
            
            file_path = data_dir + "/" + j + "/phot_" + i + "_WFC3_IR_" + Filter + ".dat" # it might be changed depending on which data you use 
            array = np.loadtxt(file_path, skiprows=2)
            ind = np.where(array[:,1] < 0)[0]  # we ignore data with f negative
            
            # Flux to magnitude
            flux = np.delete(array[:,1],ind)
            flux_err = np.delete(array[:,2],ind)

            mag = -2.5*np.log10(flux)
            mag_err = np.abs(-2.5 * np.log10(flux + flux_err) - mag)
            mag_err[mag_err < 10e-4] = 10e-4  # Avoid too close to 0 errors
            
            mjhd = np.delete(array[:,0],ind)
            
            l = pycs3.gen.lc_func.factory(np.array(mjhd), np.array(mag), np.array(mag_err), object=i)
            lcs1.append(l)

        pycs3.gen.mrg.colourise(lcs1) # Gives each curve a different colour
        
        am = j[-3:]
        if am[1] == "_":
            am = am[2]          
        elif am[0] == "_":
            am = am[1:]
        
        
        pycs3.gen.util.writepickle(lcs1, path_dir + Filter + "/lcs_" + Filter + "_no_opt_" + am)

else :
    print("The file already exists, you probably don't want to add files into it !")





































