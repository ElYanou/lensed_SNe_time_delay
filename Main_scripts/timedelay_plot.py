#!/usr/bin/env python3
import matplotlib.pyplot as plt 
import pickle as pkl
import numpy as np 

from param_timedelay_plot import *


def delays(file_name, threshold = 30):
    """
    Filters out the delay that are considered as a critical failure and returns the their proportions
    
    :param file_name: name of the file in which the time delay differences are stored. It must be a 
    pkl containing an array with 4 elements which are arrays containing the time delay differences.
    
    :param threshold: threshold to be considered as a critical failure     
    """
    
    file = open(file_name,"rb")
    array = pkl.load(file)

    delay_s2 = array[0]
    delay_s3 = array[1]
    delay_s4 = array[2]
    delay_sx = array[3]

    
    lentot2, lentot3, lentot4, lentotx = np.size(delay_s2), np.size(delay_s3), np.size(delay_s4), np.size(delay_sx)
    
    ind = np.where(delay_s2 < -threshold)[0]
    delay_s2 = np.delete(delay_s2,ind) 
    
    ind2 = np.where(delay_s2 > threshold)[0]
    delay_s2 = np.delete(delay_s2,ind2)
 
    s2deleted = np.size(ind) + np.size(ind2)
    
    ind = np.where(delay_s3 < -threshold)[0]
    delay_s3 = np.delete(delay_s3,ind)   
    ind2 = np.where(delay_s3 > threshold)[0]
    delay_s3 = np.delete(delay_s3,ind2)
    
    s3deleted = np.size(ind) + np.size(ind2)
    
    ind = np.where(delay_s4 < -threshold)[0]
    delay_s4 = np.delete(delay_s4,ind)   
    ind2 = np.where(delay_s4 > threshold)[0]
    delay_s4 = np.delete(delay_s4,ind2)
    
    s4deleted = np.size(ind) + np.size(ind2)
    
    ind = np.where(delay_sx < -threshold)[0]
    delay_sx = np.delete(delay_sx,ind)   
    ind2 = np.where(delay_sx > threshold)[0]
    delay_sx = np.delete(delay_sx,ind2)
    
    sxdeleted = np.size(ind) + np.size(ind2)
    
    return delay_s2, delay_s3, delay_s4, delay_sx, s2deleted/lentot2, s3deleted/lentot3, s4deleted/lentot4, sxdeleted/lentotx
    

    
def plot_multiple(array, leg, color, lin, bins = 500):
    """
    Returns a 4 pannel-histogram of the distributions for each SNe image
    
    :param array: array containing the time delay differences
    
    :param leg: legends to attach in the plot
    
    :param color: array of the color of the curves
    
    :param lin: array of the linestyle of the curves
    
    :param bin: bin of the histogram
    """
    
    if len(array) != len(leg) or len(array) != len(color) or len(array) != len(lin):
        raise RuntimeError("Your arrays must be of the same lengths !")
    
    fig, axs = plt.subplots(2, 2,figsize=(12,8))
    
    for i in range(np.size(array,0)):
        if len(leg[i]) != 4:
            raise RuntimeError("Your legend array must have elements of size 4 !")

        delay_s2 = array[i][0]
        delay_s3 = array[i][1]
        delay_s4 = array[i][2]
        delay_sx = array[i][3]
        
        
        axs[0, 0].hist(delay_s2,bins=bins,label = leg[i][0] + r"$%.2f^{+ %.2f}_{- %.2f} $" %(round(np.median(delay_s2),2), round(abs(np.median(delay_s2)-np.percentile(delay_s2,84)),2),round(abs(np.percentile(delay_s2,16)-np.median(delay_s2)),2)),density=True,histtype='step',color=color[i],linestyle=lin[i],linewidth=1.4)
        
        axs[0, 1].hist(delay_s3,bins=bins,label = leg[i][1] + r"$%.2f^{+ %.2f}_{- %.2f} $" %(round(np.median(delay_s3),2), round(abs(np.median(delay_s3)-np.percentile(delay_s3,84)),2),round(abs(np.percentile(delay_s3,16)-np.median(delay_s3)),2)),density=True,histtype='step',color=color[i],linestyle=lin[i],linewidth=1.4)
        
        axs[1, 0].hist(delay_s4,bins=bins,label = leg[i][2] + r"$%.2f^{+ %.2f}_{- %.2f} $" %(round(np.median(delay_s4),2), round(abs(np.median(delay_s4)-np.percentile(delay_s4,84)),2),round(abs(np.percentile(delay_s4,16)-np.median(delay_s4)),2)),density=True,histtype='step',color=color[i],linestyle=lin[i],linewidth=1.4)
        
        axs[1, 1].hist(delay_sx,bins=bins,label = leg[i][3] + r"$%.2f^{+ %.2f}_{- %.2f} $" %(round(np.median(delay_sx),2), round(abs(np.median(delay_sx)-np.percentile(delay_sx,84)),2),round(abs(np.percentile(delay_sx,16)-np.median(delay_sx)),2)),density=True,histtype='step',color=color[i],linestyle=lin[i],linewidth=1.4)
        
        
    axs[0, 0].set_xlabel('True - Estimated, S2-S1 delay', fontsize=14)
    axs[0, 0].set_ylabel('Probability Density', fontsize=14)
    axs[0, 0].grid(True)
    axs[0, 0].legend()  
    
    axs[0, 1].set_xlabel('True - Estimated, S3-S1 delay', fontsize=14)
    axs[0, 1].set_ylabel('Probability Density', fontsize=14)
    axs[0, 1].grid(True)
    axs[0, 1].legend()    
       
    axs[1, 0].set_xlabel('True - Estimated, S4-S1 delay', fontsize=14)
    axs[1, 0].set_ylabel('Probability Density', fontsize=14)
    axs[1, 0].grid(True)
    axs[1, 0].legend()  
    
    axs[1, 1].set_xlabel('True - Estimated, Sx-S1 delay', fontsize=14)
    axs[1, 1].set_ylabel('Probability Density', fontsize=14)
    axs[1, 1].grid(True)
    axs[1, 1].legend() 
    
    
    
    




###############################################################################################


               


arr = []
colors = []
legs = []

for name, thresh, c, legend in param_set:
    
    delay_s2, delay_s3, delay_s4, delay_sx, s2deleted, s3deleted, s4deleted, sxdeleted = delays(name, threshold = thresh)
    
    t = np.array([delay_s2,delay_s3,delay_s4,delay_sx],dtype=object)
    plot_multiple(np.array([[delay_s2, delay_s3, delay_s4, delay_sx]],dtype=object), [legend], [c],['solid'], bins = bin_p)   
    arr.append(t)
    colors.append(c)
    legs.append(legend)
    
    print('Critical failures are : S2', round(100*s2deleted,3), '%, S3', round(100*s3deleted,3), '%, S4', round(100*s4deleted,3), '%, Sx', round(100*sxdeleted,3), '%')
    if len(param_set) >= 19:
        raise RuntimeError("I won't let you continue, you have too much figures to show (more than 19), your pc might explode !")



#plot_multiple(arr,[r'5 cycles D$_2$',r'10 cycles D$_2$',r'15 cycles D$_2$'],['r','b','g'],bins = 100)

plot_multiple(arr,legs,colors,final_plot_style,bins = bin_p)

plt.show()














