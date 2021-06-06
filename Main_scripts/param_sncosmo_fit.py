from module_sncosmo_fit import *

"""
The code implements 5 tasks: 
- a multirun process to compute time-delay differences
- an optimization over a set of light curves and a plot as PyCS3 with time and mag shifts
- a visual template research, plot the reference light curve with optimized given model for a few templates (<19) 
- a visual time axis scale factor research, plot the reference light curve with optimized given model for a few scale factor (<19) 
- a more complete time axis scale factor research, plot of the chi squares (you can use more at values)

"""
# The booleans to do or not the corresponding tasks
mul_processing = False
single_run = False
template_search = False
at_visual_search = True
at_search = True

# Parameters that might be used by on or several tasks, numbers are tasks that use them, order is the same as above
"""
: param template_p : template to use / 1,2,4,5

: param folder_p : path to data folder / 1

: param z_p : redshift / 1,2,3,4,5

: param Filter_p : filter name as implemented in sncosmo / 1,2,3,4,5

: param filter_dat : filter name of your data, might be irrelevant depending on how your dataset is built, you should them modify the code / 1

: param path_stockage_p : path to a folder in which you store temporarily your time-delay differences, the folder must not exist yet / 1

: param startpts_p : starting points, time shifts and flux scale factors / 1,2

: param optimizer_p : optimizer to use (one of the 4) / 1,2

: param t_in_p : initial guess for {at}, must be given accurately if solu is None / 1,2

: param solu_p : true {at} to use, if not None at isn't fitted / 1,2 and must not be None for 3,4,5

: param nit_p : number of run of optimizers (some of them don't use it) / 1,2

: param th1_p : threshold for time delay guess (bounds are (guess-th1_p,guess+th1_p)) / 1,2

: param th2_p : threshold for flux scale guess (bounds are (guess-th2_p,guess+th2_p)). I would recommend to always avoid to have 0 in the interval / 1,2

: param bounds_t_p : bounds of {at}, must be tight or the fit of {at} will be irrelevant / 1,2

: param final_path : path to store the whole set of time-delay differences / 1

: param sim_path : path to the folder where the light curves as pkl are stored / 1

: param name_0 : name of light curve to delete / 1

: param file_list_p : list of light curve files, should be .dat, might need a change depending on how your dataset is built. First file is reference one. / 2,3,4,5

: param t_p : time range of the light curves, should be a linspace of at least 1000 points  / 2

: param template_list_p : list of template. should be small (<19) / 3

: param at_list_p : list of {at}, should be small (<19) / 4

: param at_big_list_p : list of {at} can be larger (>1000) / 5

"""
################### PATHS ######################
folder_p = './../Chabrier/Chabrier'  
final_path = "./sim/sim_results/tst3.pkl"
sim_path = './sim/sim_F125W'
path_stockage = './sim/sim_results/txt_timedelays_sncosmo'
name_0 = 'phot_dolphot_2_1.0000_True_1_0'
file_list_p = ['../Chabrier/Chabrier/phot_dolphot_2_1.0000_True_1_60/phot_s1_WFC3_IR_F125W.dat', '../Chabrier/Chabrier/phot_dolphot_2_1.0000_True_1_60/phot_s2_WFC3_IR_F125W.dat', '../Chabrier/Chabrier/phot_dolphot_2_1.0000_True_1_60/phot_s3_WFC3_IR_F125W.dat', '../Chabrier/Chabrier/phot_dolphot_2_1.0000_True_1_60/phot_s4_WFC3_IR_F125W.dat', '../Chabrier/Chabrier/phot_dolphot_2_1.0000_True_1_60/phot_sx_WFC3_IR_F125W.dat']

################## MODEL AND OPTIMIZER #############
template_p = 'snf-2011fe'  
template_list_p = ['salt2', 'v19-1987a','v19-2006aa']
z_p = 1.49     
Filter_p = 'f125w'
filt_dat = 'F125W'
startpts_p = [31., 1, -26., 1, 20., 0.5, 390., 0.5]
optimizer_p = optimize_simple
nit_p = 1
th1_p = 30
th2_p = 0.4

################## TIME AXIS SCALE #################
t_in_p = 1
solu_p = 1/(z_p+1)**2
bounds_t_p = None

################## PLOT #####################
t_p = np.linspace(50000,60000,10000)
at_list_p = [0.16,0.17,0.18,0.19,0.2,0.21]
at_biglist_p = np.linspace(0.16,0.2,10)










