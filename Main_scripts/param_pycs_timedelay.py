"""
:param sim_path: path to the set of light curves as pkl
    
:param file0_name: the name of the light curve 0 for which we have not the true time delays
   
:param data_folder: path to the folder of data
    
:param numbml_p: polynomial microlensing = 1, else microlensing is full spline
    
:param degree_p: degree of polynomial microlensing  
    
:param nb_cycle_p: numbre of times the optimizer is run
    
:param do_you_want_ml_p: boolean, do you want to add microlensing to your light curves ?
    
:param do_you_display_p: boolean, do you want to display the optimization ? If True multiprocessing is off

:param kn_p: spline knotstep
    
:param writer: folder in which we store the txt files, if the run is complete the folder is erased if not you should erase it manually. It allows to control the time delay differences during the run.
    
:param sttpos_p: time delay starting points
    
:param magstt_p: magnitude shift starting points
    
:param final_name: the path in which the results are stored
"""


sim_path = './sim/sim_F125W'
file0_name = 'lcs_F125W_no_opt_0'
data_folder = './../Chabrier/Chabrier'
numbml_p = 1
degree_p = 0
nb_cycle_p = 10
do_you_want_ml_p = True   
do_you_display_p = False   #Should we display each light curves optimization
kn_p = 110
writer = './sim/sim_results/txt_timedelays'
sttpos_p = [0., -31., 26., -20., -390.]
magstt_p = [0., 0., 0., -1.2, -1.3]
final_name = "./sim/sim_results/tst.pkl"
