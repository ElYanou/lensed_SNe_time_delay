# define here some set of distribtuion that you want to study
param = [("./sim/sim_results/tst.pkl", 30, "orangered",['','','','']),
         ("./sim/sim_results/tst2.pkl",30,"royalblue",['','','','']),
         ("./sim/sim_results/tst3.pkl",30,"seagreen",['','','','']),
         ]     

######################################################################################################
# define here the parameters
"""
:param bin_p: bin of the histogram

:param param_set: which set do you want to use, each element is composed of a path, a threshold a color and a legend

:param final_plot_style: linestyle of each distribution, useful when they superpose 
"""

bin_p = 100
param_set = param  
final_plot_style = ['solid', 'dashed', 'solid']
         
