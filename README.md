# lensed_SNe_time_delay
TP IVb project at EPFL

We provide here a couple of scripts to work on a dataset that we had, you should have a similar dataset to make the scripts work well. 

lc_maker.py : transforms the dataset to PyCS3 light curves object.

pycs_timedelay.py : fit the time delays from the light curves, create a pkl file with all the time delay differences stored in it.

timedelay_plot.py : plot the distributions of the time delay differences.

sncosmo_fit.py : fit time delays without microlensing using template of SNe from SNCosmo.

module_sncosmo_fit.py : basis to adapt PyCS3 to lensed SNe.
We have 2 strategies to fit :
- Merge the light curves and fit this new superlight curve
- Adapt the parameters of the template used on each light curve
