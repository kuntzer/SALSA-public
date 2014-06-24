''' 20-plot_SL_coverage_restriction.py
===============================================
AIM:	Plots the Stray Light coverage restriction in %.

INPUT:	files: 	- <orbit_id>_misc/ : files from 12-<...>.py
	variables: see section PARAMETERS (below)

OUTPUT:	in <orbit_id>_<SL_angle>figures/ : <orbit_id>_<threshold_obs_time>_<max_mag><_SAA?>_SL_coverage.png/pdf/eps

CMD:	python 20-plot_SL_coverage_restriction.py

ISSUES:	<NONE KNOWN>

REQUIRES:- standard python libraries, specific libraries in resources/ (+ SciPy)
	 - Structure of the root folder:
	   * <orbit_id>_flux/ --> flux files
	   * <orbit_id>_figures/ --> figures
	   * <orbit_id>_misc/ --> storages of data

REMARKS: <NONE>
'''
###########################################################################
### INCLUDES
import numpy as np
import pylab as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm

from resources.routines import *
from resources.TimeStepping import *
import resources.constants as const
import resources.figures as figures

import time
from matplotlib import dates

from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter

###########################################################################
### PARAMETERS
# orbit_id
orbit_id = 702

# Show plots
show = True

# Save the picture ?
save = True

# Fancy plots ?
fancy = True

threshold_obs_time = 50

# Magnitudes to plot
mag = np.array([10.1,11.1,12.,12.1,12.2])
labels = [r'$10$',r'$11$',r'$12\ \mathrm{processed}$',r'$12$',r'$12\ \mathrm{No\ SAA}$']

###########################################################################
### INITIALISATION

if fancy: figures.set_fancy()

# Formatted folders definitions
folder_flux, folder_figures, folder_misc = init_folders(orbit_id)

############################################################################
### LOADS AND PLOTS

fig = plt.figure()
ax = plt.subplot(111)

for mag_max,label in zip(mag,labels):
	input_fname = 'cumultative_SL_forbidden_%d_mag_%02.1f.dat' % (threshold_obs_time, mag_max)
	print 'Loading %s' % input_fname
	data = np.loadtxt(folder_misc+input_fname)

	plt.plot(data[1:,0],data[1:,1],label=label)

	# convert epoch to matplotlib float format
t_ini, junk, minute_ini, junk = orbit2times(data[1,0],orbit_id)
junk, junk, junk, minute_end = orbit2times(data[-1,0],orbit_id)

labels = np.linspace(minute_ini, minute_end, 13) * 60. + const.timestamp_2018_01_01
labels = labels[1:]


plt.xlim([data[1,0], data[-1,0]])
ax.xaxis.set_major_locator(MultipleLocator((data[-1,0]-data[1,0]+1)/12))

plt.legend(loc=9,prop={'size':14}, mode="expand", ncol=5)
plt.ylabel(r'$\%\mathrm{\ of\ observable\ sky\ for\ which }\frac{F_\star}{F_{SL}} > T$')

# to human readable date
pre = map (time.gmtime, labels)
labels = map(figures.format_day, pre)

ax.set_xticklabels(labels)
fig.autofmt_xdate()

plt.grid(True)

plt.show()



# Save plot
if save:
	fname =  'cumultative_SL_forbidden_%d.dat' % (threshold_obs_time)
	figures.savefig(folder_figures+fname, fig, fancy)
	print 'saved as %s' % folder_figures+fname

if show: plt.show()
