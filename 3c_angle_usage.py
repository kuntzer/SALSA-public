''' 3c_angle_usage.py
=========================
AIM:	Plots the diagonistic angle usage of the PST in SALSA. Requires the monitor_angle_usage=True in 1_compute_<p>.py and log_all_data = .true. in straylight_<orbit_id>_<p>/CODE/parameter. 

INPUT:	files: 	- <orbit_id>_misc/orbits.dat
				- <orbit_id>_flux/angles_<orbit_number>.dat
	variables: see section PARAMETERS (below)

OUTPUT:	in <orbit_id>_misc/ : file one stat file
		in <orbit_id>_figures/ : step distribution, step in function of time

CMD:	python 3b_angle_usage.py

ISSUES:	<none known>

REQUIRES:- LATEX, epstopdf, pdfcrop, standard python libraries, specific libraries in resources/
	 - Structure of the root folder:
	   * <orbit_id>_flux/ --> flux files
	   * <orbit_id>_figures/ --> figures
	   * <orbit_id>_misc/ --> storages of data
	   * all_figures/ --> comparison figures

REMARKS: This is a better version than the 3b_angle_usage.py
'''

###########################################################################
### INCLUDES
import numpy as np
import pylab as plt

from resources.routines import *
from resources.TimeStepping import *
import resources.figures as figures
from matplotlib import cm

###########################################################################
orbit_id = 704
sl_angle = 35

fancy = True
show = True
save = True

# Bins and their legends
orbit_ini = [1,441,891,1331,1771,2221,2661,3111,3551,3991,4441,4881,1]
orbit_end = [441,891,1331,1771,2221,2661,3111,3551,3991,4441,4881,5322,5322]
legends = ['Jan','Feb','Mar','Avr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Year']

###########################################################################

if fancy: figures.set_fancy()

# Formatted folders definitions
folder_flux, folder_figures, folder_misc = init_folders(orbit_id)

fig, ax = plt.subplots(1)

ii = 0.
size = len(orbit_ini)
minv=100
maxv=0
for ini, end,label in zip(orbit_ini,orbit_end,legends):
	print ini, end, label
	c = cm.rainbow(ii/float(size))
	fname = '%sangle_usage_%d_%d_%d-%d.dat' % (folder_misc,orbit_id,sl_angle,ini,end)
	values = np.loadtxt(fname)

	plt.plot(values[:,0], values[:,1],label=label, lw=2, c=c)

	if np.min(values[:,1]) < minv: minv = np.min(values[:,1])
	if np.max(values[:,1]) > maxv: maxv = np.max(values[:,1])

	ii += 1
plt.ylim( [ np.floor(minv), np.ceil(maxv) ] )


plt.xlabel(r'$\theta\ \mathrm{Angular}\ \mathrm{distance}\ \mathrm{to}\ \mathrm{limb}\ \mathrm{[deg]}$')
plt.ylabel('\% of calls')
plt.legend(loc=2,prop={'size':14}, ncol=2)

plt.grid()

if show: plt.show()

# Saves the figure
if save: 
	fname = '%stot_angle_usage_%d_%d' % (folder_figures,orbit_id,sl_angle)
	figures.savefig(fname,fig,fancy)
