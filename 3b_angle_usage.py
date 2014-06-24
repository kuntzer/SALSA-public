''' 3b_angle_usage.py
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

REMARKS: <none>
'''

###########################################################################
### INCLUDES
import numpy as np
import pylab as plt

from resources.routines import *
from resources.TimeStepping import *
import parameters as param
import resources.constants as const
import resources.figures as figures

###########################################################################

orbit_id = 1001
sl_angle = 35

fancy = True
show = True
save = True

orbit_ini = 1
orbit_end = 101

n_sampling = 55

# File name for the computed orbit file
orbits_file = 'orbits.dat'
###########################################################################

if fancy: figures.set_fancy()

# Formatted folders definitions
folder_flux, folder_figures, folder_misc = init_folders(orbit_id)

orbits = np.loadtxt(folder_misc+orbits_file,dtype='i4')

orbits = orbits[:,0]
orbits = orbits[np.where(orbits>=orbit_ini)]
orbits = orbits[np.where(orbits<=orbit_end)]

orbit_current = orbit_ini

orbit_calls = np.zeros([np.size(orbits),n_sampling])
total_calls = np.zeros(n_sampling)

for ii, orbit_current in enumerate(orbits):
	filename = ('%sangles_%d.dat') % (folder_flux, orbit_current)
	data=np.loadtxt(filename)

#	print data.sum(axis=0)/data.sum()

	orbit_calls[ii] = data.sum(axis=0)
	total_calls += data.sum(axis=0)

tot = total_calls.sum()/100

stddev = 3.*orbit_calls.std(axis=0)/tot
#print stddev
#print '*'*23
#print total_calls
#print total_calls.sum()

angles = np.linspace(sl_angle,89,n_sampling)

fname = '%sangle_usage_%d_%d_%d-%d.dat' % (folder_misc,orbit_id,sl_angle,orbit_ini,orbit_end)

tosave = np.zeros([np.size(angles),2])
ii = 0
for a,c in zip(angles,total_calls/tot):
	tosave[ii,0] = a
	tosave[ii,1] = c
	ii += 1
np.savetxt(fname,tosave,fmt='%d %1.6f')

fig, ax = plt.subplots(1)

plt.xlabel(r'$\theta\ \mathrm{Angular}\ \mathrm{distance}\ \mathrm{to}\ \mathrm{limb}\ \mathrm{[deg]}$')
plt.ylabel('\% of calls')
plt.plot(angles, total_calls/tot)

ax.fill_between(angles, total_calls/tot-stddev, total_calls/tot+stddev, facecolor='yellow', alpha=0.5,
                label='3 sigma range')

plt.grid()

if show: plt.show()

# Saves the figure
if save: 
	fname = '%sangle_usage_%d_%d_%d-%d' % (folder_figures,orbit_id,sl_angle,orbit_ini,orbit_end)
	figures.savefig(fname,fig,fancy)
