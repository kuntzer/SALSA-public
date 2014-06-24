''' 5-statistics-error.py
=========================
AIM:	Perform basic statistics on the data and gets the maximal stray light flux for one orbit

INPUT:	files: 	- <orbit_id>_misc/orbits.dat
		variables: see section PARAMETERS (below)

OUTPUT:	in <orbit_id>_misc/ : file one stat file
	in <orbit_id>_figures/ : error evolution, max. stray light evolution

CMD:	python 5-statistics-error.py

ISSUES:	<none known>

REQUIRES:- standard python libraries, specific libraries in resources/
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
import os

from resources.routines import *
from resources.TimeStepping import *
import parameters as param
import resources.figures as figures

from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter

###########################################################################
### PARAMETERS
# Orbit id
orbit_id = 1001

# Error threshold
p = 0.1

# Flux limitation [ph/(px s)]
rqmt_flux = 1

# File name for the output data file (same as in 2-statistics-step.py)
data_file = 'statistics-error.dat'

# Show plots and detailled analysis ?
show = True

# Fancy plots ?
fancy = True

###########################################################################
### INITIALISATION
# File name for the computed orbit file
error_file = 'error_evolution.dat'

# Formatted folders definitions
folder_flux, folder_figures, folder_misc = init_folders(orbit_id)

if fancy: figures.set_fancy()

if os.path.isfile(folder_misc+data_file):
	os.remove(folder_misc+data_file)

f = open(folder_misc+data_file,'w')
###########################################################################
### Load which orbits were computed
data = np.loadtxt(folder_misc+error_file, delimiter=',')
# Data type:
# ref,val,step,error,max_sl,shift

### Error evolution
print >> f, '# ERRORS'
print >> f, '# ! All errors are normalised to 1'
print >> f, '# ! ie 1.0 = 100%'
print >> f, 'error_max:', np.amax(data[:,3])
print >> f, 'error_min:', np.amin(data[:,3])
print >> f, 'error_mean:', np.mean(data[:,3])
print >> f, 'error_std:', np.std(data[:,3])

fig=plt.figure()
ax=plt.subplot(111)

ax.yaxis.set_major_locator(MultipleLocator(5))
ax.yaxis.set_minor_locator(MultipleLocator(1))

ax.xaxis.grid(True,'minor')
ax.yaxis.grid(True,'minor')
ax.xaxis.grid(True,'major',linewidth=2)
ax.yaxis.grid(True,'major',linewidth=2)

xx = data[:,1]/param.last_orbits[orbit_id]*365.
xx = figures.convert_date(xx)
plt.plot(xx, data[:,3]*100, linewidth=1.5)
plt.plot([xx[0],xx[-1]], [p*100., p*100.], color='r', lw=3)
fig.autofmt_xdate()

plt.ylim([0, 15])

plt.ylabel(r'$\mathrm{Error\ to\ previous\ step\ [\%]}$')

# Saves the figure
fname = '%serror_evolution_%d_%d' % (folder_figures,orbit_id,sl_angle)
figures.savefig(fname,fig,fancy)

############ STRAY LIGHT
print >> f, '# STRAY LIGHT'

# Get the direction of minimum stray light
id_min = find_nearest(data[:,4],np.amin(data[np.where(data[:,4]>0)]))
orbit_max = data[id_min, 1]

time_min, ra_min, dec_min, sl_min = find_direction_flux(orbit_max, orbit_id,find='min', folder=folder_flux)
print >> f, 'min:', sl_min
print >> f, 'minute_min:', time_min
print >> f, 'RA_min:', ra_min
print >> f, 'DEC_min:', dec_min

print >> f, 'mean:', np.mean(data[:,4])
print >> f, 'stddev:', np.std(data[:,4])

# Get the direction of maximum stray light
id_max = find_nearest(data[:,4],np.amax(data[:,4]))
orbit_max = data[id_max, 1]

time_max, ra_max, dec_max, sl_max = find_direction_flux(orbit_max, orbit_id, folder=folder_flux)
print >> f, 'max:', np.amax(sl_max)
print >> f, 'minute_max:', time_max
print >> f, 'RA_max:', ra_max
print >> f, 'DEC_max:', dec_max

print >> f, 'mean:', np.mean(data[:,4])
print >> f, 'stddev:', np.std(data[:,4])

print >> f, 'orbit_above_rqmt:', np.shape(data[np.where(data[:,4]>rqmt_flux)])[0]
print >> f, 'total_orbits:', np.shape(data)[0]

### Maximal sl
fig=plt.figure()
ax=plt.subplot(111)

ax.yaxis.set_major_locator(MultipleLocator(0.2))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))

ax.xaxis.grid(True,'minor')
ax.yaxis.grid(True,'minor')
ax.xaxis.grid(True,'major',linewidth=2)
ax.yaxis.grid(True,'major',linewidth=2)

plt.plot(xx, data[:,4], linewidth=3)
plt.plot([xx[0],xx[-1]], [rqmt_flux, rqmt_flux], color='r', lw=3)
fig.autofmt_xdate()

plt.ylabel(r'$\mathrm{Maximum\ stray\ light\ flux\ }\left[\frac{\mathrm{ph}}{\mathrm{px}\cdot\mathrm{s}}\right]$')

# Saves the figure
fname = '%sstray_light_flux_%d_%d' % (folder_figures,orbit_id,sl_angle)
figures.savefig(fname,fig,fancy)
####################################################################
fig=plt.figure()
ax=plt.subplot(111)
# zooms
ax.yaxis.set_major_locator(MultipleLocator(0.1))
ax.yaxis.set_minor_locator(MultipleLocator(0.02))

#ax.xaxis.set_major_locator(MultipleLocator(20.))

ax.xaxis.grid(True,'minor')
ax.yaxis.grid(True,'minor')
ax.xaxis.grid(True,'major',linewidth=2)
ax.yaxis.grid(True,'major',linewidth=2)

plt.plot(xx, data[:,4], linewidth=3)
fig.autofmt_xdate()

plt.ylim([0, 0.2])
plt.ylabel(r'$\mathrm{Maximum\ stray\ light\ flux\ }\left[\frac{\mathrm{ph}}{\mathrm{px}\cdot\mathrm{s}}\right]$')
# Saves the figure
fname = '%sstray_light_flux_zoom_%d_%d' % (folder_figures,orbit_id,sl_angle)
figures.savefig(fname,fig,fancy)

if show: plt.show()

f.close()
