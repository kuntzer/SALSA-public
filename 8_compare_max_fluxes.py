''' 8-compare_max_fluxes.py
==============================
AIM:	Plots the comparison of maximal fluxes.

INPUT:	files: 	- <orbit_id>_misc/error_evolution.dat
	variables: see section PARAMETERS (below)

OUTPUT:	in all_figures/ : comparison for every orbit_iditude

CMD:	python 8-compare_max_fluxes.py

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

from resources.routines import *
from resources.TimeStepping import *
import parameters as param
import resources.figures as figures

from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter

###########################################################################
### PARAMETERS

# Error threshold
p = 0.1

# Flux limitation [ph/(px s)]
rqmt_flux = 1

# Show plots and detailled analysis ?
show = True

# Fancy plots ?
fancy = True

###########################################################################
### INITIALISATION
# File name fot the computed orbit file
error_file = 'error_evolution.dat'

# Formatted folders definitions
orbit_id = 1001
folder_flux, folder_figures, folder_misc = init_folders(orbit_id)
folder_figures= 'all_figures/'

if fancy: figures.set_fancy()

###########################################################################
### Load which orbits were computed
data_800 = np.loadtxt('%d_%ddeg_misc/%s' % (orbit_id,error_file), delimiter=',')
data_700 = np.loadtxt('%d_%ddeg_misc/%s' % (orbit_id,error_file), delimiter=',')
data_620 = np.loadtxt('%d_%ddeg_misc/%s' % (orbit_id,error_file), delimiter=',')

corrected_700 = np.zeros([np.shape(data_800)[0], 2])
corrected_620 = np.zeros([np.shape(data_800)[0], 2])

k = 0
for junk, o, junk, junk, sl, junk in data_800:
	if np.shape(data_700[data_700[:,0]==o,1])[0] > 0:
		corrected_700[k] = o, data_700[data_700[:,0]==o,4][0]/sl
	else:
		corrected_700[k] = o, data_700[data_700[:,0]<o,4][-1]/sl

	if np.shape(data_620[data_620[:,1]==o,1])[0] > 0:
		corrected_620[k] = o, data_620[data_620[:,0]==o,4][0]/sl
	else:
		corrected_620[k] = o, data_620[data_620[:,0]<o,4][-1]/sl

	k += 1


xx = data_800[:,1]/param.last_orbits[orbit_id]*365.
xx = figures.convert_date(xx)

fig=plt.figure()
ax=plt.subplot(111)
# zooms
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))

#pax.yaxis.set_major_locator(MultipleLocator(1.))
#pax.yaxis.set_minor_locator(MultipleLocator(0.5))

#ax.xaxis.set_major_locator(MultipleLocator(20.))

ax.xaxis.grid(True,'minor')
ax.yaxis.grid(True,'minor')
ax.xaxis.grid(True,'major',linewidth=2)
ax.yaxis.grid(True,'major',linewidth=2)

plt.plot(xx, corrected_620[:,1], 'b' , linewidth=2, label='620 km')
#pplt.plot(xx, corrected_620[:,1], 'Darkorange' , linewidth=2, label='Worst case')
#pplt.plot(xx, corrected_700[:,1], 'g', linewidth=3, label='Current INAF')
plt.plot(xx, corrected_700[:,1], 'r' , linewidth=2, label='700 km')

fig.autofmt_xdate()
plt.legend(loc=2)
#plt.ylim([0, 0.022])
plt.ylabel(r'$\mathrm{Relative\ maximum\ stray\ light\ flux\ to\ 800\ km}$')
# Saves the figure
fname = '%srelative_flux_%d' % (folder_figures,sl_angle)
figures.savefig(fname,fig,fancy)

fig=plt.figure()
ax=plt.subplot(111)
# zooms
ax.yaxis.set_major_locator(MultipleLocator(0.2))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))

#pax.yaxis.set_major_locator(MultipleLocator(1.))
#pax.yaxis.set_minor_locator(MultipleLocator(0.5))

#ax.xaxis.set_major_locator(MultipleLocator(20.))

ax.xaxis.grid(True,'minor')
ax.yaxis.grid(True,'minor')
ax.xaxis.grid(True,'major',linewidth=2)
ax.yaxis.grid(True,'major',linewidth=2)

xx = data_620[:,1]/param.last_orbits[620]*365.
xx = data_620[:,1]/param.last_orbits[800]*365.
xx = figures.convert_date(xx)
#pplt.plot(xx, data_620[:,4], 'b' , linewidth=2, label='620 km')
#pplt.plot(xx, data_620[:,4], 'Darkorange' , linewidth=2, label='Worst case')
plt.plot(xx, data_620[:,4], 'Indigo' , linewidth=2, label=r'$28^\circ\mathrm{\ INAF}$')
#plt.plot(xx, data_620[:,4], 'Darkorange' , linewidth=2, label='Worst case')

#pxx = data_700[:,1]/param.last_orbits[700]*365.
xx = data_700[:,1]/param.last_orbits[800]*365.
xx = figures.convert_date(xx)
#pplt.plot(xx, data_700[:,4], 'r', linewidth=3, label='700 km')
plt.plot(xx, data_700[:,4], 'g', linewidth=3, label='Current INAF')


xx = data_800[:,1]/param.last_orbits[800]*365.
xx = figures.convert_date(xx)
plt.plot(xx, data_800[:,4], 'k' , linewidth=2, label='800 km')
#pplt.plot(xx, data_800[:,4], 'k' , linewidth=2, label='RUAG')

xx = data_800[:,1]/param.last_orbits[800]*365.
xx = figures.convert_date(xx)

plt.xlim([xx[0],xx[-1]])
plt.plot([xx[0],xx[-1]], [rqmt_flux, rqmt_flux], color='r', lw=3)

fig.autofmt_xdate()
plt.legend(loc=9)

locs, labels = plt.yticks()
#plt.yticks(locs, map(lambda x: r"$%g$" % (float(x) * 1e2), locs))
#fig.text(0.12, 0.91,  r'$\times 10^{-2}$')
# r'$\times 10^{-2}$'


plt.ylabel(r'$\mathrm{Maximum\ stray\ light\ flux\ }\left[\frac{\mathrm{ph}}{\mathrm{px}\cdot\mathrm{s}}\right]$')
# Saves the figure
fname = '%sall_fluxes_800_com%d' % (folder_figures,sl_angle)
figures.savefig(fname,fig,fancy)

if show: plt.show()
