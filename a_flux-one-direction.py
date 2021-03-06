''' a-flux-one-direction.py
=========================
AIM:	Perform basic statistics on the data and gets the maximal stray light flux for one orbit

INPUT:	files: 	- <height>_flux/flux_*.dat : data files
	variables: see section PARAMETERS (below)

OUTPUT:	<height>_figures/ : evolution of stray light

CMD:	python a-flux-one-direction.py

ISSUES:	<none known>

REQUIRES:- standard python libraries, specific libraries in resources/
	 - Structure of the root folder:
	   * <height>_flux/ --> flux files
	   * <height>_figures/ --> figures
	   * <height>_misc/ --> storages of data
	   * all_figures/ --> comparison figures

REMARKS: <none>
'''
################################################################################
import numpy as np
import pylab as plt

from resources.routines import *
from resources.TimeStepping import *
import parameters as param
import resources.constants as const
import resources.figures as figures

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino'],'size':14})
rc('text', usetex=True)

# Orbit in which to search
o = 41

orbitID = 1

t_ini, t_end, a_ini, a_end = orbit2times(o,orbitID)

sl = np.zeros(t_end-t_ini)

for minute in range(a_ini,a_end):
	ra, dec, S_sl = load_flux_file(minute, 'flux_', folder='%d_flux/' % orbitID)

	idr = np.where(np.abs(ra-1.806) < 0.1)[0]
	idd = np.where(np.abs(dec+0.55) < 0.1)[0]
	id_ = np.intersect1d(idr, idd)

	sl[minute-a_ini] = S_sl[id_]
	print S_sl[id_]


fig=plt.figure()
ax = plt.subplot(111)
fig.text(0.12, 0.91,  r'$\times 10^{-2}$')
fig.text(0.2, 0.8,  r'$\alpha=103^\circ,\ \delta=-31.5^\circ$')
plt.xlabel('Minutes in orbit %d [min]' % o)
plt.ylabel(r'$\mathrm{Mean\ stray\ light\ flux\ }\left[\frac{\mathrm{ph}}{\mathrm{px}\cdot\mathrm{s}}\right]$')
plt.grid()
plt.plot(sl*100,lw=2)
from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.yaxis.set_minor_locator(MultipleLocator(0.25))

ax.xaxis.set_major_locator(MultipleLocator(20.))
ax.xaxis.set_minor_locator(MultipleLocator(10.))

ax.xaxis.grid(True,'minor')
ax.yaxis.grid(True,'minor')
ax.xaxis.grid(True,'major',linewidth=2)
ax.yaxis.grid(True,'major',linewidth=2)
figures.savefig('all_figures/flux_fixed_direction',fig,True)
plt.show()
