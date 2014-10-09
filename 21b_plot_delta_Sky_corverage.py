''' 21-plot_Sky-coverage.py
===============================================
AIM:	Plots the Stray Light coverage difference between two 21-plot-Sky-coverage.py in terms of period of observation and accumulated observation time.

INPUT:	files: 	- <orbit_id>_misc/ : files from 21-plot-Sky-coverage.py
	variables: see section PARAMETERS (below)

OUTPUT:	in <orbit_id>_<SL_angle>figures/ : <orbit_id>_<threshold_obs_time>_<max_mag><_SAA?>_SL_coverage-other-<ID1-ID2>.png/pdf/eps

CMD:	21b-plot_Sky-coverage.py

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
import os
import matplotlib.cm as cm
import time

from resources.routines import *
from resources.TimeStepping import *

import parameters as param
import resources.constants as const
import resources.figures as figures
from resources.targets import *

from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter
from mpl_toolkits.basemap import Basemap
###########################################################################
### PARAMETERS
orbit_id_ref='700_25_AKTAR'
orbit_id_other='700_25_AKTAR'

# Give fname without the extension !
#sky_coverage_map_fname_ref = '301-sky_map-78-mag12-days10-to-40-accumulated' 
#sky_coverage_map_fname_other = '1000-sky_map-78-mag12-days10-to-40-accumulated'

#sky_coverage_map_fname_ref = '301-sky_map-50-mag9-days10-to-90-accumulated'
#sky_coverage_map_fname_other = '1000-sky_map-50-mag9-days10-to-90-accumulated'

sky_coverage_map_fname_ref = '700_25_AKTAR-sky_map-78-mag12-days10-to-16_noSAA_0.0pst-obs_time'
sky_coverage_map_fname_other = '700_25_AKTAR-sky_map-78-mag12-days10-to-16_SAA_0.0pst-obs_time'

# Nice plots?
fancy=True

# Save plots?
save = True

# Show figures ?
show = True

###########################################################################

# Formatted folders definitions
_, folder_figures_ref, folder_misc_ref = init_folders(orbit_id_ref)
_, _, folder_misc_other = init_folders(orbit_id_other)

data_ref = np.load(folder_misc_ref+sky_coverage_map_fname_ref+'.npz')
data_other = np.load(folder_misc_other+sky_coverage_map_fname_other+'.npz')

ra_grid = data_ref['ra_grid']
dec_grid = data_ref['dec_grid']
ticks = data_ref['dec_grid']

ref = data_ref['data_grid']
other = data_other['data_grid']
whereAreNaNs = np.isnan(ref);
ref[whereAreNaNs] = 0;
whereAreNaNs = np.isnan(other);
other[whereAreNaNs] = 0;
delta = ref-other
#delta[delta == 0]=np.nan

### Plotting
# transform 0 into no plotting in the data matrix

if fancy: figures.set_fancy()
fig = plt.figure()
ax=plt.subplot(111)
ax.set_aspect(2.)

plt.grid() 

CS = plt.contour( ra_grid,dec_grid,delta,colors='k')

plt.clabel(CS, inline=1,fmt='%d',colors='k', fontsize=12)

CS = plt.contourf( ra_grid ,dec_grid,delta,200,cmap=plt.cm.Paired)

plt.yticks(np.arange(-80, 100, 20.))

v = np.linspace(np.floor(np.nanmin(delta)),np.ceil(np.nanmax(delta)), 10, endpoint=True)
 
cbar = plt.colorbar(CS)
cbar.set_ticklabels(np.round(v,1))
cbar.set_label(r'$\mathrm{Days}$')

plt.xlabel('RA [deg]')
plt.ylabel('Dec [deg]')

if show: plt.show()

if save:
	fname = sky_coverage_map_fname_ref+'---'+sky_coverage_map_fname_other
	figures.savefig(folder_figures_ref+fname, fig, fancy)
	print 'saved as %s' % folder_figures_ref+fname

