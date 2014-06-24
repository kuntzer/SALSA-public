''' 4-error_evolution.py
=========================
AIM:	Similarly to 1-orbits_computed.py, checks the error evolution of the computed orbits

INPUT:	files: 	
		- all flux_*.dat files in <orbit_id>_flux/
		- <orbit_id>_misc/orbits.dat
		variables: see section PARAMETERS (below)

OUTPUT:	in <orbit_id>_misc/ : one file 'orbits.dat' containing the list

CMD:	python 4-error_evolution.py

ISSUES:	<none known>

REQUIRES:- standard python libraries, specific libraries in resources/
	 - Structure of the root folder:
	   * <orbit_id>_flux/ --> flux files
	   * <orbit_id>_figures/ --> figures
	   * <orbit_id>_misc/ --> storages of data
	   * all_figures/ --> comparison figures

REMARKS: Takes a long time to compute, see resources/routines.py for more expanations on compare_two_orbits().
'''

###########################################################################
### INCLUDES
import numpy as np
import pylab as plt
import time 

from resources.routines import *
from resources.TimeStepping import *
###########################################################################
### PARAMETERS
# Orbit id
orbit_id = 1001

# Error threshold
p = 0.1

# Show plots and detailled analysis ?
show = False

###########################################################################
### INITIALISATION
# File name for the computed orbit file
orbits_file = 'orbits.dat'
index_file = 'index.dat'

# Formatted folders definitions
folder_flux, folder_figures, folder_misc = init_folders(orbit_id)

# setup a few variables
shift = 0



###########################################################################
### Load which orbits were computed
start = time.time()

orbits = np.loadtxt(folder_misc+orbits_file,dtype='i4')

data = np.zeros([np.shape(orbits)[0]-1,6])

previous_part = -1
for ii, orbit_current in enumerate(orbits[:,0]):
	if ii == 0: continue

	preceeding = orbits[ii-1,0]

	shift = 1
	pp, sl = compare_two_orbits(preceeding, orbit_current, orbit_id, folder=folder_flux,shift=shift,return_max=True)
	pp_old = pp
	if pp > p : 
		shift = 2
		pp, sl = compare_two_orbits(preceeding, orbit_current, orbit_id, folder=folder_flux, shift = shift,return_max=True)
	if pp > p :
		shift = 0
		pp, sl = compare_two_orbits(preceeding, orbit_current, orbit_id, folder=folder_flux, shift = shift,return_max=True)
	if pp > p :
		pp = pp_old
		shift = 1

	data[ii-1] = preceeding, orbit_current, orbits[ii,1], pp, sl, shift

	percentage_done = round(float(ii)/float(np.shape(orbits)[0])*100,1)
	if show: print percentage_done, preceeding, orbit_current, pp, sl, shift
	else:
		sys.stdout.write( '\r%3.1f%%' % (percentage_done) )
		sys.stdout.flush()

header = 'ref,val,step,error,max_sl,shift'
fname = 'error_evolution.dat'
np.savetxt(folder_misc+fname,data,header=header, fmt='%4d,%4d,%2d,%1.8f,%2.8f,%d')

end = time.time()
elapsed_time = round((end-start)/60.,1)
print 'Done. Time needed : %3.1f minutes,' % elapsed_time
