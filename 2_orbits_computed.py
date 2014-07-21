''' 2-orbits_computed.py
=========================
AIM:	Verify which orbits were actually computed by the pipeline

INPUT:	files: 	- all flux_*.dat files in <orbit_id>_flux/
	variables: see section PARAMETERS (below)

OUTPUT:	in <orbit_id>_misc/ : one file 'orbits.dat' containing the list

CMD:	python 2-orbits_computed.py	

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
import os.path

from resources.routines import *
from resources.TimeStepping import *
###########################################################################
### PARAMETERS

# Orbit id same as in 0-*.py
orbit_id = '620_35_AKTAR'

# First orbit in data set
orbit_ini = 1

# Last orbit to look for
orbit_end = minute2orbit(1440*365+1,orbit_id)

# File name for the output data file
data_file = 'orbits.dat'

###########################################################################
### INITIALISATION
# Formatted folders definitions
folder_flux, folder_figures, folder_misc = init_folders(orbit_id)

predictor = int((orbit_end - orbit_ini) / 10. * 4.)

orbits = np.zeros(predictor*2).reshape(predictor, 2)

kk = 0
orbit_old = 0
###########################################################################
### Look for computed orbits

for orbit in range(orbit_ini, orbit_end+1):
	# lookup the start and end time for the orbit
	t_ini, t_end, a_ini, a_end = orbit2times(orbit,orbit_id)

	# check that minute 0 and 60 exists (they should!)
	fname = '%sflux_%d.dat' % (folder_flux, a_ini)
	if not os.path.isfile(fname): continue

	fname = '%sflux_%d.dat' % (folder_flux, a_ini+60)
	if not os.path.isfile(fname): continue

	if kk > 0:
		step = orbit-orbit_old
	else:
		step = 1

	try:
		orbits[kk,0] = orbit
		orbits[kk,1] = step
	except IndexError:
		orbits = np.vstack([orbits, [orbit, step]])

	orbit_old = orbit

	print 'Orbit %4d was computed for orbit ID %s, step %d' % (orbit, orbit_id, step)

	kk += 1

# Remove extra line in the array
orbits = orbits[orbits[:,1]>0]

np.savetxt(folder_misc+data_file, orbits, fmt='%d')
