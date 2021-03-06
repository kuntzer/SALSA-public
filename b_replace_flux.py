''' Replace files in <id>deg_flux/, see parameters '''

###########################################################################
### INCLUDES
import numpy as np
import os
import shutil

from resources.TimeStepping import *

###########################################################################
### PARAMETERS

# Orbit id
alt = 701

# Boundary orbits
orbit_ini = 4445
orbit_end = 4475

# name of directory to move from
path_of_new_flux = 'orbitID_<p>/'

# name of directory to move to
path_to_trash = 'orbitID_trash/'

###########################################################################
### PARAMETERS

message = 'Are you sure to replace orbits %d through %d ? (y/N) ' % (orbit_ini, orbit_end)
var = raw_input(message)
if not var == 'Y' and not var == 'y': 
	print 'Aboarting...'
	exit()




file_flux = 'flux_'

junk, junk, a_ini, junk = orbit2times(orbit_ini,alt)
junk, junk, junk, a_end = orbit2times(orbit_end,alt)

folder_flux = '%d_flux/' % (alt)

i = 0
k = 0
for minute in range(a_ini, a_end+1):
	to_be_replaced = '%s%s%d.dat' % (folder_flux, file_flux, minute)
	if os.path.exists(to_be_replaced) :
		k += 1
		os.system('mv %s %s' % (to_be_replaced, path_to_trash))

	replacement = '%s%s%d.dat' % (path_of_new_flux, file_flux, minute)
	if os.path.exists(replacement) :
		i += 1
		os.system('cp %s %s' % (replacement, folder_flux))
#	os.system('cp %s %s' % (replacement, folder_flux))

print k, 'files moved to trash and', i, 'copied to replace'
