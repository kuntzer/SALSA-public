''' Plot a Fourrier transform of the SAA interruption '''

###########################################################################
### INCLUDES
import numpy as np
import pylab as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm

from resources.routines import *
from resources.TimeStepping import *
import parameters as param
import resources.constants as const
import resources.figures as figures

import datetime
import time
from matplotlib import dates
import os

from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter

from numpy.fft import fft, fftfreq

###########################################################################
### PARAMETERS
# orbit_id
orbit_id = 800

# Show plots
show = True

# Save the picture ?
save = False

# Fancy plots ?
fancy = True

# Line of sight (LOS) to Earth's limb angle
sl_angle = 35

# First minute in data set
minute_ini = 0

# last minute in data set
minute_end = 1440*365

# First minute in data set
orbit_ini = 1

# last minute in data set
orbit_end = param.last_orbits[orbit_id]

# File name for the list of orbit file
orbits_file = 'orbits.dat'

###########################################################################
### INITIALISATION

if fancy: figures.set_fancy()

# Formatted folders definitions
folder_flux, folder_figures, folder_misc = init_folders(orbit_id)

############################################################################

SAA_data = np.loadtxt('resources/SAA_table_%d.dat' % orbit_id, delimiter=',')


# Formatted folders definitions
folder_flux = '%d_%ddeg_flux/' % (orbit_id, sl_angle)
folder_figures= '%d_%ddeg_figures/' % (orbit_id, sl_angle)
folder_misc	= '%d_%ddeg_misc/' % (orbit_id, sl_angle)

if not os.path.isdir(folder_figures):
	print '\tError: figure folder %s does not exists.' % (folder_figures)
#	exit()


times = np.loadtxt('resources/minute_table_%d.dat' % orbit_id, delimiter=',',dtype='Int32')
list_orbit = np.zeros( [orbit_end-orbit_ini+1, 4] )
id_o = 0
for orbit_current in range(orbit_ini, orbit_end+1):
	t_ini, t_end, a_ini, a_end = fast_orbit2times(times,orbit_current,orbit_id) 
	list_orbit[id_o] = orbit_current, a_ini, a_end, 0
	id_o += 1

count = 0
SAA = np.zeros( [np.shape(SAA_data)[0]/2, 2] )
for off, on in zip(SAA_data[::2,0], SAA_data[1::2,0]):
	off_orbit = list_orbit[list_orbit[:,1]<off][-1]
	on_orbit =  list_orbit[list_orbit[:,2]>on][0]

	if off_orbit[0] == on_orbit[0]:
		list_orbit[list_orbit[:,0]==on_orbit[0],3] += on-off+1
	else:
		try: 
			list_orbit[list_orbit[:,0]==off_orbit[0],3] += off_orbit[2]-off+1
		except: pass
		try: 
			list_orbit[list_orbit[:,0]==on_orbit[0],3] += on-on_orbit[1]+1
		except: pass


list_orbit[:,3]/=(list_orbit[:,2]-list_orbit[:,1])/100.

############################################################################
### PLOTS

signal = list_orbit[:,3]

import scipy
import scipy.fftpack
FFT = abs(scipy.fft(signal))
freqs = scipy.fftpack.fftfreq(signal.size, (list_orbit[1,0]-list_orbit[0,0]))

plt.subplot(211)
plt.plot(list_orbit[:,0], signal)
plt.subplot(212)
plt.plot(freqs,20*scipy.log10(FFT),'x')







fig = plt.figure()
ax = plt.subplot(111)

plt.plot(list_orbit[:,0],list_orbit[:,3])

print 'Mean', np.mean(list_orbit[:,3])
print 'Max', np.max(list_orbit[:,3])

plt.show()
#
exit()

