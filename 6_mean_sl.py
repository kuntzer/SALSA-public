''' 6-mean_sl.py
=========================
AIM:	Compute the mean of the stray light in the sense of:
	- Mean across all points then every minutes
	- Mean of max value
	- Mean of direction of maximum (in the orbit)

INPUT:	files: 	- <orbit_id>_misc/orbits.dat
		- <orbit_id>_flux/flux_*.dat
	variables: see section PARAMETERS (below)

OUTPUT:	in <orbit_id>_misc/ : stat files

CMD:	python 6-mean_sl.py

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

from resources.routines import *
from resources.TimeStepping import *
###########################################################################
### PARAMETERS
# Orbit id
orbit_id = 1001

###########################################################################
### INITIALISATION
# File name fot the computed orbit file
orbits_file = 'orbits.dat'
index_file = 'index.dat'

# Formatted folders definitions
folder_flux, folder_figures, folder_misc = init_folders(orbit_id)


###########################################################################
### Load which orbits were computed
start = time.time()

orbits = np.loadtxt(folder_misc+orbits_file,dtype='i4')

data = np.zeros([np.shape(orbits)[0]-1,6])
mean = np.zeros([np.shape(orbits)[0],2])
mean_max = np.zeros([np.shape(orbits)[0],2])
mean_maxdir = np.zeros([np.shape(orbits)[0],2])

previous_part = -1
for ii, orbit_current in enumerate(orbits[:,0]):
#	print orbit, orbits[ii+1]

	t_ini, t_end, a_ini, a_end = orbit2times(orbit_current,orbit_id)
	mean_orbit = 0
	mean_max_orbit = 0
	mean_maxdir_orbit = 0
	nb_points = 0	
	sl_max = 0.0

	t_ini = np.ceil(t_ini) # not sure it's useful.
	minute = t_ini
	###########################################################################
	# Iterate on every time.
	while (minute <= t_end):
		# initialise the array for the current minute (ie. make sure nothing is left in the table from last time step.

		try:
		# Try to load the fluxes for a given minute (minute goes from 0 to period whereas a_ini is the absolute time from 0:0:0.0 1/1/2018 in min
			ra, dec, S_sl = load_flux_file(int(minute+a_ini), 'flux_', folder=folder_flux)

		except IOError:
		# If not found it means that the satellite is (propably) over the SAA. Skip this step.
			minute +=1
			continue	

		mean_orbit += np.mean(S_sl)
		mean_max_orbit += np.amax(S_sl)
		if (sl_max < np.amax(S_sl)): 
			sl_max = np.amax(S_sl)
			index_max = S_sl.argmax()
		nb_points += 1

		minute += 1
	mean[ii] = orbit_current, mean_orbit/nb_points
	mean_max[ii] = orbit_current, mean_max_orbit/(float(t_end)-float(t_ini))
	minute = t_ini
	# second run to get max in axe 
	skipped_minute = 0
	while (minute <= t_end):
		# initialise the array for the current minute (ie. make sure nothing is left in the table from last time step.

		try:
		# Try to load the fluxes for a given minute (minute goes from 0 to period whereas a_ini is the absolute time from 0:0:0.0 1/1/2018 in min
			ra, dec, S_sl = load_flux_file(int(minute+a_ini), 'flux_', folder=folder_flux)

		except IOError:
		# If not found it means that the satellite is (propably) over the SAA. Skip this step.
			minute +=1
			continue	

		try: mean_maxdir_orbit += S_sl[index_max]
		except : skipped_minute += 1

		minute += 1

	mean_maxdir[ii] = orbit_current, mean_maxdir_orbit/(float(t_end)-float(t_ini)-skipped_minute)
	#TODO: take direction of max flux for this orbit
	print orbit_current

header = 'orbit,mean'
fname = 'mean_sl.dat'
np.savetxt(folder_misc+fname,mean,header=header, fmt='%4d,%g')

header = 'orbit,mean_max'
fname = 'mean_max_sl.dat'
np.savetxt(folder_misc+fname,mean_max,header=header, fmt='%4d,%g')

header = 'orbit,mean_maxdir'
fname = 'mean_maxdir_sl.dat'
np.savetxt(folder_misc+fname,mean_maxdir,header=header, fmt='%4d,%g')

end = time.time()
elapsed_time = round((end-start)/60.,1)
print 'Done. Time needed : %3.1f minutes,' % elapsed_time
