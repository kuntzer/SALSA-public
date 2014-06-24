''' 17b-observation_regions.py
=========================
AIM:	Determines the time spent in one region of the sky.

INPUT:	files: 	- <orbit_id>_<SL_angle>misc/ephemerids_obs<transit_duration>h_<max_interruptions>inter_V<mag_max><_SAA?>.npz (from 17....py
	variables: see section PARAMETERS (below)

OUTPUT: <none>

CMD:	python 17b-observation_regions.py

ISSUES:	<none known>

REQUIRES:- standard python libraries, specific libraries in resources/ (+ SciPy)
	 - Structure of the root folder:
	   * <orbit_id>_flux/ --> flux files
	   * <orbit_id>_figures/maps/ --> figures
	   * <orbit_id>_misc/ --> storages of data

REMARKS: <none>
'''

###########################################################################
### INCLUDES
import numpy as np
import os
import matplotlib.cm as cm

from resources.routines import *
from resources.TimeStepping import *

import parameters as param
import resources.figures as figures
from resources.targets import *
import resources.geometry as geometry

from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter

###########################################################################
### PARAMETERS
# Orbit id
orbit_id = 301
apogee=700
perigee=700

# File name for the list of orbit file
orbits_file = 'orbits.dat'

# Minimum observable time for plots [h]
transit_duration = 6

# Maximum interruption time tolerated [min]
max_interruptions = 97

# Maximum visible magnitude
mag_max = 12.

# Take SAA into account?
SAA = True

# Print much information ?
verbose = False

# If set to True, then it will be observations of at least (period - max_interruptions)
# If set to False, then it is minimum (period - max_interruptions) minutes per orbit, 
# not necesseraly consecutive.
consecutive = False

# Factor in the SL post treatment correction ?
SL_post_treat = True

# Stop before last reshaping of data and saving (if consecutive == False)
early_stop = False

# Minimal # of days of obs (if consecutive == False), must be a list
nb_obs_days = [13]#range(10,110,10)#range(5,45,5)#

# Minimal minutes to be observed per orbit (if consecutive == False)
min_t_obs_per_orbit = 79

# This is a way to vary the results by multiplying the whole pst by a number.
# This is very easy as if the pst is multiplied by a constant, it can be taken out of the
# integral and only multplying the flux is equivalent to re-running all the simulations
pst_factor=1.

"""
Examples:
observation_region = geometry.Polygon(points=[(25, 2), (15, 3), (15, 7), (45, 7), (45, 2)])
observation_region = geometry.Interval(axis='delta',min_val=-30,max_val=30])
observation_region = geometry.Interval(axis='alpha',min_val=3.14,max_val=6.28],unit='rad')
The points are in degree by default.
"""
observation_region=geometry.Interval(axis="delta",max_val=0) # Southern sky


# File name for the input file (in a compressed binary Python format)
if SAA: note = '_SAA'
else: note = ''
if not pst_factor == 1.: note += '_%1.1fpst' % pst_factor
if SL_post_treat: note+= '_%4.3fSLreduction' % param.SL_post_treat_reduction

if not consecutive: note += '_cumul_'
skycoverage_fname = 'skycoverage_region_%dmin_V%3.1f%s.txt' % (min_t_obs_per_orbit,mag_max,note)

### INITIALISATION
# Formatted folders definitions
folder_flux, folder_figures, folder_misc = init_folders(orbit_id)
output=open(os.path.join(folder_misc,skycoverage_fname),"a") 
print >> output, observation_region
output.close()
for nb_obs_day in nb_obs_days:
	# File name for the input file (in a compressed binary Python format)
	if consecutive:
		input_fname = 'ephemerids_obs%dh_%dinter_V%3.1f%s.npz' % (transit_duration,max_interruptions,mag_max,note)
	else: 
		input_fname = 'ephemerids_%ddays_%dmin_V%3.1f%s.npz' % (nb_obs_day,min_t_obs_per_orbit,mag_max,note)
	#####################################################################################################################
	# CONSTANTS AND PHYSICAL PARAMETERS
	period = altitude2period(apogee, perigee)
	###########################################################################
	### INITIALISATION

	output=open(os.path.join(folder_misc,skycoverage_fname),"a") 

	print 'ORBIT ID:\t\t%d\nPST factor:\t\t%d\nMin Days of Coverage:\t%d\nmin_t_obs_per_orbit\t%d\nMAGNITIUDE:\t\t%02.1f' % (orbit_id,pst_factor,nb_obs_day,min_t_obs_per_orbit, mag_max)

	# loading data
	sys.stdout.write("Loading worthy targets...\t")
	sys.stdout.flush()

	worthy_targets = np.load(folder_misc+input_fname)
	obs_tot = worthy_targets['obs_tot']
	worthy_targets = worthy_targets['worthy_targets']

	count_T = 0
	count_all=0
	for tgt, obs in zip(worthy_targets, obs_tot):
		if obs == 0.: continue
		count_all+=1
		if observation_region.is_inside(tgt.Coordinates()):
			count_T+=1
	print
	print observation_region
	print 'Percentage of observation in region: %2.1f' % (float(count_T)/float((count_all))*100.), '%'
	
	print >> output, nb_obs_day,'\t',float(count_T)/float((count_all))*100.,'\t\%'

	output.close()