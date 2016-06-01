''' 17b-observation_regions.py
=========================
AIM:	Determines the time spent in one region of the sky.

INPUT:	files: 	- <orbit_id>_<SL_angle>misc/ephemerids_obs<transit_duration>h_<max_interruptions>inter_V<mag_max><_SAA?>.npz (from 17a...py
	variables: see section PARAMETERS (below)

OUTPUT: 'skycoverage_region_%dmin_V%3.1f%s.txt' % (min_t_obs_per_orbit,mag_max,note)

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
alt = 700
orbit_id = '6am_%d_5_conf4e' % alt
apogee=alt
perigee=alt

# File name for the list of orbit file
orbits_file = 'orbits.dat'

# Minimum observable time for plots [h]
transit_duration = None

# Maximum interruption time tolerated [min]
max_interruptions = None

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
early_stop = True

# Minimal # of days of obs (if consecutive == False), must be a list
nb_obs_days = [13]#range(10,17,1)#[13]#range(10,110,10)#

# Minimal minutes to be observed per orbit (if consecutive == False)
min_t_obs_per_orbit = 79

# This is a way to vary the results by multiplying the whole pst by a number.
# This is very easy as if the pst is multiplied by a constant, it can be taken out of the
# integral and only multplying the flux is equivalent to re-running all the simulations
pst_factor=1.

#Examples:
#observation_region = geometry.Polygon(points=[(25, 2), (15, 3), (15, 7), (45, 7), (45, 2)])
#observation_region = geometry.Interval(axis='delta',min_val=-30,max_val=30])
#observation_region = geometry.Interval(axis='alpha',min_val=3.14,max_val=6.28],unit='rad')
#The points are in degree by default.
observation_region=geometry.Interval(axis="delta",max_val=0) # Southern sky


# File name for the input file (in a compressed binary Python format)
if SAA: note = '_SAA'
else: note = ''
if not pst_factor == 1.: note += '_%1.1fpst' % pst_factor
if SL_post_treat: note+= '_%4.3fSLreduction' % param.SL_post_treat_reduction

if not consecutive: note += '_cumul_'
skycoverage_fname = 'skycoverage_region_%dmin_V%3.1f%s.txt' % (min_t_obs_per_orbit,mag_max,note)

### INITIALISATION
## Prepare grid
n_alpha = param.resx
n_delta = param.resy

ra_i = -np.pi
ra_f = np.pi

dec_i = -np.pi/2.
dec_f = np.pi/2.

ra_step = (ra_f-ra_i)/n_alpha
dec_step = (dec_f-dec_i)/n_delta

iterable = (ra_i + ra_step/2+ i*ra_step for i in range(n_alpha))
ras = np.fromiter(iterable, np.float)

iterable = (dec_i + dec_step/2+ i*dec_step for i in range(n_delta))
decs = np.fromiter(iterable, np.float)

ra_grid, dec_grid = np.meshgrid(ras, decs)

data_grid = np.zeros(np.shape(ra_grid))*np.nan

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
	data_grid*=np.nan

	output=open(os.path.join(folder_misc,skycoverage_fname),"a") 

	print 
	print 'ORBIT ID:\t\t%s\nPST factor:\t\t%d\nMin Days of Coverage:\t%d\nmin_t_obs_per_orbit\t%d (%1.4g%%)\nMAGNITIUDE:\t\t%02.1f' % (orbit_id,pst_factor,nb_obs_day,min_t_obs_per_orbit, min_t_obs_per_orbit/period*100., mag_max)

	print "Loadind from %s" % input_fname

	# loading data
	sys.stdout.write("Loading worthy targets...\t")
	sys.stdout.flush()

	worthy_targets = np.load(folder_misc+input_fname)
	obs_tot = worthy_targets['obs_tot']
	worthy_targets = worthy_targets['worthy_targets']

	count_T = 0
	count_all=0
	coords=[]
	obs_sky_region = 0
	obs_sky = 0
#	cc=[]
	for tgt, obs in zip(worthy_targets, obs_tot):
		if obs == 0.: continue
		count_all+=1
		rat, dect = tgt.Coordinates()
		obs_sky +=0.5/param.resx/param.resy*np.pi*np.cos(dect)
		if observation_region.is_inside(tgt.Coordinates()):
			count_T+=1

			obs_sky_region +=0.5/param.resx/param.resy*np.pi*np.cos(dect)

		id_ra, id_dec = find_nearest(ras+np.pi, tgt.Coordinates()[0]), find_nearest(decs, tgt.Coordinates()[1])
		data_grid[id_dec,id_ra]=obs
#		coords.append(np.asarray(tgt.Coordinates()))


	print observation_region
	print '%d cells on %d tot' % (count_T,count_all)
	print 'Percentage of observation in region: %2.1f' % ((obs_sky_region/obs_sky)*100.), '%'
	
	print >> output, nb_obs_day,'\t',(obs_sky_region/obs_sky)*100.,'\t\%'

	output.close()

	import pylab as plt	
	plt.figure()

	plt.contourf(np.degrees(ra_grid), np.degrees(dec_grid),data_grid)
	plt.grid()
	plt.xlim([-180,180])
	plt.ylim([-90,90])
	plt.title(input_fname)

#	coords=np.asarray(coords)
#	plt.scatter(coords[:,0],coords[:,1])

plt.show()

