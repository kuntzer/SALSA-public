''' 15-observation_fixed_direction
===============================================
AIM:	Similar to 14-<...>.py, but for only one traget.

INPUT:	files: 	- <orbit_id>_misc/orbits.dat
		- <orbit_id>_flux/flux_*.dat
	variables: see section PARAMETERS (below)

OUTPUT:	in <orbit_id>_figures/ : (see below for file name definition)

CMD:	python 15-observation_fixed_direction

ISSUES: ! DOES NOT WORK !

REQUIRES:- standard python libraries, specific libraries in resources/ (+ SciPy)
	 - BaseMap --> http://matplotlib.org/basemap/
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
import matplotlib.cm as cm
import time

from resources.routines import *
from resources.TimeStepping import *
from resources.targets import *

import parameters as param
import resources.constants as const
import resources.figures as figures
import time
from matplotlib import dates

from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter

###########################################################################
### PARAMETERS
# Name of object of interest OI:
OI = 'BD-082823'

# orbit_iditude of the orbit in km
orbit_id = 701
apogee=700
perigee=700

# First minute analysis
minute_ini = 30.*1440.

# Last minute to look for
minute_end = 50.*1440.

# Include SAA ?
SAA = False

# Show plots
show = True

# Save the picture ?
save = False

# Fancy plots ?
fancy = True

# Take into account the stray light?
straylight = False

# Minimum observable time for plots
threshold_obs_time = 50

# Time to acquire a target
t_acquisition = 6

# Catalogue name (in resources/)
catalogue = 'cheops_target_list_v0.1.dat'

# Maximum magnitude that can be seen by CHEOPS, only for cosmetics purposes
CHEOPS_mag_max = 12.5

# File name for the list of orbit file
orbits_file = 'orbits.dat'

# Factor in the SL post treatment correction ?
SL_post_treat = True

# Factor in mirror efficiency for the equivalent star magnitude ?
mirror_correction = True
#####################################################################################################################
# CONSTANTS AND PHYSICAL PARAMETERS
period = altitude2period(apogee,perigee)

###########################################################################
### INITIALISATION
file_flux = 'flux_'

# changes the threshold by addition the acquisition time:
threshold_obs_time += t_acquisition

# Formatted folders definitions
folder_flux, folder_figures, folder_misc = init_folders(orbit_id)

## Prepare grid
n_alpha = param.resx
n_delta = param.resy

ra_i = 0
ra_f = 2.*np.pi

dec_i = -np.pi/2.
dec_f = np.pi/2.

ra_step = (ra_f-ra_i)/n_alpha
dec_step = (dec_f-dec_i)/n_delta

iterable = (ra_i + ra_step/2+ i*ra_step for i in range(n_alpha))
ras = np.fromiter(iterable, np.float)

iterable = (dec_i + dec_step/2+ i*dec_step for i in range(n_delta))
decs = np.fromiter(iterable, np.float)

ra_grid, dec_grid = np.meshgrid(ras, decs)

if SAA:
	SAA_data = np.loadtxt('resources/SAA_table_%d.dat' % orbit_id, delimiter=',')
	SAA_data = SAA_data[SAA_data[:,0]>= minute_ini]
	SAA_data = SAA_data[SAA_data[:,0]<= minute_end]

computed_orbits = np.loadtxt(folder_misc+orbits_file)[:,0]
############################################################################
### Load catalogue and assign them to the nearest grid point

name_cat, ra_cat, dec_cat, mag_cat = load_catalogue(catalogue)



index_ra_cat = np.zeros(np.shape(ra_cat))
index_dec_cat= np.zeros(np.shape(ra_cat))

ii = 0
for name in name_cat:
	if name == OI: 
		break
	ii += 1

print 'Target is >>>',  name_cat[ii]

name_cat= name_cat[ii]
ra=ra_cat[ii]
dec=dec_cat[ii]
mag=mag_cat[ii]
id_ra = find_nearest(ras, ra/const.RAD)
id_dec = find_nearest(decs, dec/const.RAD)

obj = target_list(name, ra/const.RAD, id_ra, dec/const.RAD, id_dec, mag, int(period+3))

# Apply the flux correction (SL post-treatment removal and the mirror efficiency)
corr_fact = 1.0
if mirror_correction: corr_fact /= param.mirror_efficiency
if SL_post_treat: corr_fact *= (1.0 - param.SL_post_treat_reduction)
############################################################################
### Start the anaylsis

start = time.time()

# Prepare the arrays
visibility = np.zeros(np.shape(ra_grid))
#observations = np.zeros(len(name_cat)*)
workspace = np.zeros(np.shape(ra_grid))
#data = np.zeros(np.shape(ra_grid))

# Load the reference times
orbits = np.loadtxt(folder_misc+orbits_file,dtype='i4')
minutes_orbit_iditude = np.loadtxt('resources/minute_table_%d.dat' % orbit_id, delimiter=',',dtype='Int32')

# Set variables for printing the advance
numberofminutes = minute_end+1 - minute_ini
lo = fast_minute2orbit(minutes_orbit_iditude,minute_end, orbit_id)
fo = fast_minute2orbit(minutes_orbit_iditude,minute_ini, orbit_id)
lp = -1

junk, junk, at_ini, junk = fast_orbit2times(minutes_orbit_iditude, fo, orbit_id)
first_computed = computed_orbits[computed_orbits<=fo][-1]

first_minute = minute_ini
last_minute = minute_end

if not fo == first_computed: 
	junk, junk, minute_ini, junk = fast_orbit2times(minutes_orbit_iditude, first_computed, orbit_id)
#	print '1st referenced orbit: %d\twanted orbit: %d' % (first_computed, fo)

try:
	for minute in range(minute_ini,int(minute_end)+1+int(period)):
		minute = int(minute)
	
		if SAA and fast_SAA(SAA_data, minute): SAA_at_minute = True
		else: SAA_at_minute = False

		orbit_current = fast_minute2orbit(minutes_orbit_iditude, minute, orbit_id)
		if orbit_current > lp: 
			lp = orbit_current
			message = "Analysing orbit %d on %d...\t" % (lp,lo)

			sys.stdout.write( '\r'*len(message) )
			sys.stdout.write(message)
			sys.stdout.flush()

		junk, len_orbit, atc_ini, junk = fast_orbit2times(minutes_orbit_iditude, orbit_current, orbit_id)
		try:
			ra, dec, S_sl = load_flux_file(minute, file_flux, folder=folder_flux)
			load = True
			minute_to_load = minute-atc_ini#+shift
		except IOError:
		# if there is nothing then well, do nothing ie we copy the past values
		# in which orbit are we ?
		
		# get the previous orbit computed and copy the stray light data of this orbit :
			#orbit_previous = orbits[orbits[:,0] < orbit_current][-1,0]
			#minute_replacement = minute - atc_ini + shift #+ at_ini
			minute_to_load = minute-atc_ini

			if SAA_at_minute:
				obj.current_visibility = 0
			else:
				obj.current_visibility = obj.visible_save[minute_to_load]
			load = False
			# populate the visbility matrix
#			for ii in range(0, targets[0].CountObjects()):

		if load:
			ra_ = obj.ra
			dec_ = obj.dec
			a = np.where(np.abs(ra_-ra)<ra_step/2)[0]			
			b = np.where(np.abs(dec_-dec)<dec_step/2)[0]	
			INT = np.intersect1d(a,b)
			if np.shape(INT)[0] == 0 or (straylight and S_sl[INT]*corr_fact > obj.maximum_flux()): 
				obj.visible_save[minute_to_load] = 0
				obj.current_visibility = 0
				continue
			else:				
				obj.visible_save[minute_to_load] = 1

				if SAA_at_minute: obj.current_visibility = 0
				else: obj.current_visibility = 1

		if minute == minute_ini:
			obj.workspace=obj.current_visibility
			continue

		obj.Next(minute,threshold_obs_time)

except KeyboardInterrupt: print hilite('\nWARNING! USER STOPPED LOADING AT MINUTE %d' % minute,False,False)

obj.Next(minute,threshold_obs_time)

print
############################################################################

end = time.time()
elapsed_time = round((end-start)/60.,2)
sys.stdout.write( '\r'*len(message) )
sys.stdout.flush()
print "Time needed: %2.2f min" % elapsed_time

### Plot a few things
if fancy: figures.set_fancy()

### Plot time line
figures.set_fancy()

minute_ini = first_minute
minute_end = last_minute

fig = plt.figure()
ax = plt.subplot(111)

ii = 0
#ax.yaxis.set_major_locator(MultipleLocator(1))
plt.grid(True)

visi = obj.Visibility()
invi = obj.Invisibility()

dist = 0
##for v, i in zip(visi, invi):
##	print v, i, i-v, v-dist
##	dist = i


timestamps = np.zeros(lo+1-fo)
obs_time = np.zeros(lo+1-fo)
for orbit in range(fo, lo+1):
	ii = orbit-fo

	junk, junk, a, e = fast_orbit2times(minutes_orbit_iditude, orbit, orbit_id) 
	timestamps[ii] = a

	visi_c = visi[(visi <= e) & (visi >= a)]
	next_inv = invi[(visi <= e) & (visi >= a)]
	invi_c = invi[(invi <= e) & (invi >= a)]

	if np.shape(visi_c)[0] == 2:
		print np.shape(visi_c)[0]
		exit()
	if np.shape(next_inv)[0] == 2:
		print np.shape(visi_c)[0]
		exit()
	if np.shape(visi_c)[0] > 0 and next_inv[0] > e:
		obs_time[ii] += e - visi_c + 1
	elif np.shape(visi_c)[0] > 0:
		print orbit
		obs_time[ii] += next_inv - visi_c
	

#2@	current_in = invi[(invi >= a) & (invi <= e)]
#2@	current_vi = visi[(visi >= a) & (visi <= e)]
	#2@shape_in = np.shape(current_in)[0]
	#2@shape_vi = np.shape(current_vi)[0]

	#2@if shape_in == 2 :
	#2@	obs_time[ii] += current_in[0]-a
	#2@	np.delete(current_in, 0)
	#2@	shape_in = np.shape(current_in)[0]

	#2@if shape_in == 1 and shape_vi == 1:
	#2@	obs_time[ii] += current_in[0] - current_vi[0]
	#2@elif shape_in == 1 and shape_vi == 0:
	#2@	obs_time[ii] += current_in[0] - a
	#2@elif shape_in == 0 and shape_vi == 1:
	#2@	obs_time[ii] += e - current_vi[0]

	if obs_time[ii] < 0:
		print a,e
		print current_in
		print current_vi
		exit()

#print timestamps
#print obs_time

plt.plot (timestamps, obs_time, lw=2)
plt.ylabel('Available Obs. Time per Orbit [min]')

# convert epoch to matplotlib float format
labels = timestamps * 60. + const.timestamp_2018_01_01
labels = np.linspace(minute_ini, minute_end+1, 12) * 60. + const.timestamp_2018_01_01
plt.xlim([minute_ini, minute_end+1])

#plt.xlim([minute_ini, minute_end+1])
#ax.xaxis.set_major_locator(MultipleLocator((minute_end-minute_ini+1)/11))

# to human readable date
pre = map(time.gmtime, labels)
labels = map(figures.format_second, pre)

ax.set_xticklabels(labels)
fig.autofmt_xdate()

if save:
	threshold_obs_time -= t_acquisition
	if SAA: note = '_SAA'
	else: note = ''
	fname = '%svisibility_%s_obs_%d_o_%d_to_%d%s' % (folder_figures, OI, threshold_obs_time,fo,lo, note)
	figures.savefig(fname,fig,fancy)

if show: plt.show()
