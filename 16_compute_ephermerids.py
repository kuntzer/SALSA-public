''' 16-compute-ephemerids.py
=========================
AIM:	Computes the ephermerids for all the cells that constitute the sky grids.
	The cell must be visible for at least (period)-(max_interruptions)+t_acquisition
	To be used by the three next scripts (17, 18, 19) to treat and plot.

INPUT:	files: 	- <orbit_id>_misc/orbits.dat
		- <orbit_id>_flux/flux_*.dat
		- resources/saa_*.dat
	variables: see section PARAMETERS (below)

OUTPUT:	- <orbit_id>_misc/ephemerids_inter_<max_interruptions>_mag_<mag_max><_SAA?>.npz

CMD:	python 16-compute-ephemerids.py

ISSUES:	<none known>

REQUIRES:- standard python libraries, specific libraries in resources/ (+ SciPy)
	 - Structure of the root folder:
	   * <orbit_id>_flux/ --> flux files
	   * <orbit_id>_figures/maps/ --> figures
	   * <orbit_id>_misc/ --> storages of data

REMARKS: Not with real catalogue.
'''

###########################################################################
### INCLUDES
import numpy as np
import os
import matplotlib.cm as cm
import time

from resources.routines import *
from resources.TimeStepping import *

import parameters as param
import resources.figures as figures
from resources.targets import *

from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter
###########################################################################
### PARAMETERS
# orbit_iditude of the orbit in km
alt = 700
orbit_id = '6am_%d_5_conf4e' % alt
apogee=alt
perigee=alt

# First minute in data set !
minute_ini = 0

# Last minute to look for
minute_end = 1440 * 365

# File name for the list of orbit file
orbits_file = 'orbits.dat'

# Maximum interruption time tolerated [min] (acquisition time not included)
max_interruptions = 0 # see below period = 

# Time to acquire a target [min]
t_acquisition = 0

# Take into account the stray light?
straylight = True

# Maximum visible magnitude
mag_max = 9.
conversion_V_class = {9: 'G', 12: 'K'}

# Include SAA ?
SAA = True

# This is a way to vary the results by multiplying the whole pst by a number.
# This is very easy as if the pst is multiplied by a constant, it can be taken out of the
# integral and only multplying the flux is equivalent to re-running all the simulations
pst_factor = 1.

# Factor in the SL post treatment correction ?
SL_post_treat = True

#####################################################################################################################
# CONSTANTS AND PHYSICAL PARAMETERS
period = altitude2period(apogee,perigee)
max_interruptions = period - 1

###########################################################################
### INITIALISATION
flux_threshold = param.ppm_thresholds[mag_max] * 1e-6

file_flux = 'flux_'

# changes the threshold by addition the acquisition time:
threshold_obs_time = period - max_interruptions + t_acquisition

print 'ORBIT ID:\t%s\nmax_interruptions:%d+%d min\nMAGNITIUDE:\t%02.1f\nPST factor\t%g\nSAA\t\t%s' % (orbit_id,max_interruptions,t_acquisition, mag_max,pst_factor,SAA)

# Formatted folders definitions
folder_flux, folder_figures, folder_misc = init_folders(orbit_id)

if not os.path.isdir(folder_figures):
	print '\tError: figure folder %s does not exists.' % (folder_figures)
	exit()

sys.stdout.write("Loading list of computed orbits...\t")
sys.stdout.flush()

orbits = np.loadtxt(folder_misc+orbits_file,dtype='i4')

#list_minutes = -1. * np.ones( ( np.shape(orbits)[0] + 2 ) * period )
list_minutes=[]

id_min = 0
times = np.loadtxt('resources/minute_table_%s.dat' % orbit_id, delimiter=',',dtype='Int32')
for ii, orbit_current in enumerate(orbits[:,0]):
	t_ini, t_end, a_ini, a_end = fast_orbit2times(times,orbit_current,orbit_id)
	for minute in range(a_ini, a_end+1):
		list_minutes.append(int(minute))
		id_min += 1
list_minutes=np.asarray(list_minutes)

list_minutes = list_minutes[list_minutes > -1]

# apply conditions
list_minutes = list_minutes[list_minutes >= minute_ini]
list_minutes = list_minutes[list_minutes <= minute_end]
minute_end = int(list_minutes[-1])
print 'Done.'

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

visibility = np.zeros(np.shape(ra_grid))
visibility_save = np.zeros([np.shape(ra_grid)[0], np.shape(ra_grid)[1], int(period+2)])

workspace = np.zeros(np.shape(ra_grid))
data = np.zeros(np.shape(ra_grid))

numberofminutes = minute_end+1 - minute_ini

if SAA:
	SAA_data = np.loadtxt('resources/SAA_table_%s.dat' % orbit_id, delimiter=',')
	SAA_data = SAA_data[SAA_data[:,0]>= minute_ini]
	SAA_data = SAA_data[SAA_data[:,0]<= minute_end]

computed_orbits = np.loadtxt(folder_misc+orbits_file)[:,0]

stellar_type = conversion_V_class[mag_max]
stellar_flux = param.stellar_fluxes[stellar_type][mag_max]

aperture_aera_in_px = np.pi * param.aperture_size * param.aperture_size


############################################################################
### Load catalogue and assign them to the nearest grid point

message = 'Preparing the target list...\t\t'
sys.stdout.write(message)
sys.stdout.flush()

targets = []

for ra, dec in zip(np.ravel(ra_grid), np.ravel(dec_grid)):
	id_ra = find_nearest(ras, ra)
	id_dec = find_nearest(decs, dec)

	targets.append(target_list('%3.1f/%2.1f' % (ra,dec), ra, id_ra, dec, id_dec, mag_max, int(period+3), flux=stellar_flux))

message = 'Done, %d targets prepared.\n' % len(targets)
sys.stdout.write(message)
sys.stdout.flush()

# Apply the flux correction (SL post-treatment removal and the mirror efficiency)
corr_fact = 1.
if SL_post_treat: corr_fact *= (1. - param.SL_post_treat_reduction)

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
minutes_altitude = np.loadtxt('resources/minute_table_%s.dat' % orbit_id, delimiter=',',dtype='Int32')

# Set variables for printing the status
numberofminutes = minute_end+1 - minute_ini
lo = fast_minute2orbit(minutes_altitude,minute_end, orbit_id)
fo = fast_minute2orbit(minutes_altitude,minute_ini, orbit_id)
lp = -1

_, _, at_ini, _ = fast_orbit2times(minutes_altitude, fo, orbit_id)
first_computed = computed_orbits[computed_orbits<=fo][-1]

first_minute = minute_ini
last_minute = minute_end

if not fo == first_computed: 
	_, _, minute_ini, _ = fast_orbit2times(minutes_altitude, first_computed, orbit_id)
#	print '1st referenced orbit: %d\twanted orbit: %d' % (first_computed, fo)
try:
	for minute in range(minute_ini,minute_end+1):#+int(period)):
		minute = int(minute)
	
		if SAA and fast_SAA(SAA_data, minute): SAA_at_minute = True
		else: SAA_at_minute = False

		orbit_current = fast_minute2orbit(minutes_altitude, minute, orbit_id)
		if orbit_current > lp: 
			lp = orbit_current
			message = "Analysing orbit %d on %d...\t" % (lp,lo)

			sys.stdout.write( '\r'*len(message) )
			sys.stdout.write(message)
			sys.stdout.flush()

		_, len_orbit, atc_ini, atc_end = fast_orbit2times(minutes_altitude, orbit_current, orbit_id)
		try:
			ra, dec, S_sl = load_flux_file(minute, file_flux, folder=folder_flux)
			S_sl *= pst_factor
			load = True
			minute_to_load = minute-atc_ini#+shift
		except IOError:
		# if there is nothing then well, do nothing ie we copy the past values
		# in which orbit are we ?
		
		# get the previous orbit computed and copy the stray light data of this orbit :
			#orbit_previous = orbits[orbits[:,0] < orbit_current][-1,0]
			#minute_replacement = minute - atc_ini + shift #+ at_ini
			minute_to_load = minute-atc_ini

			
			for obj in targets:
				try:
					obj.current_visibility = obj.visible_save[minute_to_load]
				except IndexError:
					print minute_to_load, minute_end, atc_end, minute
					raise IndexError()

				if SAA_at_minute and obj.current_visibility==1:
					obj.current_SAA_interruption+=1


			load = False
			# populate the visbility matrix
#			for ii in range(0, targets[0].CountObjects()):

		if load:
			for obj in targets:
				ra_ = obj.ra
				dec_ = obj.dec
				a = np.where(np.abs(ra_-ra)<ra_step/10)[0]
				b = np.where(np.abs(dec_-dec)<dec_step/10)[0]
				INT = np.intersect1d(a,b)
				assert np.size(INT)<2

				F_sl_for_obj = S_sl[INT] * corr_fact * param.SL_QE * aperture_aera_in_px
				F_star = obj.get_flux() * flux_threshold # global throughput already included in the stellar flux!
				
				
				if np.shape(INT)[0] == 0 or (straylight and F_sl_for_obj > F_star): 
					obj.visible_save[minute_to_load] = 0
					obj.current_visibility = 0
					continue
				else:		
					#print S_sl[INT] * corr_fact * param.SL_QE * aperture_aera_in_px, S_sl[INT], corr_fact, param.SL_QE, aperture_aera_in_px, F_sl_for_obj, F_star		
					obj.visible_save[minute_to_load] = 1

				if SAA_at_minute: 
					obj.current_SAA_interruption+=1
					obj.current_visibility = 1
				else: obj.current_visibility = 1

		if minute == minute_ini:
			for obj in targets:
				obj.obs_time=obj.current_visibility
			continue

		for obj in targets: obj.Next(minute,threshold_obs_time)

except KeyboardInterrupt: 
	print hilite('\nWARNING! USER STOPPED LOADING AT MINUTE %d' % minute,False,False)
	raise KeyboardInterrupt()

for ii in range(0, targets[0].CountObjects()): targets[ii].Next(minute,threshold_obs_time)


print
worthy_targets = []
for obj in targets: obj.PrepareSave()

for ii in range(0, targets[0].CountObjects()):
	if np.shape(targets[ii].visible)[0] > 0:
		worthy_targets.append(targets[ii])


############################################################################

end = time.time()
elapsed_time = round((end-start)/60.,2)
sys.stdout.write( '\r'*len(message) )
sys.stdout.flush()
print "Time needed: %2.2f min" % elapsed_time

threshold_obs_time -= t_acquisition
if SAA: note = '_SAA'
else: note = ''

if not pst_factor == 1.: note += '_%1.1fpst' % pst_factor
if SL_post_treat: note+= '_%4.3fSLreduction' % param.SL_post_treat_reduction
fname = 'ephemerids_inter_%d_mag_%3.1f%s' % (max_interruptions,mag_max,note)#,threshold_obs_time,fo,lo, note)
print 'Filed saved as %s' % fname
np.savez_compressed(folder_misc+fname, worthy_targets=worthy_targets)


