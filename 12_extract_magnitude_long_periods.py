''' 12-extract_magnitude_long_periods.py
===============================================
AIM:	Prepare cumulative plots (THIS SCRIPT IS with STRAY LIGHT)

INPUT:	files: 	- <orbit_id>_misc/orbits.dat
		- <orbit_id>_flux/flux_*.dat
	variables: see section PARAMETERS (below)

OUTPUT:	in <orbit_id>_misc/ : complicated name files depending on the case (handled by 13-<...>.py)

CMD:	python 12-extract_magnitude_long_periods.py

ISSUES:	<none known>

REQUIRES:- standard python libraries, specific libraries in resources/ (+ SciPy)
	 - Structure of the root folder:
	   * <orbit_id>_flux/ --> flux files
	   * <orbit_id>_figures/ --> figures
	   * <orbit_id>_misc/ --> storages of data
	   * all_figures/ --> comparison figures

REMARKS: THIS SCRIPT IS with STRAY LIGHT
'''
###########################################################################
### INCLUDES
import numpy as np
import pylab as plt
import os
import time

from resources.routines import *
from resources.TimeStepping import *

import parameters as param
import resources.figures as figures

###########################################################################
### PARAMETERS
# orbit_id
orbit_id = 301
apogee=700
perigee=700

# First minute in data set
minute_ini = 0

# Last minute to look for
minute_end = 1440*365

# File name for the list of orbit file
orbits_file = 'orbits.dat'

# Minimum consecutive observable time for plots
threshold_obs_time = 78

# Take a minimum observation time per orbit and Minimum observable time per orbit (NON-CONSECUITIVE)
min_obs_per_orbit = True
threshold_obs_time_per_orbit = 78

# Time to acquire a target
t_acquisition = 3

# Maximum visible magnitude
mag_max = 12.

# File name for the output file 
output_fname = 'mag_over_year_%d_mag_%02.1f' % (threshold_obs_time, mag_max)
extension = '.dat'

# Show preview ?
show = True

# Include SAA ?
SAA = False

# File name that contains the orbit used and the percentage of the sky unviewable because of SL
namefile = 'cumultative_SL_forbidden_%d_mag_%02.1f.dat' % (threshold_obs_time, mag_max)

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

print 'ORBIT ID:\t\t%d\nTHRESHOLD OBS TIME:\t%d+%d min' % (orbit_id,threshold_obs_time,t_acquisition)
if min_obs_per_orbit:
	print 'obs time per orbit\t%d min' % threshold_obs_time_per_orbit
print 'MAGNITIUDE:\t\t%02.1f\nN/S max:\t\t%d ppm\nIncluding SAA?\t\t%g' % (mag_max,param.ppm_threshold,SAA)

# changes the threshold by addition the acquisition time:
threshold_obs_time += t_acquisition

# Formatted folders definitions
folder_flux, folder_figures, folder_misc = init_folders(orbit_id)

if not os.path.isdir(folder_figures):
	print '\tError: figure folder %s does not exists.' % (folder_figures)
	exit()

sys.stdout.write("Loading list of computed orbits...\t\t")
sys.stdout.flush()

orbits = np.loadtxt(folder_misc+orbits_file,dtype='i4')

list_minutes = -1. * np.ones( ( np.shape(orbits)[0] + 2 ) * period )

id_min = 0
times = np.loadtxt('resources/minute_table_%d.dat' % orbit_id, delimiter=',',dtype='Int32')
for ii, orbit_current in enumerate(orbits[:,0]):
	t_ini, t_end, a_ini, a_end = fast_orbit2times(times,orbit_current,orbit_id)
	for minute in range(a_ini, a_end+1):
		list_minutes[id_min] = int(minute)
		id_min += 1

list_minutes = list_minutes[list_minutes > -1]

# apply conditions
list_minutes = list_minutes[list_minutes >= minute_ini]
list_minutes = list_minutes[list_minutes <= minute_end]
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
if min_obs_per_orbit: data_orbit = np.zeros(np.shape(ra_grid))

numberofminutes = minute_end+1 - minute_ini

minutes_orbit_iditude = np.loadtxt('resources/minute_table_%d.dat' % orbit_id, delimiter=',',dtype='Int32')

maximum_sl_flux = mag2flux(mag_max)
if mirror_correction: maximum_sl_flux *= param.mirror_efficiency

if SAA:
	SAA_data = np.loadtxt('resources/SAA_table_%d.dat' % orbit_id, delimiter=',')
	SAA_data = SAA_data[SAA_data[:,0]>= minute_ini]
	SAA_data = SAA_data[SAA_data[:,0]<= minute_end]

if os.path.isfile(folder_misc+namefile):
	os.remove(folder_misc+namefile)

f = open(folder_misc+namefile,'w')
###########################################################################
### LOAD AND COMPUTE LARGEST OBSERVATION PERIOD
start = time.time()

lp = -1
previous_part = -1
sl_rejection_orbit = 0.
shutdown_time = 0.
previous_period_rel = 1.
orbit_previous = 0.
do_print = False
load = True

try:
	for minute in range(minute_ini,minute_end+1):
		minute = int(minute)

		if SAA and fast_SAA(SAA_data, minute): SAA_at_minute = True
		else:  SAA_at_minute = False

		orbit_current = fast_minute2orbit(minutes_orbit_iditude, minute, orbit_id)
		junk, period_rel, atc_ini, junk = fast_orbit2times(minutes_orbit_iditude, orbit_current, orbit_id)

		if orbit_current > lp: 
			lp = orbit_current
			message = "Loading stray light data orbit %d on %d...\t" % (lp, minutes_orbit_iditude[-1,0])

			sys.stdout.write( '\r'*len(message) )
			sys.stdout.write(message)
			sys.stdout.flush()

			if load: print >> f, orbit_previous, sl_rejection_orbit/float(previous_period_rel-shutdown_time+1)*100.

			shutdown_time = 0
			sl_rejection_orbit = 0.
			sl_rejection_orbit_save = 0.
			previous_period_rel = period_rel
			orbit_previous = orbit_current

			if min_obs_per_orbit:
				data[ (data_orbit>threshold_obs_time_per_orbit-1)] += \
				data_orbit[(data_orbit>threshold_obs_time_per_orbit-1)]

				data_orbit = np.zeros(np.shape(ra_grid))
		try:
			ra, dec, S_sl = load_flux_file(minute, file_flux, folder=folder_flux)
			load=True
			# Apply the flux correction (SL post-treatment removal)
			if SL_post_treat: S_sl *= (1.0 - param.SL_post_treat_reduction)
			nb_targets = np.size(S_sl)
		except IOError:
		# if there is nothing then well, do nothing ie we copy the past values
		# in which orbit are we ?
		
		# get the previous orbit computed and copy the stray light data of this orbit :
			load = False
			
			minute_replacement = minute - atc_ini# + at_ini
			
		# populate the visbility matrix
		if SAA_at_minute: 
			visibility = np.zeros(np.shape(ra_grid))
			shutdown_time += 1

		elif load:
			sl_rejection_minute = 0.
			visibility_save[...,minute-atc_ini] = 0
			for ra_, dec_, sl in zip(ra,dec,S_sl):
				if sl > maximum_sl_flux: 
					sl_rejection_minute += 1.
					continue
				id_ra = find_nearest(ras,ra_)
				id_dec = find_nearest(decs,dec_)
				visibility[id_dec,id_ra] = 1
				visibility_save[id_dec,id_ra,minute-atc_ini] = 1
			sl_rejection_orbit += sl_rejection_minute/nb_targets

		else: visibility = visibility_save[...,minute_replacement]


		if minute == minute_ini: workspace=visibility.copy()
		else : 
		# if there is an interruption then, reset the value in workspace
		# but before saves the value if it is larger than "threshold_obs_time" minutes
			if min_obs_per_orbit:
				data_orbit[ (workspace>threshold_obs_time-1) & (visibility < 1) ] += \
				workspace[(workspace>threshold_obs_time-1)&(visibility< 1)]
			else:
				data[ (workspace>threshold_obs_time-1) & (visibility < 1) ] += \
				workspace[(workspace>threshold_obs_time-1)&(visibility< 1)]
			

			workspace[visibility < 1] = 0

		# if the point existed already, then add one minute
			workspace[visibility > 0] += 1

		# reset visibility without taking a chance of a wrong something
		del visibility
		visibility = np.zeros(np.shape(ra_grid))

except KeyboardInterrupt: print hilite('\nWARNING! USER STOPPED LOADING AT MINUTE %d' % minute,False,False)


# Check that we did not left anything behind (in a try structure to avoid weird things...)
try:
	data[ (workspace>threshold_obs_time-1) ] += \
	workspace[(workspace>threshold_obs_time-1)&(visibility< 1)]
except ValueError: pass
del workspace

end = time.time()
elapsed_time = round((end-start)/60.,1)
sys.stdout.write( '\r'*len(message) )
sys.stdout.flush()
print
print "Loaded stray light data\tTime needed: %2.2f min" % elapsed_time

if SAA: note = '_SAA'
else: note = ''

np.savetxt(folder_misc+output_fname+note+extension,data)

print "Data saved in %s%s" % (folder_misc,output_fname+note+extension)

if not show : exit()

plt.figure()
ax = plt.subplot(111)

extent = (-np.pi,np.pi,-np.pi/2.,np.pi/2.)
CS = ax.contour((ra_grid-np.pi)*180. / np.pi,dec_grid*180. / np.pi,data,colors='k',extent=extent)
CS = ax.contourf((ra_grid-np.pi)*180. / np.pi,dec_grid*180. / np.pi,data,cmap=plt.cm.jet,extent=extent)

plt.xlim([-180, 180])
plt.ylim([-90, 90])
plt.colorbar(CS)

ax.grid(True)

ax.set_xlabel(r'$\alpha$')
ax.set_ylabel(r'$\delta$')
plt.title('PREVIEW OF THE DATA [MINUTES]')

plt.show()

f.close()

