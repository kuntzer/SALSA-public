''' 11-extract_long_periods.py
===============================================
AIM:	Prepare cumulative plots (THIS SCRIPT IS WITHOUT STRAY LIGHT)

INPUT:	files: 	- <orbit_id>_misc/orbits.dat
		- <orbit_id>_flux/flux_*.dat
	variables: see section PARAMETERS (below)

OUTPUT:	in <orbit_id>_misc/ : complicated name files depending on the case (handled by 12-<...>.py)

CMD:	python 11-extract_long_periods.py

ISSUES:	<none known>

REQUIRES:- standard python libraries, specific libraries in resources/ (+ SciPy)
	 - Structure of the root folder:
	   * <orbit_id>_flux/ --> flux files
	   * <orbit_id>_figures/ --> figures
	   * <orbit_id>_misc/ --> storages of data
	   * all_figures/ --> comparison figures

REMARKS: THIS SCRIPT IS WITHOUT STRAY LIGHT
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
import resources.constants as const
import resources.figures as figures

###########################################################################
### PARAMETERS
# altitude of the orbit in km
apogee=700
perigee=700
orbit_id = 301

# First minute in data set !
minute_ini = 0

# Last minute to look for
minute_end = 1440*365

# File name for the list of orbit file
orbits_file = 'orbits.dat'

# Minimum observable time for plots
threshold_obs_time = 50

# Time to acquire a target
t_acquisition = 3

# Show preview ?
show = True

# Include SAA ?
SAA = True

# File name for the output file
output_fname = 'TEST-data_%d' % (threshold_obs_time)
extension = '.dat'

# Factor in the SL post treatment correction ?
SL_post_treat = True
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

numberofminutes = minute_end+1 - minute_ini

minutes_orbit_iditude = np.loadtxt('resources/minute_table_%d.dat' % orbit_id, delimiter=',',dtype='Int32')

if SAA:
	SAA_data = np.loadtxt('resources/SAA_table_%d.dat' % orbit_id, delimiter=',')
	SAA_data = SAA_data[SAA_data[:,0]>= minute_ini]
	SAA_data = SAA_data[SAA_data[:,0]<= minute_end]
###########################################################################
### LOAD AND COMPUTE LARGEST OBSERVATION PERIOD
start = time.time()

lp = -1

try:
	for minute in range(minute_ini,minute_end+1):
		minute = int(minute)

		if SAA and fast_SAA(SAA_data, minute): SAA_at_minute = True
		else:  SAA_at_minute = False

		orbit_current = fast_minute2orbit(minutes_orbit_iditude, minute, orbit_id)
		junk, junk, atc_ini, junk = fast_orbit2times(minutes_orbit_iditude, orbit_current, orbit_id)
		if orbit_current > lp: 
			lp = orbit_current
			message = "Loading stray light data orbit %d on %d...\t" % (lp, minutes_orbit_iditude[-1,0])

			sys.stdout.write( '\r'*len(message) )
			sys.stdout.write(message)
			sys.stdout.flush()
		try:
			ra, dec, S_sl = load_flux_file(minute, file_flux, folder=folder_flux)
			# Apply the flux correction (SL post-treatment removal)
			if SL_post_treat: S_sl *= (1.0 - param.SL_post_treat_reduction)
			load = True
		except IOError:
		# if there is nothing then well, do nothing ie we copy the past values
		# in which orbit are we ?
		
		# get the previous orbit computed and copy the stray light data of this orbit :
#			orbit_previous = orbits[orbits[:,0] < orbit_current][-1,0]
#			junk, junk, at_ini, junk = fast_orbit2times(minutes_orbit_iditude, orbit_previous, orbit_id)
			
			minute_replacement = minute - atc_ini# + at_ini
			load = False
#			try: ra, dec, S_sl = fast_load_flux_file(minute_replacement, file_flux, folder=folder_flux)
#			except IOError: ra, dec, S_sl = fast_load_flux_file(minute_replacement-1, file_flux, folder=folder_flux)
#			pass

		# populate the visbility matrix
		if SAA_at_minute: visibility = np.zeros(np.shape(ra_grid))
		elif load:
			visibility_save[...,minute-atc_ini] = 0
			for index, ra_ in enumerate(ra):
				id_ra = find_nearest(ras,ra_)
				id_dec = find_nearest(decs,dec[index])
				visibility[id_dec,id_ra] = 1
				visibility_save[id_dec,id_ra,minute-atc_ini] = 1
		else: visibility = visibility_save[...,minute_replacement]

		if minute == minute_ini: workspace=visibility.copy()
		else : 
		# if there is an interruption then, reset the value in workspace
		# but before saves the value if it is larger than "threshold_obs_time" minutes
			data[ (workspace>threshold_obs_time-1) & (visibility < 1) ] += \
			workspace[(workspace>threshold_obs_time-1)&(visibility< 1)]

			workspace[visibility < 1] = 0

		# if the point existed already, then add one minute
			workspace[visibility > 0] += 1

		# reset visibility without taking a chance of wrong something
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

