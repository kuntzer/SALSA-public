#!/bin/env python
######################################################################
''' COMPUTE_MASTER.py
=========================
AIM:	Cycle through the observability maps and execute stray light computations
	Temporal resolution is 60 seconds and orbit step size is given by the max. error

INPUT:	files: 	- observability maps in form of minutes,ra,dec [rad] in a file orbit_<#>.dat in raw_maps_<orbit_id>/
		- in resources/ : minute_tables__<orbit_id>, moon__<orbit_id>, sun_<orbit_id>
		- The complete package of the stray_light.f code (including the orbit file)
	variables: see section PARAMETERS (see below)

OUTPUT:	in height_part/ (see PARAMETERS): one file per minute with ra,dec,flux [ra,ra, ph/(s.px)]

CMD:	python (name_file).py	

ISSUES:	<NONE KNOWN>

REQUIRES: standard python libraries, specific libraries in resources/

REMARKS: <none>
'''
######################################################################
# DEFINITIONS AND INCLUDES
import numpy as np
import subprocess
import os
import sys
import time

import resources.constants as const
from resources.routines import *
from resources.TimeStepping import *
######################################################################
# PARAMETERS
# km (only required to compute orbit's period)
apogee = 800
perigee= 800

# orbit id 
orbit_id = '620_35_AKTAR'
 
# 1, 2 or 3 (or, ...) into which folder ? Scheme is orbit_id_part/
part = 1

# Orbit to start by:
orbit_ini = 1

# Last orbit to compute (see # of orbits in one year in parameters.py)
orbit_end = 16

# Defintions of the steps
orbit_step = 10

# If after a given computation, the precision is not good enough, the step of the orbit is adapted.
adaptative_timestep = True

# Minimum orbital step
min_step = 5

# must be in straylight_xxx/ORBIT
file_orbit = 'orbit_%s.dat' % orbit_id
file_sun = 'sun_%s.dat' % orbit_id
folder = 'absolute_map_%s' % orbit_id

# Standard values are:
# folder_flux = '%d_%d' % (orbit_id,part)
# file_flux = 'flux_'
folder_flux = '%s_%d' % (orbit_id,part)
file_flux = 'flux_'

# Tolerance of the maximal magnitude difference (5% = 0.05)
p = 0.1
# Recalculate previously computed minutes ? (Recommendation: False -- time consuming)
override = False
# Save the date to a file for straylight.f (used only if the log_non_nominal is set to .true.)
save_date_JD = True
# Recompiles the Fortran code to ensure that the correct parameters and constants are loaded.
# Make sure "compile" file CHMOD is set to 777. (Recommendation: True -- very fast)
recompile = True

monitor_angle_usage = False
######################################################################
# INITIALISATION
# path to stray light folder (containting CODE, OUTPUT, INPUT, ORBIT)
path = 'straylight_%s_%d' % (orbit_id,part)

# Says hi
print '\nSALSA v%s' % const.version
print 'Computing max. every %d orbits for orbit ID %s' % (orbit_step,orbit_id)
print '------------------------------------------------------'

start = time.time()

# Loads the heavy tables (orbit, sun)
sys.stdout.write("Loading trajectory file...\t\t")
sys.stdout.flush()
 
try:
	orbit = np.loadtxt(path+'/ORBIT/'+file_orbit, delimiter='\t')
except ValueError:
	orbit = np.loadtxt(path+'/ORBIT/'+file_orbit, delimiter=' ')
print "Done."

sys.stdout.write("Loading Sun position file...\t\t")
sys.stdout.flush()

try:
	sun = np.loadtxt(path+'/ORBIT/'+file_sun, delimiter='\t')
except ValueError:
	sun = np.loadtxt(path+'/ORBIT/'+file_sun, delimiter=' ')
print "Done."

# initialise a few variables
total_targets = 0
map_obstot = np.empty(3)
is_first_run = True
orig_dir = os.getcwd()

# Computes the period of the orbit given the altitude
period = altitude2period(apogee,perigee)


# Recompiles the code
if recompile:
	sys.stdout.write("Compilation of the code...\t\t")
	sys.stdout.flush()
	os.chdir(os.path.join(os.path.abspath(sys.path[0]), '%s/CODE/' % path))
	os.system("./compile")
	os.chdir(orig_dir)
	print "Done."
######################################################################
# LOOPING OVER ALL GIVEN MAPS
orbit_current = orbit_ini
# Loop on the orbits
preceeding = orbit_current
former_step = orbit_step
while (orbit_current <= orbit_end):
	print '\n---------------- ORBIT %d --- ID %s ----------------------------' % (orbit_current, orbit_id)
	# Get the initial and final time
	start_minute = time.time()
	t_ini, t_end, a_ini, a_end = orbit2times(orbit_current,orbit_id)
	minute = a_ini

	# Load the observability map for the orbit
	try:
		map_obstot = load_map(orbit_current,folder)
	except IOError:
		print 'Error, could not load orbit %d' % orbit_current
		exit()

	while (minute <= a_end):
		if not override:
			try:
			# Try to load the fluxes for a given minute (minute goes from 0 to period whereas a_ini is the absolute time from 0:0:0.0 1/1/2018 in min
				ra, dec, S_sl = load_flux_file(minute, file_flux,folder=folder_flux)
				minute += 1
				continue
			# If there is an error while loading the file, give up and go to the next minute
			except IOError: pass

		sys.stdout.write( '\rComputing minute: %3d\tAbsolute: %6d\t' % ((minute-a_ini),minute) )
		sys.stdout.flush()

		if save_date_JD:
			# compute the julian date for a given minute and saves it to a file read by stray_light.f
			JD = minutes_2_JD(minute)
			f = open('%s/INPUT/date_JD.dat' % path,'w')
			f.write(str(JD)+'d0\n')
			f.close()


		# find the position of the satellite in the orbit and compute RA, DEC of the sat with respect to the Earth
		id_sat = find_nearest(orbit[:,0],minute)
		x = orbit[id_sat,1]
		y = orbit[id_sat,2]
		z = orbit[id_sat,3]

		r = R_3d(x,y,z)
		ra_sat = rev( right_ascension(x,y) )
		dec_sat= declination(x,y,z)


		np.savetxt('%s/INPUT/coord_sat.dat' % path,[x,y,z,ra_sat,dec_sat,r],delimiter='\n')

		# find the position of the sun in the orbit and compute RA, DEC of the sat with respect to the Earth
		id_sun = find_nearest(sun[:,0],minute)
		xs = sun[id_sun,1]
		ys = sun[id_sun,2]
		zs = sun[id_sun,3]

		rs = R_3d(xs,ys,zs)
		ra_sun = rev( right_ascension(xs,ys) )
		dec_sun= declination(xs,ys,zs)


		np.savetxt('%s/INPUT/coord_sun.dat' % path,[xs,ys,zs,ra_sun,dec_sun,rs],delimiter='\n')

		# select only none zero value with resilience to rounding errors
		# See resources/routines.py
		map_obs = slice_map(map_obstot, minute)

		# Count the number of points for that particular minute
		total_targets += np.shape(map_obs)[0]

		if np.shape(map_obs)[0] > 0 :
		# Execute the stray light code only if there is more than 0 target !
			sys.stdout.write(str(np.shape(map_obs)[0])+' points\t')
			sys.stdout.flush()

		# Save the targets and the number of line to two separate files to optimise Fortran read.
			np.savetxt('%s/INPUT/coord_targets.dat' % path, map_obs[:,1:3], delimiter=' ', fmt='%3.3f')
			f = open('%s/INPUT/number_of_targets.dat' % path, 'w')
			f.write(str(np.shape(map_obs)[0])+'\n')
			f.close()
			os.chdir(os.path.join(os.path.abspath(sys.path[0]), '%s/CODE/' % path))
			subprocess.call(["./stray_light"])
		        # Move the files to the right ouput folder
			os.chdir(orig_dir)
			subprocess.call(["mv", "%s/OUTPUT/straylight.out" % path,'%s/%s%d.dat' % (folder_flux,file_flux,minute)])
		else:
		      sys.stdout.write('No points. ')
		      sys.stdout.flush()


		minute += 1
	end_minute = time.time()
	elapsed_time = round( end_minute - start_minute , 1)

	# end of orbit
	message = '\r%3.1f minutes to compute %d minutes of the orbit no %d / %d       \n' % (elapsed_time/60, t_end+1, orbit_current, orbit_end)
	sys.stdout.write(message)
	sys.stdout.flush()

	if monitor_angle_usage:
		subprocess.call(["mv", "%s/OUTPUT/angle_usage.out" % path,'%s/angles_%d.dat' % (folder_flux,orbit_current)])

	# Computes the differences to the reference orbit
	if adaptative_timestep and orbit_current>orbit_ini:
		# See details of compare_two_orbits in resources/routines
		# tries with minute 0 of reference is same as current
		pp = compare_two_orbits(preceeding, orbit_current, orbit_id, p=p, file_flux=file_flux, folder=folder_flux)
		pp_old = pp
		# tries with minute 0 of reference is same as current+1
		if pp > p :
			pp = compare_two_orbits(preceeding, orbit_current, orbit_id, p=p, file_flux=file_flux, folder=folder_flux, shift = 1)
		# tries with minute 0 of reference is same as current+2
		if pp > p : 
			pp = compare_two_orbits(preceeding, orbit_current, orbit_id, p=p, file_flux=file_flux, folder=folder_flux, shift = 2)
		if pp > p :
			# Tried, still bad
			status = False
			pp = pp_old
		else : status = True
		print 'Precision former step:', hilite( str(np.round(pp*100,2)), status, False),'%'

		# if the test worked, continue with optimum step size or increase it
		if pp < p:
			preceeding = orbit_current
			if former_step == orbit_step:
				current_step = orbit_step
			else:
				current_step = int(former_step*2)
				if current_step > orbit_step : current_step = orbit_step
				else : print 'Next adaptative step of :', current_step, 'orbit(s)'
			orbit_current += current_step
			former_step = current_step
		# if the test failed, reduce the step size
		else :
			current_step = int(former_step/2)
			if current_step < min_step: current_step = min_step 
			if current_step == former_step and current_step ==min_step  :
				message = 'Could not reach precision (reference orbit=%s, current=%s)' % (preceeding, orbit_current)
				print hilite(message, False, True)
				preceeding = orbit_current
				orbit_current += min_step 
			else:
				orbit_current += current_step - former_step
				former_step = current_step
			print 'Next adaptative step of :', current_step, 'orbit(s)'

	# Ensure that we do not get stuck somehow.
	else: orbit_current += orbit_step

######################################################################
# OUTPUTS THE FINAL REMARK
end = time.time()
elapsed_time = round((end-start)/60.,1)
print 'Stray light calculation carried out in '+ str(elapsed_time) +' min for '+ str(total_targets)+ ' points. Have nice day!'
