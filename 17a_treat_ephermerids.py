''' 17-treat-ephemerids.py
=========================
AIM:	Using the ephemerids computed by 16-compute-ephemerids.py and observational constraints
	(period of the planet, transit time) calculates observations period.
	To be used by the two next scripts (18, 19) to treat and plot.

INPUT:	files: 	- <orbit_id>_misc/ephemerids_inter_<max_interruptions>_mag_<mag_max><_SAA?>.npz
	variables: see section PARAMETERS (below)

OUTPUT:	<orbit_id>_<SL_angle>misc/ephemerids_obs<transit_duration>h_<max_interruptions>inter_V<mag_max><_SAA?>.npz

CMD:	python 17-treat-ephemerids.py

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
# Orbit id
orbit_id = '800_35_AKTAR'
apogee=800
perigee=800

# File name for the list of orbit file
orbits_file = 'orbits.dat'

# Minimum observable time for plots [h] (Only used for consecutive observation time)
transit_duration = 10

# Maximum interruption time tolerated [min]
max_interruptions = 99

# Maximum visible magnitude
mag_max = 12.

# Take SAA into account?
SAA = False

# Print much information ?
verbose = False

# If set to True, then it will be observations of at least (period - max_interruptions)
# If set to False, then it is minimum (period - max_interruptions) minutes per orbit, 
# not necesseraly consecutive.
consecutive = False

# Factor in the SL post treatment correction ?
SL_post_treat = True

# Stop before saving results to file.
early_stop = False

# Minimal # of days of obs (if consecutive == False), must be a list
nb_obs_days = range(0,10,5)#[13]#range(20,45,5)#[13]#range(5,45,5)#[0,10,20,30,40]#range(10,17,1)##range(10,110,10)#

# Minimal minutes to be observed per orbit (if consecutive == False)
mins_t_obs_per_orbit = [80]#[78]#range(68,78,1)

# This is a way to vary the results by multiplying the whole pst by a number.
# This is very easy as if the pst is multiplied by a constant, it can be taken out of the
# integral and only multplying the flux is equivalent to re-running all the simulations
pst_factor=0.

# File name for the input file (in a compressed binary Python format)
if SAA: note = '_SAA'
else: note = ''
if not pst_factor == 1.: note += '_%1.1fpst' % pst_factor
if SL_post_treat: note+= '_%4.3fSLreduction' % param.SL_post_treat_reduction
input_fname = 'ephemerids_inter_%d_mag_%3.1f%s.npz' % (max_interruptions,mag_max,note)


if not consecutive: note += '_cumul_'

for min_t_obs_per_orbit in mins_t_obs_per_orbit:
	print '*'*30, 'min_t_obs_per_orbit %1.1f' % min_t_obs_per_orbit

	skycoverage_fname = 'skycoverage_%dmin_V%3.1f%s.txt' % (min_t_obs_per_orbit,mag_max,note)
	for nb_obs_day in nb_obs_days:
		# File name for the input file (in a compressed binary Python format)
		if consecutive:
			output_fname = 'ephemerids_obs%dh_%dinter_V%3.1f%s.npz' % (transit_duration,max_interruptions,mag_max,note)
		else: 
			output_fname = 'ephemerids_%ddays_%dmin_V%3.1f%s.npz' % (nb_obs_day,min_t_obs_per_orbit,mag_max,note)
		#####################################################################################################################
		# CONSTANTS AND PHYSICAL PARAMETERS
		period = altitude2period(apogee, perigee)
		###########################################################################
		### INITIALISATION
		# Formatted folders definitions
		folder_flux, folder_figures, folder_misc = init_folders(orbit_id)

		sky_coverage=0.
	
		output=open(os.path.join(folder_misc,skycoverage_fname),"a") 
	
		print 'ORBIT ID:\t\t%s\nPST factor:\t\t%d\nMin Days of Coverage:\t%d\nmin_t_obs_per_orbit\t%d (%.1f%%)\nMAGNITIUDE:\t\t%02.1f\nSAA :\t%g' % (orbit_id,pst_factor,nb_obs_day,min_t_obs_per_orbit,min_t_obs_per_orbit/period*100., mag_max, SAA)
	
		# loading data
		sys.stdout.write("Loading worthy targets from %s ...\t" % input_fname)
		sys.stdout.flush()
	
		worthy_targets = np.load(folder_misc+input_fname)
		worthy_targets = worthy_targets['worthy_targets']

		check=np.zeros_like(worthy_targets)
	
		max_len = 0
		for k in range(0, len(worthy_targets)):
			if max_len < np.shape(worthy_targets[k].Visibility())[0]:
				max_len = np.shape(worthy_targets[k].Visibility())[0] 
	
		# too optimistic
		max_len = int(max_len)
	
		start_obs = np.empty([len(worthy_targets),max_len])
		stop_obs = np.empty([len(worthy_targets),max_len])
		interruptions_obs = np.empty([len(worthy_targets),max_len])
	
		print 'Done\n.%d targets loaded' % len(worthy_targets)
		
		###########################################################################
		### COMPUTATIONS
		###########################################################################
	
		###########################################################################
		# consecutive
		if consecutive:
		# Loop on all worthy targets
			for ii in range (0, len(worthy_targets)):
				y = float(ii)
				visi = worthy_targets[ii].Visibility()
				invi = worthy_targets[ii].Invisibility()
	
			# for every region in the sky/worthy target:
			# >> Find when you can look with transit_duration [h] with maximal max_interruptions [min]
			# >>>> return start and end time of observations with duration of interruptions [min]
	
				# Initialise all variables
				k = 0
				j = 0
				total_interruptions = 0
				start_observation_time = 0
				count_observation_time = 0
		
				do_observe = False
				has_observed=False
	
			# iterate on the visibility (i.e. time when the target becomes visible)	
				for k in range(0, len(visi)):
					# shorthand notations
					vis = visi[k]
					ini = invi[k]
	
					# Try to compute the interruption time with the next observability window
					try: 
						time_to_next_vis = visi[k+1]-ini
						stop_to_observe = False
					except IndexError: 
						stop_to_observe = True
		
					# if the time to next visiblity is larger than the max interruption time or no next window --> can't observe anymore
					if stop_to_observe or max_interruptions < time_to_next_vis: 
						# if the observation time is larger than the transit duration, then it can be observed.
						if do_observe and count_observation_time >= transit_duration*60. :
						# if you have been observing for longer than transit_duration [h], then remember when and remember the interruptions
							start_obs[ii,j] = start_observation_time
							stop_obs[ii,j] = ini
							interruptions_obs[ii,j] = total_interruptions
							j+=1
							has_observed = True
			
						do_observe = False
						total_interruptions = 0		
						start_observation_time = 0
						count_observation_time = 0
			
						k+=1
						if stop_to_observe: break
						else: continue
		
					# if the time to next visiblity is smaller than the max interruption time --> save the interruption time
					else: total_interruptions += time_to_next_vis
		
					# if you were not observing, you can now.
					if max_interruptions > time_to_next_vis and not do_observe: 
						do_observe=True
						start_observation_time = vis
		
					# count the time you can observe
					count_observation_time += ini-vis + time_to_next_vis
					
					k+=1
					if stop_to_observe: break
	
				# Debugging infos
				has_observed = False
				if has_observed and verbose: print start_obs[ii,0], stop_obs[ii,0], interruptions_obs[ii,0]
	
	###########################################################################
	# non-consecutive
		count = 0
		count_day=0
		nogo_count=0
		if not consecutive:
			times = np.loadtxt('resources/minute_table_'+str(orbit_id)+'.dat', delimiter=',',dtype='Int32')
			a_end_orbits = np.linspace(1,param.last_orbits[orbit_id],param.last_orbits[orbit_id])
			tmps = a_end_orbits
			a_end_orbits = fast_orbit2a_end_vect(times, a_end_orbits, orbit_id)
		
			obs_tot = np.zeros(len(worthy_targets))
			obs_orbit_tot = np.zeros(len(worthy_targets))
			min_obs_time = nb_obs_day*24.*60.*min_t_obs_per_orbit/period
		
			nb_orbit_per_day = 24.*60. / period
	
			for ii in range (0, len(worthy_targets)):
				y = float(ii)
				message = '\r%3.1f %%' % (y/float(len(worthy_targets))*100.)
				sys.stdout.write(message)
				sys.stdout.flush()
	
				visi = worthy_targets[ii].Visibility()
				invi = worthy_targets[ii].Invisibility()
	
	
				# Initialise all variables
				k = 0
				t_obs_per_orbit = np.zeros(param.last_orbits[orbit_id]+1)
	
				observations = invi-visi
				observations=observations[observations>=min_t_obs_per_orbit]
				if np.size(observations)>0:
					check[ii] += observations.sum()
				else: continue

			sky_coverage=0.
			#print np.size(check[check>nb_obs_day*24.*60.])
			for ii in range(len(worthy_targets)):
				if check[ii]<nb_obs_day*24.*60.: continue
				rat, dect = worthy_targets[ii].Coordinates()
				sky_coverage+=0.5/param.resx/param.resy*np.pi*np.cos(dect)
		
			message = '\rComputations done.' 
			sys.stdout.write(message)
			sys.stdout.flush()
	
			print '\nSky coverage for %d days' % nb_obs_day

			rslts = sky_coverage
	
			print nb_obs_day,'\t***', round(rslts*100.,3), ' % ***'
			print >> output, nb_obs_day,'\t', round(rslts*100.,3)
	
			if early_stop: continue
	
			if verbose: print obs_tot/24./60.
			np.savez_compressed(folder_misc+output_fname, worthy_targets=worthy_targets, obs_tot=check)
			print 'Filed saved as %s' % output_fname
			output.close()

