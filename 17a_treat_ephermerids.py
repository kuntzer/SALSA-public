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
transit_duration = 13

# Maximum interruption time tolerated [min]
max_interruptions = 99

# Maximum visible magnitude
mag_max = 9.

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

# Stop before saving results to file.
early_stop = False

# Minimal # of days of obs (if consecutive == False), must be a list
nb_obs_days = range(10,110,10)#[5,10,13,20,30,40]#[0,10,20,30,40]#range(10,17,1)#range(5,45,5)#

# Minimal minutes to be observed per orbit (if consecutive == False)
mins_t_obs_per_orbit = [50]#[78]#range(68,78,1)

# This is a way to vary the results by multiplying the whole pst by a number.
# This is very easy as if the pst is multiplied by a constant, it can be taken out of the
# integral and only multplying the flux is equivalent to re-running all the simulations
pst_factor=1.

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
	
	
		output=open(os.path.join(folder_misc,skycoverage_fname),"a") 
	
		print 'ORBIT ID:\t\t%s\nPST factor:\t\t%d\nMin Days of Coverage:\t%d\nmin_t_obs_per_orbit\t%d (%.1f%%)\nMAGNITIUDE:\t\t%02.1f\nSAA :\t%g' % (orbit_id,pst_factor,nb_obs_day,min_t_obs_per_orbit,min_t_obs_per_orbit/period*100., mag_max, SAA)
	
		# loading data
		sys.stdout.write("Loading worthy targets from %s ...\t" % input_fname)
		sys.stdout.flush()
	
		worthy_targets = np.load(folder_misc+input_fname)
		worthy_targets = worthy_targets['worthy_targets']
	
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
	
			# iterate on the visibility (i.e. time when the target becomes visible)	
				for k in range(0, len(visi)):
					# shorthand notations
					vis = visi[k]
					ini = invi[k]
	
	
					orbit = fast_minute2orbit(times, vis,orbit_id)
					if orbit > param.last_orbits[orbit_id]: continue
	
					#print orbit, fast_minute2orbit(times, ini,orbit_id)
	
					a_end = a_end_orbits[orbit-1]
					if a_end >= ini:
						t_obs_per_orbit[orbit-1] += ini-vis
					else:	
						t_obs_per_orbit[orbit-1] += a_end - vis
	
						orbit_currently_in=orbit
						while True:
							if fast_minute2orbit(times, ini,orbit_id)>orbit_currently_in:
								#print fast_minute2orbit(times, ini,orbit_id),orbit_currently_in
								#else: 
								t_obs_per_orbit[orbit_currently_in]+=period
							else:
								#print orbit_currently_in, fast_minute2orbit(times, ini,orbit_id)
								try: t_obs_per_orbit[orbit_currently_in]+=ini-a_end_orbits[orbit_currently_in-1]	
								except:
									pass
								break
							orbit_currently_in+=1
							if orbit_currently_in>=param.last_orbits[orbit_id]-1:
								t_obs_per_orbit[orbit_currently_in-2]+=ini-a_end_orbits[orbit_currently_in-2]	
								
							
						'''
						if orbit<param.last_orbits[orbit_id]:
							if ini>a_end_orbits[orbit]:
								t_obs_per_orbit[orbit] += a_end_orbits[orbit]-a_end
								print a_end_orbits[orbit], ini, , 
								t_obs_per_orbit[orbit+1] += ini-a_end_orbits[orbit]
							else:
								t_obs_per_orbit[orbit] += ini-a_end
						else: t_obs_per_orbit[orbit-1] += ini-a_end'''
		
				sufficent_obs = np.where(t_obs_per_orbit>=min_t_obs_per_orbit)[0]
				t_obs_per_orbit[t_obs_per_orbit<min_t_obs_per_orbit] = 0
	
				if np.size(sufficent_obs) > 0:
					t_obs_per_orbit=t_obs_per_orbit[sufficent_obs]
					obs_tot[ii] = np.sum(t_obs_per_orbit)
					obs_orbit_tot[ii] = np.size(sufficent_obs)
	
					sufficent_days = fast_orbit2day(times, sufficent_obs, orbit_id)
	#				sufficent_days=np.size(np.unique(sufficent_days))
	
	
					try:
						yy = np.bincount(np.int64(sufficent_days))
						keys = np.nonzero(yy)[0]
						yy = yy[keys]
					except ValueError:
						yy = np.int64(sufficent_days)
		
				
					sufficent_days = np.float64(yy)/nb_orbit_per_day
					sufficent_days = sufficent_days.sum()
	
					if verbose: print '\nThere are %g independant days of observations of %g %% of obs time' % (sufficent_days, 100.*obs_tot[ii]/24./60. / (nb_obs_day*(min_t_obs_per_orbit/period))) 
	
					#print obs_orbit_tot[ii], nb_obs_day*(min_t_obs_per_orbit/period)
	
					if obs_tot[ii]/24./60. >= nb_obs_day*(min_t_obs_per_orbit/period) and sufficent_days >= nb_obs_day :
						count += 1
						#print count
						aini = fast_orbit2a_ini_vect(times, sufficent_obs+1, orbit_id)
						start_obs[ii,:np.size(aini)] = aini
						stop_obs[ii,:np.size(aini)] = aini + t_obs_per_orbit
						interruptions_obs[ii,:] = np.nan
						rat, dect = worthy_targets[ii].Coordinates()
						if verbose: print rat*180./np.pi, dect*180./np.pi, '>>', obs_tot[ii]/24./60., 'days\t', sufficent_days, nb_obs_day
					else:
						obs_tot[ii] = 0.
						obs_orbit_tot[ii] = 0.
	
	#				if np.any(obs_orbit_tot>0.):
	#					print obs_orbit_tot/nb_orbit_per_day
						
			message = '\rComputations done.' 
			sys.stdout.write(message)
			sys.stdout.flush()
	
			print '\nSky coverage for %d days' % nb_obs_day
			"""for min_percentage in range(0, 110, 10):
	
				condition = nb_obs_day*24.*60.*min_percentage/period
				indexes = np.where(obs_tot>condition)
				rslts = np.size(obs_tot[indexes])"""
	
			#condition = nb_obs_day*nb_orbit_per_day
			#indexes = np.where(obs_orbit_tot>=condition)
			rslts = np.shape(obs_orbit_tot[obs_orbit_tot>0])[0]
			#rslts = np.size(obs_orbit_tot[indexes])
			#print rslts
	
			#indexes = np.where(obs_orbit_tot<condition)
			#obs_tot[indexes] = 0.
	
				#print min_percentage, '%\t', rslts, ' targets\t', round(float(rslts)/param.total_nb_targets*100.,3), '%'
	
			print nb_obs_day,'\t',rslts, ' targets\t***', round(float(rslts)/param.total_nb_targets*100.,3), ' % ***'
			print >> output, nb_obs_day,'\t',rslts, ' targets\t', round(float(rslts)/param.total_nb_targets*100.,3), ' %'
	
			if early_stop: continue
	
			if verbose: print obs_tot/24./60.
			np.savez_compressed(folder_misc+output_fname, worthy_targets=worthy_targets, obs_tot=obs_tot)
			print 'Filed saved as %s' % output_fname
			output.close()
			
	"""		exit()
	# Checks which targets are actually visible
		truth_table=[False]*len(worthy_targets)
		iii=0
	
		for ii in range(0, len(worthy_targets)):
			if np.abs(stop_obs[ii,0]) > 1e-5 :
				iii+=1
				truth_table[ii]=True
		truth_table=np.asarray(truth_table)
		
		interruptions_obs=interruptions_obs[truth_table,:]
		start_obs=start_obs[truth_table,:]
		stop_obs=stop_obs[truth_table,:]
		worthy_targets=worthy_targets[truth_table]
		
	###########################################################################
	###########################################################################
	###########################################################################
	
	
	
		#print np.shape(interruptions_obs), np.shape(start_obs), np.shape(stop_obs), np.shape(worthy_targets)
		print '\r%d targets can be observed' % (len(worthy_targets))
	
	
	
		# Saves stuff
		np.savez_compressed(folder_misc+output_fname, worthy_targets=worthy_targets, start_obs=start_obs, stop_obs=stop_obs, interruptions_obs=interruptions_obs)
		print 'Filed saved as %s' % output_fname
	"""
