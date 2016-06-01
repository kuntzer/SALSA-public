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
import pylab as plt

from resources.routines import *
from resources.TimeStepping import *

import resources.constants as const
import parameters as param
from resources.coordinates import ecliptic2equatorial

import resources.figures as figures
from resources.targets import *

from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter
###########################################################################
### PARAMETERS
# Orbit id
orbit_id = '6am_700_25_conf4'
apogee=700
perigee=700

ddate = [{'name':'1st January', 'dayofyear':1}, {'name':'5th April', 'dayofyear':95}, {'name':'21st June', 'dayofyear':172}, {'name':'5th September', 'dayofyear':248}, {'name':'21st December', 'dayofyear':355}]
ddate = [{'name':'21 January 2018', 'dayofyear':20},
	 {'name':'21 February 2018', 'dayofyear':51},
	 {'name':'21 March 2018', 'dayofyear':79},
	 {'name':'21 April 2018', 'dayofyear':110},
	 {'name':'21 May 2018', 'dayofyear':140},
	 {'name':'21 June 2018', 'dayofyear':171},
	 {'name':'21 July 2018', 'dayofyear':201},
	 {'name':'19 August 2018', 'dayofyear':230},
	 {'name':'18 September 2018', 'dayofyear':260},
	 {'name':'17 October 2018', 'dayofyear':289},
	 {'name':'16 November 2018', 'dayofyear':309},
	 {'name':'15 December 2018', 'dayofyear':343},
]

# Maximum interruption time tolerated [min]
max_interruptions = 97

# Maximum visible magnitude
mag_max = 9.

# Take SAA into account?
SAA = True

# Show plot ?
show = False

# Save plot ?
save = True

fancy = True

# min of scale
min_val=0
# max of scale
max_val=24
#
step_scale=2

# Print much information ?
verbose = False

# Factor in the SL post treatment correction ?
SL_post_treat = True

# Stop before saving results to file.
early_stop = False

# Minimal minutes to be observed per orbit (if consecutive == False), must be a list
min_t_obs_per_orbit = 0 # 79

# This is a way to vary the results by multiplying the whole pst by a number.
# This is very easy as if the pst is multiplied by a constant, it can be taken out of the
# integral and only multplying the flux is equivalent to re-running all the simulations
pst_factor=1.



for iddate in range(len(ddate)):

	date_looked_for = ddate[iddate]['name']

	# First minute in data set !
	minute_ini = (ddate[iddate]['dayofyear'] - 1) * 1440

	# Last minute to look for
	minute_end = (ddate[iddate]['dayofyear']) * 1440

	print '*'*30, 'min_t_obs_per_orbit %1.1f' % min_t_obs_per_orbit

	
	# File name for the input file (in a compressed binary Python format)
	if SAA: note = '_SAA'
	else: note = ''
	if not pst_factor == 1.: note += '_%1.1fpst' % pst_factor
	if SL_post_treat: note+= '_%4.3fSLreduction' % param.SL_post_treat_reduction
	input_fname = 'ephemerids_inter_%d_mag_%3.1f%s.npz' % (max_interruptions,mag_max,note)


	skycoverage_fname = 'skycoverage_%dmin_V%3.1f%s.txt' % (min_t_obs_per_orbit,mag_max,note)
	#####################################################################################################################
	# CONSTANTS AND PHYSICAL PARAMETERS
	period = altitude2period(apogee, perigee)
	###########################################################################
	### INITIALISATION
	# Formatted folders definitions
	folder_flux, folder_figures, folder_misc = init_folders(orbit_id)

	sky_coverage=0.
	
	print 'ORBIT ID:\t\t%s\nPST factor:\t\t%d\nmin_t_obs_per_orbit\t%d (%.1f%%)\nMAGNITIUDE:\t\t%02.1f\nSAA :\t%g' % (orbit_id,pst_factor,min_t_obs_per_orbit,min_t_obs_per_orbit/period*100., mag_max, SAA)
	
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
	
	count = 0
	period = altitude2period(apogee,perigee)
	check=np.zeros(len(worthy_targets))
	
	sky_coverage=0.
	
	for ii in range(len(worthy_targets)):
		y = float(ii)
		message = '\r%3.1f %%' % (y/float(len(worthy_targets))*100.)
		sys.stdout.write(message)
		sys.stdout.flush()
	
		visi = worthy_targets[ii].Visibility()
		invi = worthy_targets[ii].Invisibility()
		inter = worthy_targets[ii].get_interruption_time()

		# *much* slower than with broadcast, but so much easier too because of the conditions
		for a, d, i in zip(visi, invi, inter) :
			if (a >= minute_ini or d >= minute_end) and a <= minute_end :
				o = d - a -i
				if o < min_t_obs_per_orbit : continue

				if a < minute_ini :				
					effective_percentage = (d - minute_end) / (d - a)
					a = minute_ini
					i = effective_percentage * i
					o = d - a -i

				# Some of the observations can continue after the end time, stop them and compensate the interruption time by reducing the amount of interruption by the fraction of time spent observing
				if d > minute_end:
					effective_percentage = (minute_end - a) / (d - a)
					d = minute_end
					i = effective_percentage * i
					o = d - a -i

				check[ii] += o
		'''exit()
		ids_in_range = np.where(np.logical_and(np.logical_or(visi >= minute_ini, visi >= minute_ini), visi <= minute_end))
		print ids_in_range
		ids_in_range = ids_in_range[0]
		if len(ids_in_range) == 0: continue

		visi = visi[ids_in_range]
		invi = invi[ids_in_range]
		inter = inter[ids_in_range]

		observations = invi - visi - inter

		validated_ids = observations>=min_t_obs_per_orbit

		validated_observations = observations[validated_ids]
		vinter = inter[validated_ids]
		vvis = visi[validated_ids]
		vinvi = invi[validated_ids]

		# Some of the observations can start before the ini time, stop them and compensate the interruption time by reducing the amount of interruption by the fraction of time spent observing
		before_end_obs_id = np.where(vvis < minute_ini)[0]
		effective_percentage = (vinvi[before_end_obs_id] - minute_ini) / (vinvi[before_end_obs_id] - vvis[before_end_obs_id])
		vinter[before_end_obs_id] = effective_percentage * vinter[before_end_obs_id]
		vvis[before_end_obs_id] = minute_ini

		validated_observations = vinvi - vvis - vinter

		# Some of the observations can continue after the end time, stop them and compensate the interruption time by reducing the amount of interruption by the fraction of time spent observing
		after_end_obs_id = np.where(vinvi > minute_end)[0]
		effective_percentage = (minute_end - vvis[after_end_obs_id]) / (vinvi[after_end_obs_id] - vvis[after_end_obs_id])
		vinter[after_end_obs_id] = effective_percentage * vinter[after_end_obs_id]
		vinvi[after_end_obs_id] = minute_end

		validated_observations = vinvi - vvis - vinter

		if np.size(validated_observations)>0:
			if validated_observations.sum() > 1440:
				for i, vo in enumerate(validated_observations):
					print 
					print vo, vvis[i], vinvi[i], vinter[i], vinter[i]/(vinvi[i] - vvis[i])
				exit()
			check[ii] += validated_observations.sum()
			"""#print validated_observations; 
			obs_efficiency_in_orbit = validated_observations/period
			time_lost = np.ceil(obs_efficiency_in_orbit) - obs_efficiency_in_orbit
			#print obs_efficiency_in_orbit
			#print validated_visi; print validated_ids
			print
			for i, vo in enumerate(validated_observations):
				if i == len(validated_observations) : continue
				print (vvis[i+1] - vinvi[i]) / period, vinter[i] / period
				print vinvi[i] - vvis[i], vinvi[i] - vvis[i]-vinter[i], vo, vo >=min_t_obs_per_orbit
				if i == len(validated_visi) - 1 : continue
				dt_between_obs = (validated_visi[i + 1] - vs)/period
				if dt_between_obs < 0.9 : 
					print
					print dt_between_obs
					print (validated_visi[i + 1] - vs)
					print vs
					print validated_visi[i + 1]
					print validated_observations[i]
			exit()"""
		'''
		rat, dect = worthy_targets[ii].Coordinates()
		sky_coverage+=0.5/param.resx/param.resy*np.pi*np.cos(dect)
	
	message = '\rComputations done.' 
	sys.stdout.write(message)
	sys.stdout.flush()
	
	print '\nSky coverage for', ddate[iddate]['name']

	print '\t***', round(sky_coverage*100.,3), ' % ***'

	worthy_targets=worthy_targets
	obs_tot=check




	
	###########################################################################
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
	
	data_grid = np.zeros(np.shape(ra_grid))
	data_grid_days = np.zeros(np.shape(ra_grid))
	
	###########################################################################
	# cycling through the targets:
	obs_time = np.zeros(len(worthy_targets))
	totxtdata = []
	for index_target, target in enumerate(worthy_targets):
	#	tar_start = start_obs[index_target,:]
	#	tar_stop = stop_obs[index_target,:]
	
	#print target.Coordinates()[0]*180./np.pi, target.Coordinates()[1]*180./np.pi
		#if verbose: print index_target, target.Coordinates()[0]*180./np.pi, target.Coordinates()[1]*180./np.pi
	
		if obs_tot[index_target]>0.:
			obs_time[index_target]=obs_tot[index_target]/60.##/1440. * 100.
	
		# Associate the density to a grid point
		if target.Coordinates()[0] < np.pi:
			id_ra = np.where(np.abs(ras-target.Coordinates()[0]) < 0.05)[0]
		else:
			id_ra = np.where(np.abs(ras-(target.Coordinates()[0]-2.*np.pi)) < 0.05)[0]
		id_dec= np.where(np.abs(decs-target.Coordinates()[1]) < 0.05)[0]
	
		if data_grid[id_dec, id_ra] == 0 and obs_tot[index_target]>0.:
			data_grid_days[id_dec, id_ra] = 1
			data_grid[id_dec, id_ra] = obs_tot[index_target]/60.#1440. * 100.
			if verbose: print target.Coordinates()[0]*180./np.pi,'\t',target.Coordinates()[1]*180./np.pi,'\t', obs_tot[index_target]/60.##/1440. * 100.
			totxtdata.append([target.Coordinates()[0], target.Coordinates()[1], obs_tot[index_target]])
	
	if verbose: print 'obs start | obs end | hours of obs'
	
	print np.amin(data_grid), np.amax(data_grid)

	print np.shape(obs_tot)
	totxtdata = np.asarray(totxtdata)
	
	
	###########################################################################
	### Plotting
	# transform 0 into no plotting in the data matrix
	mag_min= np.amin(data_grid[data_grid>0])
	data_grid[data_grid < mag_min] = np.nan
	
	mag_min= np.amin(data_grid_days[data_grid_days>0])
	data_grid_days[data_grid_days < mag_min] = np.nan
	
	if fancy: figures.set_fancy()
	fig = plt.figure()
	ax=plt.subplot(111)
	ax.set_aspect(2.)
	
	min_nb_obs_day = np.nanmin(data_grid)
	max_nb_obs_day = np.nanmax(data_grid)
	
	plt.grid()
	
	ra_grid *= const.RAD
	dec_grid *= const.RAD
	data_grid[data_grid<min_nb_obs_day]=0
	
	v = np.arange(min_val,max_val+step_scale, step_scale)
	
	CS = plt.contour(ra_grid,dec_grid,data_grid,colors='k',levels=v)
	
	plt.clabel(CS, inline=1,fmt='%d',colors='red', fontsize=12)
	
	CS = plt.contourf(ra_grid, dec_grid, data_grid, levels=v, cmap=plt.cm.winter)
	
	plt.yticks(np.arange(-80, 100, 20.))
	
	
	#print v
	#print np.nanmin(data_grid)
	#print np.nanmax(data_grid)
	
	#v = np.arange(0,1440, 60)
	
	cbar = plt.colorbar(CS, ticks=v)
	#cbar.ax.set_yticklabels([r'$%g\%%$' % tv for tv in v])
	
	cbar.set_label(r'$\mathrm{Observation\ time\ [h]}$')
	
	plt.xlabel('RA [deg]')
	plt.ylabel('Dec [deg]')
	plt.title(date_looked_for)	
	
	###########################################################################
	if not SAA: note = '_noSAA'
	else: note = '_SAA'
	if not pst_factor == 1.: note += '_%1.1fpst' % pst_factor
	# Save plot
	if save:
		fname = '%s-sky_map-%d-mag%d_onday%d%s' % (orbit_id,min_t_obs_per_orbit,mag_max,ddate[iddate]['dayofyear'],note)
		figures.savefig(folder_figures+fname, fig, fancy)
	
		print 'figure saved as %s' % fname

		np.savetxt(folder_misc+fname+'.dat', totxtdata, fmt='%1.3f, %1.3f, %d')#, header='ra [rad], dec [rad], obstime [min]')
		print 'ASCII file saved as %s' % fname


if show: plt.show()



