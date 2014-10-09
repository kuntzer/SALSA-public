''' 21-plot_Sky-coverage.py
===============================================
AIM:	Plots the Stray Light coverage in % in terms of period of observation and accumulated observation time.

INPUT:	files: 	- <orbit_id>_misc/ : files from 17-<...>.py
	variables: see section PARAMETERS (below)

OUTPUT:	in <orbit_id>_<SL_angle>figures/ : <orbit_id>_<threshold_obs_time>_<max_mag><_SAA?>_SL_coverage.png/pdf/eps

CMD:	21-plot_Sky-coverage.py

ISSUES:	<NONE KNOWN>

REQUIRES:- standard python libraries, specific libraries in resources/ (+ SciPy)
	 - Structure of the root folder:
	   * <orbit_id>_flux/ --> flux files
	   * <orbit_id>_figures/ --> figures
	   * <orbit_id>_misc/ --> storages of data

REMARKS: <NONE>
'''
###########################################################################
### INCLUDES
import numpy as np
import pylab as plt
import os
import matplotlib.cm as cm
import time

from resources.routines import *
from resources.TimeStepping import *

import parameters as param
import resources.constants as const
import resources.figures as figures
from resources.targets import *
from resources.coordinates import ecliptic2equatorial

from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter
from mpl_toolkits.basemap import Basemap
###########################################################################
### PARAMETERS
# orbit_id
orbit_id = '700_35_AKTAR'
apogee = 700
perigee = 700

# First minute in data set !
minute_ini = 0

# Last minute to look for
minute_end = 1440*365/12

# File name for the list of orbit file
orbits_file = 'orbits.dat'

# Maximum visible magnitude
mag_max = 9.

# Min nb_obs_day
min_nb_obs_day = 10#10
max_nb_obs_day = 90#16#90
step_nb_obs_day = 10#10#1#10

# Take SAA into account?
SAA = True

# Print much information ?
verbose = False


# Factor in the SL post treatment correction ?
SL_post_treat = True

# This is a way to vary the results by multiplying the whole pst by a number.
# This is very easy as if the pst is multiplied by a constant, it can be taken out of the
# integral and only multplying the flux is equivalent to re-running all the simulations
pst_factor = 1.

# If set to True, then it will be observations of at least (period - max_interruptions)
# If set to False, then it is minimum (period - max_interruptions) minutes per orbit, 
# not necesseraly consecutive.
consecutive = False

# Minimal minutes to be observed per orbit (if consecutive = False)
min_t_obs_per_orbit = 49

# Nice plots?
fancy=True

# Save plots?
save = True

# Show figures ?
show = True

# Show ecliptic ?
show_ecliptic=True

# File name for the input file (in a compressed binary Python format)
if SAA: note = '_SAA'
else: note = ''
if not pst_factor == 1.: note += '_%1.1fpst' % pst_factor
if SL_post_treat: note+= '_%4.3fSLreduction' % param.SL_post_treat_reduction
if not consecutive: note += '_cumul_'

#####################################################################################################################

	# for every region in the sky/worthy target:
	# >> Find when you can look with transit_duration [h] with maximal max_interruptions [min]
	# >>>> return start and end time of observations with duration of interruptions [min]
	# >> can we observe a transit ?
	# >>>> Vary the start of transit time by transit_duration [h] until exoplanet_period [h]

#####################################################################################################################
# CONSTANTS AND PHYSICAL PARAMETERS
period = altitude2period(apogee, perigee)

###########################################################################
### INITIALISATION

# Formatted folders definitions
folder_flux, folder_figures, folder_misc = init_folders(orbit_id)


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


for nb_obs_day in np.arange(max_nb_obs_day, min_nb_obs_day-step_nb_obs_day, -1.*step_nb_obs_day):
	# File name for the input file (in a compressed binary Python format)
	input_fname = 'ephemerids_%ddays_%dmin_V%3.1f%s.npz' % (nb_obs_day,min_t_obs_per_orbit,mag_max,note)

	# loading data
	print 'loading %s' % input_fname
	sys.stdout.write("Loading worthy targets...\t")
	sys.stdout.flush()

	data = np.load(folder_misc+input_fname)
	worthy_targets = data['worthy_targets']
	obs_tot=data['obs_tot']
#	start_obs=data['start_obs']
#	stop_obs=data['stop_obs']
#	interruptions_obs=data['interruptions_obs']

	print 'Done, %d targets loaded for nb_obs_day %3.1f' % (len(worthy_targets), nb_obs_day)

	###########################################################################
	# cycling through the targets:
	obs_time = np.zeros(len(worthy_targets))
	for index_target, target in enumerate(worthy_targets):
#		tar_start = start_obs[index_target,:]
#		tar_stop = stop_obs[index_target,:]

#	print target.Coordinates()[0]*180./np.pi, target.Coordinates()[1]*180./np.pi
		if verbose: print index_target, target.Coordinates()[0]*180./np.pi, target.Coordinates()[1]*180./np.pi

		if obs_tot[index_target]>0.:
			obs_time[index_target]=obs_tot[index_target]/60./24.

		'''
		for i, f in zip(tar_start, tar_stop):
			if i >= 0 and f > 0:
				if verbose: print i, f, f-i
				obs_time[index_target]+=f-i
		if verbose: print '-'*30
		'''
		# Associate the density to a grid point
		if target.Coordinates()[0] < np.pi:
			id_ra = np.where(np.abs(ras-target.Coordinates()[0]) < 0.05)[0]
		else:
			id_ra = np.where(np.abs(ras-(target.Coordinates()[0]-2.*np.pi)) < 0.05)[0]
		id_dec= np.where(np.abs(decs-target.Coordinates()[1]) < 0.05)[0]
	# Transform density in prob of transit:
		if data_grid[id_dec, id_ra] == 0 and obs_tot[index_target]>0.:
			data_grid_days[id_dec, id_ra] = nb_obs_day
			data_grid[id_dec, id_ra] = obs_tot[index_target]/60./24.
			if verbose: print target.Coordinates()[0]*180./np.pi,'\t',target.Coordinates()[1]*180./np.pi,'\t', obs_tot[index_target]/24./60.

	if verbose: print 'obs start | obs end | hours of obs'

print np.amin(data_grid), np.amax(data_grid)

co = np.size(data_grid[np.where(data_grid>0)])
print 'coverage', float(co)/float(np.size(data_grid))*100, '%'

#plt.figure()

#for index_target, target in enumerate(worthy_targets):
#	c = density[index_target]
#	plt.scatter(target.Coordinates()[0]*180./np.pi,target.Coordinates()[1]*180./np.pi,c=c, cmap=cm.jet, vmin=np.amin(density), vmax=np.amax(density), edgecolor='none', s=50)

#plt.xlim([0,360])
#plt.ylim([-90,90])
#plt.grid()
#cb=plt.colorbar()
#cb.set_label('Probabilty of transit of min. %d hours' % transit_duration)

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

plt.grid()

ra_grid *= const.RAD
dec_grid *= const.RAD

CS = plt.contour( ra_grid,dec_grid,data_grid,int((max_nb_obs_day-min_nb_obs_day)/step_nb_obs_day),colors='k')

plt.clabel(CS, inline=1,fmt='%d',colors='red', fontsize=12)

CS = plt.contourf( ra_grid ,dec_grid,data_grid,200,cmap=plt.cm.winter)

plt.yticks(np.arange(-80, 100, 20.))

v = np.linspace(min_nb_obs_day,max_nb_obs_day, (max_nb_obs_day-min_nb_obs_day)/step_nb_obs_day+1, endpoint=True)
 
cbar = plt.colorbar(CS)#, ticks=v)
#cbar.set_ticklabels(v)
cbar.set_label(r'$\mathrm{Days}$')

plt.xlabel('RA [deg]')
plt.ylabel('Dec [deg]')

if show_ecliptic:
	a=np.linspace(-np.pi, np.pi)
	b=np.zeros_like(a)
	res=np.rad2deg(ecliptic2equatorial(a,b))
	plt.plot(res[:,0],res[:,1],lw=2,color="red")

	# Sun in june
	plt.plot([-90.],[23.433],'o',color="yellow", markersize=8, zorder=5)
	plt.text(-85.,23.433, r"$\odot_\mathrm{June}$", color='k', size='small', weight='black')
	# Sun in december
	plt.plot([90.],[-23.433],'o',color="yellow", markersize=8, zorder=5)
	plt.text(95.,-23.433, r"$\odot_\mathrm{Dec}$", color='k', size='small', weight='black')

###########################################################################
fig2 = plt.figure()
ax=plt.subplot(111)
ax.set_aspect(2.)

plt.grid()
CS = plt.contour( ra_grid,dec_grid,data_grid_days,int((max_nb_obs_day-min_nb_obs_day)/step_nb_obs_day),colors='k')
plt.clabel(CS, inline=1,fmt='%d',colors='red', fontsize=12)

CS = plt.contourf( ra_grid ,dec_grid,data_grid_days,200,cmap=plt.cm.winter)

plt.yticks(np.arange(-80, 100, 20.))

v = np.linspace(min_nb_obs_day,max_nb_obs_day, (max_nb_obs_day-min_nb_obs_day)/step_nb_obs_day+1, endpoint=True)
 
cbar = plt.colorbar(CS, ticks=v)
cbar.set_ticklabels(v)
cbar.set_label(r'$\mathrm{Days}$')

plt.xlabel('RA [deg]')
plt.ylabel('Dec [deg]')

if show_ecliptic:
	plt.plot(res[:,0],res[:,1],lw=2,color="red")

	# Sun in june
	plt.plot([-90.],[23.433],'o',color="yellow", markersize=8, zorder=5)
	plt.text(-85.,23.433, r"$\odot_\mathrm{June}$", color='white', size='small', weight='black')
	# Sun in december
	plt.plot([90.],[-23.433],'o',color="yellow", markersize=8, zorder=5)
	plt.text(95.,-23.433, r"$\odot_\mathrm{Dec}$", color='white', size='small', weight='black')

###########################################################################
if not SAA: note = '_noSAA'
else: note = '_SAA'
if not pst_factor == 1.: note += '_%1.1fpst' % pst_factor
# Save plot
if save:
	fname = '%s-sky_map-%d-mag%d-days%d-to-%d%s-accumulated' % (orbit_id,min_t_obs_per_orbit,mag_max,min_nb_obs_day,max_nb_obs_day,note)
	figures.savefig(folder_figures+fname, fig, fancy)
	np.savez_compressed(folder_misc+fname, ra_grid=ra_grid, dec_grid=dec_grid, data_grid=data_grid, ticks=int((max_nb_obs_day-min_nb_obs_day)/step_nb_obs_day))
	print 'saved as %s' % fname


	fname = '%s-sky_map-%d-mag%d-days%d-to-%d%s-period_obs' % (orbit_id,min_t_obs_per_orbit,mag_max,min_nb_obs_day,max_nb_obs_day,note)
	figures.savefig(folder_figures+fname, fig2, fancy)
	np.savez_compressed(folder_misc+fname, ra_grid=ra_grid, dec_grid=dec_grid, data_grid=data_grid_days, ticks=int((max_nb_obs_day-min_nb_obs_day)/step_nb_obs_day))
	print 'saved as %s' % fname


if show: plt.show()

