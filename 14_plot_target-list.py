''' 14-plot_target-list.py
===============================================
AIM:	Given a catalogue of objects, plots when the targets are visible according to their magnitude for a given period of time.

INPUT:	files: 	- <orbit_id>_misc/orbits.dat
		- <orbit_id>_flux/flux_*.dat
	variables: see section PARAMETERS (below)

OUTPUT:	in <orbit_id>_figures/ : (see below for file name definition)

CMD:	python 14-plot_target-list.py

ISSUES: <none known>

REQUIRES:- standard python libraries, specific libraries in resources/ (+ SciPy)
	 - BaseMap --> http://matplotlib.org/basemap/
	 - Structure of the root folder:
	   * <orbit_id>_flux/ --> flux files
	   * <orbit_id>_figures/ --> figures
	   * <orbit_id>_misc/ --> storages of data
	   * all_figures/ --> comparison figures

REMARKS: based on 11-<...>.py, but has a better way of saving appearance and disapperance of the targets, using the class in resources/targets.py
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
# orbit_id
orbit_id = 701
apogee=700
perigee=700

# First minute analysis
minute_ini = 0

# Last minute to look for
minute_end = 1440

# Include SAA ?
SAA = False

# Show plots
show = True

# Save the picture ?
save = True

# Fancy plots ?
fancy = True

# Take into account the stray light?
straylight = True

# Minimum observable time for plots
threshold_obs_time = 50

# Time to acquire a target
t_acquisition = 6

# Catalogue name (in resources/)
catalogue = 'cheops_target_list_v0.1.dat'

# Maximum magnitude that can be seen by CHEOPS, only for cosmetics purposes
CHEOPS_mag_max = 12

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

targets = []

for name, ra, dec, mag in zip(name_cat, ra_cat, dec_cat, mag_cat):
	id_ra = find_nearest(ras, ra/const.RAD)
	id_dec = find_nearest(decs, dec/const.RAD)

	targets.append(target_list(name, ra/const.RAD, id_ra, dec/const.RAD, id_dec, mag, int(period+3)))


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
	for minute in range(minute_ini,minute_end+1+int(period)):
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

			for obj in targets:
				if SAA_at_minute:
					obj.current_visibility = 0
				else:
					obj.current_visibility = obj.visible_save[minute_to_load]
			load = False
			# populate the visbility matrix
#			for ii in range(0, targets[0].CountObjects()):

		if load:
			for obj in targets:
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
			for obj in targets:
				obj.workspace=obj.current_visibility
			continue

		for obj in targets: obj.Next(minute,threshold_obs_time)

except KeyboardInterrupt: print hilite('\nWARNING! USER STOPPED LOADING AT MINUTE %d' % minute,False,False)

for ii in range(0, targets[0].CountObjects()): targets[ii].Next(minute,threshold_obs_time)

### #TODO if first minute look for past orbits anyways
print
worthy_targets = []
for ii in range(0, targets[0].CountObjects()):
	if np.shape(targets[ii].visible)[0] > 0:
		worthy_targets.append(targets[ii])
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

maxy = len(worthy_targets)
print 'Number of star visible in period selected: %d' % maxy

size = 2 + maxy/3
figsize = (17.,size) # fig size in inches (width,height)
fig = plt.figure(figsize=figsize)
ax = plt.subplot(111)

ii = 0
ax.yaxis.set_major_locator(MultipleLocator(1))
plt.grid(True)

for ii in range (0, len(worthy_targets)):
	y = float(ii)

	visi = worthy_targets[ii].Visibility()
	invi = worthy_targets[ii].Invisibility()
	
	for vis, ini in zip(visi, invi):
		plt.hlines(y, vis, ini, lw=3, color=cm.Dark2(y/(maxy+5)))

	if ii > maxy: break
	else: ii+=1

labels = ['%s (%2.1f)' % (wt.name, wt.mag) for wt in worthy_targets[0:maxy]]

ax.set_yticklabels(labels)

ax.set_ylim(-0.5,maxy-0.5)

# convert epoch to matplotlib float format
labels = np.linspace(minute_ini, minute_end+1, 12) * 60. + const.timestamp_2018_01_01

plt.xlim([minute_ini, minute_end+1])
ax.xaxis.set_major_locator(MultipleLocator((minute_end-minute_ini+1)/11))

# to human readable date
pre = map (time.gmtime, labels)
labels = map(figures.format_second, pre)

ax.set_xticklabels(labels)
fig.autofmt_xdate()

if save:
	threshold_obs_time -= t_acquisition
	if SAA: note = '_SAA'
	else: note = ''
	fname = '%svisibility_stars_obs_%d_o_%d_to_%d%s' % (folder_figures,threshold_obs_time,fo,lo, note)
	figures.savefig(fname,fig,fancy)


### A spatial plot of the targets
fig = plt.figure()
ax = plt.subplot(111, projection='mollweide')

plt.scatter((ra_cat-180)/const.RAD,dec_cat/const.RAD, c=mag_cat, marker='*', s=50, edgecolor='none', vmin=param.magnitude_min,vmax=param.magnitude_max+0.2)

v = np.linspace(param.magnitude_min,param.magnitude_max, (param.magnitude_max-param.magnitude_min+1), endpoint=True)
t = map(figures.format_mag, v)
cbar = plt.colorbar(ticks=v, orientation='horizontal',shrink=.8)
cbar.set_ticklabels(t)
l,b,w,h = plt.gca().get_position().bounds
ll,bb,ww,hh = cbar.ax.get_position().bounds
cbar.ax.set_position([ll, bb+0.1, ww, hh])

ax.grid(True)
ax.set_xticklabels([r'$30^{\circ}$',r'$60^{\circ}$',r'$90^{\circ}$',r'$120^{\circ}$',r'$150^{\circ}$',r'$180^{\circ}$',r'$210^{\circ}$',r'$240^{\circ}$',r'$270^{\circ}$',r'$300^{\circ}$',r'$330^{\circ}$']) #,r'$360^{\circ}$'
ax.set_xlabel(r'$\alpha$')
ax.set_ylabel(r'$\delta$')

if save:
	fname = '%stargets_distribution' % folder_figures
	figures.savefig(fname,fig,fancy)

### A histogram of the magnitudes
fig = plt.figure(dpi=100)
ax = fig.add_subplot(111)

bins=np.linspace(np.amin(mag_cat),np.amax(mag_cat), 50)

n, bins, patches = plt.hist(mag_cat,bins=bins)
plt.setp(patches, 'edgecolor', 'black', 'linewidth', 2, 'facecolor','blue','alpha',1)

ax.xaxis.set_major_locator(MultipleLocator(2))
ax.xaxis.set_minor_locator(MultipleLocator(1))

ax.yaxis.set_major_locator(MultipleLocator(5))
ax.yaxis.set_minor_locator(MultipleLocator(1))

ax.xaxis.grid(True,'minor')
ax.yaxis.grid(True,'minor')
ax.xaxis.grid(True,'major',linewidth=2)
ax.yaxis.grid(True,'major',linewidth=2)

plt.xlim([np.amin(mag_cat)*0.95, 1.05*np.amax(mag_cat)])

plt.xlabel(r'$m_V$')
plt.ylabel(r'$\mathrm{distribution}$')

x1,x2,y1,y2 = plt.axis()
plt.axvline(CHEOPS_mag_max, lw=2, color='r')

if save:
	fname = '%stargets_hist_mag' % folder_figures
	figures.savefig(fname,fig,fancy)

if show: plt.show()
