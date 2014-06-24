''' 19-plot-transit-proba-mag-limits.py
=========================
AIM:	Plots deepest achiveable magnitude according to 17-treat-ephemerids.py.
	A minimum detection capability of 100% means that a whole orbit must be observable

INPUT:	files: 	- <orbit_id>_misc/ephemerids_obs<transit_duration>h_<max_interruptions>inter_V<mag_max><_SAA?>.npz
	variables: see section PARAMETERS (below)

OUTPUT:	in <orbit_id>_figures/ with the following format:
			<orbit_id>-inter<max_interruptions>-mag-<typep>-<detection_rate>.png/.eps/.pdf

CMD:	python 19-plot-transit-proba-mag-limits.py

ISSUES:	<none known>

REQUIRES:
	 - Latex
	 - standard python libraries, specific libraries in resources/ (+ SciPy)
	 - Structure of the root folder:
	   * <orbit_id>_flux/ --> flux files
	   * <orbit_id>_figures/maps/ --> figures
	   * <orbit_id>_misc/ --> storages of data

REMARKS: Not with real catalogue.
'''

###########################################################################
### INCLUDES
import numpy as np
import pylab as plt
import matplotlib.cm as cm

from resources.routines import *
from resources.TimeStepping import *

import parameters as param
import resources.constants as const
import resources.figures as figures
from resources.targets import *

from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter
from mpl_toolkits.basemap import Basemap
###########################################################################
### PARAMETERS
# orbit_id
orbit_id = 702
apogee = 700
perigee = 700

# First minute in data set !
minute_ini = 0

# Last minute to look for
minute_end = 1440*365/12

# File name for the list of orbit file
orbits_file = 'orbits.dat'


# Maximum interruption time tolerated [min]
max_interruptions = 97

# Planet type
typep = 'Super-Earths'

if typep == 'coverage':
	# Orbit period [d]
	exoplanet_period = 50
	# Minimum observable time for plots [h]
	transit_duration = 0.0001
	min_detection_rate = 100
elif typep == 'Super-Earths':
	# Orbit period [d]
	exoplanet_period = 50
	# Minimum observable time for plots [h]
	transit_duration = 6
	min_detection_rate = 50
elif typep == 'Neptunes':
	# Orbit period [d]
	exoplanet_period = 13
	# Minimum observable time for plots [h]
	transit_duration = 3
	min_detection_rate = 200

# Maximum visible magnitude
mag_max = 6.
# Minimum visible magnitude
mag_min = 6.
# Magnitude step
mag_sep = 1.



# Plot a few stars as well ?
stars= False
targets_exo=False

# Take SAA into account?
SAA = True

# Print much information ?
verbose = False

# Nice plots?
fancy=True

# Save plots?
save = True

# Show figures ?
show = True

# If set to True, then it will be observations of at least (period - max_interruptions)
# If set to False, then it is minimum (period - max_interruptions) minutes per orbit, 
# not necesseraly consecutive.
consecutive = False

# Minimal # of days of obs (if consecutive = False)
nb_obs_day = 50

# Minimal minutes to be observed per orbit (if consecutive = False)
min_t_obs_per_orbit = 50

if SAA: note = '_SAA'
else: note = ''
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
number_of_transit = exoplanet_period * 24. / transit_duration

###########################################################################
### INITIALISATION

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

data_grid = np.zeros(np.shape(ra_grid))



if stars:
	ra_stars=[101.2833, 95.9875, 213.9167, 219.9, 279.2333, 78.6333, 114.8250, 88.7917]
	dec_stars=[-16.7161, -52.6956, 19.1822, -60.8339, 38.7836, -8.2014, 5.2250, 7.4069]
	y_offset=[0.5e6,0.5e6,-0.8e6,0.5e6,0.5e6,0.5e6,-0.8e6,0.5e6]
	labels = ['Sirius','Canopus','Arcturus',r'$\alpha\mathrm{Centauri}$','Vega','Rigel','Procyon','Betelgeuse']

if targets_exo: ra_tar, dec_tar, magn = np.loadtxt('resources/defined-exo.csv', delimiter=';', unpack=True)


for current_mag in np.arange(mag_max, mag_min-mag_sep, -1.*mag_sep):
	# File name for the input file (in a compressed binary Python format)
	if consecutive:
		input_fname = 'ephemerids_obs%dh_%dinter_V%3.1f%s.npz' % (transit_duration,max_interruptions,current_mag,note)
	else: 
		input_fname = 'ephemerids_%ddays_%dmin_V%3.1f%s.npz' % (nb_obs_day,min_t_obs_per_orbit,current_mag,note)

	# loading data
	print 'loading %s' % input_fname
	sys.stdout.write("Loading worthy targets...\t")
	sys.stdout.flush()

	data = np.load(folder_misc+input_fname)
	worthy_targets = data['worthy_targets']
	start_obs=data['start_obs']
	stop_obs=data['stop_obs']
	interruptions_obs=data['interruptions_obs']

	print 'Done, %d targets loaded for mag %3.1f' % (len(worthy_targets), current_mag)

	###########################################################################
	# cycling through the targets:
	density = np.zeros(len(worthy_targets))
	for index_target, target in enumerate(worthy_targets):
		tar_start = start_obs[index_target,:]
		tar_stop = stop_obs[index_target,:]

#	print target.Coordinates()[0]*180./np.pi, target.Coordinates()[1]*180./np.pi
		if verbose: print index_target, target.Coordinates()[0]*180./np.pi, target.Coordinates()[1]*180./np.pi

		for i, f in zip(tar_start, tar_stop):
			if i >= 0 and f > 0:
				if verbose: print i, f, (f-i)/60
				density[index_target]+=np.floor((f-i)/60 / transit_duration)
		if verbose: print '-'*30

		density[index_target]=float(density[index_target]) / number_of_transit * 1.e2

		if density[index_target] < min_detection_rate : density[index_target] = 0

		# Associate the density to a grid point
		id_ra = np.where(np.abs(ras-target.Coordinates()[0]) < 0.05)[0]
		id_dec= np.where(np.abs(decs-target.Coordinates()[1]) < 0.05)[0]
	# Transform density in prob of transit:
		if data_grid[id_dec, id_ra] == 0 and density[index_target] > 0:
			data_grid[id_dec, id_ra] = current_mag

	if verbose: print 'obs start | obs end | hours of obs'

co = np.size(data_grid[np.where(data_grid>0)])
print 'coverage', float(co)/float(np.size(data_grid))*100, '%'
if typep == 'coverage': exit()

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

if fancy: figures.set_fancy()
fig = plt.figure()
axes=fig.add_subplot(1,1,1, axisbg='black')
m = Basemap(projection='moll',lon_0=180,resolution='c')

extent = (-np.pi,np.pi,-np.pi/2.,np.pi/2.)

ra_grid *= const.RAD
#ra_grid -= 180.
#ra_grid = ra_grid - 180 #= (ra_grid-np.pi)  #*180. / np.pi
dec_grid *= const.RAD

m.contour( ra_grid,dec_grid,data_grid,10,colors='k',latlon=True)
CS = m.contourf( ra_grid ,dec_grid,data_grid,int((mag_max-mag_min)/mag_sep+1),cmap=plt.cm.gist_rainbow,latlon=True)
#m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
m.drawparallels(np.arange(-60.,90.,30.),labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,30.))


ra__ = np.arange(0., 360., 30.)
#print ra__
x, y = m(ra__,ra__*0)
for x,y,ra in zip(x,y,ra__):
	plt.text(x, y, figures.format_degree(ra), color='black', ha='center', weight='black', size='small') ##93c6ed
v = np.linspace(mag_min,mag_max, (mag_max-mag_min+1), endpoint=True)
t = map(figures.format_mag, v)

cbar = plt.colorbar(CS, ticks=v, orientation='horizontal',shrink=.8)
cbar.set_ticklabels(t)

#cbar = plt.colorbar(CS, orientation='horizontal',shrink=.8, ticks=t)
#cbar.ax.set_xticklabels(labels)
l,b,w,h = plt.gca().get_position().bounds
ll,bb,ww,hh = cbar.ax.get_position().bounds
cbar.ax.set_position([ll, bb+0.1, ww, hh])
cbar.set_label(r'$\mathrm{faintest}\ V\ \mathrm{magnitude\ for\ %s\ (%d\%%\ detection)}$' % (typep,min_detection_rate))

if stars:
	x,y = m(ra_stars, dec_stars)
	m.plot(x,y, 'w*', markersize=10)

	for label, xpt, ypt, y_offset in zip(labels, x, y,y_offset):
		plt.text(xpt, ypt+y_offset, label, color='white', size='x-small', ha='center', weight='black') # #93a4ed

if targets_exo:
	x,y = m(ra_tar*180./np.pi, dec_tar*180./np.pi)
	x,y = m(ra_tar, dec_tar)
	m.scatter(x,y, c='white', edgecolor='k', marker='+', s=20,zorder=10, lw=0.5)


# Save plot
if save:
	fname = '%d-inter%d-mag-%s-%d%s' % (orbit_id,max_interruptions,typep,min_detection_rate,note)
	figures.savefig(folder_figures+fname, fig, fancy)
	print 'saved as %s' % folder_figures+fname

if show: plt.show()

