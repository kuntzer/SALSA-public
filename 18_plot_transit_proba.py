''' 18-plot-transit-proba.py
=========================
AIM:	Plots transit probabilities according to 17-treat-ephemerids.py.
	A probability of 100% corresponds to being able to observe the target for
	its whole period.

INPUT:	files: 	- <orbit_id>_misc/ephemerids_obs<transit_duration>h_<max_interruptions>inter_V<mag_max><_SAA?>.npz
	variables: see section PARAMETERS (below)

OUTPUT:	in <orbit_id>_figures/ with the following format:
			proba_<orbit_id>_<exoplanet_period>obs_<max_interruptions>inter_V%3.1f.png/.eps/.pdf

CMD:	python 18-plot-transit-proba.py

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
# orbit_iditude of the orbit in km
orbit_id = 702
apogee = 700
perigee = 700

# First minute in data set !
minute_ini = 0

# Last minute to look for
minute_end = 1440*365/12

# File name for the list of orbit file
orbits_file = 'orbits.dat'

typep = 'Neptunes'

if typep == 'Super-Earths':
	# Orbit period [d]
	exoplanet_period = 50
	# Minimum observable time for plots [h]
	transit_duration = 6
elif typep == 'Neptunes':
	# Orbit period [d]
	exoplanet_period = 13
	# Minimum observable time for plots [h]
	transit_duration = 3

# Maximum visible magnitude
mag_max = 10.

# Plot a few stars as well ?
stars= False
targets_exo=False

# Maximum interruption time tolerated [min]
max_interruptions = 20

# Take SAA into account?
SAA = True

# File name for the input file (in a compressed binary Python format)
if SAA: note = '_SAA'
else: note = ''

# File name for the input file (in a compressed binary Python format)
input_fname = 'ephemerids_obs%dh_%dinter_V%3.1f%s.npz' % (transit_duration,max_interruptions,mag_max,note)

# Print much information ?
verbose = False

# Nice plots?
fancy=True

# Save plots?
save = True

# Show figures ?
show = True

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




# Formatted folders definitions
folder_flux, folder_figures, folder_misc = init_folders(orbit_id)

# loading data
print 'loading %s' % input_fname
sys.stdout.write("Loading worthy targets...\t")
sys.stdout.flush()

data = np.load(folder_misc+input_fname)
worthy_targets = data['worthy_targets']
start_obs=data['start_obs']
stop_obs=data['stop_obs']
interruptions_obs=data['interruptions_obs']

print 'Done, %d targets loaded' % len(worthy_targets)
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
			if verbose: print i/60/24, f/60/24, (f-i)/60
			density[index_target]+=np.floor((f-i)/60 / transit_duration)
	if verbose: print '-'*30

	density[index_target]=float(density[index_target]) / number_of_transit * 100.

	# Associate the density to a grid point
	id_ra = np.where(np.abs(ras-target.Coordinates()[0]) < 0.05)[0]
	id_dec= np.where(np.abs(decs-target.Coordinates()[1]) < 0.05)[0]
	# Transform density in prob of transit:
	data_grid[id_dec, id_ra] = density[index_target]

if verbose: print 'obs start | obs end | hours of obs'

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

if fancy: figures.set_fancy()
fig = plt.figure()
m = Basemap(projection='moll',lon_0=180,resolution='c')

extent = (-np.pi,np.pi,-np.pi/2.,np.pi/2.)

ra_grid *= const.RAD
#ra_grid -= 180.
#ra_grid = ra_grid - 180 #= (ra_grid-np.pi)  #*180. / np.pi
dec_grid *= const.RAD
m.contour( ra_grid,dec_grid,data_grid,10,colors='k',latlon=True)
CS = m.contourf( ra_grid ,dec_grid,data_grid,100,cmap=plt.cm.gist_stern,latlon=True,vmin=0)
#m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
m.drawparallels(np.arange(-60.,90.,30.),labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,30.))


ra__ = np.arange(0., 360., 30.)
#print ra__
x, y = m(ra__,ra__*0)
for x,y,ra in zip(x,y,ra__):
	plt.text(x, y, figures.format_degree(ra), color='black', ha='center', weight='black', size='small') ##93c6ed


t = np.linspace(0., np.amax(density),5)
labels = ['%3.1f\%%' % a for a in t]
cbar = plt.colorbar(CS, orientation='horizontal',shrink=.8, ticks=t)
cbar.ax.set_xticklabels(labels)
l,b,w,h = plt.gca().get_position().bounds
ll,bb,ww,hh = cbar.ax.get_position().bounds
cbar.ax.set_position([ll, bb+0.1, ww, hh])
cbar.set_label('Probabilty of seeing a transit of %d hours for V=%3.1f' % (transit_duration,mag_max))

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
	fname = 'proba_%d_%dobs_%dinter_V%3.1f' % (orbit_id, transit_duration, max_interruptions, mag_max)
	figures.savefig(folder_figures+fname, fig, fancy)
	print 'saved as %s' % folder_figures+fname

if show: plt.show()

