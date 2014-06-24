''' 13-plot_long_periods.py
===============================================
AIM:	Prepare cumulative plots (THIS SCRIPT IS with STRAY LIGHT)

INPUT:	files: 	- <orbit_id>_misc/ : files from 12-<...>.py or 11-<...>.py
	variables: see section PARAMETERS (below)

OUTPUT:	in <orbit_id>_figures/ : <height>_<threshold_obs_time>_<max_mag><_SAA?>.png/pdf/eps

CMD:	python 13-plot_long_periods.py

ISSUES:	<none known>

REQUIRES:- standard python libraries, specific libraries in resources/ (+ SciPy)
	 - BaseMap --> http://matplotlib.org/basemap/
	 - Structure of the root folder:
	   * <orbit_id>_flux/ --> flux files
	   * <orbit_id>_figures/ --> figures
	   * <orbit_id>_misc/ --> storages of data
	   * all_figures/ --> comparison figures

REMARKS: see typeplot to know which map to plot
'''
###########################################################################
### INCLUDES
import numpy as np
import pylab as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm

from resources.routines import *
from resources.TimeStepping import *
import parameters as param
import resources.constants as const
import resources.figures as figures

###########################################################################
### PARAMETERS
# orbit_id
orbit_id = 301

# Show plots ?
show = True

# Save the picture ?
save = True

# max_mag
max_mag = 12.

# Plot a few stars as well ?
stars= False
targets=False

# Fancy plots ?
fancy = True

# Scale of the plot il log form ?
log_plot = False

# SAA ?
SAA = True

threshold_obs_time = 50

# what to plot ?
# mag -> magnitude
# raw -> using raws maps data (unused)
# anything else: without magnitude
typeplot = 'mag'

## Do not touch
if SAA: note = '_SAA'
else: note = ''
## end do not touch

# File name for the input data file (WITHOUT EXTENSION)
input_fname_wo_mag = 'TEST-data_%d%s' % (threshold_obs_time, note)
input_fname_wi_mag = 'mag_over_year_%d_mag_%02.1f%s' % (threshold_obs_time, max_mag, note)
input_fname_raw    = 'data_raw_%d%s' % (threshold_obs_time, note)
# Extension (example: .dat)
ext = '.dat'

###########################################################################
### INITIALISATION

if typeplot == 'mag': input_fname = input_fname_wi_mag
elif typeplot == 'raw': input_fname = input_fname_raw
else: input_fname =  input_fname_wo_mag
input_fname += ext

print 'Loading %s' % input_fname

if fancy: figures.set_fancy()

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

visibility = np.zeros(np.shape(ra_grid))
workspace = np.zeros(np.shape(ra_grid))
data = np.zeros(np.shape(ra_grid))

if stars:
	ra_stars=[101.2833, 95.9875, 213.9167, 219.9, 279.2333, 78.6333, 114.8250, 88.7917]
	dec_stars=[-16.7161, -52.6956, 19.1822, -60.8339, 38.7836, -8.2014, 5.2250, 7.4069]
	y_offset=[0.5e6,0.5e6,-0.8e6,0.5e6,0.5e6,0.5e6,-0.8e6,0.5e6]
	labels = ['Sirius','Canopus','Arcturus',r'$\alpha\mathrm{Centauri}$','Vega','Rigel','Procyon','Betelgeuse']

if targets: ra_tar, dec_tar = np.loadtxt('resources/targets.dat', unpack=True)
if targets: ra_tar, dec_tar, magn = np.loadtxt('resources/defined-exo.csv', delimiter=';', unpack=True)

############################################################################
### LOADS AND PLOTS
data = np.loadtxt(folder_misc+input_fname)

sky = np.size(data)
seeable_points = 0.
for i in range(0,np.shape(data)[0]):
	seeable_points += np.size(data[i,data[i,:]>0])
seeable_points = seeable_points / sky * 100.

# mask point in the map for which observation time is zero
#data[data<1] = np.nan

fig = plt.figure()
m = Basemap(projection='moll',lon_0=180,resolution='c')

extent = (-np.pi,np.pi,-np.pi/2.,np.pi/2.)

if log_plot:
	from matplotlib.colors import LogNorm
	max_level = np.ceil(np.max(data/60.))

	levels = np.logspace(np.log10(1), np.log10(max_level),100)
	levels_lines_cb = np.logspace(np.log10(1), np.log10(max_level),10, endpoint=True)

	levels_lines = np.linspace(1, max_level,5)

	fmt={}
	for l in levels_lines:
		fmt[l] = '%3.0f' % l

	fmt_cb={}
	for l in levels_lines_cb:
		fmt_cb[l] = '%3.0f' % l
else:
	levels_lines = 10
	levels = 100

ra_grid *= const.RAD
#ra_grid -= 180.
#ra_grid = ra_grid - 180 #= (ra_grid-np.pi)  #*180. / np.pi
dec_grid *= const.RAD
CS1=m.contour( ra_grid,dec_grid,data/60.,levels_lines,colors='k',latlon=True)
if log_plot:
	CS = m.contourf( ra_grid ,dec_grid,data/60.,levels,norm=LogNorm(), cmap=plt.cm.jet,latlon=True)
else: 
	CS = m.contourf( ra_grid ,dec_grid,data/60.,levels, cmap=plt.cm.jet,latlon=True)
#m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
m.drawparallels(np.arange(-60.,90.,30.),labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,30.))


if stars:
	x,y = m(ra_stars, dec_stars)
	m.plot(x,y, 'w*', markersize=10)

	for label, xpt, ypt, y_offset in zip(labels, x, y,y_offset):
		plt.text(xpt, ypt+y_offset, label, color='white', size='x-small', ha='center', weight='black') # #93a4ed

ra__ = np.arange(0., 360., 30.)
#print ra__
x, y = m(ra__,ra__*0)
for x,y,ra in zip(x,y,ra__):
	plt.text(x, y, figures.format_degree(ra), color='black', ha='center', weight='black', size='small') ##93c6ed

if targets:
	x,y = m(ra_tar*180./np.pi, dec_tar*180./np.pi)
	x,y = m(ra_tar, dec_tar)
	m.scatter(x,y, c='white', edgecolor='k', marker='+', s=20,zorder=10, lw=0.5)

if log_plot:
	plt.clabel(CS1, inline=1, fontsize=10, fmt=fmt)
	cbar = plt.colorbar(CS,ticks=levels_lines_cb, orientation='horizontal',shrink=.8, format='%.f', spacing='proportional')
else:
	cbar = plt.colorbar(CS, orientation='horizontal',shrink=.8)
l,b,w,h = plt.gca().get_position().bounds
ll,bb,ww,hh = cbar.ax.get_position().bounds
cbar.ax.set_position([ll, bb+0.1, ww, hh])
cbar.set_label(r'$\mathrm{Observation\ time\ [Hours]}$')



# Save plot
if save:
	if SAA: note='_SAA' 
	else: note=''
	if log_plot: note = '%s_log' % note

	if typeplot == 'mag': fname = '%d_%d_mag_%3.1f%s' % (orbit_id, threshold_obs_time, max_mag, note)
	elif typeplot == 'raw': fname = '%d_%d%s' % (orbit_id, threshold_obs_time, note)
	else: fname =  '%d_%d%s' % (orbit_id, threshold_obs_time, note)
	figures.savefig(folder_figures+fname, fig, fancy)
	print 'saved as %s' % folder_figures+fname

if show: plt.show()

print 'Percentage of the sky which is visible: %3.1f%%' % seeable_points

