''' 9-plot_flux.py
=========================
AIM:	Plot maps of the stray light flux or equivalent magnitude given a particular date

INPUT:	files: 	- <orbit_id>_misc/orbits.dat
		- <orbit_id>_flux/flux_*.dat
		- resources/moon_*.dat, sun_*.dat, orbits_*.dat
	variables: see section PARAMETERS (below)

OUTPUT:	in <orbit_id>_figures/maps/ : map with the following name: flux_%07d.png

CMD:	python 9-plot_flux.py

ISSUES:	- Boundary zone when observing zone is centred on 360 deg

REQUIRES:- standard python libraries, specific libraries in resources/ (+ SciPy)
	 - Structure of the root folder:
	   * <orbit_id>_flux/ --> flux files
	   * <orbit_id>_figures/maps/ --> figures
	   * <orbit_id>_misc/ --> storages of data
	   * all_figures/ --> comparison figures

REMARKS: <none>
'''

#####################################################################################################################
# DEFINITIONS AND INCLUDES
import numpy as np
import pylab as plt

import os
import random
import time

import copy
from scipy.interpolate import griddata
from matplotlib.patches import Rectangle, Circle

from resources.routines import *
from resources.TimeStepping import *
import parameters as param
import resources.constants as const
import resources.figures as figures

from matplotlib.collections import PatchCollection
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

#####################################################################################################################
# PARAMETERS
# Orbital elements
apogee=650
perigee=650
orbit_id = '6am_650_5_conf4e'

# First minute in data set !
minute_ini = 75#1440 * 172
minute_end = 81#1440 * 172 + 100

# File name for the output data file
orbits_file = 'orbits.dat'

# is true outputs .eps and .pdf for every step. Much slower and heavier
fancy = False 

# compare the data to the flux of different stars
magnitudes = False 

# Show stray light contour map ?
straylight = True

# Draw boundaries ?
boundaries = False

# Speeds up by not loading the different position files, does not save, but shows --> for new implementation or debugging
future = True
save = False
# TODO:
# - For minute (orbit_ID 800) 363500, contour of zone better

file_orbit      = 'orbit_%s.dat' % orbit_id
file_sun        = 'sun_%s.dat' % orbit_id
file_moon       = 'moon_%s.dat' % orbit_id

# Factor in the SL post treatment correction ?
SL_post_treat = True
# Factor in mirror efficiency for the equivalent star magnitude ?
mirror_correction = False
#####################################################################################################################
# CONSTANTS AND PHYSICAL PARAMETERS
period = altitude2period(apogee,perigee)

sl_min = 1e-9
sl_max = 0.1
#####################################################################################################################
# INITIALISATION
file_flux = 'flux_'
# Formatted folders definitions
folder_flux, folder_figures, folder_misc = init_folders(orbit_id)
folder_figures= '%s_figures/maps/' % (orbit_id)

#params = {'backend': 'ps','axes.labelsize': 14,'text.fontsize': 18,'legend.fontsize': 18,'xtick.labelsize': 14,'ytick.labelsize': 14,'text.usetex': True}
#plt.rcParams.update(params)

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino'],'size':14})
rc('text', usetex=True)

print '\nObservability Map Plotting'
print 'ORBIT ID:', orbit_id
print '-------------------------------------------'

if not os.path.isdir(folder_figures):
	print '\tError: figure folder %s does not exists.' % (folder_figures)
	exit()

ra = dec = S_sl=np.zeros(1)

# load the postions of the Moon, Sun and Earth
if not future :

	sys.stdout.write("Loading orbit file...\t\t\t")
	sys.stdout.flush()

	try:
		sat = np.loadtxt('resources/'+file_orbit, delimiter='\t')
	except ValueError:
		sat = np.loadtxt('resources/'+file_orbit, delimiter=' ')

	# apply time conditions
	sat = sat[sat[:,0] >= minute_ini]
	sat = sat[sat[:,0] <= minute_end]
	print "Done."

	sys.stdout.write("Loading Sun file...\t\t\t")
	sys.stdout.flush()

	try:
		sun = np.loadtxt('resources/'+file_sun, delimiter=' ')
	except ValueError:
		sun = np.loadtxt('resources/'+file_sun, delimiter='\t')

	sun = sun[sun[:,0] >= minute_ini]
	sun = sun[sun[:,0] <= minute_end]
	print "Done."

	sys.stdout.write("Loading Moon file...\t\t\t")
	sys.stdout.flush()

	try:
		moon = np.loadtxt('resources/'+file_moon, delimiter=',')
	except ValueError:
		moon = np.loadtxt('resources/'+file_moon, delimiter=' ')

	moon = moon[moon[:,0] >= minute_ini]
	moon = moon[moon[:,0] <= minute_end]
	print "Done."
################################################################################
# Prepare the grid
n_alpha = param.resx
n_delta = param.resy

ra_i = 0
ra_f = 2.*np.pi

dec_i = -np.pi/2.
dec_f = np.pi/2.

ra_step = (ra_f-ra_i)/n_alpha
dec_step = (dec_f-dec_i)/n_delta

iterable = (ra_i + i*ra_step for i in range(n_alpha))
ras = np.fromiter(iterable, np.float)-np.pi

iterable = (dec_i + i*dec_step for i in range(n_delta))
decs = np.fromiter(iterable, np.float)

ra_grid, dec_grid = np.meshgrid(ras, decs)

grid_points0 = np.zeros(np.shape(ra_grid))

iterable = (ra_i + ra_step/2 + i*ra_step for i in range(n_alpha))
ras2 = np.fromiter(iterable, np.float)-np.pi

iterable = (dec_i + dec_step/2 + i*dec_step for i in range(n_delta))
decs2 = np.fromiter(iterable, np.float)

#####################################################################################################################
# Prepares the list of minutes

sys.stdout.write("Loading computed orbits...\t\t")
sys.stdout.flush()

orbits = np.loadtxt(folder_misc+orbits_file,dtype='i4')

list_minutes = -1. * np.ones( ( np.shape(orbits)[0] + 2 ) * period )

id_min = 0
times = np.loadtxt('resources/minute_table_%s.dat' % orbit_id, delimiter=',',dtype='Int32')
for ii, orbit_current in enumerate(orbits[:,0]):
	t_ini, t_end, a_ini, a_end = fast_orbit2times(times,orbit_current,orbit_id)
	for minute in range(a_ini, a_end+1):
		list_minutes[id_min] = int(minute)
		id_min += 1

list_minutes = list_minutes[list_minutes > -1]

# apply time conditions
list_minutes = list_minutes[list_minutes >= minute_ini]
list_minutes = list_minutes[list_minutes <= minute_end]

print 'Done.'


#####
def draw_boundaries(ax,xx,yy,ra_step,dec_step, force=False):
	# sort the outer points rights to draw the limits correctly
	xx, yy, ox, oy = sort_boundary(xx,yy,ra_step,dec_step, force)

	# the last point is the first to close the path
	xx = np.append(xx,xx[0])
	yy = np.append(yy,yy[0])

	# Format the vertices of the path correctly ie [(x1,y1), (x2,y2), ..., (x1, y1)]
	verts = np.array(zip(xx,yy))
	# Create the path and say which points are the 1st, verticies or last point
	Path = mpath.Path
	codes = [Path.MOVETO]
	for ii in range(1, len(xx)-1): codes.append(Path.LINETO)
	codes.append(Path.CLOSEPOLY)

	path = mpath.Path(verts+[0,0], codes)
	patch = mpatches.PathPatch(path, facecolor='none', edgecolor='black',zorder=100)
	patch = ax.add_patch(patch)
	
	return ox, oy

####

if not magnitudes: from matplotlib import colors, ticker, cm

#####################################################################################################################
# Loops on every time step to get the data

for id_min, minute in enumerate(list_minutes):
	minute = int(minute)
	sys.stdout.write("Plotting stray light map "+str(minute)+'...\t')
	sys.stdout.flush()
	#####################################################################################################################
	# Loads the data
	# As it's produced by fortran, it may be a bit sketchy. We read line by line (which is slower than the native numpy method of reading text btw and filter out all space af$
	try:
		ra, dec, S_sl = load_flux_file(minute, file_flux, folder=folder_flux)

		# Apply the flux correction (SL post-treatment removal and the mirror efficiency)
		if mirror_correction: S_sl /= param.mirror_efficiency
		if SL_post_treat: S_sl *= (1.0 - param.SL_post_treat_reduction)
		S_sl *= param.SL_QE

		# Now we compare the flux with the magnitude of the star.
		if magnitudes: S_sl = flux2mag(S_sl,param.ppm_threshold)

		# Data is loaded and processed. --> We do a plot
		doplot = True
	except IOError: # There is no file. That means that there is nothing to plot
		raise 
		S_sl=np.zeros(1)
		ra = dec = S_sl

		doplot = False
		sys.stdout.write("Warning: No point\t")
		sys.stdout.flush()
	#####################################################################################################################
	# Start ploting for minute min

	plt.figure()
	ax = plt.subplot(111, projection="mollweide")
	ax.grid(True)
	ax.set_xticklabels([r'$30^{\circ}$',r'$60^{\circ}$',r'$90^{\circ}$',r'$120^{\circ}$',r'$150^{\circ}$',r'$180^{\circ}$',r'$210^{\circ}$',r'$240^{\circ}$',r'$270^{\circ}$',r'$300^{\circ}$',r'$330^{\circ}$']) #,r'$360^{\circ}$'
	ax.set_xlabel(r'$\alpha$')
	ax.set_ylabel(r'$\delta$')
	#####################################################################################################################
	# Prepare plots...
	extent = (-np.pi,np.pi,-np.pi/2.,np.pi/2.)
	if doplot:


		if magnitudes: 
			vf = np.linspace(param.magnitude_min,param.magnitude_max+0.2, (param.magnitude_max-param.magnitude_min+1)*2, endpoint=True)
			v = np.linspace(param.magnitude_min,param.magnitude_max, (param.magnitude_max-param.magnitude_min+1), endpoint=True)
			t = map(figures.format_mag, v)
		else:
			vf = np.logspace(np.log10(sl_min), np.log10(sl_max), np.log10(sl_max)-np.log10(sl_min)*2, endpoint=True)
			v = np.logspace(np.log10(sl_min), np.log10(sl_max), np.log10(sl_max)-np.log10(sl_min)+1, endpoint=True)
			t = map(figures.format_log10, v)


		# prepare the surface plot
		xi = np.linspace(-np.pi,np.pi,n_alpha*2)
		yi = np.linspace(-np.pi/2,np.pi/2,n_delta*2)

		# grid the data.

# Sometimes if the points are too well aligned, griddate hangs up as it uses 2d interpolation. Introducing therefore a small random noise on the position.
# or http://stackoverflow.com/questions/10886971/orbit_idernatives-to-scipy-interpolate-griddata-that-dont-hang-on-aligned-points
		ra = ra + np.random.random(ra.shape[0]) * 1e-6 - np.pi
		dec = dec + np.random.random(ra.shape[0]) * 1e-6
		zi = griddata((ra, dec), S_sl, (xi[None,:], yi[:,None]), method='cubic')

#		print np.log10(sl_max)-np.log10(sl_min)

		# if not close enough, mask value
		# cosmetics
		for ii, rag in enumerate(xi) :
			a = np.where(abs(ra-rag)<0.1)[0]
			if np.shape(a) == 0: continue
			# if no a in array, continue
			for jj, decg in enumerate(yi) :
				b = np.where(abs(dec-decg)<0.1)[0]
				# if no dec in array, continue
				if np.shape(b)[0] == 0: continue

				# if not close enough, mask value
				if np.shape(np.intersect1d(a,b))[0] == 0: zi[jj,ii] = np.nan

		if straylight:
		# complete the holes in the interpolation
			cmap = plt.cm.jet
			w = ra_step
			h = dec_step

			for x, y, c_sl in zip(ra, dec, S_sl):
				if magnitudes: cc = (c_sl-param.magnitude_min)/(param.magnitude_max+0.2-param.magnitude_min)
				else: cc = (np.log10(c_sl)-np.log10(sl_min))/(np.log10(sl_max)-np.log10(sl_min))
				ax.add_artist(Rectangle(xy=(x-w/2,y-h/2), color=cmap(cc), width=w, height=h, zorder=0))#,


		# make a plot of the centers of the cells
		scat=plt.plot(ra,dec,'o',c='k', markersize=2)

		if straylight:
			if magnitudes:
				CS = plt.contour(xi,yi,zi,vf,linewidths=0.5,colors='k',extent=extent)
				CS = plt.contourf(xi,yi,zi,vf,cmap=plt.cm.jet,extent=extent)
			else:
				CS = plt.contour(xi,yi,zi,vf,linewidths=0.5,colors='k',extent=extent,locator=ticker.LogLocator())
				CS = plt.contourf(xi,yi,zi,vf,cmap=plt.cm.jet,extent=extent,locator=ticker.LogLocator())

			cbar = plt.colorbar(CS, ticks=v, orientation='horizontal',shrink=.8)
			cbar.set_ticklabels(t)
			l,b,w,h = plt.gca().get_position().bounds
			ll,bb,ww,hh = cbar.ax.get_position().bounds
			cbar.ax.set_position([ll, bb+0.1, ww, hh])

		
			if magnitudes: cbar.set_label(r'$m_V$')
			else: cbar.set_label(r'$\mathrm{Stray\ light\ flux}\ [\frac{\mathrm{ph}}{\mathrm{px}\cdot s}]$')


##########################################################################################
		# Get the outer limit of the data
		# make deep copy to avoid pointers issues from one image to the next.
		grid_points = copy.deepcopy(grid_points0)
		grid_points2 = copy.deepcopy(grid_points0)

		# For each point, find out what are the neighbouring points on the (half) grid
		for ii, ra_ in enumerate(ra):
			dec_ = dec[ii]

			id_ra = find_nearest(ras,ra_+ra_step/2)
			id_dec = find_nearest(decs,dec_+dec_step/2)
			grid_points[id_dec,id_ra] = 1

			id_ra = find_nearest(ras,ra_+ra_step/2)
			id_dec = find_nearest(decs,dec_-dec_step/2)
			grid_points[id_dec,id_ra] = 1

			id_ra = find_nearest(ras,ra_-ra_step/2)
			id_dec = find_nearest(decs,dec_+dec_step/2)
			grid_points[id_dec,id_ra] = 1

			id_ra = find_nearest(ras,ra_-ra_step/2)
			id_dec = find_nearest(decs,dec_-dec_step/2)
			grid_points[id_dec,id_ra] = 1


		if boundaries:
		# Size of "board"
			X = n_alpha-1
			Y = n_delta-1

		# Create the list of all points interesting
			neighbours = lambda x, y : [(x2, y2) for x2 in range(x-1, x+2) for y2 in range(y-1, y+2) if -1 < x <= X and -1 < y <= Y and (x != x2 or y != y2)]
		
		# restrict to only outer points
			for ii in range(0, n_alpha):
				for jj in range(0, n_delta):
					if grid_points[jj, ii]>0. :
						v = [grid_points[neighbours(ii,jj)[kk][1], neighbours(ii,jj)[kk][0]] for kk in range( 0, len( neighbours(ii,jj) ) )]
						# If there is point grid point which is not neighbour to a stray light calculation point, then it is a outer cell point.
						if 0 in v: grid_points2[jj,ii] = 1 

		# restrict the list of the point in the cell grid to only outer points
			ra_grid2 = ra_grid[np.where(grid_points2>0.)]
			dec_grid2 = dec_grid[np.where(grid_points2>0.)]

		# change the name of the variable for clarity
			xx = ra_grid2
			yy = dec_grid2
			del ra_grid2, dec_grid2, grid_points, grid_points2

			ox, oy = draw_boundaries(ax,xx,yy,ra_step,dec_step)	
			if np.shape(ox)[0] > 2 and future: ox, oy = draw_boundaries(ax,ox,oy,ra_step,dec_step,True)

			del xx, yy, ox, oy


		if not future:
			# find the position of the satellite in the orbit and to compute RA, DEC of the Earth and the Sun with respect to the sat's center.

			# EARTH
			id_sat = find_nearest(sat[:,0],minute)		

			x = -1*sat[id_sat,1]
			y = -1*sat[id_sat,2]
			z = -1*sat[id_sat,3]

			ra_sat = rev( right_ascension(x,y) )
			dec_sat= declination(x,y,z)

			circle = Circle((ra_sat-np.pi,dec_sat), 0.07, facecolor='blue',
	                edgecolor='none', alpha=1, zorder=5)
			ax.add_patch(circle)

			# SUN
			id_sat = find_nearest(sun[:,0],minute)
			x = sun[id_sat,1]
			y = sun[id_sat,2]
			z = sun[id_sat,3]

			ra_sun = rev( right_ascension(x,y) )
			dec_sun= declination(x,y,z)

			plt.plot(ra_sun - np.pi,dec_sun,'o',color="yellow", markersize=8, zorder=5)

			# MOON
			id_sat = find_nearest(moon[:,0],minute)
			x = moon[id_sat,1]
			y = moon[id_sat,2]
			z = moon[id_sat,3]

			ra_moon = rev( right_ascension(x,y) )
			dec_moon= declination(x,y,z)
			circle = Circle((ra_moon - np.pi,dec_moon), 0.03, facecolor='white',
	                edgecolor='black', linewidth=0.5, alpha=1, zorder=5)
			ax.add_patch(circle)

			# reduces the size of the arrays speed up the plots
			sat = sat[sat[:,0] >= minute]
			sun = sun[sun[:,0] >= minute]
			moon = moon[moon[:,0] >= minute]
	else:
		v = np.linspace(param.magnitude_min,param.magnitude_max+0.2, (param.magnitude_max-param.magnitude_min+1)*2, endpoint=True)

		scat=plt.scatter(ra,dec,c=S_sl, s=2,vmin=param.magnitude_min, vmax=v[-1])
		

		t = np.linspace(param.magnitude_min,param.magnitude_max, (param.magnitude_max-param.magnitude_min+1), endpoint=True)
		plt.colorbar(ticks=t)
		scat.remove()

	plt.grid(True)

	# add the time, the orbit number and the stray light angle.

	# convert epoch to matplotlib float format
	labels = minute * 60. + const.timestamp_2018_01_01
	# to human readable date
	pre = time.gmtime(labels)
	labels = figures.format_second(pre)

	orbit_current = fast_minute2orbit(times,minute,orbit_id)

	plt.text(-0.1, 1.0,'%s' % labels, transform = ax.transAxes)
	plt.text(-0.1, 0.9,r'$\mathrm{orbit}\ %d$' % orbit_current, transform = ax.transAxes)
#	plt.text(-0.1, 0.8,r'$\mathrm{id}\ %d$' % orbit_id, transform = ax.transAxes)

	if future:
		plt.show()
		exit()
	if not save:
		plt.show()

	if magnitudes: fname = '%s/flux_%07d' % (folder_figures, minute)
	else: fname = '%s/straylight_%07d' % (folder_figures, minute)
	if (fancy and save):
		plt.savefig(fname+'.eps')
		os.system("epstopdf "+fname+".eps")
		os.system('pdfcrop '+fname+'.pdf')
		os.system('mv '+fname+'-crop.pdf '+fname+'.pdf')
		os.system('pdftocairo -png '+fname+'.pdf'+' '+fname)

	if save: plt.savefig(fname+'.png', dpi=param.dpi)

	plt.close()
	del ra, dec, S_sl
	print "Done."

