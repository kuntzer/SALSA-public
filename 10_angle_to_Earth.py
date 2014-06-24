''' 10-angle_to_earth.py
=========================
AIM:	Similar to 9-plot_flux.py, but shows angle to Earth limb instead of sl flux

INPUT:	files: 	- <orbit_id>_misc/orbits.dat
		- <orbit_id>_flux/flux_*.dat
		- resources/moon_*.dat, sun_*.dat, orbits_*.dat
	variables: see section PARAMETERS (below)

OUTPUT:	in <orbit_id>_figures/maps/ : map with the following name: dist_earth_%07d.png

CMD:	python 10-angle_to_earth.py

ISSUES:	<none known>

REQUIRES:- standard python libraries, specific libraries in resources/
	 - Structure of the root folder:
	   * <orbit_id>_flux/ --> flux files
	   * <orbit_id>_figures/maps/ --> figures
	   * <orbit_id>_misc/ --> storages of data
	   * all_figures/ --> comparison figures

REMARKS: <none>
'''

######################################################################
import numpy as np
import pylab as plt
import time

from resources.constants import *
from resources.routines import *
from resources.TimeStepping import *
import parameters as param
import matplotlib.cm as cm
import resources.figures as figures

from matplotlib.ticker import MaxNLocator

#p1 = np.array([-86.67, 36.12])
#p2 = np.array([-118.40, 33.94])
#p1 *= np.pi / 180.
#p2 *= np.pi / 180.
#print vSphericalDistance(p1[0],p1[1],p2[0],p2[1])
######################################################################
# orbit_id
orbit_id = 1000

# Line of sight (LOS) to Earth's limb angle
sl_angle = 35
 
file_orbit = 'orbit_'+str(orbit_id)+'.dat'
minute_ini = 4036
minute_ini = 0
minute_end = 4136
minute_end = 10

n_alpha = param.resx
n_delta = param.resy

fancy = False
######################################################################
# Initialisation
file_flux = 'flux_'
# Formatted folders definitions
folder_flux, folder_figures, folder_misc = init_folders(orbit_id)
folder_figures= '%d_figures/maps/' % (orbit_id)

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino'],'size':14})
rc('text', usetex=True)

ra_i = 0
ra_f = 2.*np.pi

dec_i = -np.pi/2.
dec_f = np.pi/2.

ra_step = (ra_f-ra_i)/n_alpha
dec_step = (dec_f-dec_i)/n_delta

iterable = (ra_i + ra_step/2 + i*ra_step for i in range(n_alpha))
ras = np.fromiter(iterable, np.float)

iterable = (dec_i + dec_step/2 + i*dec_step for i in range(n_delta))
decs = np.fromiter(iterable, np.float)



distance = np.zeros([n_delta,n_alpha])
ra_grid, dec_grid = np.meshgrid(ras,decs)
######################################################################

sys.stdout.write("Loading orbit file...\t\t\t")
sys.stdout.flush()
try:
	sat = np.loadtxt('resources/'+file_orbit, delimiter='\t')
except ValueError:
	sat = np.loadtxt('resources/'+file_orbit, delimiter=' ')
print "Done."

minute = minute_ini
while ( minute < minute_end+1 ):
	sys.stdout.write("Ploting minute "+str(minute)+'...\t\t')
	sys.stdout.flush()

	ra, dec, S_sl = load_flux_file(minute, file_flux, folder=folder_flux)

	id_sat = find_nearest(sat[:,0],minute)
	x = -1*sat[id_sat,1]
	y = -1*sat[id_sat,2]
	z = -1*sat[id_sat,3]

	ra_sat = rev( right_ascension(x,y) )
	dec_sat= declination(x,y,z)

	points = np.empty([np.shape(ra)[0],2])
	points[:,0] = vrev(ra)
	points[:,1] = dec

	r = R_3d(x,y,z)
	d_earthsat = R_3d(x,y,z)
	psi0 = np.arcsin((atmosphere + R_Earth)/1e5/d_earthsat)
	psi = sl_angle / 180. * np.pi

	for ii in range(0, n_alpha):
		for jj in range(0, n_delta):
			distance[jj, ii] = vSphericalDistance(ra_sat,dec_sat,ras[ii],decs[jj])


	alpha= np.arcsin((R_Earth+atmosphere)/r/1e5)
	beta = np.pi/2. - alpha
	gamma= np.abs(dec_sat) + beta
	distance = (distance-beta)*180/np.pi
	distance_t = (vSphericalDistance(ra_sat,dec_sat,ra,dec) - beta)*180/np.pi
	LIMIT_e = limit_Earth(gamma,ra_sat,dec_sat)

	beta = np.pi/2. - alpha + psi
	gamma= np.abs(dec_sat) + beta
	LIMIT = limit_Earth(gamma,ra_sat,dec_sat)

	plt.figure(dpi=param.dpi)
	ax = plt.subplot(111)

	if sl_angle == 35 : vmax = sl_angle+120
	else : vmax = sl_angle+145
	step = 13
	v = np.linspace(sl_angle, vmax,step, endpoint=True)

	CS=plt.contour(ra_grid*RAD,dec_grid*RAD,distance,v,colors='k')
	plt.contourf(ra_grid*RAD,dec_grid*RAD,distance,v,cmap=cm.RdYlGn)
	plt.clabel(CS, inline=1, fontsize=10)
	cbar = plt.colorbar()
	cbar.set_label(r'$\theta\ \mathrm{Angular}\ \mathrm{distance}\ \mathrm{to}\ \mathrm{limb}\ \mathrm{[deg]}$')

	plt.scatter(LIMIT_e[:,0]*RAD,(LIMIT_e[:,1])*RAD,color="blue", s=8, edgecolor='none')
	plt.scatter(LIMIT[:,0]*RAD,(LIMIT[:,1])*RAD,color="red", s=8, edgecolor='none')
	plt.scatter(ra*RAD,dec*RAD,c='k')

	plt.scatter(ra_sat*RAD,dec_sat*RAD,s=80)
	plt.grid(True)
	plt.xlabel(r'$\alpha$')
	plt.ylabel(r'$\delta$')

	plt.xlim([0, 360])
	plt.ylim([-90,90])

	plt.gca().xaxis.set_major_locator( MaxNLocator(nbins = 13) )
	plt.gca().yaxis.set_major_locator( MaxNLocator(nbins = 7) )

	current_date = minute * 60. + timestamp_2018_01_01
	current_date = figures.format_second(time.gmtime(current_date))

	plt.text(0.01, 0.96,current_date, transform = ax.transAxes)

	fname = '%sdist_earthlimb_%07d' % (folder_figures, minute)
	plt.savefig(fname+'.png', dpi=param.dpi)
	if (fancy):
		plt.savefig(fname+'.eps')
		os.system("epstopdf "+fname+".eps")
		os.system('pdfcrop '+fname+'.pdf')
		os.system('mv '+fname+'-crop.pdf '+fname+'.pdf')
		os.system('pdftocairo -png '+fname+'.pdf'+' '+fname)

	minute += 1
	plt.close()
	del LIMIT, ra, dec, S_sl, points, LIMIT_e
	print 'Done.'

print 'finished', minute, minute_end
