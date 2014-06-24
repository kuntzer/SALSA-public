''' 7-plot-means.py
=========================
AIM:	Plots the means of 6-mean_sl.py ! Comparison for an arbitrary number of cases

INPUT:	files: 	- <orbit_id>_misc/orbits.dat
		- <orbit_id>_flux/flux_*.dat
	variables: see section PARAMETERS (below)

OUTPUT:	in all_figures/ : comparison for every orbit_id

CMD:	python 7-plot-means.py

ISSUES:	<None known>

REQUIRES:- standard python libraries, specific libraries in resources/
	 - Structure of the root folder:
	   * <orbit_id>_flux/ --> flux files
	   * <orbit_id>_figures/ --> figures
	   * <orbit_id>_misc/ --> storages of data
	   * all_figures/ --> comparison figures

REMARKS: <none>
'''

###########################################################################
### INCLUDES
import numpy as np
import pylab as plt
import os

from resources.routines import *
from resources.TimeStepping import *
import parameters as param
import resources.figures as figures

from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter

###########################################################################
### PARAMETERS
# Show plots and detailled analysis ?
show = True

# Fancy plots ?
fancy = True

save = False

# what kind of average to plot ? nothing --> average of average
# '_max' average of max
# '_maxdir' maxdirection for each orbit
average = ''

orbit_ids = [1001]
legends = ['Clean']

###########################################################################
### INITIALISATION
# File name fot the computed orbit file
error_file = 'mean%s_sl.dat' % average

if average=='': model=True
else: model=False

if fancy: figures.set_fancy()

###########################################################################
fig=plt.figure()
ax=plt.subplot(111)

maxy = 0.
miny = 0.

for orbit_id, legend in zip(orbit_ids,legends):
	folder_flux, folder_figures, folder_misc = init_folders(orbit_id)
	data = np.loadtxt('%d_misc/%s' % (orbit_id, error_file), delimiter=',')
	dates = data[:,0]/param.last_orbits[orbit_id]*365.
	values= data[:,1]

	if np.amax(values) > maxy: maxy = np.amax(values)
	if np.min(values) < miny: miny = np.amin(values)

	dates = figures.convert_date(dates)
	ax.plot(dates, values,label=legend, linewidth=2)


fig.autofmt_xdate()
plt.ylim([miny*0.95,maxy*1.05])

plt.grid()
ax.xaxis.grid(True,'major',linewidth=2)
ax.yaxis.grid(True,'major',linewidth=2)
plt.legend(loc='upper center')

folder_figures= 'all_figures/'
if average == '' : plt.ylabel(r'$\mathrm{Mean\ stray\ light\ flux\ }\left[\frac{\mathrm{ph}}{\mathrm{px}\cdot\mathrm{s}}\right]$')
if average == '_max' : plt.ylabel(r'$\mathrm{Mean\ maximum\ stray\ light\ flux\ }\left[\frac{\mathrm{ph}}{\mathrm{px}\cdot\mathrm{s}}\right]$')
if average == '_maxdir' : plt.ylabel(r'$\mathrm{Flux\ in\ worst\ direction}\left[\frac{\mathrm{ph}}{\mathrm{px}\cdot\mathrm{s}}\right]$')

###########################################################################
###########################################################################
###########################################################################

### Loading the value of n(alpha)
values = np.loadtxt('704_misc/angle_usage_704_35_1-5322.dat')
xo = values[:-1,0]
yo = values[:-1,1]
from scipy import stats
slope, intercept, r_value, p_value, std_err = stats.linregress(xo,yo)
s = slope
c1 = intercept

def ni(alpha,s=s,c1=c1):
	return s * alpha + c1
n=np.vectorize(ni)

fig=plt.figure()
ax=plt.subplot(111)

maxy = 0.
miny = 0.

n_sampling = 55
angles_s = np.linspace(35,89,n_sampling)

from scipy.optimize import leastsq

orbit_ini = [1,441,891,1331,1771,2221,2661,3111,3551,3991,4441,4881]
orbit_end = [441,891,1331,1771,2221,2661,3111,3551,3991,4441,4881,5322]

def residuals(p, y, x):
	s,pw,c,s2,pw2 = p
	err = y - (s*np.power(x,pw) + s2*np.power(x,pw2) + c)
	return err

def peval(x, p):
	return p[0]*np.power(x,p[1])+p[2] + p[3]*np.power(x,p[4]) 

p0 = [0.2, 1.0, -2.,0.,1.]
'''
for oi, of in zip(orbit_ini, orbit_end):
	values = np.loadtxt('704_misc/angle_usage_704_35_%d-%d.dat' % (oi, of))

	xo = values[:-1,0]
	yo = values[:-1,1]

	plsq = leastsq(residuals, p0, args=(yo, xo), maxfev=2500)
'''
for orbit_id, legend in zip(orbit_ids,legends):
	### Loading the PST
	values = np.loadtxt('pst/%d.dat' % orbit_id)
	angles = values[:,0]
	pst_val= values[:,1]

	yinterp = np.interp(angles_s, angles, pst_val)

	n_pst = n(angles_s)

	n_pst = yinterp * n_pst
	tot = n_pst.sum()
	print tot, legend

	folder_flux, folder_figures, folder_misc = init_folders(orbit_id)
	data = np.loadtxt('%d_misc/%s' % (orbit_id, error_file), delimiter=',')
	dates = data[:,0]/param.last_orbits[orbit_id]*365.
	values= data[:,1]/tot

	if np.amax(values) > maxy: maxy = np.amax(values)
	if np.min(values) < miny: miny = np.amin(values)

	dates = figures.convert_date(dates)
	ax.semilogy(dates, values,label=legend, linewidth=2)


fig.autofmt_xdate()
#plt.ylim([miny*0.95,maxy*1.05])

plt.grid()
ax.xaxis.grid(True,'major',linewidth=2)
ax.yaxis.grid(True,'major',linewidth=2)
plt.legend(loc='upper center')
plt.title('divided by n*PST')

folder_figures= 'all_figures/'
if average == '' : plt.ylabel(r'$\mathrm{Mean\ stray\ light\ flux\ }\left[\frac{\mathrm{ph}}{\mathrm{px}\cdot\mathrm{s}}\right]$')
if average == '_max' : plt.ylabel(r'$\mathrm{Mean\ maximum\ stray\ light\ flux\ }\left[\frac{\mathrm{ph}}{\mathrm{px}\cdot\mathrm{s}}\right]$')
if average == '_maxdir' : plt.ylabel(r'$\mathrm{Flux\ in\ worst\ direction}\left[\frac{\mathrm{ph}}{\mathrm{px}\cdot\mathrm{s}}\right]$')

###########################################################################
###########################################################################
###########################################################################

if show: plt.show()

# Saves the figure
fname = '%sPST_mean_fluxes' % (folder_figures)
fname = '%smean%s_fluxes' % (folder_figures,average)

fname = '%smean%s_fluxes' % (folder_figures,average)
if save: 
	figures.savefig(fname,fig,fancy)
	print 'saved as', fname
exit()

