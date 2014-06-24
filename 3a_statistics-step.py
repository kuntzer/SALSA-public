''' 3-statistics-step.py
=========================
AIM:	Perform basic statistics on the data

INPUT:	files: 	- <orbit_id>_misc/orbits.dat
	variables: see section PARAMETERS (below)

OUTPUT:	in <orbit_id>_misc/ : file one stat file
		in <orbit_id>_figures/ : step distribution, step in function of time

CMD:	python 3-statistics-step.py

ISSUES:	<none known>

REQUIRES:- LATEX, epstopdf, pdfcrop, standard python libraries, specific libraries in resources/
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
import os.path

from resources.routines import *
from resources.TimeStepping import *
import parameters as param
import resources.figures as figures

###########################################################################
### PARAMETERS

# Orbit id
orbit_id = 1001

# First orbit in data set
orbit_ini = 1

# Last orbit to look for
orbit_end = minute2orbit(1440*365+1,orbit_id)

# File name for the output data file
data_file = 'statistics-step.dat'

# Show plots ?
show = True

# Fancy plots ? (ie. eps, pdf, png) otherwise png with less nice fonts
fancy = True

# To compute optimum value, what was the step max value ? (around 10)
max_step_allowed = 10

###########################################################################
### INITIALISATION
# File name fot the computed orbit file
orbits_file = 'orbits.dat'

# Formatted folders definitions
folder_flux, folder_figures, folder_misc = init_folders(orbit_id)

if fancy: figures.set_fancy()

if os.path.isfile(folder_misc+data_file):
	os.remove(folder_misc+data_file)

f = open(folder_misc+data_file,'w')

### Look for the computed orbits
orbits = np.loadtxt(folder_misc+orbits_file)

###########################################################################

############ ORBITS
print >> f, '# ORBITS'

print >> f, 'Number of orbits:', np.shape(orbits[:,1])[0]
print >> f, 'Optimum was:', int(np.floor((orbit_end - orbit_ini) / max_step_allowed))

############ STEPS
print >> f, '# STEPS'
# 1-Average step, mean, min, max
steps = orbits[1:,1]
step_max = np.amax(steps)
print >> f, 'Mean:', np.mean(steps)
print >> f, 'min:', np.amin(steps)
print >> f, 'max:', step_max
print >> f, 'stddev:', np.std(steps)

# 2-Histogramme
bin_edges = np.arange(1,step_max+2)
bins = np.arange(1,step_max+1)
hist, bin_edges = np.histogram(steps,bins=bin_edges)

print >> f, 'bins: ', bins
print >> f, 'histogram: ', hist

fig = plt.figure()
plt.grid(True)

# Bar plot
plt.bar(bins, hist, align='center')

# Cosemetics and labels
plt.xlim([0.5,step_max+0.5])
plt.xticks(range(1, int(step_max+1)))
plt.xlabel(r'$\mathrm{Distribution\ of\ steps\ [orbits]}$')
plt.ylabel(r'$\mathrm{Occurrence}$')

# Saves the figure
fname = '%sdistrib_steps_%d' % (folder_figures,orbit_id)
figures.savefig(fname,fig,fancy)

# 3-Step evolution with time
fig, ax=plt.subplots()
ax.margins(0.05)
plt.subplot(111)

xx = orbits[1:,0]/param.last_orbits[orbit_id]*365.
xx = figures.convert_date(xx)
plt.plot(xx, steps,linewidth=1.5)
fig.autofmt_xdate()

# Cosemetics and labels
plt.grid(True)

plt.ylabel(r'$\mathrm{Step\ [orbits]}$')

# Saves the figure
fname = '%sstep_evolution_%d' % (folder_figures,orbit_id)
figures.savefig(fname,fig,fancy)

f.close()
if show: plt.show()
