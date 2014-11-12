import numpy as np
import pylab as plt
import os

import resources.figures as figures

###############################################################################

# A list of dictionaries containing all data to plot.
# Mandatory keys:
# 'folder', 'fname'
# Optional keys:
# "label",'lw', 'color'
skycoverage=[
	{'folder':"800_35_AKTAR_100x50_misc",
	 'fname':"skycoverage_81min_V12.0_SAA_0.0pst_0.996SLreduction_cumul_.txt",
	 'color':"k","label":"Nominal case",'lw':2},
	{'folder':"800_35_AKTAR_100x50_misc",
	 'fname':"skycoverage_81min_V12.0_0.0pst_0.996SLreduction_cumul_.txt",
	 'color':"b","label":"No SAA",'lw':2}
	]

# What is the target ? If None, nothing is shown
requirement_days=13
requirement_obs=25

# Force the axis, if set to None, automatic from pylab
xlim=[0,120]
ylim=[0,50]

# Filename for the figure, if None, the figure is not saved. Without extension!
out_fname=None#"800_35_AKTAR_100x50_figures/skycoverage_stat_81min"

# Nicer plot + exportation in pdf, eps
fancy=True

###############################################################################

if fancy:figures.set_fancy()
legend=False
f=plt.figure()
for el in skycoverage:
	data=np.loadtxt(os.path.join(el['folder'],el['fname']))
	if 'lw' in el:
		lw=el['lw']
	else:
		lw=2
	if 'label' in el:
		label=el['label']
		legend=True
	else:
		label=None

	line,=plt.plot(data[:,0],data[:,1],lw=lw,label=label)
	if 'color' in el:
		line.set_color=el['color']


if legend:
	plt.legend(loc='best')

if not requirement_days is None:
	plt.axvline(x=requirement_days, c="gray", lw=2)
if not requirement_obs is None:
	plt.axhline(y=requirement_obs, c="gray", lw=2)

plt.ylabel(r"$\mathrm{Fraction\ of\ full\ sky\ [\%]}$")
plt.xlabel(r"$\mathrm{Required\ observation\ time\ [days]}$")

if not xlim is None:
	plt.xlim(xlim)
if not ylim is None:
	plt.ylim(ylim)
plt.grid()

if not out_fname is None:
	figures.savefig(out_fname, f, fancy)
	print 'saved as %s' % out_fname
plt.show()
