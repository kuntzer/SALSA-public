''' figures.py
=========================
AIM:	Provide several specific functions to save beautiful figures

INPUT:	function depend

OUTPUT:	function depend

CMD:	To include: import resources.figures as figures

ISSUES:	<none known>

REQUIRES: standard python libraries, specific libraries in resources/

REMARKS: in general fancy means latex interpreter (font is serif, Palatino) and generates *.eps and *.pdf
'''
######################################################################

import numpy as np

def savefig(fname,fig,fancy=False):
	import os
	import subprocess

	directory=os.path.dirname(fname)
	if not os.path.exists(directory):
		os.makedirs(directory)

	fig.savefig(fname+'.png',dpi=300)

	if fancy: 
		fig.savefig(fname+'.eps',transparent=True)
		os.system("epstopdf "+fname+".eps")
		command = 'pdfcrop %s.pdf' % fname
		subprocess.check_output(command, shell=True)
		os.system('mv '+fname+'-crop.pdf '+fname+'.pdf')
	

def set_fancy():
	from matplotlib import rc
	rc('font',**{'family':'serif','serif':['Palatino'],'size':16})
	rc('text', usetex=True)

def cd(xx, year=2018, day=1, month=1):
	''' converts day in 2018 to usual date '''
	import datetime
	dd = datetime.date.today()
	dd = dd.replace(year=year, month=month, day=day)
	first_ordinal = dd.toordinal()
	return datetime.date.fromordinal(first_ordinal+int(round(xx)))
convert_date = np.vectorize(cd)

def format_log10(value):
	return r'$10^{%d}$' % np.log10(value)

def format_mag(value):
	return r'$%d$' % value

def format_degree(value):
	return r'$%d^\circ$' % value

def format_second(xx):
	import time
	return time.strftime('%d %b %H:%M', xx)

def format_day(xx):
	import time
	return time.strftime('%d %b', xx)
