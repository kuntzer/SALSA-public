''' TimeStepping.py
=========================
AIM:	Library to provide function when computing the step size and start/end times of orbits

INPUT:	function depend

OUTPUT:	function depend

CMD:	To include: from resources.targets import *

ISSUES:	<none known>

REQUIRES: standard python libraries, specific libraries in resources/

REMARKS: <none>
'''


def DeltaMagnitude(value,ref):
	''' Returns the relative difference in magnitudes '''
	import numpy as np
	import resources.routines as rr
	import parameters as param

	value = rr.flux2mag(value,param.ppm_threshold)
	ref = rr.flux2mag(ref,param.ppm_threshold)
	
	diff = (value - ref)/ref

	return np.abs(diff)

def orbit2times(orbit,orbit_id):
	''' Returns starting and ending times of the orbit in the format
		0, orbit lenght, absolute time of start, abs. time of finish
	    This function requires no reference and loads the reference by itself --> slow
	'''
	import numpy as np

	times = np.loadtxt('resources/minute_table_%s.dat' % orbit_id, delimiter=',',dtype='Int32')
	return fast_orbit2times(times,orbit,orbit_id)

def fast_orbit2day(times, orbit, orbit_id):
	''' converts day in 2018 to usual date '''
	import datetime
	import numpy as np

	a_ini = fast_orbit2a_ini_vect(times, orbit, orbit_id)

	return np.round(a_ini/60./24.)

def fast_orbit2a_ini_vect(times, orbit, orbit_id):
	''' converts day in 2018 to usual date '''
	import datetime
	import numpy as np

	ids = np.in1d(times[:,0], orbit)
	times = times[ids,:]
	a_ini = times[::2,2]

	return a_ini

def fast_orbit2a_end_vect(times, orbit, orbit_id):
	''' converts day in 2018 to usual date '''
	import datetime
	import numpy as np

	ids = np.in1d(times[:,0], orbit)
	times = times[ids,:]
	a_end = times[1::2,2]
	return a_end

def fast_orbit2times(times,orbit,orbit_id):
	import numpy as np

	times = times[times[:,0] == orbit]
	t_ini = times[0,1]
	t_end = times[-1,1]

	a_ini = times[0,2]
	a_end = times[-1,2]
	del times

	return t_ini, t_end, a_ini, a_end

def minute2orbit(minute,orbit_id):
	from numpy import loadtxt

	times = loadtxt('resources/minute_table_'+str(orbit_id)+'.dat', delimiter=',',dtype='Int32')

	return fast_minute2orbit(times, minute,orbit_id)

def fast_minute2orbit(times, minute,orbit_id):
	from resources.routines import find_nearest
	id_orbit = find_nearest(times[:,2],minute)
	orbit = times[id_orbit,0]

	return orbit

def print_rule():
	''' <NOT USED> '''
	print '\t\t\t\t1',
	for i in range(1,11):
		print '      ',
		print i*10,
	print

def print_flags():
	''' <NOT USED> '''
	print '\tFlags used in the representations of the dataset:'
	print '\t1\tArrays were compared, difference < p'
	print '\t0\tArrays were compared, equal or all values masked'
	print '\ts\tNo data (SAA?)'
	print '\te\tTried second chance and now error < p'
	print '\tp\tValue is > p'
	print


def LoadReferences(minute,t_end,ras,decs,n_alpha=None,n_delta=None,file_flux='flux_',folder='flux'):
	''' Prepares the reference orbit '''
	import numpy as np
	import resources.routines as rr

	import parameters as param

	# Prepare the grid just like for the observability maps.
	if n_alpha == None: n_alpha = param.resx
	if n_delta == None: n_delta = param.resy

	Period = t_end - minute+1

	t_ini=minute
	references = np.zeros(n_alpha*n_delta*(Period)).reshape(Period,n_alpha,n_delta)
	while (minute <= t_end):
		try:
			ra, dec, S_sl = rr.load_flux_file(minute, file_flux, folder=folder)
			for i, target in enumerate(S_sl):
				id_ra = rr.find_nearest(ras,ra[i])
				id_dec = rr.find_nearest(decs,dec[i])
				references[minute-int(t_ini), id_ra, id_dec] = S_sl[i]

		except IOError:
			
			# In case of SAA in the file.
#			print 'No file for this reference at minute '+str(minute)
			try:
				ra, dec, S_sl = rr.load_flux_file(minute-Period, file_flux, folder=folder)
				for i, target in enumerate(S_sl):
					id_ra = rr.find_nearest(ras,ra[i])
					id_dec = rr.find_nearest(decs,dec[i])
					references[minute-int(t_ini), id_ra, id_dec] = S_sl[i]
			except:
				try:
					ra, dec, S_sl = rr.load_flux_file(minute+Period, file_flux, folder=folder)
					for i, target in enumerate(S_sl):
						id_ra = rr.find_nearest(ras,ra[i])
						id_dec = rr.find_nearest(decs,dec[i])
						references[minute-int(t_ini), id_ra, id_dec] = S_sl[i]
				except:
#					print 'Assuming 0'
					pass
		minute += 1
	return references





def find_best_fit(time,value,references):
	''' <NOT USED> '''
	import numpy as np
	min = 10
	id_min=0
	nb_masked_min = 10000

	for minute in range(time-3,time+3):
		try:
			ref=references[minute,:,:]
		except:
			return -1, .5

	# Mask all targets that are not in each array.
		nb_masked = np.shape(np.where(value-ref==value))
		value[np.where(value-ref==value)] = 0
		ref[np.where(value==0)] = 0

	# Proceed with computation of mag, then their relative difference
		diff = DeltaMagnitude(value,ref)
	if min > np.amin(diff[np.where(diff>0)]):	
		min = np.amin(diff[np.where(diff>0)])
		id_min = minute

	return id_min, min



