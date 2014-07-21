''' routines.py
=========================
AIM:	Provide several specific functions for the other scripts

INPUT:	function depend

OUTPUT:	function depend

CMD:	To include: from resources.routines import *

ISSUES:	- sort_boundaries() does not work well when the observation zone spans high values of ra and small (ie there are two zones: [ra_min,360] + [0, ra_max])

REQUIRES: standard python libraries, specific libraries in resources/

REMARKS: Beware! If there are changes the version in analysis/resources must be the most uptodate.
'''
######################################################################
import numpy as np
import subprocess
import csv
import sys
import time

def init_folders(orbit_id):
	folder_flux = '%s_flux/' % (orbit_id)
	folder_figures= '%s_figures/' % (orbit_id)
	folder_misc	= '%s_misc/' % (orbit_id)
	return folder_flux, folder_figures, folder_misc

def find_nearest(array,value):
    ''' Find nearest value is an array '''
    idx = (np.abs(array-value)).argmin()
    return idx

def JulianDay(date):
	'''
	Convert a date with the format as below to JD
	format : array with:
	0->yr
	1->month
	2->day
	3->hour
	4->min
	5->sec
	'''
	fday = date[3] / 24.00 + date[4] / 1440.00 + date[5] / 86400. 
	if (date[1] <= 2.):
		year = date[0] - 1.
		m = date[1] + 12.
	else:
		year = date[0]
		m = date[1]

	A = int(year/100)
	B = 2 - A + int(A/4)
	JD= int(365.25*year) + int(30.6001*(m+1)) + date[2] + fday + 1720994.5 + B
	return JD

def minutes_from_2018(JD):
	import resources.constants as const
	return (JD - const.JD_2018) * 24. * 60.

def minutes_2_JD(minutes):
	import resources.constants as const
	return minutes / 60. / 24. + const.JD_2018

def CartesianCoordinates(a,d,r):
	from numpy import cos, sin, zeros
	pos = zeros(3)
	pos[0] = r * cos(a)*cos(d)
	pos[1] = r * sin(a)*cos(d)
	pos[2] = r * sin(d)
	return pos

def declination(x,y,z):
	return np.arctan2(z, np.sqrt(x*x+y*y) )

def right_ascension(x,y):
	return np.arctan2(y, x)

def R_3d(x,y,z):
	return np.sqrt(x*x+y*y+z*z)

def rev(angle):
	twopi = 2.*np.pi
	angle = angle - int(angle/twopi)*twopi
	if (angle < 0):
		angle = angle + twopi
	return angle
''' Vectorisation of the function (for numpy.arrays) See http://docs.scipy.org/doc/numpy/reference/generated/numpy.vectorize.html '''
vrev = np.vectorize(rev)


def flux2mag(noise, threshold=None):
	''' Returns the maximum visible magnitude for a given flux ("noise") and with a threshold in ppm.
	UNITS : in:  [ph/(s px)]
		out: [mag]
	REMARS: if no specific thershold the default one is defined in parameters
	'''
	import parameters as param
	import resources.constants as const

	if threshold == None : threshold = param.ppm_threshold

	threshold *= 1e-6

	A = const.Jy * const.Fv * param.dlambda_over_lambda
	B = (param.radius_tele * param.radius_tele) / (param.radius_psf * param.radius_psf)

	mag = -2.5 * np.log10(noise.clip(min=1e-40)  / (A * B * threshold))

	mag[np.where(mag>param.magnitude_max)]=param.magnitude_max

	return mag

def mag2flux(mag, threshold=None):
	''' Returns the stray light to got a certain magnitude with a given noise("noise") and with a threshold in ppm.
 UNITS : in:  [mag]
	  out: [ph/(s px)]
	'''
	import parameters as param
	import resources.constants as const

	if threshold == None : threshold = param.ppm_threshold

	threshold *= 1e-6

	A = const.Jy * const.Fv * param.dlambda_over_lambda
	B = (param.radius_tele * param.radius_tele) / (param.radius_psf * param.radius_psf)

	noise = A * B * 10.**(mag*-0.4)

	return noise * threshold

def fast_SAA(SAA, minute):
	''' Returns whether the satellite is in the SAA (True) or not (False)
		in numpy SAA array, minute. Requires SAA data. '''
	lower = SAA[SAA[:,0] <= minute]
	upper = SAA[SAA[:,0] >= minute]

	try:
		minute_SAA = lower[-1,0]
		lower = lower[-1,1]
		if lower == 0 and minute_SAA == minute: return False
		upper = upper[0,1]
	except IndexError: return False
	if lower == 0 and upper == 1 : return False
	else: return True


def altitude2period(apogee,perigee):
	import resources.constants as const
	import parameters as param

	a = const.R_Earth / 100. + (apogee+perigee)/2. * 1e3
	return 2.*np.pi * np.sqrt(a*a*a/const.mu_Earth)  / 60.

def load_map(orbit_id,folder='absolute_map'):
	import parameters as param
	import os
	return np.loadtxt(os.path.join(folder,"%s%s.dat" % (param.file_map,orbit_id)))	

def slice_map(map_obstot, minute, threshold=1e-5):
	return map_obstot[abs(map_obstot[:,0]-minute) < threshold]

def load_flux_file(minute,fname,folder='flux'):
	''' Loads flux files as output by Fortran, making sure the data is retrieved correctly'''
	ra_sl = list()
	dec_sl= list()
	observable = list()

	with open(folder+"/"+fname+str(minute)+".dat", "rb") as observability_file:
		observability_data = csv.reader(observability_file, delimiter="\t")
		for row in observability_data:
		# if line is empty, skip otherwise filter out the blank
			if len(row) > 0:
				line = row[0].split(' ')
				ra_sl.append(float(filter(None, line)[0]))
				dec_sl.append(float(filter(None, line)[1]))
				observable.append(float(filter(None, line)[2]))
	S_sl = np.asarray(observable)
	ra = np.asarray(ra_sl)
	dec= np.asarray(dec_sl)

	return np.asarray(ra_sl), np.asarray(dec_sl), S_sl

def fast_load_flux_file(minute,fname,folder='flux'):
	''' <NOT USED> '''
	ra_sl = list()
	dec_sl= list()
	observable = list()

	with open(folder+"/"+fname+str(minute)+".dat", "rb") as observability_file:
		observability_data = csv.reader(observability_file, delimiter="\t")
		for row in observability_data:
		# if line is empty, skip otherwise filter out the blank
			if len(row) > 0:
				line = row[0].split(' ')
				ra_sl.append(float(filter(None, line)[0]))
				dec_sl.append(float(filter(None, line)[1]))
				observable.append(float(filter(None, line)[2]))
	S_sl = np.asarray(observable)
	ra = np.asarray(ra_sl)
	dec= np.asarray(dec_sl)

	return ra_sl, dec_sl, S_sl

def ROT_MATRIX_X(theta, x, y, z):
	import numpy as np
#     For the inverse matrix replace theta by -theta

	xx = x
	yy = y
	zz = z

	x = xx
	y = yy*np.cos(theta) - zz*np.sin(theta)
	z = yy*np.sin(theta) + zz*np.cos(theta)

	return x,y,z

def ROT_MATRIX_Y(theta, x, y, z):
	import numpy as np
#     For the inverse matrix replace theta by -theta
	xx = x
	yy = y
	zz = z

	x = xx*np.cos(theta) + zz*np.sin(theta)
	y = yy
	z = -xx*np.sin(theta) + zz*np.cos(theta)

	return x,y,z


def ROT_MATRIX_Z(theta, x, y, z):
	import numpy as np
#     For the inverse matrix replace theta by -theta
	xx = x
	yy = y
	zz = z
	   
	x = xx*np.cos(theta) - yy*np.sin(theta)
	y = xx*np.sin(theta) + yy*np.cos(theta)
	z = zz

	return x,y,z

def limit_Earth(gamma,RA_SAT,DEC_SAT,n=100):
	import numpy as np
	def sign(x,y):
		if (x > 0 and y > 0) or (x < 0 and y < 0): return 1
		else: return -1

	LIMIT = np.zeros([n,2])

	x = np.cos(gamma)*np.cos(RA_SAT)
	y = np.cos(gamma)*np.sin(RA_SAT)
	z = sign(1.,DEC_SAT)*np.sin(gamma)

	theta = -RA_SAT
	x, y, z = ROT_MATRIX_Z(theta, x, y, z)

	theta = DEC_SAT
	x, y, z = ROT_MATRIX_Y(theta, x, y, z)

	xx = x
	yy = y
	zz = z

	for i in range(1, n+1):
#		print i, x, y, z
		x = xx
		y = yy
		z = zz

		theta = float(i-1)*2.*np.pi/float(n-1)

		x,y,z = ROT_MATRIX_X(theta, x, y, z) 
#		print i, x, y, z
		theta = -DEC_SAT
		x,y,z = ROT_MATRIX_Y(theta, x, y, z)
#		print i, x, y, z
		theta = RA_SAT
		x,y,z = ROT_MATRIX_Z(theta, x, y, z)
#		print i, x, y, z
		LIMIT[i-1,0]= rev(np.arctan2(y,x))
		LIMIT[i-1,1]= np.arctan2(z, np.sqrt(x*x+y*y))
	return LIMIT

def limit_theta2(RA_SAT,DEC_SAT,angle,n=300):
	import constants as cst
	import numpy as np

	theta = np.linspace(0,2*np.pi,n)
	LIMIT = np.zeros([n,2])

	LIMIT[:,0] = vrev((np.cos(theta))*angle+RA_SAT)
	LIMIT[:,1] = np.sin(theta)*angle+DEC_SAT

	import pylab as plt
#	plt.figure()
#	plt.plot(LIMIT[:,0],LIMIT[:,1])
#	plt.show()

	return LIMIT

def SphericalDistance(ra1,dec1,ra2,dec2):
	'''Compute the distance on the great circle. Computations are based upon http://en.wikipedia.org/wiki/Great-circle_distance#Formulas'''
	from numpy import sin, cos, arctan2, arccos
	ls = ra1
	lf = ra2

	ps = dec1
	pf = dec2

	dl = lf-ls
	dp = pf-ps

	num1 = cos(pf)*sin(dl)
	num2 = (cos(ps)*sin(pf) - sin(ps)*cos(pf)*cos(dl))

	den1 = sin(ps)*sin(pf)
	den2 = cos(ps)*cos(pf)*cos(dl)

	num = np.sqrt(num1 * num1 + num2*num2)
	den = den1 + den2

	return arctan2(num,den)
vSphericalDistance = np.vectorize(SphericalDistance)

def rev_d(angle):
	import numpy as np
	if angle < -np.pi/2: angle += np.pi
	elif angle > np.pi/2: angle -= np.pi
	return angle
vrev_d = np.vectorize(rev_d)


def compare_two_orbits(orbit_ref, orbit_val, orbit_id, p=0.05, file_flux='flux_', folder='flux',shift=0, return_max=False, return_direction=False):
	'''compare_two_orbits returns the maximum difference between them in relative terms. Unit : - (ie. 0.05 means 5%)'''
	import numpy as np
	import TimeStepping as ts
	import parameters as param

	# Prepare the grid just like for the observability maps.
	resx = param.resx
	resy = param.resy

	ra_i = 0
	ra_f = 2.*np.pi

	dec_i = -np.pi/2.
	dec_f = np.pi/2.

	ra_step = (ra_f-ra_i)/resx
	dec_step = (dec_f-dec_i)/resy

	iterable = (ra_i + ra_step/2 + i*ra_step for i in range(resx))
	ras = np.fromiter(iterable, np.float)

	iterable = (dec_i + dec_step/2 + i*dec_step for i in range(resy))
	decs = np.fromiter(iterable, np.float)
	
	t_ini2, t_end2, a_ini2, a_end2 = ts.orbit2times(orbit_ref,orbit_id)
	references=ts.LoadReferences(a_ini2, a_end2,ras,decs,file_flux=file_flux,folder=folder)
	if t_end2 == 99: shift+=1

	t_ini, t_end, a_ini, a_end = ts.orbit2times(orbit_val,orbit_id)
	error_max = 0
	second_chance= False
	count_error  = True

	t_ini = np.ceil(t_ini) # not sure it's useful.
	minute = t_ini

	sl_max = 0

	###########################################################################
	# Iterate on every time.
	while (minute <= t_end):
		# initialise the array for the current minute (ie. make sure nothing is left in the table from last time step.

		value = np.zeros([resx,resy])

		try:		
			ref = references[minute,:,:]
		except IndexError:
			minute = t_end+1
			continue

		try:
		# Try to load the fluxes for a given minute (minute goes from 0 to period whereas a_ini is the absolute time from 0:0:0.0 1/1/2018 in min
			ra, dec, S_sl = load_flux_file(int(minute+a_ini+shift), file_flux, folder=folder)

		except IOError:
		# If not found it means that the satellite is (propably) over the SAA. Skip this step.
			minute +=1
			continue	

		# Maps the data to the full grid.
		for i, target in enumerate(S_sl):
			id_ra = find_nearest(ras,ra[i])
			id_dec = find_nearest(decs,dec[i])
			value[id_ra, id_dec] = S_sl[i]

		if np.amax(value) > sl_max: 
			sl_max = np.amax(value)

		# Mask all targets that are not in each array.
		value[np.where(value-ref==value)] = 0
		ref[np.where(value==0)] = 0

		# Proceed with computation of mag, then their relative difference
		diff = ts.DeltaMagnitude(value,ref)

		# Select all points that are on both arrays and if we output the flags, do it ;)
		tostats = diff[np.where(diff>0)]

		err = np.amax(diff)
		# Verify that no more than number_illegal points are different than p % of the reference
		if (diff[np.where(diff>p)]).shape[0] > 0:
			# Well, that's illegal. Is it a problem of indicies ?
			# We give a second chance to the code :
			if not (second_chance):
				second_chance = True
				minute -= 1
				shift += 1
			else:
				second_chance = False
				count_error = False
				# print the flag i for illegal as we already had a second chance.
				shift -= 1
		else:
			# Do not take into account if there are values > p and their # is less then the number_illegal.
			diff[np.where(diff>p)]=0
			second_chance = False

		# Record the maximum error, but if there is another run, dismiss.
		if error_max < np.amax(diff) and count_error : 
			error_max = err
		del value
		if not count_error : count_error = True
		minute += 1

	###########################################################################
	# Treats the whole orbit

	if return_max: return error_max, sl_max
	else:	return error_max

def find_direction_flux(orbit, orbit_id, find='max', file_flux='flux_', folder='flux'):
	'''Find out the direction in which the flux is the maximum or minimum for a given orbit'''
	import numpy as np
	import TimeStepping as ts
	import parameters as param

	# Prepare the grid just like for the observability maps.
	resx = param.resx
	resy = param.resy

	ra_i = 0
	ra_f = 2.*np.pi

	dec_i = -np.pi/2.
	dec_f = np.pi/2.

	ra_step = (ra_f-ra_i)/resx
	dec_step = (dec_f-dec_i)/resy

	iterable = (ra_i + ra_step/2 + i*ra_step for i in range(resx))
	ras = np.fromiter(iterable, np.float)

	iterable = (dec_i + dec_step/2 + i*dec_step for i in range(resy))
	decs = np.fromiter(iterable, np.float)

	t_ini, t_end, a_ini, a_end = ts.orbit2times(orbit,orbit_id)

	t_ini = np.ceil(t_ini) # not sure it's useful.
	minute = t_ini

	if find == 'max' : sl = 0
	else: sl = 1e9
	
	dd = np.nan
	
	###########################################################################
	# Iterate on every time.
	while (minute <= t_end):
		# initialise the array for the current minute (ie. make sure nothing is left in the table from last time step.

#		value = np.zeros(resx*resy)

		try:
		# Try to load the fluxes for a given minute (minute goes from 0 to period whereas a_ini is the absolute time from 0:0:0.0 1/1/2018 in min
			ra, dec, S_sl = load_flux_file(int(minute+a_ini), file_flux, folder=folder)

		except IOError:
		# If not found it means that the satellite is over the SAA. Skip this step.
			minute +=1
			continue	

		# Maps the data to the full grid.
#		for i, target in enumerate(S_sl):
#			id_ra = find_nearest(ras,ra[i])
#			id_dec = find_nearest(decs,dec[i])
#			value[id_ra, id_dec] = S_sl[i]

		if find == 'max' and np.amax(S_sl) > sl: 
			sl = np.amax(S_sl)
			id_max = find_nearest(S_sl,sl)

			ra_f = ra[id_max]
			dec_f= dec[id_max]
			dd= minute+a_ini

		elif find == 'min' and np.amin(S_sl[np.where(S_sl>0)]) < sl: 
			sl = np.amin(S_sl[np.where(S_sl>0)])
			id_max = find_nearest(S_sl,sl)

			ra_f = ra[id_max]
			dec_f = dec[id_max]
			dd = minute+a_ini

		minute += 1
	return dd, ra_f, dec_f, sl


def hilite(string, status, bold):
	'''Graphism: colors and bold in the terminal'''
	import sys

	if not sys.stdout.isatty() : return '*'+string+'*'
	
	attr = []
	if status:
		# green
		attr.append('32')
	else:
		# red
		attr.append('31')
	if bold:
		attr.append('1')
	return '\x1b[%sm%s\x1b[0m' % (';'.join(attr), string)



def sort_boundary(x,y, dra, ddec, force=False):
	''' Sort the targets in such a way that the boundaries are drawn correctly'''
	# prepare the array from the sorted coordinates
	new_x = np.array([x[0]])
	new_y = np.array([y[0]])

	# remove element in list:
	x = np.delete(x, 0, axis=0)
	y = np.delete(y, 0, axis=0)

	goto = len(y)
	for ii in range(0, goto):
		id_np = -1
		
		# is there one next to it (from left to right) :
		a = np.where(np.abs(new_x[ii]+dra-x)<1e-3)[0]
		b = np.where(np.abs(y-new_y[ii])<1e-3)[0]
		if not np.shape(np.intersect1d(a,b))[0] == 0: 
			id_np = np.intersect1d(a,b)[0]
			nx = x[id_np]
			ny = y[id_np]


		# looking directly up
		if np.shape(np.intersect1d(a,b))[0] == 0:
			# first look if there is one immediately above
			a = np.where(np.abs(new_y[ii]+ddec-y)<1e-3)[0]			
			b = np.where(np.abs(new_x[ii]-x)<1e-3)[0]			
			if not np.shape(np.intersect1d(a,b))[0] == 0: 
				id_np = np.intersect1d(a,b)[0]
				nx = x[id_np]
				ny = y[id_np]


		# is there one next to it (from right to left) :
		if np.shape(np.intersect1d(a,b))[0] == 0:
			a = np.where(np.abs(new_x[ii]-dra-x)<1e-3)[0]
			b = np.where(np.abs(y-new_y[ii])<1e-3)[0]
			if not np.shape(np.intersect1d(a,b))[0] == 0: 
				id_np = np.intersect1d(a,b)[0]
				nx = x[id_np]
				ny = y[id_np]


		# looking directly down
		if np.shape(np.intersect1d(a,b))[0] == 0:
			a = np.where(np.abs(new_y[ii]-ddec-y)<1e-3)[0]			
			b = np.where(np.abs(new_x[ii]-x)<1e-3)[0]			
			if not np.shape(np.intersect1d(a,b))[0] == 0: 
				id_np = np.intersect1d(a,b)[0]
				nx = x[id_np]
				ny = y[id_np]

		# down, close
		if np.shape(np.intersect1d(a,b))[0] == 0:
			# first look if there is one immediately below
			a = np.where(np.abs(new_y[ii]-ddec-y)<1e-3)[0]			
			b = np.where(np.abs(new_x[ii]-x)<1.5)[0]	
			xx = x[b]
			try :
				b = np.array([find_nearest(new_x[ii], xx)])
			except ValueError:
				#return new_x, new_y
				pass
			if not np.shape(np.intersect1d(a,b))[0] == 0: 
				id_np = np.intersect1d(a,b)[0]
				nx = xx[id_np]
				ny = y[id_np]
				id_np = find_nearest(xx[id_np], x)

		# further down, close
		if np.shape(np.intersect1d(a,b))[0] == 0:
			# first look if there is one immediately below
			a = np.where(np.abs(new_y[ii]-2.*ddec-y)<1e-3)[0]			
			b = np.where(np.abs(new_x[ii]-x)<1.5)[0]	
			xx = x[b]
			try :
				b = np.array([find_nearest(new_x[ii], xx)])
			except ValueError:
				#return new_x, new_y
				pass
			if not np.shape(np.intersect1d(a,b))[0] == 0: 
				id_np = np.intersect1d(a,b)[0]
				nx = xx[id_np]
				ny = y[id_np]
				id_np = find_nearest(xx[id_np], x)

		# way up, close
		if np.shape(np.intersect1d(a,b))[0] == 0 and force:
			# first look if there is one immediately below
			a = np.where(np.abs(new_y[ii]-y)<2.)[0]			
			b = np.where(np.abs(new_x[ii]-x)<1.5)[0]	
			xx = x[b]
			try :
				b = np.array([find_nearest(new_x[ii], xx)])
			except ValueError:
				#return new_x, new_y
				pass
			if not np.shape(np.intersect1d(a,b))[0] == 0: 
				id_np = np.intersect1d(a,b)[0]
				nx = xx[id_np]
				ny = y[id_np]
				id_np = find_nearest(xx[id_np], x)
		
		# if nothing return the empty sorted list to avoid errors
		if id_np == -1: return new_x, new_y, x, y

		# otherwise, add the points to the list
		new_x = np.append(new_x,nx)
		new_y = np.append(new_y,ny)

		# remove element in list of available points:
		x = np.delete(x, id_np, axis=0)
		y = np.delete(y, id_np, axis=0)
	

	return new_x, new_y, x, y

