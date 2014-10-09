import numpy as np

def equatorial2ecliptic(a,b,obliquity=np.deg2rad(23.4333)):
	# http://en.wikipedia.org/wiki/Celestial_coordinate_system
	e=obliquity

	up=np.sin(a)*np.cos(e)+np.tan(b)*np.sin(e)
	dw=np.cos(a)
	ra=np.arctan2(up,dw)

	dec = np.sin(b)*np.cos(e)-np.cos(b)*np.sin(e)*np.sin(a)
	dec = np.arcsin(dec)

	return np.asarray([ra,dec]).T


def ecliptic2equatorial(a,b,obliquity=np.deg2rad(23.4333)):
	# http://en.wikipedia.org/wiki/Celestial_coordinate_system
	e=obliquity

	up=np.sin(a)*np.cos(e)-np.tan(b)*np.sin(e)
	dw=np.cos(a)
	lamb=np.arctan2(up,dw)

	beta = np.sin(b)*np.cos(e)+np.cos(b)*np.sin(e)*np.sin(a)
	beta = np.arcsin(beta)

	return np.asarray([lamb,beta]).T
	
