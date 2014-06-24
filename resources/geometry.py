
import numpy as np
try: import shapely.geometry
except:
	print "I can't import shapely."
	print "To install, try"
	print "pip install --user Shapely"

def check_cood(coord):
	try:
		assert (coord[0]>=0. and coord[0]<=2.*np.pi)
		assert (coord[0]>=-np.pi and coord[1]<=np.pi)
	except: raise AssertionError("Coordinate [%3.2f,%3.2f] out of bounds..." % (coord[0],coord[1]))
	return True

class Polygon:
	def __init__(self, points, unit='deg'):
		"""
		points is a list of tuples (x, y) with the corners of the polygons.
		This list will be closed automatically
		For full info, see
		http://toblerity.org/shapely/manual.html#polygons
		
		"""
		if unit=='deg':self.points=[(point[0]/180.*np.pi, point[1]/180.*np.pi) for point in points]
		elif unit=='rad': self.points = points
		else: raise ValueError("The unit must be either deg or rad")

		self.poly = shapely.geometry.Polygon(self.points)

	def __str__(self):
		return 'Polygon with edges (in rad): '+str(self.points)

	def is_inside(self, coord):
		"""
		is_insides a direction. Returns True if the direction is inside of the polygon, False if it is outside.
		"""
		check_cood(coord)
		p = shapely.geometry.Point(coord[0], coord[1])
		
		if not self.poly.intersects(p):
			return False
		else:
			return True

class Interval:

	def __init__(self, axis, min_val=None, max_val=None, unit='deg'):
		"""
		direction must be alpha or delta
		min_val and max_val are optional. Specify only one of them to is_inside only upper or lower bounds...
		"""
		if axis == 'alpha': self.axis=0
		elif axis == 'delta': self.axis=1
		else: raise ValueError("Axis must be alpha or delta")

		self.min_val = min_val
		self.max_val = max_val
		if unit=='deg':
			if not min_val==None: self.min_val = min_val/180.*np.pi
			if not max_val==None: self.max_val = max_val/180.*np.pi
		elif not unit=='rad': raise ValueError("The unit must be either deg or rad")
	

	def __str__(self):
		names={0:'alpha',1:'delta'}
		if self.min_val==None: min_val="None"
		else: min_val=self.min_val
		if self.max_val==None: max_val="None"
		else: max_val=self.max_val
		return 'Interval on %s with edges (in rad): [%s,%s]' % (names[self.axis],str(min_val),str(max_val))
	
	def is_inside(self, coord):
		"""
		is_insides a direction. Returns True if the direction is in the interval, False if it is outside.
		"""
		check_cood(coord)
		if self.min_val != None:
			if coord[self.axis] < self.min_val:
				return False
				
		if self.max_val != None:
			if coord[self.axis] > self.max_val:
				return False
				
		return True
