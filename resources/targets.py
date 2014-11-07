''' targets.py
=========================
AIM:	Class headers for the target_list class and load a catalogue of objects of interest

INPUT:	function depend

OUTPUT:	function depend

CMD:	To include: from resources.targets import *

ISSUES:	<none known>

REQUIRES: standard python libraries, specific libraries in resources/

REMARKS: <none>
'''


class target_list(object):
	''' CLASS CHEOPS OBJECT OF INTEREST
	Handle the targets and their visibility'''
	import numpy as np
	nb_targets = 0
	obs_max=2500

	def __init__(self, name, ra, id_ra, dec, id_dec, mag, Int_period):
		from numpy import zeros
		from resources.routines import mag2flux
		self.name = name
		self.ra = ra
		self.id_ra = id_ra
		self.id_dec= id_dec
		self.dec= dec
		self.mag= mag
		self.max_mag=mag2flux(mag)

		self.visible = zeros([0])
		self.invisible = zeros([0])
		self.workspace = 0
		self.current_visibility = 0
		self.current_SAA_interruption = 0
		self.visible_save = zeros(Int_period)
		target_list.nb_targets+=1
		self.ii = 0

	def CountObjects(self):
		return target_list.nb_targets

	def Coordinates(self):
		return self.ra, self.dec

	def maximum_flux(self):
		return self.max_mag

	def Visibility(self):
		return self.visible

	def Invisibility(self):
		return self.invisible

	def PrepareSave(self):
		if self.ii < target_list.obs_max: 
			self.visible = self.visible[:self.ii]
			self.invisible = self.invisible[:self.ii]
		return 1

	def In_cell(self, id_ra_up, id_dec_up):
		if id_ra_up == self.id_ra and id_dec_up == self.id_dec: return True
		else: return False

	def Appear(self,minute):
		from numpy import append	
		try:
			self.visible[self.ii] = minute
		except IndexError:
			self.visible = append(self.visible, minute)

	def Disappear(self,minute):
		from numpy import append
		try:
			self.invisible[self.ii] = minute
		except IndexError:
			self.invisible = append(self.invisible, minute)
		self.workspace = 0
		self.current_SAA_interruption=0

	def Is_visible(self,minute):
		from numpy import shape, where
		if shape(self.visible)[0] == 0 : return False

		index = where(self.visible <= minute)

		try: dis = self.invisible[index][-1]
		except IndexError: return True
		
		if dis > minute : return True
		else: return False

	def Interruption(self, minimum):
		if (not self.current_visibility == 1) and self.workspace >= minimum: return True
		else : return False

	def Continuity(self):
		if (self.current_visibility == 1): return True
		else : return False

	def Next(self,minute,minimum):
		if self.Interruption(minimum):
#			print self.workspace
			self.Appear(minute-self.workspace)
			self.Disappear(minute-self.current_SAA_interruption)
			if (self.workspace<=self.current_SAA_interruption and self.workspace>0):
				# This almost could be an assert, but we want a nice error message
				print
				print "obj", self.ra, self.dec
				print minute
				print 'Appeared at %d' % (minute-self.workspace)
				print 'Disappeared at %d' % (minute)
				print 'Time spent in SAA is %d' % (self.current_SAA_interruption)
				print 'Thus >> obs time = %d' % (self.invisible[-1]-self.visible[-1])
				raise RuntimeError("More time in SAA than total visit time!")
			#if self.current_SAA_interruption>10:
			#	exit()
			
#			print self.workspace
#			print '*'*13
			self.ii += 1
		elif self.Continuity(): self.workspace += 1
		else: self.workspace = 0
			
		self.current_visibility = 0





def load_catalogue(fname,folder='resources', remove_duplicates = True, verbose=True):
	import numpy as np
	import csv

	name = list()
	ra = np.empty(0)
	dec = np.empty(0)
	mag = np.empty(0)

	with open(folder+"/"+fname, "rb") as catalogue_file:
		catalogue_data = csv.reader(catalogue_file, delimiter="\t")
		for row in catalogue_data:
		# if line is empty, skip otherwise filter out the blank
			if len(row) > 0:
				line = row[0].split(' ')
				if filter(None, line)[0] == '#': continue

				name.append( row[0][0:12].replace(' ', '') )

				ra = np.append(ra, float(row[0][29:38].replace(' ', '')) )
				dec= np.append(dec, float(row[0][55:64].replace(' ', '')) )
				mag= np.append(mag, float(row[0][127:133].replace(' ', '')) )

	nbtargets = len(name)

	if remove_duplicates:
		def idfun(x): return x
		seen = {}
		result = []
		ra_ = np.empty(0)
		dec_= np.empty(0)
		mag_= np.empty(0)
		for index, item in enumerate(name):
			marker = idfun(item)
       # in old Python versions:
       # if seen.has_key(marker)
       # but in new ones:
			if marker in seen: continue
			seen[marker] = 1
			result.append(item)
			ra_=np.append(ra_,ra[index])
			dec_=np.append(dec_,dec[index])
			mag_=np.append(mag_,mag[index])
		name = result
		ra = ra_
		dec= dec_
		mag= mag_

	if verbose:
		message = 'Loaded %d targets' % (len(name))
		if remove_duplicates: message = '%s (%d duplicates were removed)' % (message, nbtargets-len(name))
		message = '%s from catalogue %s' % (message, fname)
		print message

	return name, ra, dec, mag

