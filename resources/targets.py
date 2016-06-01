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

import numpy as np

class target_list(object):
	''' CLASS CHEOPS OBJECT OF INTEREST
	Handle the targets and their visibility'''
	nb_targets = 0
	obs_max=2500

	def __init__(self, name, ra, id_ra, dec, id_dec, mag, Int_period, flux=None):
		from numpy import zeros

		self.name = name
		self.ra = ra
		self.id_ra = id_ra
		self.id_dec= id_dec
		self.dec= dec
		self.mag= mag

		if flux is None:
			from resources.routines import mag2flux
			self.flux=mag2flux(mag)
		else:
			self.flux = flux

		self.visible = zeros([0])
		self.invisible = zeros([0])
		self.interruption_time = zeros([0])
		self.obs_time = 0
		self.current_visibility = 0
		self.current_SAA_interruption = 0
		self.visible_save = zeros(Int_period)
		target_list.nb_targets+=1
		self.ii = 0

	def CountObjects(self):
		return target_list.nb_targets

	def Coordinates(self):
		return self.ra, self.dec

	def get_flux(self):
		return self.flux

	def Visibility(self):
		return self.visible

	def Invisibility(self):
		return self.invisible

	def get_interruption_time(self):
		return self.interruption_time

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
		try:
			self.invisible[self.ii] = minute
		except IndexError:
			self.invisible = np.append(self.invisible, minute)
		self.obs_time = 0
		self.Count_interruption_time()
		
	def Count_interruption_time(self):
		try:
			self.interruption_time[self.ii] = self.current_SAA_interruption
		except IndexError:
			self.interruption_time = np.append(self.interruption_time, self.current_SAA_interruption)
		self.current_SAA_interruption=0

	def Is_visible(self,minute):
		if np.shape(self.visible)[0] == 0 : return False

		index = np.where(self.visible <= minute)

		try: dis = self.invisible[index][-1]
		except IndexError: return True
		
		if dis > minute : return True
		else: return False

	def Interruption(self, minimum):
		if (not self.current_visibility == 1) and self.obs_time >= minimum: return True
		else : return False

	def Continuity(self):
		if (self.current_visibility == 1): return True
		else : return False

	def Next(self,minute,minimum):
		if self.Interruption(minimum):
#			print self.obs_time
			self.Appear(minute-self.obs_time)
			self.Disappear(minute)
			if (self.obs_time<=self.current_SAA_interruption and self.obs_time>0):
				# This almost could be an assert, but we want a nice error message
				print
				print "obj", self.ra, self.dec
				print minute
				print 'Appeared at %d' % (minute-self.obs_time)
				print 'Disappeared at %d' % (minute)
				print 'Time spent in SAA is %d' % (self.current_SAA_interruption)
				print 'Thus >> obs time = %d' % (self.invisible[-1]-self.visible[-1])
				raise RuntimeError("More time in SAA than total visit time!")
			#if self.current_SAA_interruption>10:
			#	exit()
			
#			print self.obs_time
#			print '*'*13
			self.ii += 1
		elif self.Continuity(): self.obs_time += 1
		else: self.obs_time = 0
			
		self.current_visibility = 0





def load_catalogue(fname,folder='resources', remove_duplicates = True, verbose=True):
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

