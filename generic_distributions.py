"""Deprecated. Use SPARCBasicDistributionManager instead."""

from distributions import *

class NonconsecutiveFreqDistManager(FrequencyDistributionManager):
	def __init__(self, path, freq_path, axisradii_path=None):
		super(NonconsecutiveFreqDistManager,self).__init__(path, freq_path, axisradii_path)
		self.type = frequency_nonconsec_disttype

class ConsecutiveFreqDistManager(FrequencyDistributionManager):
	def __init__(self, path, freq_path, axisradii_path=None):
		super(ConsecutiveFreqDistManager,self).__init__(path, freq_path, axisradii_path)
		self.type = frequency_consec_disttype
		self.blocks_secondary_structures = False

	def _update_score(self, protein, data, score, prior=2):
		nearby = []
		if aa.tag > 0 and prior != False: nearby.append(protein.aminoacids[aa.tag - 1])
		if aa.tag < len(protein.aminoacids) - 1 and prior != True: nearby.append(protein.aminoacids[aa.tag + 1])
		localscore = 0.0
		for aa2 in nearby:
			hypo = next((x for x in data if x.tag == aa2.tag), None)
			if hypo is not None: aa2 = hypo
			tag1 = aacode(aa.type)
			tag2 = aacode(aa2.type)
			if tag1 >= AMINO_ACID_COUNT: tag1 = 0
			if tag2 >= AMINO_ACID_COUNT: tag2 = 0
			# Determine the score based on the reference amino acid's perspectives and its neighbors'.
			# Formula: S = -ln(F/F0) if F > 0, +5 otherwise
			
			if aa2.acarbon.distanceto(aa.acarbon) > 5.0:
				score += 10.0
				continue
			zone = aa2.tolocal(aa.acarbon).floor()
			if aa.tag < aa2.tag and zone.x > 0.0:
				score += 10.0
				continue
			freq = self.alpha_frequency(tag2, tag1, zone)
			if freq > 0:
				subscore = -math.log(freq / self.median_frequencies[tag2][tag1])
			else:
				subscore = 10.0
			score += subscore
			
			zone = aa.tolocal(aa2.acarbon).floor()
			freq = self.alpha_frequency(tag1, tag2, zone)
			if freq > 0:
				subscore = -math.log(freq / self.median_frequencies[tag1][tag2])
			else:
				subscore = 10.0
			score += subscore
			localscore += subscore
		#aa.localscore = localscore
		return score

	def score(self, protein, data, prior=2):
		"""This subclassed implementation takes into account consecutive amino acids even when they are not within 10 angstroms. Also, pass True or False for prior to use only the preceding/succeeding amino acids of each aa in data."""
		score = 0.0
		consec = 2
		if self.type == frequency_consec_disttype: consec = 1
		elif self.type == frequency_nonconsec_disttype: consec = 0
		for aa in data:
			if protein.secondary_structure_aa(aa.tag) and self.blocks_secondary_structures: continue
			score = self._update_score(protein, data, score, prior)
		return score * self.weight

class SecondaryFreqDistManager(ConsecutiveFreqDistManager):
	"""This subclass of the generic consecutive frequency distribution manager handles only those interactions which occur within a secondary structure."""
	
	def __init__(self, path, freq_path, axisradii_path=None):
		super(SecondaryFreqDistManager,self).__init__(path, freq_path, axisradii_path)
		self.type = frequency_consec_disttype

	def score(self, protein, data, prior=2):
		"""This subclassed implementation takes into account consecutive amino acids even when they are not within 10 angstroms, but only when they are within a secondary structure. Also, pass True or False for prior to use only the preceding/succeeding amino acids of each aa in data."""
		score = 0.0
		consec = 2
		if self.type == frequency_consec_disttype: consec = 1
		elif self.type == frequency_nonconsec_disttype: consec = 0
		for aa in data:
			if not protein.secondary_structure_aa(aa.tag): continue
			score = self._update_score(protein, data, score, prior)
		return score * self.weight
