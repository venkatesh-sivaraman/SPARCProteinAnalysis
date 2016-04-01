from distributions import *
from secondary_structure import *
import reference_state
from loading_indicator import *

#MARK: - Distributions

sec_structs = [	secondary_struct_helix + "1", secondary_struct_helix + "2",
			   secondary_struct_helix + "3", secondary_struct_helix + "4",
			   secondary_struct_helix + "5", secondary_struct_helix + "6",
			   secondary_struct_helix + "7", secondary_struct_helix + "8",
			   secondary_struct_helix + "9", secondary_struct_helix + "10",
			   secondary_struct_sheet + "0", secondary_struct_sheet + "1",
			   secondary_struct_sheet + "-1"]

def _secondary_structures_dict(inner_value={}):
	global sec_structs
	sec_dict = {}
	for ss in sec_structs:
		sec_dict[ss] = inner_value
	return sec_dict

def _is_valid_secondary_structure(struct_type):
	global sec_structs
	if not struct_type:
		return False
	if isinstance(struct_type, basestring):
		return struct_type in sec_structs
	if len(struct_type) > 1:
		return (struct_type[0].type + str(struct_type[1].identifiers[0])) in sec_structs

class SPARCBasicDistributionManager (FrequencyDistributionManager):
	"""The basic distribution handles just frequencies. It represents a noncommutative potential function for EITHER nonconsecutive or consecutive data."""
	
	def __init__(self, frequencies_path, isconsec=False, blocks_sec_struct=False, short_range=2, references=None):
		"""frequencies_path should be a path to a directory of alpha zones paired with frequencies for individual amino acid pairs. If blocks_sec_struct is True, then secondary structures will be ignored. To treat secondary structures exclusively, use a SPARCSecondaryDistributionManager. The default value for short_range is 2, which means that all nonconsecutive interactions are captured. If you pass True or False for short_range, those interactions will be either considered exclusively or excluded."""
		self.alpha_frequencies = {}
		self.reference_frequencies = None #{}
		self.total_interactions = [[0 for i in xrange(AMINO_ACID_COUNT)] for j in xrange(AMINO_ACID_COUNT)]
		self.median_frequencies = [[0 for i in xrange(AMINO_ACID_COUNT)] for j in xrange(AMINO_ACID_COUNT)]
		self.blocks_secondary_structures = blocks_sec_struct
		self.short_range = short_range
		self.total_median = 0
		self.identifier = os.path.basename(frequencies_path)
		self.load_frequencies(frequencies_path)
		if references:
			self.load_references(references)
		self.weight = 1.0
		if isconsec:
			self.type = frequency_consec_disttype
		else:
			self.type = frequency_nonconsec_disttype
		self.defaultvalue = 0
		self.refstate = False

	def __repr__(self):
		return "<Distribution Manager for '{}' data>".format(self.identifier)
	
	def alpha_frequency(self, type1, type2, zone):
		"""This helper function retrieves the frequency of 'zone' in the loaded frequency data. type1 refers to the type (code) of the amino acid that serves as the origin of the zone space, while type2 is the amino acid to which the zone refers."""
		if zone in self.alpha_frequencies:
			'''if self.reference_frequencies:
				ref_zone = zone
				if ref_zone.x < 0.0: ref_zone.x = -ref_zone.x - 1.0
				if ref_zone.y < 0.0: ref_zone.y = -ref_zone.y - 1.0
				if ref_zone.z < 0.0: ref_zone.z = -ref_zone.z - 1.0
				if ref_zone in self.reference_frequencies:
					return (self.alpha_frequencies[zone][type1][type2], self.reference_frequencies[ref_zone][type1][type2])'''
			return self.alpha_frequencies[zone][type1][type2]
		else:
			return 0

	def subscore(self, protein, aa, aa2, onlyone=False, consec=2, zero_value=0.01):
		tag1 = aacode(aa.type)
		tag2 = aacode(aa2.type)
		if tag1 >= AMINO_ACID_COUNT: tag1 = 0
		if tag2 >= AMINO_ACID_COUNT: tag2 = 0
		
		if self.refstate:
			zone = aa2.tolocal(aa.acarbon).floor()
			subscore = self.alpha_frequency(tag2, tag1, zone)
			if subscore == 0:
				subscore = zero_value
			zone = aa.tolocal(aa2.acarbon).floor()
			subscore2 = self.alpha_frequency(tag1, tag2, zone)
			if subscore2 == 0:
				subscore2 = zero_value
			if onlyone:
				print subscore, subscore2
				return (subscore * self.weight, subscore2 * self.weight, self.total_interactions[tag1][tag2], self.total_interactions[tag2][tag1])
			vol = protein.volume()
			contact_type = 2
			if self.short_range: contact_type = 1
			elif consec == 1: contact_type = 0
			return -math.log(((subscore * subscore2) / (self.total_interactions[tag1][tag2] ** 2)) * reference_state.contact_probability(tag1, tag2, vol, contact_type) / (4000.0 * math.pi / (3 * vol) * reference_state.position_ref(zone)))
			#return -math.log((subscore * subscore2 / (self.total_interactions[tag2][tag1] ** 2)) / reference_state.position_ref(zone))
		else:
			zone = aa2.tolocal(aa.acarbon).floor()
			subscore = self.alpha_frequency(tag2, tag1, zone)
			if subscore == 0:
				subscore = zero_value
			subscore = -math.log(subscore / self.median_frequencies[tag2][tag1] * self.total_interactions[tag2][tag1] / self.total_median)
			zone2 = aa.tolocal(aa2.acarbon).floor()
			subscore2 = self.alpha_frequency(tag1, tag2, zone2)
			if subscore2 == 0:
				subscore2 = zero_value
			subscore2 = -math.log(subscore2 / self.median_frequencies[tag1][tag2] * self.total_interactions[tag1][tag2] / self.total_median)
			if onlyone:
				print subscore, subscore2
				return (subscore * self.weight, subscore2 * self.weight, self.total_interactions[tag1][tag2], self.total_interactions[tag2][tag1])
			return subscore + subscore2

	def score(self, protein, data, system=None, isolate=False, onlyone=False, prior=2, zero_value=0.01):
		"""For frequency distributions, pass in an array of hypothetical aminoacids. This implementation returns the product of the frequencies of each pairwise interaction. If isolate=True, only the amino acids in data will be considered for the energy calculation.
			Pass prior to consider ONLY the amino acid before (True) or after (False) each amino acid in data. This works best for consecutive distributions."""
		score = 0.0
		consec = 2
		taglist = {}
		if self.type == frequency_consec_disttype: consec = 1
		elif self.type == frequency_nonconsec_disttype: consec = 0
		for aa in data:
			if not aa: continue
			sec_struct = protein.secondary_structure_aa(aa.tag)
			if self.blocks_secondary_structures and _is_valid_secondary_structure(sec_struct): continue
			if prior != 2:
				nearby = []
				if aa.tag > 0 and prior != False: nearby.append(protein.aminoacids[aa.tag - 1])
				if aa.tag < len(protein.aminoacids) - 1 and prior != True: nearby.append(protein.aminoacids[aa.tag + 1])
			else:
				if system and not consec and self.short_range != 0:
					nearby = system.nearby_aa(aa, protein, 10.0, consec=consec)
				else:
					nearby = protein.nearby_aa(aa, 10.0, consec=consec)
			for aa2 in nearby:
				if not aa2: continue
				if aa2.tag - aa.tag == 1 and aa.has_break: continue
				elif math.fabs(aa2.tag - aa.tag) > 5 and self.short_range == 1: continue
				elif math.fabs(aa2.tag - aa.tag) <= 5 and self.short_range == 0: continue
				if (aa.tag in taglist and aa2.tag in taglist[aa.tag]) or (aa2.tag in taglist and aa.tag in taglist[aa2.tag]):
					continue
				hypo = next((x for x in data if x and x.tag == aa2.tag), None)
				if hypo is not None: aa2 = hypo
				elif isolate: continue
				
				try:
					subscore = self.subscore(protein, aa, aa2, onlyone, consec, zero_value)
				except ZeroDivisionError:
					subscore = 0.0
				
				if onlyone: return subscore
				else: score += subscore
				
				if aa.tag in taglist:
					taglist[aa.tag].append(aa2.tag)
				else:
					taglist[aa.tag] = [aa2.tag]
		return score * self.weight
					
	def read_frequency_line(self, line, tag1, tag2, all):
		"""Helper method for load_frequencies, intended for subclasses to easily modify the reading procedure."""
		ptcomps, freq = line.strip().split(";")
		alpha = Point3D(*ptcomps.split(","))
		if all:
			self.default_alpha_frequencies[alpha] = float(freq.strip().split(" ")[0])
		else:
			if alpha not in self.alpha_frequencies:
				self.alpha_frequencies[alpha] = [[0 for i in xrange(AMINO_ACID_COUNT)] for k in xrange(AMINO_ACID_COUNT)]
			self.alpha_frequencies[alpha][tag1][tag2] = float(freq)
			self.total_interactions[tag1][tag2] += float(freq)
	
	def load_frequencies(self, path):
		files = os.listdir(path)
		percentage = 0
		self.total_interactions = [[0 for i in xrange(AMINO_ACID_COUNT)] for j in xrange(AMINO_ACID_COUNT)]
		loading_indicator.add_loading_data(len(files))
		for n, indfile in enumerate(files):
			loading_indicator.update_progress(1)
			if indfile.find(".txt") == -1 or indfile[0] == ".": continue
			all = ("all" in indfile)
			if not all:	tag1, tag2 = indfile[0:-4].split('-')
			tag1 = int(tag1)
			tag2 = int(tag2)
			with open(join(path, indfile), 'r') as file:
				for line in file:
					if ";" not in line and len(line.strip()) > 0:
						#if not self.median_frequencies[tag1][tag2]:
						self.median_frequencies[tag1][tag2] = float(line)
						continue
					self.read_frequency_line(line, tag1, tag2, all)
		#Compute median total frequency
		s = sorted([x for list1 in self.total_interactions for x in list1])
		self.total_median = sum(s) / float(len(s)) #s[int(len(s) / 2.0)]

	def load_references(self, path):
		files = os.listdir(path)
		percentage = 0
		self.reference_frequencies = {}
		loading_indicator.add_loading_data(len(files))
		for n, indfile in enumerate(files):
			loading_indicator.update_progress(1)
			if indfile.find(".txt") == -1: continue
			tag1, tag2 = indfile[0:-4].split('-')
			tag1 = int(tag1)
			tag2 = int(tag2)
			with open(join(path, indfile), 'r') as file:
				for line in file:
					if ";" not in line and len(line.strip()) > 0:
						continue
					ptcomps, freq = line.strip().split(";")
					alpha = Point3D(*ptcomps.split(","))
					if alpha not in self.reference_frequencies:
						self.reference_frequencies[alpha] = [[0 for i in xrange(AMINO_ACID_COUNT)] for k in xrange(AMINO_ACID_COUNT)]
					self.reference_frequencies[alpha][tag1][tag2] = float(freq)
		print "Loaded references"


class SPARCSecondaryDistributionManager (SPARCBasicDistributionManager):
	"""This specialized distribution manager handles frequencies for secondary structures."""
	
	def __init__(self, frequencies_path):
		"""frequencies_path should be a path to a directory of alpha zones paired with frequencies for individual secondary structure types."""
		self.alpha_frequencies = {}
		self.total_interactions = _secondary_structures_dict(inner_value=0)
		self.median_frequencies = _secondary_structures_dict(inner_value=0)
		self.blocks_secondary_structures = False
		self.identifier = os.path.basename(frequencies_path)
		self.total_median = 0
		self.load_frequencies(frequencies_path)
		self.weight = 1.0
		self.type = frequency_consec_disttype
		self.defaultvalue = 0
		self.refstate = False
	
	def alpha_frequency(self, sec_struct_id, zone):
		"""This helper function retrieves the frequency of 'zone' in the loaded frequency data. sec_struct_id is the string, which is also a key in the data (e.g. "helix1", "sheet0")."""
		if zone in self.alpha_frequencies:
			return self.alpha_frequencies[zone][sec_struct_id]
		else:
			return 0

	def subscore(self, protein, aa, aa2, sec_struct_type, onlyone=False):
		tag1 = aacode(aa.type)
		tag2 = aacode(aa2.type)
		if tag1 >= AMINO_ACID_COUNT: tag1 = 0
		if tag2 >= AMINO_ACID_COUNT: tag2 = 0
		
		if self.refstate == True:
			zone = aa2.tolocal(aa.acarbon).floor()
			subscore = self.alpha_frequency(sec_struct_type, zone)
			if subscore == 0:
				subscore = 1e-2
			zone2 = aa.tolocal(aa2.acarbon).floor()
			subscore2 = self.alpha_frequency(sec_struct_type, zone2)
			if subscore2 == 0:
				subscore2 = 1e-2
			vol = protein.volume()
			contact_type = 2
			if self.short_range: contact_type = 1
			elif consec == 1: contact_type = 0
			subscore = -math.log(((subscore * subscore2) / (self.total_interactions[sec_struct_type] ** 2)) * reference_state.contact_probability(tag1, tag2, vol, contact_type) / (4000.0 * math.pi / (3 * vol) * reference_state.position_ref(zone)))
			if onlyone:
				print subscore
				return (subscore * self.weight, subscore * self.weight, self.total_interactions[sec_struct_type], self.total_interactions[sec_struct_type])
			return subscore
		else:
			zone = aa2.tolocal(aa.acarbon).floor()
			subscore = self.alpha_frequency(sec_struct_type, zone)
			if subscore == 0:
				subscore = 1e-2
			subscore = -math.log(subscore / self.median_frequencies[sec_struct_type] * self.total_interactions[sec_struct_type] / self.total_median)
			zone = aa.tolocal(aa2.acarbon).floor()
			subscore2 = self.alpha_frequency(sec_struct_type, zone)
			if subscore2 == 0:
				subscore2 = 1e-2
			subscore2 = -math.log(subscore2 / self.median_frequencies[sec_struct_type] * self.total_interactions[sec_struct_type] / self.total_median)
			if onlyone:
				print subscore, subscore2
				return (subscore * self.weight, subscore2 * self.weight, self.total_interactions[sec_struct_type], self.total_interactions[sec_struct_type])
			return subscore + subscore2

	def score(self, protein, data, system=None, isolate=False, onlyone=False, prior=2):
		"""For frequency distributions, pass in an array of hypothetical aminoacids. This implementation returns the product of the frequencies of each pairwise interaction. If isolate=True, only the amino acids in data will be considered for the energy calculation.
			Pass prior to consider ONLY the amino acid before (True) or after (False) each amino acid in data. This works best for consecutive distributions."""
		score = 0.0
		consec = 2
		taglist = {}
		if self.type == frequency_consec_disttype: consec = 1
		elif self.type == frequency_nonconsec_disttype: consec = 0
		for aa in data:
			if not aa: continue
			sec_struct = protein.secondary_structure_aa(aa.tag)
			if not _is_valid_secondary_structure(sec_struct): continue
			sec_struct_type = sec_struct[0].type + str(sec_struct[1].identifiers[0])
			if prior != 2:
				nearby = []
				if aa.tag > 0 and prior != False: nearby.append(protein.aminoacids[aa.tag - 1])
				if aa.tag < len(protein.aminoacids) - 1 and prior != True: nearby.append(protein.aminoacids[aa.tag + 1])
			else:
				if system:
					nearby = system.nearby_aa(aa, protein, 10.0, consec=consec)
				else:
					nearby = protein.nearby_aa(aa, 10.0, consec=consec)
			for aa2 in nearby:
				if (aa.tag in taglist and aa2.tag in taglist[aa.tag]) or (aa2.tag in taglist and aa.tag in taglist[aa2.tag]):
					continue
				hypo = next((x for x in data if x and x.tag == aa2.tag), None)
				if hypo is not None: aa2 = hypo
				elif isolate: continue
				if aa2.tag - aa.tag == 1 and aa.has_break: continue

				try:
					subscore = self.subscore(protein, aa, aa2, sec_struct_type, onlyone)
				except ZeroDivisionError:
					subscore = 0.0

				if onlyone: return subscore
				else: score += subscore

				if aa.tag in taglist:
					taglist[aa.tag].append(aa2.tag)
				else:
					taglist[aa.tag] = [aa2.tag]
		return score * self.weight

	def read_frequency_line(self, line, sec_struct_type, all):
		ptcomps, freq = line.strip().split(";")
		alpha = Point3D(*ptcomps.split(","))
		if alpha not in self.alpha_frequencies:
			self.alpha_frequencies[alpha] = _secondary_structures_dict(inner_value=0)
		self.alpha_frequencies[alpha][sec_struct_type] = float(freq)
		self.total_interactions[sec_struct_type] += float(freq)
	
	def load_frequencies(self, path):
		files = os.listdir(path)
		percentage = 0
		self.total_interactions = _secondary_structures_dict(inner_value=0)
		loading_indicator.add_loading_data(len(files))
		for n, indfile in enumerate(files):
			loading_indicator.update_progress(1)
			if indfile.find(".txt") == -1: continue
			if "all" in indfile:
				continue
			sec_struct_type = indfile[0:-4]
			with open(join(path, indfile), 'r') as file:
				for line in file:
					if ";" not in line and len(line.strip()) > 0:
						self.median_frequencies[sec_struct_type] = float(line)
						continue
					self.read_frequency_line(line, sec_struct_type, False)
		#Compute median total frequency
		s = sorted(self.total_interactions.values())
		self.total_median = sum(s) / float(len(s)) #s[int(len(s) / 2.0)]

#MARK: - Experimental Both-Orientation Potential

class SPARCBothOrientationDistributionManager (SPARCBasicDistributionManager):
	"""This class manages distributions in which both amino acids' orientations are specified together."""

	def __init__(self, frequencies_path, isconsec=False, blocks_sec_struct=False, short_range=2, references=None):
		"""frequencies_path should be a path to a directory of alpha zones paired with frequencies for individual amino acid pairs. If blocks_sec_struct is True, then secondary structures will be ignored. To treat secondary structures exclusively, use a SPARCSecondaryDistributionManager. The default value for short_range is 2, which means that all nonconsecutive interactions are captured. If you pass True or False for short_range, those interactions will be either considered exclusively or excluded."""
		self.alpha_frequencies = [[{} for i in xrange(AMINO_ACID_COUNT)] for k in xrange(AMINO_ACID_COUNT)]
		self.reference_frequencies = None #{}
		self.total_interactions = [[0 for i in xrange(AMINO_ACID_COUNT)] for j in xrange(AMINO_ACID_COUNT)]
		self.median_frequencies = [[0 for i in xrange(AMINO_ACID_COUNT)] for j in xrange(AMINO_ACID_COUNT)]
		self.blocks_secondary_structures = blocks_sec_struct
		self.short_range = short_range
		self.total_median = 0
		self.identifier = os.path.basename(frequencies_path)
		self.load_frequencies(frequencies_path)
		if references:
			self.load_references(references)
		self.weight = 1.0
		if isconsec:
			self.type = frequency_consec_disttype
		else:
			self.type = frequency_nonconsec_disttype
		self.defaultvalue = 0
		self.refstate = True

	def read_frequency_line(self, line, tag1, tag2, all):
		ptcomps, freq = line.strip().split(";")
		coords = ptcomps.split(",")
		alpha = int(sum((int(coords[i]) + 10) * (20 ** i) for i in xrange(len(coords))))
		if all:
			self.default_alpha_frequencies[alpha] = float(freq.strip().split(" ")[0])
		else:
			self.alpha_frequencies[tag1][tag2][alpha] = float(freq)
			self.total_interactions[tag1][tag2] += float(freq)

	def alpha_frequency(self, type1, type2, zone):
		"""This helper function retrieves the frequency of 'zone' in the loaded frequency data. For this subclass, zone is a tuple (zone1, zone2). type1 refers to the type (code) of the amino acid that serves as the origin of zone1, while type2 is the amino acid for zone2."""
		coords = (zone[0].x, zone[0].y, zone[0].z, zone[1].x, zone[1].y, zone[1].z)
		alpha = int(sum((coords[i] + 10) * (20 ** i) for i in xrange(len(coords))))
		if alpha in self.alpha_frequencies[type1][type2]:
			return self.alpha_frequencies[type1][type2][alpha]
		else:
			return 0

	def subscore(self, protein, aa, aa2, onlyone=False, consec=2, zero_value=0.01):
		tag1 = aacode(aa.type)
		tag2 = aacode(aa2.type)
		if tag1 >= AMINO_ACID_COUNT: tag1 = 0
		if tag2 >= AMINO_ACID_COUNT: tag2 = 0
		
		#Switch the amino acids around
		if tag1 > tag2 or (tag1 == tag2 and aa.tag > aa2.tag):
			cache = aa
			aa = aa2
			aa2 = cache
			cache = tag1
			tag1 = tag2
			tag2 = cache
		
		if self.refstate:
			zone1 = aa.tolocal(aa2.acarbon).floor()
			zone2 = aa2.tolocal(aa.acarbon).floor()
			subscore = self.alpha_frequency(tag1, tag2, (zone1, zone2))
			if subscore == 0:
				subscore = zero_value
			if onlyone:
				print subscore
				return (subscore * self.weight, subscore * self.weight, self.total_interactions[tag1][tag2], self.total_interactions[tag2][tag1])
			vol = protein.volume()
			contact_type = 2
			if self.short_range: contact_type = 1
			elif consec == 1: contact_type = 0
			return -math.log((subscore / self.total_interactions[tag1][tag2]) * reference_state.contact_probability(tag1, tag2, vol, contact_type) / (4000.0 * math.pi / (3 * vol) * reference_state.position_ref(zone2)))
			#return -math.log((subscore / (self.total_interactions[tag1][tag2])) / reference_state.position_ref(zone2))
		else:
			zone1 = aa.tolocal(aa2.acarbon).floor()
			zone2 = aa2.tolocal(aa.acarbon).floor()
			subscore = self.alpha_frequency(tag2, tag1, (zone1, zone2))
			if subscore == 0:
				subscore = zero_value
			subscore = -math.log(subscore / self.median_frequencies[tag2][tag1] * self.total_interactions[tag2][tag1] / self.total_median)
			if onlyone:
				print subscore
				return (subscore * self.weight, self.total_interactions[tag1][tag2], self.total_interactions[tag2][tag1])
			return subscore

class SPARCSecondaryBothDistributionManager (SPARCSecondaryDistributionManager):
	"""Subclass, similar to SPARCBothOrientationDistributionManager, that handles secondary structure information where both orientations are given together."""
	
	def read_frequency_line(self, line, sec_struct_type, all):
		ptcomps, freq = line.strip().split(";")
		coords = ptcomps.split(",")
		alpha = int(sum((int(coords[i]) + 10) * (20 ** i) for i in xrange(len(coords))))
		if alpha not in self.alpha_frequencies:
			self.alpha_frequencies[alpha] = _secondary_structures_dict(inner_value=0)
		self.alpha_frequencies[alpha][sec_struct_type] = float(freq)
		self.total_interactions[sec_struct_type] += float(freq)

	def alpha_frequency(self, sec_struct_id, zone):
		"""This helper function retrieves the frequency of 'zone' in the loaded frequency data, where zone is a tuple (zone1, zone2). sec_struct_id is the string, which is also a key in the data (e.g. "helix1", "sheet0")."""
		coords = (zone[0].x, zone[0].y, zone[0].z, zone[1].x, zone[1].y, zone[1].z)
		alpha = int(sum((coords[i] + 10) * (20 ** i) for i in xrange(len(coords))))
		coords = (zone[1].x, zone[1].y, zone[1].z, zone[0].x, zone[0].y, zone[0].z)
		alpha2 = int(sum((coords[i] + 10) * (20 ** i) for i in xrange(len(coords))))
		if alpha in self.alpha_frequencies:
			return self.alpha_frequencies[alpha][sec_struct_id]
		elif alpha2 in self.alpha_frequencies:
			return self.alpha_frequencies[alpha2][sec_struct_id]
		else:
			return 0

	def subscore(self, protein, aa, aa2, sec_struct_type, onlyone=False):
		tag1 = aacode(aa.type)
		tag2 = aacode(aa2.type)
		if tag1 >= AMINO_ACID_COUNT: tag1 = 0
		if tag2 >= AMINO_ACID_COUNT: tag2 = 0
		
		if self.refstate == True:
			zone1 = aa.tolocal(aa2.acarbon).floor()
			zone2 = aa2.tolocal(aa.acarbon).floor()
			subscore = self.alpha_frequency(sec_struct_type, (zone1, zone2))
			if subscore == 0:
				if consec == 0:	subscore = 1e-3
				else:			subscore = 1e-2
			subscore = -math.log((subscore / self.total_interactions[sec_struct_type]) / reference_state.position_ref(zone2))
			if onlyone:
				print subscore
				return (subscore * self.weight, subscore * self.weight, self.total_interactions[sec_struct_type], self.total_interactions[sec_struct_type])
			return subscore
		else:
			zone1 = aa.tolocal(aa2.acarbon).floor()
			zone2 = aa2.tolocal(aa.acarbon).floor()
			subscore = self.alpha_frequency(sec_struct_type, (zone1, zone2))
			if subscore == 0:
				if consec == 0:	subscore = 1e-3
				else:			subscore = 1e-2
			subscore = -math.log(subscore / self.median_frequencies[sec_struct_type] * self.total_interactions[sec_struct_type] / self.total_median)
			if onlyone:
				print subscore
				return (subscore * self.weight, subscore * self.weight, self.total_interactions[sec_struct_type], self.total_interactions[sec_struct_type])
			return subscore
