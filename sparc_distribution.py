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
		"""frequencies_path should be a path to a directory of alpha zones paired with frequencies for individual amino acid pairs. If blocks_sec_struct is True, then secondary structures will be ignored. To treat secondary structures exclusively, use a SPARCSecondaryDistributionManager."""
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

	def __repr__(self):
		return "<Distribution Manager for '{}' data>".format(self.identifier)
	
	def alpha_frequency(self, type1, type2, zone):
		"""This helper function retrieves the frequency of 'zone' in the loaded frequency data. type1 refers to the type (code) of the amino acid that serves as the origin of the zone space, while type2 is the amino acid to which the zone refers."""
		if zone in self.alpha_frequencies:
			if self.reference_frequencies:
				ref_zone = zone
				if ref_zone.x < 0.0: ref_zone.x = -ref_zone.x - 1.0
				if ref_zone.y < 0.0: ref_zone.y = -ref_zone.y - 1.0
				if ref_zone.z < 0.0: ref_zone.z = -ref_zone.z - 1.0
				if ref_zone in self.reference_frequencies:
					return (self.alpha_frequencies[zone][type1][type2], self.reference_frequencies[ref_zone][type1][type2])
			return self.alpha_frequencies[zone][type1][type2]
		else:
			return 0
	
	def score(self, protein, data, isolate=False, onlyone=False, prior=2):
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
				tag1 = aacode(aa.type)
				tag2 = aacode(aa2.type)
				if tag1 >= AMINO_ACID_COUNT: tag1 = 0
				if tag2 >= AMINO_ACID_COUNT: tag2 = 0
				
				zone = aa2.tolocal(aa.acarbon).floor()
				subscore = self.alpha_frequency(tag2, tag1, zone)
				if subscore == 0:
					if consec == 0:	subscore = 1e-3
					else:			subscore = 1e-2
				subscore = -math.log(subscore / self.median_frequencies[tag2][tag1] * self.total_interactions[tag2][tag1] / self.total_median)
				zone2 = aa.tolocal(aa2.acarbon).floor()
				subscore2 = self.alpha_frequency(tag1, tag2, zone2)
				if subscore2 == 0:
					if consec == 0:	subscore2 = 1e-3
					else:			subscore2 = 1e-2
				subscore2 = -math.log(subscore2 / self.median_frequencies[tag1][tag2] * self.total_interactions[tag1][tag2] / self.total_median)
				'''if self.short_range == 1:
					if self.alpha_frequency(tag2, tag1, zone) == 0 or self.alpha_frequency(tag1, tag2, zone2) == 0:
						print "SR1:", aa2.tag - aa.tag, zone, self.alpha_frequency(tag2, tag1, zone), subscore
						print "SR2:", aa2.tag - aa.tag, zone2, self.alpha_frequency(tag1, tag2, zone2), subscore2'''
				score += subscore + subscore2
				if onlyone:
					print subscore, subscore2
					return (subscore * self.weight, subscore2 * self.weight, self.total_interactions[tag1][tag2], self.total_interactions[tag2][tag1])
				if aa.tag in taglist:
					taglist[aa.tag].append(aa2.tag)
				else:
					taglist[aa.tag] = [aa2.tag]
					
				'''zone = aa2.tolocal(aa.acarbon).floor()
				subscore = self.alpha_frequency(tag2, tag1, zone)
				if subscore == 0:
					subscore = 1e-10
				zone2 = aa.tolocal(aa2.acarbon).floor()
				subscore2 = self.alpha_frequency(tag1, tag2, zone2)
				if subscore2 == 0:
					subscore2 = 1e-10
				subscore = -math.log(subscore / self.median_frequencies[tag2][tag1] * self.total_interactions[tag2][tag1] / self.total_median)
				subscore2 = -math.log(subscore2 / self.median_frequencies[tag1][tag2] * self.total_interactions[tag1][tag2] / self.total_median)
				#if not consec:
				#	print subscore * subscore2, self.total_interactions[tag2][tag1] ** 2, reference_state.position_ref(zone), -math.log((subscore * subscore2 / self.total_interactions[tag2][tag1] ** 2) / reference_state.position_ref(zone))
				#subscore = -math.log((subscore * subscore2 / (self.total_interactions[tag2][tag1] ** 2)) / reference_state.position_ref(zone))
				score += subscore + subscore2
				'''
		return score * self.weight
	
	def load_frequencies(self, path):
		files = os.listdir(path)
		percentage = 0
		self.total_interactions = [[0 for i in xrange(AMINO_ACID_COUNT)] for j in xrange(AMINO_ACID_COUNT)]
		loading_indicator.add_loading_data(len(files))
		for n, indfile in enumerate(files):
			loading_indicator.update_progress(1)
			if indfile.find(".txt") == -1: continue
			all = ("all" in indfile)
			if not all:	tag1, tag2 = indfile[0:-4].split('-')
			tag1 = int(tag1)
			tag2 = int(tag2)
			with open(join(path, indfile), 'r') as file:
				for line in file:
					if ";" not in line and len(line.strip()) > 0:
						self.median_frequencies[tag1][tag2] = float(line)
						continue
					ptcomps, freq = line.strip().split(";")
					alpha = Point3D(*ptcomps.split(","))
					if all:
						self.default_alpha_frequencies[alpha] = float(freq.strip().split(" ")[0])
					else:
						if alpha not in self.alpha_frequencies:
							self.alpha_frequencies[alpha] = [[0 for i in xrange(AMINO_ACID_COUNT)] for k in xrange(AMINO_ACID_COUNT)]
						self.alpha_frequencies[alpha][tag1][tag2] = float(freq)
						self.total_interactions[tag1][tag2] += float(freq)
		#Compute median total frequency
		s = sorted([x for list1 in self.total_interactions for x in list1])
		self.total_median = s[int(len(s) / 2.0)] #sum(s) / float(len(s))

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
	
	def alpha_frequency(self, sec_struct_id, zone):
		"""This helper function retrieves the frequency of 'zone' in the loaded frequency data. sec_struct_id is the string, which is also a key in the data (e.g. "helix1", "sheet0")."""
		if zone in self.alpha_frequencies:
			return self.alpha_frequencies[zone][sec_struct_id]
		else:
			return 0
	
	def score(self, protein, data, isolate=False, onlyone=False, prior=2):
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
				nearby = protein.nearby_aa(aa, 10.0, consec=consec)
			for aa2 in nearby:
				if (aa.tag in taglist and aa2.tag in taglist[aa.tag]) or (aa2.tag in taglist and aa.tag in taglist[aa2.tag]):
					continue
				hypo = next((x for x in data if x and x.tag == aa2.tag), None)
				if hypo is not None: aa2 = hypo
				elif isolate: continue
				if aa2.tag - aa.tag == 1 and aa.has_break: continue
				tag1 = aacode(aa.type)
				tag2 = aacode(aa2.type)
				if tag1 >= AMINO_ACID_COUNT: tag1 = 0
				if tag2 >= AMINO_ACID_COUNT: tag2 = 0
				
				zone = aa2.tolocal(aa.acarbon).floor()
				subscore = self.alpha_frequency(sec_struct_type, zone)
				if subscore == 0:
					subscore = 0.001
				subscore = -math.log(subscore / self.median_frequencies[sec_struct_type] * self.total_interactions[sec_struct_type] / self.total_median)
				zone = aa.tolocal(aa2.acarbon).floor()
				subscore2 = self.alpha_frequency(sec_struct_type, zone)
				if subscore2 == 0:
					subscore2 = 0.001
				subscore2 = -math.log(subscore2 / self.median_frequencies[sec_struct_type] * self.total_interactions[sec_struct_type] / self.total_median)
				score += subscore + subscore2
				if onlyone:
					print subscore, subscore2
					return (subscore * self.weight, subscore2 * self.weight, self.total_interactions[sec_struct_type], self.total_interactions[sec_struct_type])
				if aa.tag in taglist:
					taglist[aa.tag].append(aa2.tag)
				else:
					taglist[aa.tag] = [aa2.tag]
		return score * self.weight
	
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
					ptcomps, freq = line.strip().split(";")
					alpha = Point3D(*ptcomps.split(","))
					if alpha not in self.alpha_frequencies:
						self.alpha_frequencies[alpha] = _secondary_structures_dict(inner_value=0)
					self.alpha_frequencies[alpha][sec_struct_type] = float(freq)
					self.total_interactions[sec_struct_type] += float(freq)
		#Compute median total frequency
		s = sorted(self.total_interactions.values())
		self.total_median = sum(s) / float(len(s)) #s[int(len(s) / 2.0)]
