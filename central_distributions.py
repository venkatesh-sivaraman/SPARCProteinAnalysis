"""This module contains the SPARCCentralDistributionManager class, which can manage all orientation-based terms of SPARC."""

from distributions import *
from secondary_structure import *
from memory_profiler import profile

sparc_consecutive_mode = 'consec'
sparc_short_range_mode = 'short_range'
sparc_long_range_mode = 'long_range'
sparc_secondary_mode = 'secondary'
sparc_consec_secondary_mode = 'consec_secondary'
sparc_default_mode = 'default'

class SPARCCentralDistributionManager(FrequencyDistributionManager):
	
	def __init__(self, frequencies_path, references=None):
		"""frequencies_path should be a path to a directory of alpha zones paired with frequencies for individual amino acid pairs."""
		self.alpha_frequencies = {}
		self.reference_frequencies = None #{}
		self.reference_totals = [[[0 for n in xrange(42)] for i in xrange(AMINO_ACID_COUNT)] for j in xrange(AMINO_ACID_COUNT)]
		self.total_interactions = [[[0 for n in xrange(42)] for i in xrange(AMINO_ACID_COUNT)] for j in xrange(AMINO_ACID_COUNT)]
		self.median_frequencies = [[[0 for n in xrange(42)] for i in xrange(AMINO_ACID_COUNT)] for j in xrange(AMINO_ACID_COUNT)]
		self.total_median = 0
		self.identifier = os.path.basename(frequencies_path)
		if references:
			self.load_references(references)
		self.load_frequencies(frequencies_path)
		self.weight = 1.0
		self.defaultvalue = 0
		self.refstate = True

	def __repr__(self):
		return "<Distribution Manager for '{}' data>".format(self.identifier)
	
	def alpha_frequency(self, aa, aa2, sec_name):
		"""This helper function retrieves the frequency of the orientation between aa and aa2 in the loaded frequency data. sec_name is the string type for the secondary structure shared by both amino acids, if any."""
		zone = aa.tolocal(aa2.acarbon).floor()
		
		alpha_freq = 0
		reference_freq = 0
		if zone in self.alpha_frequencies:
			data_dict = self.alpha_frequencies[zone][aacode(aa.type)][aacode(aa2.type)]
			separation = int(min(math.fabs(aa.tag - aa2.tag), 6) - 1) * 7
			if sec_name:
				struct_idx = next((i for i, ss in enumerate([None, "helix1", "helix5", "helix7", "sheet0", "sheet1", "sheet-1"]) if ss == sec_name), 0)
				separation += struct_idx
			if separation in data_dict:
				alpha_freq = float(data_dict[separation])
		if zone in self.reference_frequencies:
			data_dict = self.reference_frequencies[zone][sec_name]
			separation = int(min(math.fabs(aa.tag - aa2.tag), 6) - 1)
			if separation in data_dict:
				reference_freq = float(data_dict[separation])
		return (alpha_freq, reference_freq)

	def subscore(self, protein, aa, aa2, onlyone=False, zero_value=0.01):
		tag1 = aacode(aa.type)
		tag2 = aacode(aa2.type)
		if tag1 >= AMINO_ACID_COUNT: tag1 = 0
		if tag2 >= AMINO_ACID_COUNT: tag2 = 0
		sec_struct = protein.secondary_structure_aa(aa.tag)
		sec_struct_2 = protein.secondary_structure_aa(aa2.tag)
		if sec_struct and sec_struct_2 and sec_struct[1].start == sec_struct_2[1].start:
			sec_name = sec_struct[0].type + str(sec_struct[1].identifiers[0])
		else:
			sec_name = "default"
		
		separation = int(min(math.fabs(aa.tag - aa2.tag), 6) - 1) * 7
		if sec_name:
			struct_idx = next((i for i, ss in enumerate([None, "helix1", "helix5", "helix7", "sheet0", "sheet1", "sheet-1"]) if ss == sec_name), 0)
			separation += struct_idx
		
		if self.refstate:
			subscore, ref = self.alpha_frequency(aa2, aa, sec_name)
			if subscore == 0:
				subscore = zero_value
			if ref == 0:
				ref = zero_value
			subscore2, ref2 = self.alpha_frequency(aa, aa2, sec_name)
			if subscore2 == 0:
				subscore2 = zero_value
			if ref2 == 0:
				ref2 = zero_value
			if onlyone:
				print subscore, subscore2
				return (subscore * self.weight, subscore2 * self.weight, self.total_interactions[tag1][tag2], self.total_interactions[tag2][tag1])
			return -math.log((subscore / self.total_interactions[tag1][tag2][separation]) / (ref / self.reference_totals[sec_name][int(min(math.fabs(aa.tag - aa2.tag), 6) - 1)])) - math.log((subscore2 / self.total_interactions[tag2][tag1][separation]) / (ref2 / self.reference_totals[sec_name][int(min(math.fabs(aa.tag - aa2.tag), 6) - 1)]))
		'''else:
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
			return subscore + subscore2'''

	def score(self, protein, data, system=None, isolate=False, onlyone=False, prior=2, zero_value=0.01, mode='default'):
		"""For frequency distributions, pass in an array of hypothetical aminoacids. This implementation returns the product of the frequencies of each pairwise interaction. If isolate=True, only the amino acids in data will be considered for the energy calculation.
			Pass prior to consider ONLY the amino acid before (True) or after (False) each amino acid in data. This works best for consecutive modes."""
		score = 0.0
		taglist = {}

		consec = 2
		use_secondary = 2
		use_short_range = 2
		if mode == sparc_consecutive_mode:
			consec = 1
			use_secondary = 0
		elif mode == sparc_secondary_mode:
			consec = 1
			use_secondary = 1
		elif mode == sparc_consec_secondary_mode:
			consec = 1
		elif mode == sparc_short_range_mode:
			consec = 0
			use_short_range = 1
		elif mode == sparc_long_range_mode:
			consec = 0
			use_short_range = 0
		
		for aa in data:
			if not aa: continue
			if prior != 2:
				nearby = []
				if aa.tag > 0 and prior != False: nearby.append(protein.aminoacids[aa.tag - 1])
				if aa.tag < len(protein.aminoacids) - 1 and prior != True: nearby.append(protein.aminoacids[aa.tag + 1])
			else:
				if system and not consec and use_short_range != 0:
					nearby = system.nearby_aa(aa, protein, 10.0, consec=consec)
				else:
					nearby = protein.nearby_aa(aa, 10.0, consec=consec)
			for aa2 in nearby:
				if not aa2: continue
				sec_struct = protein.secondary_structure_aa(aa.tag)
				sec_struct_2 = protein.secondary_structure_aa(aa2.tag)
				if aa2.tag - aa.tag == 1 and aa.has_break: continue
				elif math.fabs(aa2.tag - aa.tag) > 5 and use_short_range == 1: continue
				elif math.fabs(aa2.tag - aa.tag) <= 5 and use_short_range == 0: continue
				elif use_secondary == 0 and sec_struct and sec_struct_2 and sec_struct[1].start == sec_struct_2[1].start: continue
				elif use_secondary == 1 and not (sec_struct and sec_struct_2 and sec_struct[1].start == sec_struct_2[1].start): continue
				if (aa.tag in taglist and aa2.tag in taglist[aa.tag]) or (aa2.tag in taglist and aa.tag in taglist[aa2.tag]):
					continue
				hypo = next((x for x in data if x and x.tag == aa2.tag), None)
				if hypo is not None: aa2 = hypo
				elif isolate: continue
				
				try:
					subscore = self.subscore(protein, aa, aa2, onlyone, zero_value)
				except ZeroDivisionError:
					subscore = 0.0
				
				if onlyone: return subscore
				else: score += subscore
				
				if aa.tag in taglist:
					taglist[aa.tag].append(aa2.tag)
				else:
					taglist[aa.tag] = [aa2.tag]
		return score * self.weight
					
	def read_frequency_line(self, line, tag1, tag2):
		"""Helper method for load_frequencies, intended for subclasses to easily modify the reading procedure."""
		ptcomps, freqs = line.strip().split(";")
		alpha = Point3D(*ptcomps.split(","))
		if alpha not in self.alpha_frequencies:
			self.alpha_frequencies[alpha] = [[{} for i in xrange(AMINO_ACID_COUNT)] for k in xrange(AMINO_ACID_COUNT)]
		freqs = [freq for freq in freqs.split(",") if len(freq)]
		for sep, freq in enumerate(freqs):
			freq = int(freq)
			if freq != 0:
				self.alpha_frequencies[alpha][tag1][tag2][sep] = freq
				self.total_interactions[tag1][tag2][sep] += freq
	
	def load_frequencies(self, path):
		files = os.listdir(path)
		self.total_interactions = [[[0 for k in xrange(42)] for i in xrange(AMINO_ACID_COUNT)] for j in xrange(AMINO_ACID_COUNT)]
		loading_indicator.add_loading_data(len(files))
		for n, indfile in enumerate(files):
			loading_indicator.update_progress(1)
			if indfile.find(".txt") == -1 or indfile[0] == ".": continue
			tag1, tag2 = indfile[0:-4].split('-')
			tag1 = int(tag1)
			tag2 = int(tag2)
			with open(join(path, indfile), 'r') as file:
				for line in file:
					if ";" not in line:
						if len(line.strip()) > 0:
							self.median_frequencies[tag1][tag2] = [float(y) for y in line.split(",")]
						continue
					self.read_frequency_line(line, tag1, tag2)
		#Compute median total frequency
		#s = sorted([x for list1 in self.total_interactions for x in list1])
		#self.total_median = sum(s) / float(len(s)) #s[int(len(s) / 2.0)]

	def load_references(self, path):
		files = os.listdir(path)
		percentage = 0
		self.reference_frequencies = {}
		self.reference_totals = secondary_structures_dict([0 for i in xrange(7)])
		self.reference_totals["default"] = [0 for i in xrange(7)]
		self.reference_totals["all"] = [0 for i in xrange(7)]
		loading_indicator.add_loading_data(len(files))
		for n, indfile in enumerate(files):
			loading_indicator.update_progress(1)
			if indfile.find(".txt") == -1: continue
			sec_name = indfile[0:-4]
			with open(join(path, indfile), 'r') as file:
				for line in file:
					if ";" not in line or len(line.strip()) == 0:
						continue
					ptcomps, freqs = line.strip().split(";")
					alpha = Point3D(*ptcomps.split(","))
					if alpha not in self.reference_frequencies:
						self.reference_frequencies[alpha] = secondary_structures_dict()
						self.reference_frequencies[alpha]["default"] = {}
						self.reference_frequencies[alpha]["all"] = {}
					freqs = [freq for freq in freqs.split(",") if len(freq)]
					for sep, freq in enumerate(freqs):
						if int(freq) != 0.0:
							self.reference_frequencies[alpha][sec_name][sep] = int(freq)
							self.reference_totals[sec_name][sep] += int(freq)
		print "Loaded references"

class SPARCCentralDistributionPuppet (object):
	"""The SPARCCentralDistributionPuppet class provides objects that can act like individual frequency managers, while all the time referring back to a single centralized distribution manager. Simply initialize the object with a manager object and a mode (see the top of the module), then use it by calling the score() method as you would any FrequencyDistributionManager."""
	
	def __init__(self, manager, mode='default', weight=1.0):
		self.manager = manager
		self.mode = mode
		self.weight = weight
		self.identifier = self.mode
		self.short_range = 2
		self.blocks_secondary_structures = 2
		if self.mode == sparc_consecutive_mode:
			self.type = frequency_consec_disttype
			self.blocks_secondary_structures = 1
		elif self.mode == sparc_secondary_mode or self.mode == sparc_consec_secondary_mode:
			self.type = frequency_consec_disttype
			self.blocks_secondary_structures = 0
		elif self.mode == sparc_long_range_mode:
			self.type = frequency_nonconsec_disttype
			self.short_range = 0
		elif self.mode == sparc_short_range_mode:
			self.type = frequency_nonconsec_disttype
			self.short_range = 1

	def score(self, protein, data, system=None, isolate=False, onlyone=False, prior=2, zero_value=0.01):
		"""This method funnels through to the puppet's original central manager, passing in the mode parameter."""
		return self.manager.score(protein, data, system, isolate, onlyone, prior, zero_value, self.mode) * self.weight