from aminoacids import *
from os.path import join
import os
from math import floor, fabs
import datetime
import loading_indicator

frequency_disttype = 0
frequency_nonconsec_disttype = 1
frequency_consec_disttype = 2
medium_disttype = 3
other_disttype = 4

class DistributionManager(object):
	'''Abstract base class whose concrete instances are passed to the folding simulator. Instances of DistributionManager's subclasses should provide a score for a particular amino acid positional relationship. Instances may be passed a list of amino acids, in which the first amino acid should be considered the subject of scoring. 
		
		Only interactions in which the second amino acid's tag is greater than the subject amino acid should be scored (preventing scores from being doubled).'''
	
	def __init__(self):
		self.type = other_disttype
		self.weights = []

	def score(self, protein, data):
		"""score() must be implemented by all instances, and is used to determine the feasibility of a protein conformation."""
		raise NotImplementedError("Instances of DistributionManager must implement score().")

#MARK: - Frequency
class FrequencyDistributionManager(DistributionManager):
	'''Semi-abstract subclass of DistributionManager that provides file-handling methods that parse the standard amino acid frequency rankings. Can be used for nonconsecutive and consecutive distribution data.'''

	def __init__(self, rankings_path, frequencies_path, axisradii_path=None):
		"""FrequencyDistributionManager maintains a series of datasets for retrieving the frequencies of various situations, including alpha carbon zones and axis configurations. 
			++ rankings_path should provide the path to a directory of top-ranked alpha zones and their most common axis configurations. 
			++ frequencies_path should be a path to a directory of alpha zones paired with frequencies for individual amino acid pairs (like rankings_path, except with all zones represented). Provide a file named "all.txt" inside the directory to define zone frequencies to be used for all amino acid pairs.
			++ axisradii_path should be a path to a directory with two files, x.txt and y.txt. These two files should contain regression information for predicting the frequency of an axis being found at a certain distance from the reference amino acid: zone.x, zone.y, zone.z; vertex.x of quadratic regression; vertex.y; leading coefficient."""
		super(FrequencyDistributionManager,self).__init__()
		self.default_alpha_frequencies = {}
		self.alpha_frequencies = {}
		self.axis_regressions = {}	#Format - key: alpha zone, value: [h1, k1, a1, h2, k2, a2] for quadratic model with 1 = x-axis, 2 = y-axis
		self.axis_frequencies = [[{} for i in xrange(AMINO_ACID_COUNT)] for i in xrange(AMINO_ACID_COUNT)]
		self.median_frequencies = [[0 for i in xrange(AMINO_ACID_COUNT)] for i in xrange(AMINO_ACID_COUNT)]
		print "Loading..."
		a = datetime.datetime.now()
		self.load_rankings(rankings_path)
		self.load_frequencies(frequencies_path)
		self.load_axis_regressions(axisradii_path)
		b = datetime.datetime.now()
		print "Loaded in {1:.4f} sec.".format(self, (b - a).total_seconds())
		self.weight = 1.0
		self.type = frequency_disttype
		self.defaultvalue = 0
	
	def alpha_frequency(self, type1, type2, zone):
		"""This helper function retrieves the frequency of 'zone' in the loaded frequency data. type1 refers to the type (code) of the amino acid that serves as the origin of the zone space, while type2 is the amino acid to which the zone refers."""
		if zone in self.default_alpha_frequencies:
			return self.default_alpha_frequencies[zone]
		elif zone in self.alpha_frequencies:
			return self.alpha_frequencies[zone][type1][type2]
		else:
			return 0
	
	def score(self, protein, data):
		"""For frequency distributions, pass in an array of hypothetical aminoacids."""
		score = 0.0
		consec = 2
		if self.type == frequency_consec_disttype: consec = 1
		elif self.type == frequency_nonconsec_disttype: consec = 0
		for aa in data:
			if not aa: continue
			nearby = protein.nearby_aa(aa, 10.0, consec=consec)
			for aa2 in nearby:
				hypo = next((x for x in data if x and x.tag == aa2.tag), None)
				if hypo is not None: aa2 = hypo
				tag1 = aacode(aa.type)
				tag2 = aacode(aa2.type)
				if tag1 >= AMINO_ACID_COUNT: tag1 = 0
				if tag2 >= AMINO_ACID_COUNT: tag2 = 0
				# Determine the score based on the reference amino acid's perspectives and its neighbors'.
				# Formula: S = -ln(F/F0) if F > 0, +5 otherwise
				
				zone = aa2.tolocal(aa.acarbon).floor()
				freq = self.alpha_frequency(tag2, tag1, zone)
				if freq > 0:
					subscore = -math.log(freq / self.median_frequencies[tag2][tag1])
				else:
					subscore = 5
				score += subscore
	
				zone = aa.tolocal(aa2.acarbon).floor()
				freq = self.alpha_frequency(tag1, tag2, zone)
				if freq > 0:
					subscore = -math.log(freq / self.median_frequencies[tag1][tag2])
				else:
					subscore = 5
				score += subscore
		return score * self.weight
	
	def load_frequencies(self, path):
		files = os.listdir(path)
		for indfile in files:
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
		print "Loaded frequencies"

	def load_rankings(self, path):
		files = os.listdir(path)
		for indfile in files:
			if indfile.find(".txt") == -1: continue
			(tag1, tag2) = indfile[0:-4].split('-')
			tag1 = int(tag1)
			tag2 = int(tag2)
			#alpha = self.alpha_frequencies[tag1][tag2]
			axis = self.axis_frequencies[tag1][tag2]
			with open(join(path, indfile), 'r') as file:
				current_alpha = None
				#total_count = 0
				for line in file:
					if "(" not in line:
						'''if len(line) > 1:
							self.total_count[tag1][tag2] = int(line)
							break
						else:'''
						continue
					comps = line.translate(None, "()").split(";")
					if "\t" not in line:		# New alpha carbon zone
						current_alpha = Point3D.list(comps[0].split(","))
						#alpha[current_alpha] = int(comps[1])
						axis[current_alpha] = {}
					else:
						x_axis = Point3D.list(comps[1].split(","))
						y_axis = Point3D.list(comps[2].split(","))
						z_axis = Point3D.list(comps[3].split(","))
						axis[current_alpha][PositionZone(current_alpha, x_axis, y_axis, z_axis)] = int(comps[4])
		print "Loaded rankings"

	def load_axis_regressions(self, path):
		if path is None:
			print "No path specified for axis radius regressions."
			return
		files = os.listdir(path)
		for indfile in files:
			if indfile == "x.txt": startidx = 0
			elif indfile == "y.txt": startidx = 3
			else: continue

			with open(join(path, indfile), 'r') as file:
				for line in file:
					alpha, h, k, a = line.strip().split(";")
					alpha_zone = Point3D(*alpha.split(","))
					if alpha_zone in self.axis_regressions:
						self.axis_regressions[alpha_zone][startidx] = float(h)
						self.axis_regressions[alpha_zone][startidx + 1] = float(k)
						self.axis_regressions[alpha_zone][startidx + 1] = float(a)
					else:
						item = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
						item[startidx] = float(h)
						item[startidx + 1] = float(k)
						item[startidx + 1] = float(a)
						self.axis_regressions[alpha_zone] = item
		print "Loaded axis regressions"

#MARK: - Medium

class MediumDistributionManager(DistributionManager):
	'''Concrete subclass of DistributionManager that provides file-handling methods that parse medium amino acid density data (number of amino acids within 10 angstroms of each type).'''
	
	def __init__(self, frequencies_path):
		"""MediumDistributionManager maintains a dataset sourced from frequencies_path that contains the (scaled) frequency of each amino acid having a particular number of other amino acids around it. The format should be a directory of files, one for each amino acid type, with the following structure in each file (coords need not be consecutive):
			1,(frequency of 1 amino acid)
			2,(frequency of 2 amino acids)
			.
			.
			.
			64,(frequency of 64 amino acids)
			67,(frequency of 67 amino acids)
			(median frequency)"""
		super(MediumDistributionManager,self).__init__()
		self.weight = 1.0
		self.type = medium_disttype
		self.defaultvalue = 0
		self.frequencies = [{} for i in xrange(AMINO_ACID_COUNT)]
		self.median_frequencies = [0.0 for i in xrange(AMINO_ACID_COUNT)]
		self.load_frequencies(frequencies_path)
		self.identifier = os.path.basename(frequencies_path)

	def load_frequencies(self, frequencies_path):
		files = os.listdir(frequencies_path)
		loading_indicator.add_loading_data(len(files))
		for indfile in files:
			loading_indicator.update_progress(1)
			if indfile.find(".txt") == -1: continue
			idx = int(indfile[0:-4])
			with open(join(frequencies_path, indfile), 'r') as file:
				for line in file:
					line = line.strip()
					if len(line) == 0: continue
					if "," not in line and " " not in line:
						if self.median_frequencies[idx] == 0.0:
							self.median_frequencies[idx] = float(line)
					else:
						if "," in line: comps = line.split(",")
						else: comps = line.split()
						coord = int(comps[0])
						self.frequencies[idx][coord] = float(comps[1])
		for idx in xrange(AMINO_ACID_COUNT):
			self.median_frequencies[idx] = sum(self.frequencies[idx])# / float(len([f for f in self.frequencies[idx] if f != 0]))

	def score(self, protein, data):
		"""Pass in an array of amino acids for data."""
		score = 0.0
		density = 1.0 / (1.410 + 0.145 * math.exp(-protein.mass / 13.0)) # 0.73
		volume = (density * 1e24) / (6.02e23) * protein.mass
		if volume <= 4000.0 / 3.0 * math.pi:
			return score
		for aa in data:
			if not aa: continue
			tag = aacode(aa.type)
			if tag >= AMINO_ACID_COUNT: tag = 0
			coord = len(protein.nearby_aa(aa, 10.0))
			freqs = self.frequencies[tag]
			if coord in freqs:
				f = freqs[coord]
			elif len(freqs) > 0:
				try: less = max(i for i in freqs if freqs[i] <= coord)
				except ValueError: less = -1
				try: more = min(i for i in freqs if freqs[i] >= coord)
				except ValueError: more = -1
				if less != -1 and more != -1:
					f = (freqs[more] - freqs[less]) / (more - less) * (coord - less) + freqs[less]
				if less == -1 and more == -1:
					f = 0.0
				elif less == -1:
					f = more
				elif more == -1:
					f = less
			else:
				f = 0.0
			#f0 = self.median_frequencies[tag]
			f = f / self.median_frequencies[tag]
			N = len(protein.aminoacids)
			f0 = ((1 - (4000 * math.pi) / (3 * volume)) ** (N - 1)) / (((3 * volume) / (4000 * math.pi) - 1) ** coord)
			f0 *= float(math.factorial(N - 1) / (math.factorial(coord) * math.factorial(N - coord - 1)))
			if f <= 0: f = 0.01
			subscore = -math.log(f / f0)
			score += subscore
		return score * self.weight


#MARK: - Helper

def place_element_sorted(array, element, sorter=None, min=0, max=0):
	"""Performs binary search to place the element in its proper place in the sorted array."""
	if sorter is None:
		def sort(a, b):
			if a > b: return 1
			elif a == b: return 0
			else: return -1
		sorter = sort
	
	if max == 0: max = len(array)
	
	if max - min == 2:
		value = sorter(array[min], element)
		if value >= 0: array.insert(min, element)
		else:
			value = sorter(array[min + 1], element)
			if value >= 0: array.insert(min + 1, element)
			else: array.insert(min + 2, element)
	elif max - min == 1:
		value = sorter(array[min], element)
		if value >= 0: array.insert(min, element)
		else: array.insert(min + 1, element)
	elif max - min == 0:
		array.append(element)
	else:
		mid = int(floor((max + min) / 2))
		value = sorter(array[mid], element)
		if value > 0: place_element_sorted(array, element, sorter=sorter, min=min, max=mid)
		elif value == 0: array.insert(mid + 1, element)
		else: place_element_sorted(array, element, sorter=sorter, min=mid, max=max)

def place_elements_sorted(array, elements, sorter=None):
	"""Places elements in array so that the array is sorted according to 'sorter'. Assumes that array is already sorted."""
	if sorter is None:
		def sort(a, b):
			if a > b: return 1
			elif a == b: return 0
			else: return -1
		sorter = sort
	for x in elements:
		place_element_sorted(array, x, sorter=sorter)