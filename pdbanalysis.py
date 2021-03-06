from polypeptide import *
from sparc_distribution import *
import os, gc, time
from os.path import join
import numpy
import random
from scipy.stats import linregress
from scipy.spatial import ConvexHull
import planeregression as plreg
from functools import partial
import multiprocessing
import shutil, subprocess

#MARK: Aggregation and helpers

def aggregate_aadata(source, dest, two_d=True):
	if not os.path.exists(dest): os.mkdir(dest)
	contents = os.listdir(source)

	def _aggregate_mid(ksource, ksources, kfile, kfileName):
		for folder in ksources:
			sourcePath = join(join(ksource, folder), kfileName)
			if not os.path.exists(sourcePath): continue
			sfile = open(sourcePath, 'r')
			for line in sfile:
				file.write(line)
			sfile.close()
	
	for i in xrange(AMINO_ACID_COUNT):
		if two_d:
			for j in xrange(AMINO_ACID_COUNT):
				destPath = "%d-%d.txt" % (i, j)
				file = open(join(dest, destPath), 'w')
				_aggregate_mid(source, contents, file, destPath)
				file.close()
				print "Saved %r to file" % destPath
		else:
			destPath = "%d.txt" % i
			file = open(join(dest, destPath), 'w')
			_aggregate_mid(source, contents, file, destPath)
			file.close()
			print "Saved %r to file" % destPath

def read_alpha_zones(file):
	pzs = []
	for line in file:
		if not len(line.strip()) or ";" not in line: continue
		(x, y, z) = line.split(";")[0].split(",")
		pzs.append(Point3D(x, y, z))
	return pzs

def read_positional_zones(file, azone=Point3D.zero()):
	pzs = []
	file.seek(0)
	for line in file:
		i = 0
		zp = PositionZone()
		for pointstr in line.split(";"):
			try:
				(x, y, z) = pointstr.split(",")
				if i == 0:
					zp.alpha_zone = Point3D(x, y, z)
					if azone != Point3D.zero() and zp.alpha_zone != azone:
						i = -1
						break
				elif i == 1: zp.x_axis = Point3D(x, y, z)
				elif i == 2: zp.y_axis = Point3D(x, y, z)
				elif i == 3: zp.z_axis = Point3D(x, y, z)
				else: assert False, "More components than necessary: %r" % line
			except ValueError:
				print "Exception for", pointstr
			i += 1
		if i != -1: pzs.append(zp); zp.hash = zp.calchash()
	return pzs

def read_axis_zones(file, freqdata, azones, axis=0):
	"""Provide 0 for x-axis, 1 for y-axis, 2 for z-axis. Returns a list of points."""
	file.seek(0)
	for z in azones:
		if z not in freqdata: freqdata[z] = {}
	for line in file:
		try:
			comps = line.split(";")
			if len(comps) < 4: continue
			(x, y, z) = comps[0].split(",")
			azone = Point3D(x, y, z)
			if azone not in azones: continue
			if axis == 0:
				(x1, y1, z1) = comps[1].split(",")
			elif axis == 1:
				(x1, y1, z1) = comps[2].split(",")
			elif axis == 2:
				(x1, y1, z1) = comps[3].split(",")
			else: assert False, "Invalid axis number (must be 0, 1, or 2): %d" % axis
			pz = Point3D(x1, y1, z1)
			if pz in freqdata[azone]: freqdata[azone][pz] += 1
			else: freqdata[azone][pz] = 1
		except ValueError:
			print "Exception for", pointstr

def read_alpha_frequencies(file):
	"""Reads out the alpha frequencies and returns them as a dictionary {zone : frequency }."""
	dist = {}
	for line in file:
		line = line.strip()
		if len(line) == 0 or "," not in line: continue
		comps = line.split(";")
		x, y, z = comps[0].split(",")
		comps[1] = comps[1].strip()
		if " " in comps[1]: comps[1] = comps[1][:comps[1].find(" ")]
		if float(comps[1]) == 0: continue
		dist[Point3D(float(x), float(y), float(z))] = float(comps[1])
	return dist


def positional_zone_rankings(input, output):
	"""This function writes to 'output' the most common position zone rankings found in the file specified by 'input'. The output of this function (most common alpha zones, followed by common axis positions) is to be the standard format for frequency distributions."""
	print os.path.basename(input)
	file = open(input, 'r')
	pzs = read_alpha_zones(file)
	
	frequencies = {}
	for pz in pzs:
		if pz in frequencies:
			frequencies[pz] += 1
		else:
			frequencies[pz] = 1
	file.close()

	rankings = sorted(frequencies.items(), key=lambda x: x[1], reverse=True)
	del frequencies
	writefile = open(output, 'w')
	for x in rankings[0:20]:
		print str(x[0]) + ";", x[1]
		writefile.write(str(x[0]) + "; " + str(x[1]) + "\n")
		pzs = read_positional_zones(file, azone=x[0])
		frequencies = {}
		for pz in pzs:
			if pz in frequencies:
				frequencies[pz] += 1
			else:
				frequencies[pz] = 1
		rankings2 = sorted(frequencies.items(), key=lambda x: x[1], reverse=True)
		for x in rankings2[0:10]:
			writefile.write("\t" + str(x[0]) + "; " + str(x[1]) + "\n")
	writefile.close()
	del rankings
	print "Done"

def list_position_zone_frequencies(input, output):
	"""This function writes to 'output' ALL position zone frequencies found in the file specified by 'input'."""
	print os.path.basename(input)
	file = open(input, 'r')
	pzs = read_alpha_zones(file)
	
	frequencies = {}
	for pz in pzs:
		if pz in frequencies:
			frequencies[pz] += 1
		else:
			frequencies[pz] = 1
	file.close()
	writefile = open(output, 'w')

	i = 0
	for x in xrange(-10, 10):
		for y in xrange(-10, 10):
			for z in xrange(-10, 10):
				pt = Point3D(x, y, z)
				'''if math.fabs(z) > 20   - math.fabs(x) - math.fabs(y):
					if pt in frequencies:
						print "Nonzero frequency for %r: %d" % (str(pt), frequencies[pt])
					continue
				print str(pt), i'''
				if pt in frequencies:
					freq = frequencies[pt]
					writefile.write(str(freq) + "\n")
					#print str(pt), freq
				else:
					writefile.write("0\n")
				i += 1
	writefile.close()
	return

	def sorting_function(a, b):
		if a.x != b.x: return int(a.x - b.x)
		else:
			if a.y != b.y: return int(a.y - b.y)
			else: return int(a.z - b.z)
	data = sorted(frequencies.items(), cmp=sorting_function, key=lambda x: x[0])
	
	for (pz, freq) in data:
		writefile.write(str(freq) + "\n")
	writefile.close()

def write_median_frequencies(input):
	"""Writes the median frequency to the end of the file for each alpha zone file in input."""
	for path in os.listdir(input):
		if ".txt" not in path: continue
		fncomps = path[:path .find(".txt")].split("-")
		print path[:path.find(".txt")]
		dist = []
		with open(os.path.join(input, path), 'r') as file:
			for line in file:
				line = line.strip()
				if len(line) == 0 or ";" not in line: continue
				comps = line.split(";")
				nums = [float(x) for x in comps[1].split(",")]
				i = 0
				while i < len(dist):
					dist[i].append(nums[i])
					i += 1
				while i < len(nums):
					dist.append([nums[i]])
					i += 1
		means = []
		interaction_medians = []
		medians = []
		for subdist in dist:
			s = sorted(subdist)
			if len(s) > 0:
				median = s[int(len(s) / 2.0)]
				interaction_median = sum(s) / 2.0
				runner = 0.0
				for x in s:
					runner += x
					if runner >= interaction_median:
						interaction_median = x
						break
				mean = sum(s) / 5111.0 #float(len(s))
			else:
				median = 0.0
				interaction_median = 0.0
				mean = 0.0
			means.append(mean)
			interaction_medians.append(interaction_median)
			medians.append(median)
		with open(os.path.join(input, path), 'r') as file:
			text = file.readlines()
		with open(os.path.join(input, path), 'w') as file:
			for line in text:
				if len(line.split(";")) > 1:
					file.write(line)
			def list_str(numbers_string):
				ret = ""
				for num in numbers_string:
					ret += str(num) + ","
				return ret[:-1]
			file.write(list_str(interaction_medians) + "\n" + list_str(medians) + "\n" + list_str(means) + "\n")

def write_mean_hydrophobicity_dist(input):
	"""Writes the mean frequency to the end of the file for each medium file in input. (Deprecated; use format_hydrophobicity_dist instead.)"""
	for path in os.listdir(input):
		if ".txt" not in path: continue
		idx1 = int(path[:path.find(".txt")])
		print idx1
		dist = []
		with open(os.path.join(input, path), 'r') as file:
			for line in file:
				line = line.strip()
				if len(line) == 0 or "," not in line: continue
				comps = line.split(",")
				dist.append(float(comps[1]))
		s = sorted(dist)
		if len(s) > 0:
			#median = s[int(len(s) / 2.0)]
			mean = sum(s) / float(len(s))
		else:
			#median = 0.0
			mean = 0.0
		#print median, mean
		print mean
		with open(os.path.join(input, path), 'a') as file:
			file.write("\n" + str(mean) + "\n")


#MARK: - Frequency function
'''The following functions were used as attempts to determine a function relating zone and frequency based on the existing frequency tables.'''

def compute_neighbor_zones(point, offset=1):
	neighbor_zones = []
	for x in xrange(-offset, offset + 1):
		for y in xrange(-offset, offset + 1):
			for z in xrange(-offset, offset + 1):
				if x > -offset and x < offset and y > -offset and y < offset and z > -offset and z < offset: continue
				else: neighbor_zones.append(Point3D(point.x + x, point.y + y, point.z + z))
	return neighbor_zones

def iter_neighbor_zones(point, offset=1):
	neighbor_zones = []
	for x in xrange(-offset, offset + 1):
		for y in xrange(-offset, offset + 1):
			for z in xrange(-offset, offset + 1):
				if x > -offset and x < offset and y > -offset and y < offset and z > -offset and z < offset: continue
				else: yield Point3D(point.x + x, point.y + y, point.z + z)

def compute_box_zones(point, offset=1):
	neighbor_zones = []
	for x in xrange(-offset, offset + 1):
		for y in xrange(-offset, offset + 1):
			for z in xrange(-offset, offset + 1):
				neighbor_zones.append(Point3D(point.x + x, point.y + y, point.z + z))
	return neighbor_zones

def coeff_of_determination(data, function):
	mean = numpy.mean([x[1] for x in data])
	totsum = sum([(x[1] - mean) ** 2 for x in data])
	ressum = sum([(x[1] - function(x[0])) ** 2 for x in data])
	return 1.0 - ressum / totsum

def average_zp_frequency_change(input, zoneoffset=1, dist=None):
	"""This function returns the average frequency percentage change that is undergone when moving 'zoneoffset' units away from a high-ranking position zone in any direction. Pass in a path to a file containing position zones and the relevant distribution."""
	print os.path.basename(input)
	file = open(input, 'r')
	pzs = read_alpha_zones(file)
	frequencies = {}
	for pz in pzs:
		if pz in frequencies:
			frequencies[pz] += 1
		else:
			frequencies[pz] = 1
	file.close()

	maxima = {}
	neighbor_zones = []
	border_zones = []
	total_average = 0.0
	for (point, frequency) in frequencies.iteritems():
		if frequency < 10: continue
		neighbor_zones = compute_neighbor_zones(point)
		maximum = True
		total_difference = 0.0
		for neighbor in neighbor_zones:
			if neighbor in frequencies:
				if frequencies[neighbor] > frequency:
					maximum = False
					break
		if maximum == True:
			border_zones = compute_neighbor_zones(point, offset=zoneoffset)
			for neighbor in border_zones:
				if neighbor in frequencies:
					total_difference += float(frequency - frequencies[neighbor]) / float(frequency) * 100.0
				else: total_difference += 100.0
			average = float(total_difference) / len(border_zones)
			maxima[point] = (frequency, average)
			total_average += average
			#This code computes linear regressions for distance data.
			'''border_zones = compute_box_zones(point, offset=2)
			distancedata = [(pz.distanceto(point), frequencies[pz]) for pz in border_zones if pz in frequencies]
			for (distance, freq) in distancedata:
				print distance, freq
			slope, intercept, r_value, p_value, std_err = linregress(distancedata)
			print "y = {}(1/d)+{}, r^2={}".format(slope, intercept, r_value ** 2)'''

		del neighbor_zones[:]
		del border_zones[:]
	maxima = sorted(maxima.items(), key=lambda x: x[1][0], reverse=True)
	for (key, value) in maxima:
		print str(key) + ", %d, %.6f" % value
	return total_average / float(len(maxima))

def positional_zone_stats(input, dist=None, quantiledata=None):
	"""This function logs the distribution of distances from the maximum frequency positional zone. It's not extremely reliable for data because the max positional zone itself is defined only approximately, with the result that apparently there are 0 cases that actually match the max zone.
		
		Input: a path to a single file containing a list of alpha carbon zones. Pass in a FrequencyDistributionManager to use the max zone as the center point; otherwise, the origin will be used. Pass in quantiledata to log the distance-quantile relationship."""
	print os.path.basename(input)
	file = open(input, 'r')
	pzs = read_alpha_zones(file)
	file.close()

	max = None
	if dist is None:
		'''frequencies = {}
		for pz in pzs:
			if pz in frequencies:
				frequencies[pz] += 1
			else:
				frequencies[pz] = 1

		rankings = sorted(frequencies.items(), key=lambda x: x[1], reverse=True)[:20]
		max = rankings[0][0]
		del frequencies'''
		max = Point3D.zero()
	else:
		tag1, tag2 = os.path.basename(input).replace(".txt", "").split("-")
		max = sorted(dist.alpha_frequencies[int(tag1)][int(tag2)].items(), key=lambda x: x[1], reverse=True)[0][0]
		print str(max)
	#for (point, frequency) in rankings[1:]:
	#	print point.distanceto(max)

	#Set up a data table with distance to the max as the IV (rounded to array indexes) and frequency as the DV.
	cutoffdist = 3
	frequencies = {"%.2f" % i: 0 for i in numpy.arange(0,cutoffdist,0.05)}
	for pz in pzs:
		distance = math.floor(pz.distanceto(max) * 20.0) / 20.0
		if quantiledata is not None:
			print pz.distanceto(max), quantiledata[pz]
		if distance < cutoffdist:
			frequencies["%.2f" % distance] += 1

	#Determine the mean, median, and mode of the frequencies based on the distance from each point to the maximum. They should be centered at 0.
	sum = 0
	count = 0
	frequencies = sorted(frequencies.items(), key=lambda x: x[0])
	for (distance, frequency) in frequencies:
		print distance + ", " + str(frequency)
		sum += float(distance) * frequency
		count += frequency
	print "Mean:", float(sum) / float(count)
	
	file.close()
	print "Done"

def median_frequency(input):
	print os.path.basename(input)
	file = open(input, 'r')
	pzs = read_alpha_zones(file)
	file.close()
	frequencies = {}
	for pz in pzs:
		if pz in frequencies:
			frequencies[pz] += 1
		else:
			frequencies[pz] = 1

	frequencies = sorted(frequencies.values())
	if len(frequencies) == 0: return -1
	return frequencies[len(frequencies) / 2]

def compare_dist_frequencies(distribution, input):
	"""This function evaluates the difference between predicted frequencies for a particular pair of amino acids and actual frequencies determined from data. Pass in a FrequencyDistributionManager and a path to a file containing positional zones."""
	print os.path.basename(input)
	file = open(input, 'r')
	pzs = read_alpha_zones(file)
	frequencies = {}
	for pz in pzs:
		if pz in frequencies:
			frequencies[pz] += 1
		else:
			frequencies[pz] = 1
	file.close()

	tag1, tag2 = os.path.basename(input).replace(".txt", "").split("-")
	total_difference = 0
	for pz, freq in frequencies.iteritems():
		total_difference += math.fabs(freq - distribution.score((int(tag1), int(tag2), pz)))
	print "Average difference: %.3f" % (float(total_difference) / len(frequencies))

def frequency_vs_zonepairs(input):
	"""This function determines the *number of zone pairs* for each given frequency. 'data' will be output as a list of tuples, with data[x][0] being the frequency and [1] being the number of zone pairs."""
	files = os.listdir(input)
	avg_frequencies = []
	total_freqcount = 0
	for path in files:
		if path.find(".txt") == -1: continue
		print os.path.basename(path)
		file = open(join(input, path), 'r')
		frequencies = read_alpha_frequencies(file)
		file.close()
		
		if len(frequencies) == 0: continue
		maximum = max(frequencies.values())
		print maximum
		data = [0 for i in xrange(maximum + 1)]
		for x in xrange(-10, 10):
			for y in xrange(-10, 10):
				for z in xrange(-10, 10):
					#if Point3D.zero().distanceto(Point3D(math.floor(math.fabs(x)), math.floor(math.fabs(y)), math.floor(math.fabs(z)))) > 10.0: continue #To stop skewing the data for samples we don't have
					if Point3D(x, y, z) in frequencies:
						freq = frequencies[Point3D(x, y, z)]
						data[freq] += 1
					else: data[0] += 1
		total_count = 0
		total_frequency = 0
		cutoff = 20
		total_freqcount += 1
		for i, zpcount in enumerate(data):
			#if zpcount != 0: print i, zpcount
			if len(avg_frequencies) > i: avg_frequencies[i] += zpcount
			elif len(avg_frequencies) == i: avg_frequencies.append(zpcount)
			else:
				while len(avg_frequencies) < i: avg_frequencies.append(0)
				avg_frequencies.append(zpcount)
			if zpcount > cutoff:
				total_count += zpcount
				total_frequency += i
		print "Average baseline frequency:", float(total_frequency) / float(total_count)
		del data[:]
	lastvalue = -1
	for i, zpcount in enumerate(avg_frequencies):
		if zpcount > 0:
			value = float(zpcount) / float(total_freqcount)
			if value != lastvalue:
				print i, value
				lastvalue = float(zpcount) / float(total_freqcount)

def _calculate_quantiles(frequencies, n=10):
	quantiles = []
	distribution = sorted(frequencies.items(), key=lambda x: x[1])
	step = int(len(distribution) / n)
	for i in xrange(1,n):
		quantile = distribution[i * step][1]
		quantiles.append(quantile)
	return quantiles

def quantiles(input, output=None, n=10, use_frequencies=True):
	"""Computes quantiles of frequency distribution into n intervals. Returns two values: a list of all zones found and their quantiles, and a list of the quantile values. Writes the list of the quantile values to 'output' if it is non-nil, otherwise prints a summary to the console. Pass use_frequencies=True to use the methods for raw alpha pz data."""
	file = open(input, 'r')
	if use_frequencies:
		frequencies = read_alpha_frequencies(file)
	else:
		pzs = read_alpha_zones(file)
		frequencies = {}
		for pz in pzs:
			if pz in frequencies:
				frequencies[pz] += 1
			else:
				frequencies[pz] = 1
	file.close()

	if len(frequencies) == 0: return {}

	#Now compute the values of each quantile by sorting the frequencies dictionary by increasing frequency
	'''quantiles = _calculate_quantiles(frequencies, n)
	if output is not None:
		file = open(output, 'w')
		for i, quantile in enumerate(quantiles):
			print i, quantile
			file.write(str(quantile) + "\n")
		file.close()
	else:
		retstr = str(len(frequencies))#str(len([x for x in frequencies if frequencies[x] > 50]))
		for quantile in quantiles:
			retstr += "\t" + str(quantile)
		print retstr'''

	#Now compute the quantile that each position zone belongs in based on its frequency.
	'''data = {}
	for x in xrange(-10, 10):
		for y in xrange(-10, 10):
			for z in xrange(-10, 10):
				zone = Point3D(x, y, z)
				if zone in frequencies:
					freq = frequencies[zone]
					try:
						data[zone] = next(key for key, value in enumerate(quantiles) if value > freq)
					except StopIteration:
						data[zone] = n - 1
					#if math.fabs(z) <= 0.01:
					#	print x, y, data[zone]'''
	return len(frequencies)

def _max_zones(flagged_data, startzone):
	flagged_data[startzone][1] = 1
	zones = [(startzone.x, startzone.y, startzone.z)]
	quantile = flagged_data[startzone][0]
	for zone in iter_neighbor_zones(startzone):
		if math.fabs(zone.x) > 10.0 or math.fabs(zone.y) > 10.0 or math.fabs(zone.z) > 10.0: continue
		if zone in flagged_data and flagged_data[zone][0] == quantile and flagged_data[zone][1] == 0:
			zones += _max_zones(flagged_data, zone)
	return zones

def same_quantile_cluster_size(data):
	"""This function determines how many consecutive zones have the same quantile. The 'data' is of the format returned by quantiles()."""
	
	flagged_data = { key: [value, 0] for key, value in data.iteritems() }
	for key, list in flagged_data.iteritems():
		if list[1] == 0:
			zones = _max_zones(flagged_data, key)
			if len(zones) > 330:
				print "%d consecutive zones for quantile %d" % (len(zones), flagged_data[key][0])
				for zone in zones: print zone[0], zone[1], zone[2]
				minx, maxx = min(x[0] for x in zones), max(x[0] for x in zones)
				miny, maxy = min(x[1] for x in zones), max(x[1] for x in zones)
				minz, maxz = min(x[2] for x in zones), max(x[2] for x in zones)
				print "Cube from (%d-%d, %d-%d, %d-%d)" % (minx, maxx, miny, maxy, minz, maxz)
				excess = math.fabs(maxx - minx + 1) * math.fabs(maxy - miny + 1) * math.fabs(maxz - minz + 1) - len(zones)
				print "Excess:", excess

def random_expected_frequency(zone, distribution=None, type1=None, type2=None):
	if distribution is not None and type1 is not None and type2 is not None:
		if zone in distribution.alpha_frequencies[type1][type2]:
			return distribution.alpha_frequencies[type1][type2][zone] * distribution.total_count[type1][type2] / 100.0
	return random.randint(0, 100)

def evaluate_random_frequency_function(input, distribution):
	print os.path.basename(input)
	tag1, tag2 = os.path.basename(input).replace(".txt", "").split("-")
	file = open(input, 'r')
	pzs = read_alpha_zones(file)
	frequencies = {}
	for pz in pzs:
		if pz in frequencies:
			frequencies[pz] += 1
		else:
			frequencies[pz] = 1
	file.close()

	total_difference = 0.0
	total_percent_difference = 0.0
	for pz, actualfreq in frequencies.iteritems():
		predfreq = random_expected_frequency(pz, distribution, int(tag1), int(tag2))
		total_difference += math.fabs(actualfreq - predfreq)
		print predfreq, actualfreq
		total_percent_difference += math.fabs(actualfreq - predfreq) / actualfreq * 100.0
	print "Average difference:", total_difference / len(frequencies)
	print "Average % difference:", total_percent_difference / len(frequencies)

def evaluate_frequency_function(input, distribution):
	print os.path.basename(input)
	tag1, tag2 = os.path.basename(input).replace(".txt", "").split("-")
	file = open(input, 'r')
	pzs = read_alpha_zones(file)
	frequencies = {}
	for pz in pzs:
		if pz in frequencies:
			frequencies[pz] += 1
		else:
			frequencies[pz] = 1
	file.close()

	distribution.load_frequencies(AminoAcid(aatype(int(tag2))), [AminoAcid(aatype(int(tag1)))])

	total_difference = 0.0
	total_percent_difference = 0.0
	for pz, actualfreq in frequencies.iteritems():
		predfreq = distribution.zone_frequency(int(tag1), int(tag2), pz)
		total_difference += math.fabs(actualfreq - predfreq)
		total_percent_difference += math.fabs(actualfreq - predfreq) / actualfreq * 100.0
	print "Average difference:", total_difference / len(frequencies)
	print "Average % difference:", total_percent_difference / len(frequencies)

def partition_zones_stdev(input, size=5, distribution=None):
	"""This function determines the viability of larger partitions of the position zone space to determine if the random frequency function can be used. Pass in distribution to exclude zones that are listed as ranks."""
	print os.path.basename(input)
	file = open(input, 'r')
	tag1, tag2 = os.path.basename(input).replace(".txt", "").split("-")
	pzs = read_alpha_zones(file)
	frequencies = {}
	for pz in pzs:
		if pz in frequencies:
			frequencies[pz] += 1
		else:
			frequencies[pz] = 1
	file.close()

	partition_freqs = {}
	for pz, freq in frequencies.iteritems():
		if distribution and pz in distribution.alpha_frequencies[int(tag1)][int(tag2)]:
			continue
		part = Point3D(math.floor(pz.x / float(size)) * size, math.floor(pz.y / float(size)) * size, math.floor(pz.z / float(size)) * size)
		if part in partition_freqs:
			partition_freqs[part].append(freq)
		else:
			partition_freqs[part] = [freq]

	for part, list in partition_freqs.iteritems():
		mean = numpy.mean(list)
		stdev = numpy.std(list)
		print str(part), mean, stdev

def similar_alpha_frequencies(input, output):
	"""Pass in the entire directory of alpha frequencies. This function will list alpha zones whose relative frequencies are fairly close to one another. Pass in the path for an alpha zone frequencies directory to write the non-similar zones to the appropriate files."""
	files = os.listdir(input)
	data = {}
	for path in files:
		if ".txt" not in path: continue
		print path
		with open(join(input, path), 'r') as file:
			total_freq = 0
			for line in file:
				freq = int(line)
				total_freq += freq
			file.seek(0)
			k = 0
			if total_freq == 0: continue
			for line in file:
				freq = float(line) / float(total_freq)
				pt = Point3D(math.floor(k / 400) - 10.0, math.floor((k % 400) / 20) - 10.0, (k % 20) - 10.0)
				if pt in data:	data[pt].append(freq)
				else:			data[pt] = [freq]
				k += 1
	for pt, freqs in data.iteritems():
		stdev = max(freqs) - min(freqs)
		mean = np.mean(freqs)
		if stdev < 0.0001:
			print str(pt), mean, stdev
			data[pt] = None
	if output is None:
		return
	#Writing the whole file
	for path in files:
		if ".txt" not in path: continue
		print path
		inddata = {}
		values = []
		with open(join(input, path), 'r') as file:
			total_freq = 0
			for line in file:
				freq = int(line)
				total_freq += freq
			file.seek(0)
			k = 0
			if total_freq == 0: continue
			for line in file:
				freq = float(line) / float(total_freq)
				pt = Point3D(math.floor(k / 400) - 10.0, math.floor((k % 400) / 20) - 10.0, (k % 20) - 10.0)
				if data[pt] is not None:
					inddata[pt] = freq
				if freq != 0.0:
					values.append(freq)
				k += 1
		with open(join(output, path), 'w') as file:
			for pt, freq in inddata.iteritems():
				file.write("{0:.0f}, {1:.0f}, {2:.0f}; {3}\n".format(pt.x, pt.y, pt.z, freq))
			median = sorted(values)[int(math.floor(len(values) / 2.0))]
			file.write(str(median) + "\n")

def consecutive_allowed_alpha_zones(input, output, cutoff=0.005):
	"""This method determines which alpha zones for every data file are "allowed" (above a certain frequency) and outputs them to the output files.
		If cutoff >= 1, an absolute frequency is assumed. If cutoff < 1, it is assumed that this means a fraction of the total number of occurrences of this interaction."""
	print os.path.basename(input)
	dist = {}
	with open(input, 'r') as file:
		for line in file:
			line = line.strip()
			if len(line) == 0 or ";" not in line: continue
			comps = line.split(";")
			(x, y, z) = comps[0].split(",")
			dist[Point3D(x, y, z)] = int(comps[1])

	s = sum(dist.values())
	min_zones = cutoff
	if min_zones < 1.0:
		min_zones *= s
	with open(output, 'w') as file:
		for pz, freq in dist.iteritems():
			if freq > min_zones:
				file.write("{}, {}, {}\n".format(pz.x, pz.y, pz.z))


#MARK: - Orientation Analysis
# Now for the orientation analysis (as opposed to alpha zone location)

def angle_frequency_dependence(input):
	'''This function logs the frequencies of various orientations with respect to the angle formed between the x-axis of the second amino acid and the line connecting the two alpha carbons.'''
	print os.path.basename(input)
	file = open(input, 'r')
	tag1, tag2 = os.path.basename(input).replace(".txt", "").split("-")
	pzs = read_positional_zones(file)
	frequencies = {}
	for pz in pzs:
		if pz in frequencies:
			frequencies[pz] += 1
		else:
			frequencies[pz] = 1
	file.close()

	for pz, freq in frequencies.iteritems():
		vector = pz.alpha_zone.multiply(-1.0)
		mags = (pz.x_axis.magnitude() * vector.magnitude())
		if mags == 0.0 or (dotproduct(pz.x_axis, vector) / mags) > 1.0 or (dotproduct(pz.x_axis, vector) / mags) < -1.0: continue
		angle = math.acos(dotproduct(pz.x_axis, vector) / mags)
		print angle, freq

def angle_frequency_vs_zonepairs(input):
	"""This function prints the *number of angle-containing position zones* for each given frequency."""
	files = os.listdir(input)
	for path in files:
		if path.find(".txt") == -1: continue
		print os.path.basename(path)
		file = open(join(input, path), 'r')
		pzs = read_positional_zones(file)
		frequencies = {}
		for pz in pzs:
			if pz in frequencies:
				frequencies[pz] += 1
			else:
				frequencies[pz] = 1
	
		if len(frequencies) == 0: continue
		maximum = max(frequencies.values())
		print maximum
		data = [0 for i in xrange(maximum + 1)]
		for pz, freq in frequencies.iteritems():
			data[freq] += 1
		total_count = 0
		total_frequency = 0
		cutoff = 20
		for i, zpcount in enumerate(data):
			#if zpcount != 0: print i, zpcount
			if zpcount > cutoff:
				total_count += zpcount
				total_frequency += i
		print "Average baseline frequency:", float(total_frequency) / float(total_count)
		del data[:]

def avg_num_angles_used(input):
	'''This function logs the frequencies of various orientations with respect to the angle formed between the x-axis of the second amino acid and the line connecting the two alpha carbons.'''
	file = open(input, 'r')
	tag1, tag2 = os.path.basename(input).replace(".txt", "").split("-")
	pzs = read_positional_zones(file)
	frequencies = {}
	for pz in pzs:
		if pz in frequencies:
			frequencies[pz] += 1
		else:
			frequencies[pz] = 1
	file.close()

	afrequencies = {}

	print os.path.basename(input)
	for pz, freq in frequencies.iteritems():
		if pz.alpha_zone in afrequencies:
			afrequencies[pz.alpha_zone][0] += 1
			afrequencies[pz.alpha_zone][1] += freq
		else:
			afrequencies[pz.alpha_zone] = [1, freq]

	for apz in afrequencies:
		values = afrequencies[apz]
		afrequencies[apz] = float(values[0]) / float(values[1])

	return float(sum(afrequencies.values())) / float(len(afrequencies))

def axes_for_mode(input):
	"""This function finds the most common alpha carbon location in the file and generates a list of all the x-axis locations and their frequencies."""
	file = open(input, 'r')
	tag1, tag2 = os.path.basename(input).replace(".txt", "").split("-")
	pzs = read_positional_zones(file)
	frequencies = {}
	for pz in pzs:
		if pz.alpha_zone in frequencies:
			frequencies[pz.alpha_zone] += 1
		else:
			frequencies[pz.alpha_zone] = 1
	file.close()

	if len(frequencies) == 0: return -1
	mode = (Point3D(-6.0, 1.0, 0.0), frequencies[Point3D(-6.0, 1.0, 0.0)]) #sorted(frequencies.items(), key=lambda x: x[1], reverse=True)[0]
	#print str(mode[0]), mode[1]
	afrequencies = {}
	
	print os.path.basename(input)
	for pz in pzs:
		if pz.alpha_zone == mode[0]:
			if pz.x_axis in afrequencies:
				afrequencies[pz.x_axis] += 1
			else:
				afrequencies[pz.x_axis] = 1

	#for x, freq in afrequencies.iteritems():
	#	print x.x, x.y, x.z, float(freq) / float(mode[1]) * 100.0

	'''y = []
	x = [ [], [] ]
	for pt in afrequencies.iterkeys():
		print pt.x, pt.y, pt.z
		for i in xrange(afrequencies[pt]):
			y.append(pt.z)
			x[0].append(pt.x)
			x[1].append(pt.y)'''
	'''s = sorted(afrequencies.items(), key=lambda x: x[1], reverse=True)[:40]
	y = map(lambda x: x[0].z, s)
	x = [map(lambda x: x[0].x, s),
		 map(lambda x: x[0].y, s)]
	plane = plreg.reg_m(y, x)
	print plane.summary()'''
	#point1, point2, point3 = (x[0] for x in sorted(afrequencies.items(), key=lambda x: x[1], reverse=True)[:3])
	'''cp = crossproduct(point2.subtract(point1), point3.subtract(point1)).normalize()
	print str(cp)
	#ax + by + cz + d = 0
	d = -cp.x * point1.x - cp.y * point1.y - cp.z * point1.z
	#z = (-ax - by - d) / c
	a = -cp.x / cp.z
	b = -cp.y / cp.z
	c = -d / cp.z'''
	'''print str(point1), str(point2), str(point3)
	vector1 = [point2.x - point1.x, point2.y - point1.y, point2.z - point1.z]
	vector2 = [point3.x - point1.x, point3.y - point1.y, point3.z - point1.z]
	cross_product = [vector1[1] * vector2[2] - vector1[2] * vector2[1], -1 * vector1[0] * vector2[2] - vector1[2] * vector2[0], vector1[0] * vector2[1] - vector1[1] * vector2[0]]
	d = cross_product[0] * point1.x - cross_product[1] * point1.y + cross_product[2] * point1.z
	a = cross_product[0]
	b = cross_product[1]
	c = cross_product[2]
	print "{}x + {}y + {}z = {}".format(a, b, c, d)'''

	'''mode = sorted(afrequencies.items(), key=lambda x: x[1], reverse=True)[0]
	print "0", mode[1]
	for xaxis, freq in afrequencies.iteritems():
		print xaxis.distanceto(mode[0]), freq'''

	return afrequencies.keys()

	'''points = afrequencies.keys()#[(p.x, p.y, p.z) for p in afrequencies.keys()]
	#a, b, c = plane_regression(points)
	point1 = point2 = point3 = Point3D.zero()
	while point1 == point2 or point2 == point3 or point1 == point3:
		point1 = random.choice(points)
		point2 = random.choice(points)
		point3 = random.choice(points)
	cp = crossproduct(point1.subtract(point2), point2.subtract(point3))
	print str(cp)
	#ax + by + cz + d = 0
	d = -cp.x * point1.x - cp.y * point1.y - cp.z * point1.z
	#z = (-ax - by - d) / c
	a = -cp.x / cp.z
	b = -cp.y / cp.z
	c = -d / cp.z
	print "z = {}x + {}y + {}".format(a, b, c)'''

def acarbon_axis_relation(input, axisattr="x_axis"):
	"""This function prints the most common axis vector for each alpha carbon location. Pass in the entire directory of data for input. For axisattr, provide a string for the attribute x_axis, y_axis, or z_axis."""
	files = os.listdir(input)
	#data is a dictionary where the keys are acarbon zones and the values are (dictionaries where the keys are x-axes and the values are frequencies).
	data = {}
	for path in files:
		if path.find(".txt") == -1: continue
		print path
		with open(join(input, path), 'r') as file:
			pzs = read_positional_zones(file)
		frequencies = {}
		temp = {}
		for pz in pzs:
			if abs(pz.alpha_zone.distanceto(Point3D.zero()) - 6.0) >= 1.0: continue
			if pz.alpha_zone in data:
				d = data[pz.alpha_zone]
				if getattr(pz, axisattr) in d:
					d[getattr(pz, axisattr)] += 1
				elif pz.alpha_zone in temp and getattr(pz, axisattr) in temp[pz.alpha_zone]:
					d[getattr(pz, axisattr)] = 2
				else:
					if pz.alpha_zone in temp:
						d = temp[pz.alpha_zone]
						if getattr(pz, axisattr) in d:
							data[pz.alpha_zone] = { getattr(pz, axisattr) : 2 }
						else:
							d[getattr(pz, axisattr)] = 1
					else:
						temp[pz.alpha_zone] = { getattr(pz, axisattr) : 1 }
			elif pz.alpha_zone in temp:
				d = temp[pz.alpha_zone]
				if getattr(pz, axisattr) in d:
					data[pz.alpha_zone] = { getattr(pz, axisattr) : 2 }
				else:
					d[getattr(pz, axisattr)] = 1
			else:
				temp[pz.alpha_zone] = { getattr(pz, axisattr) : 1 }
		frequencies.clear()
		temp.clear()
		del frequencies
		del pzs[:]
		del temp
		gc.collect()

	for alphazone, axes in data.iteritems():
		maximum = max(axes.iteritems(), key=lambda x: x[1])
		print str(alphazone), str(maximum[0]), "(%d)" % maximum[1]

def optimum_radius(input, axisattr="x_axis"):
	"""This function determines the radius of the sphere centered at the reference amino acid that passes through the most axis position zones (which axis is used depends on the string value of axisattr - e.g., "x_axis" or "y_axis" or "z_axis"). Only alpha zones with more than 25 occurrences will be used. Input is a single pair of amino acid's path. The output is a list of distances from reference acarbon to subject acarbon paired with the optimum radii."""
	frequencies = {}
	with open(input, 'r') as file:
		pzs = read_positional_zones(file)
		for i, pz in enumerate(pzs):
			if pz.alpha_zone in frequencies:
				frequencies[pz.alpha_zone][0] += 1
				frequencies[pz.alpha_zone][1].append(i)
			else:
				frequencies[pz.alpha_zone] = [ 1, [i] ]

	print os.path.basename(input)
	pairs = []
	for alpha_zone, [freq, idxs] in frequencies.iteritems():
		if freq <= 100: continue
		axes = {}
		for i in idxs:
			ax = getattr(pzs[i], axisattr)
			if ax in axes:
				axes[ax] += 1
			else:
				axes[ax] = 1
					
		#Compile a list of minimum and maximum radii for each axis zone
		radii = {}
		for ax in axes:
			dist1 = Point3D(math.floor(alpha_zone.x + ax.x), math.floor(alpha_zone.y + ax.y), math.floor(alpha_zone.z + ax.z)).distanceto(Point3D.zero())
			dist2 = Point3D(math.ceil(alpha_zone.x + ax.x), math.ceil(alpha_zone.y + ax.y), math.ceil(alpha_zone.z + ax.z)).distanceto(Point3D.zero())
			radii[ax] = ( min(dist1, dist2), max(dist1, dist2) )

		#Test out which radius will pass through the most position zones
		alphadistance = alpha_zone.distanceto(Point3D.zero())
		max_satisfied = 0
		max_radius = 0.0
		for radius in np.arange(alphadistance - 1.0, alphadistance + 1.01, 0.1):
			satisfied = 0
			for ax, bounds in radii.iteritems():
				if radius >= bounds[0] and radius <= bounds[1]:
					satisfied += 1
			if satisfied > max_satisfied:
				max_satisfied = satisfied
				max_radius = radius
		pairs.append([alphadistance, max_radius])
		print alphadistance, max_radius
		axes.clear()
		radii.clear()
		del axes
		del radii
		gc.collect()
	return pairs

def processfile(input, zoneslist, axis, path, q):
	with open(join(input, path), 'r') as file:
		dict = {}
		zones = [z for z in zoneslist if z.distanceto(Point3D.zero()) <= 10.5]
		if len(zones) == 0:
			q.put(dict)
			return
		print path
		read_axis_zones(file, dict, zones, axis=axis)
		q.put(dict)

def update_allaxes(allaxes, result):
	#Add the asynchronously-read results into the allaxes dictionary.
	for key, value in result.iteritems():
		if key in allaxes:
			for key1, value1 in value.iteritems():
				if key1 in allaxes[key]:	allaxes[key][key1] += value1
				else:						allaxes[key][key1] = value1
		else: allaxes[key] = value

def _build_processfile(path, axis, q):
	print os.path.basename(path)
	all_data = [{} for i in xrange(8000)]
	with open(path, 'r') as file:
		for line in file:
			try:
				comps = line.split(";")
				if len(comps) < 4: continue
				(x, y, z) = comps[0].split(",")
				azone = Point3D(x, y, z)
				if axis == 0:
					(x1, y1, z1) = comps[1].split(",")
				elif axis == 1:
					(x1, y1, z1) = comps[2].split(",")
				elif axis == 2:
					(x1, y1, z1) = comps[3].split(",")
				else: assert False, "Invalid axis number (must be 0, 1, or 2): %d" % axis
				pz = Point3D(x1, y1, z1)
				idx = int((azone.x + 10) * 20 * 20 + (azone.y + 10) * 20 + azone.z + 10)
				if pz in all_data[idx]: all_data[idx][pz] += 1
				else: all_data[idx][pz] = 1
			except ValueError:
				print "Exception for", pointstr
	q.put(all_data)

def build_axis_frequencies_file(input, output, axisattr="x_axis", start=0):
	axiscode = { "x_axis" : 0, "y_axis" : 1, "z_axis" : 2 }[axisattr]
	files = os.listdir(input)
	
	seg_len = 100
	worker_count = 15
	for i in xrange(start, 8000, seg_len):
		all_data = [{} for x in xrange(seg_len)]
		zoneslist = [Point3D(math.floor(k / 400) - 10.0, math.floor((k % 400) / 20) - 10.0, (k % 20) - 10.0) for k in xrange(i, min(i + seg_len, 8001))]
		print "Processing", str(zoneslist[0]), "through", str(zoneslist[-1])
		for paths in (files[pos:pos + worker_count] for pos in xrange(0, len(files), worker_count)):
			processes = []
			queue = multiprocessing.Queue()
			for path in paths:
				if ".txt" not in path: continue
				p = multiprocessing.Process(target=processfile, args=(input, zoneslist, axiscode, path, queue))
				processes.append(p)
				p.start()
			for p in processes:
				dict = queue.get()
				for azone, value in dict.iteritems():
					idx = zoneslist.index(azone)#int((azone.x + 10) * 20 * 20 + (azone.y + 10) * 20 + azone.z + 10)
					for key1, value1 in value.iteritems():
						if key1 in all_data[idx]:	all_data[idx][key1] += value1
						else:						all_data[idx][key1] = value1
			for p in processes:
				p.join()
			del processes[:]
		print "Writing"
		with open(join(output, str(int(math.floor(i / 800)))) + ".txt", 'a') as file:
			for j, axes in enumerate(all_data):
				k = i + j
				file.write(str(Point3D(math.floor(k / 400) - 10.0, math.floor((k % 400) / 20) - 10.0, (k % 20) - 10.0)) + "\n")
				for axis, freq in axes.iteritems():
					file.write("\t" + str(axis) + " " + str(freq) + "\n")
		del processes[:]
		del all_data[:]
		gc.collect()
	print "Done"

def percent_satisfied_by_optimum_radius(input, axisattr="x_axis", startpt=Point3D.zero()):
	"""This function finds the optimum radius for each alpha zone over ALL pairs of amino acids and determines what percentage of all occurrences of the alpha zone have axes that fall along that optimum radius."""

	axismapping = { "x_axis" : 0, "y_axis" : 1, "z_axis" : 2 }
	files = os.listdir(input)
	#total_sat = 0.0
	#total_satct = 0
	allaxes = {}
	for zoneslist in Point3D.zero().iteroffsetx(10.0, 25, startpt=startpt, test=lambda x: x.distanceto(Point3D.zero()) <= 10.5):
		#if alpha_zone.x < -7.0 or (alpha_zone.x == -7.0 and (alpha_zone.y < 3.0 or (alpha_zone.y == 3.0 and alpha_zone.z <= -6.0))): continue
		#zoneslist = [x for x in tzoneslist if x.distanceto(Point3D.zero()) <= 10.5]
		if len(zoneslist) == 0: continue
		print "Reading"
		allaxes = {}
		for paths in (files[pos:pos + 10] for pos in xrange(0, len(files), 10)):
			processes = []
			queue = multiprocessing.Queue()
			for path in paths:
				if path.find(".txt") == -1: continue
				p = multiprocessing.Process(target=processfile, args=(input, zoneslist, axismapping[axisattr], path, queue))
				processes.append(p)
				p.start()
			for p in processes:
				update_allaxes(allaxes, queue.get())
			for p in processes:
				p.join()
			del processes[:]

		#while num_completed < num: time.sleep(1.0)
	
		radii = {}
		while len(allaxes) > 0:
			alpha_zone, azaxes = allaxes.popitem()
			print str(alpha_zone)
			#Compile a list of minimum and maximum radii for each axis zone (simultaneously find out how many occurrences total there are)
			total_count = 0
			for ax, freq in azaxes.iteritems():
				dist1 = Point3D(math.floor((alpha_zone.x + ax.x) * 10.0) / 10.0, math.floor((alpha_zone.y + ax.y) * 10.0) / 10.0, math.floor((alpha_zone.z + ax.z) * 10.0) / 10.0).distanceto(Point3D.zero())
				dist2 = Point3D(math.ceil((alpha_zone.x + ax.x) * 10.0) / 10.0, math.ceil((alpha_zone.y + ax.y) * 10.0) / 10.0, math.ceil((alpha_zone.z + ax.z) * 10.0) / 10.0).distanceto(Point3D.zero())
				radii[ax] = ( min(dist1, dist2), max(dist1, dist2) )
				total_count += freq
			
			#Test out which radius will pass through the most position zones
			alphadistance = alpha_zone.distanceto(Point3D.zero())
			max_satisfied = 0
			max_radius = 0.0
			for radius in np.arange(alphadistance - 1.0, alphadistance + 1.01, 0.1):
				satisfied = 0
				for ax, bounds in radii.iteritems():
					if radius >= bounds[0] and radius <= bounds[1]:
						satisfied += azaxes[ax]
				if satisfied > max_satisfied:
					max_satisfied = satisfied
					max_radius = radius
			if total_count > 0:
				print "Finished reading: {} axis zones\n{} {} {}".format(len(azaxes), max_radius, float(max_satisfied) / float(total_count) * 100.0, total_count)
			else:
				print "Finished reading: {} axis zones\n{} {} {}".format(len(azaxes), max_radius, 0.0, 0.0)
			#total_sat += float(max_satisfied) / float(total_count)
			#total_satct += 1
			radii.clear()

		#Cleanup
		allaxes.clear()
		gc.collect()
	#print "Average: %.4f" % (float(total_sat) / float(total_satct))

def calc_optimum_radii(input, axisattr="x_axis"):
	pool = multiprocessing.Pool(processes=3, maxtasksperchild=10)
	#print zipped
	partial_calc = partial(percent_satisfied_by_optimum_radius, input, axisattr=axisattr)
	pool.map(partial_calc, Point3D.zero().iteroffsets(10.0, startpt=Point3D(-6.0, -8.0, 0.0)))
	pool.close()
	pool.join()

def read_axis_frequency_file(path):
	entries = {}
	with open(path, 'r') as file:
		current_pt = None
		for line in file:
			if "eading" in line: continue
			elif "(" in line:
				x, y, z = line.strip("()\n ").split(",")
				current_pt = Point3D(float(x), float(y), float(z))
			else:
				rad, percent = line.replace("\n", "").split(" ")
				entries[current_pt] = (float(rad), float(percent))
	return entries

def test_axis_frequency_function(input, axispath, axisattr="x_axis"):
	freqdata = read_axis_frequency_file(axispath)
	def freq_func(alpha_zone, axis):
		distance = abs(axis.add(alpha_zone).distanceto(Point3D.zero()) - freqdata[alpha_zone][0])
		prob = 1.0 / distance
		return prob
	#Choose random position zones to test
	zoneslist = [Point3D(random.randint(0.0, 10.0), random.uniform(0.0, 2 * math.pi), random.uniform(0.0, math.pi)).tocartesian().floor() for x in xrange(10)]
	for x in zoneslist: print str(x)
	print "Reading"
	files = os.listdir(input)
	axismapping = { "x_axis" : 0, "y_axis" : 1, "z_axis" : 2 }
	allaxes = {}
	for paths in (files[pos:pos + 10] for pos in xrange(0, len(files), 10)):
		processes = []
		queue = multiprocessing.Queue()
		for path in paths:
			if path.find(".txt") == -1: continue
			p = multiprocessing.Process(target=processfile, args=(input, zoneslist, axismapping[axisattr], path, queue))
			processes.append(p)
			p.start()
		for p in processes:
			update_allaxes(allaxes, queue.get())
		for p in processes:
			p.join()
		del processes[:]
	for alpha_zone in allaxes:
		print str(alpha_zone)
		for ax, freq in allaxes[alpha_zone].iteritems():
			print "\t" + str(ax) + " " + str(freq)
	while len(allaxes) > 0:
		alpha_zone, azaxes = allaxes.popitem()
		print str(alpha_zone)
		#Compile a list of minimum and maximum radii for each axis zone (simultaneously find out how many occurrences total there are)
		total_count = 0
		for ax, freq in azaxes.iteritems():
			total_count += freq
		for ax, freq in azaxes.iteritems():
			print freq_func(alpha_zone, ax), float(freq) / float(total_count)

def iter_read_axis_frequencies(path):
	with open(path, 'r') as file:
		current_alpha = None
		freqdata = {}
		for line in file:
			if '\t' in line:			#Axis zone
				pt, freq = line.strip().split(")")
				(x1, y1, z1) = pt.strip("( )").split(",")
				pz = Point3D(x1, y1, z1)
				freqdata[pz] = int(freq)
			else:						#Alpha zone
				if current_alpha is not None:
					yield (current_alpha, freqdata)
				(x1, y1, z1) = line.strip("( )\n").split(",")
				current_alpha = Point3D(x1, y1, z1)
				freqdata.clear()
		if current_alpha is not None:
			yield (current_alpha, freqdata)

def axis_radius_distributions(input, output):
	"""This function prints the distributions of various radii intervals. Pass in the path to the folder of axis data."""
	#files = os.listdir(input)
	#path = files[4]
	total_coeff = 0.0
	total_count = 0
	#file = open(output, 'w')
	for alpha_zone, freqdata in iter_read_axis_frequencies(input):
		#Compile a list of minimum and maximum radii for each axis zone (simultaneously find out how many occurrences total there are)
		print str(alpha_zone)
		if len(freqdata) == 0: continue
		total_count = 0
		radii = {}
		for ax, freq in freqdata.iteritems():
			dist1 = Point3D(alpha_zone.x + ax.x, alpha_zone.y + ax.y, alpha_zone.z + ax.z).distanceto(Point3D.zero())
			dist2 = Point3D(alpha_zone.x + ax.x + 0.1, alpha_zone.y + ax.y + 0.1, alpha_zone.z + ax.z + 0.1).distanceto(Point3D.zero())
			radii[ax] = ( min(dist1, dist2), max(dist1, dist2) )
			total_count += freq
	
		if total_count < 1000: continue
		#Test out which radius will pass through the most position zones
		alphadistance = alpha_zone.distanceto(Point3D.zero())
		distribution = {}
		for ax, bounds in radii.iteritems():
			for r in np.arange(math.ceil(bounds[0] * 10.0) / 10.0, math.floor(bounds[1] * 10.0) / 10.0 + 0.05, 0.1):
				key = "%.1f" % r
				if key in distribution: distribution[key] += freqdata[ax]
				else: distribution[key] = freqdata[ax]

		distarray = sorted(distribution.items(), key=lambda x: float(x[0]))
		#Calculate a quadratic regression
		xs = map(lambda x: float(x[0]), distarray)
		ys = map(lambda x: x[1], distarray)
		maximum = max(ys)
		ys = [float(y) / float(maximum) for y in ys]
		for radius, freq in distarray:
			print radius, float(freq) / float(maximum)
		coeffs = np.polyfit(xs, ys, 2)
		vertex = -coeffs[1] / (2 * coeffs[0])
		print "y={0:.4f}(x - {1:.4f})^2+({2:.4f})".format(coeffs[0], vertex, coeffs[0] * vertex ** 2 + coeffs[1] * vertex + coeffs[2])
		#file.write("{0:.0f}, {1:.0f}, {2:.0f}; {3:.6f}; {4:.6f}; {5:.6f}\n".format(alpha_zone.x, alpha_zone.y, alpha_zone.z, vertex, coeffs[0] * vertex ** 2 + coeffs[1] * vertex + coeffs[2], coeffs[0]))
		total_coeff += coeffs[0]
		total_count += 1
	#file.close()
	print "Average: %.4f" % (total_coeff / float(total_count))

def permissible_triples(permissibility, output):
	allowed_prezones = set()
	allowed_postzones = set()
	for aatype1 in xrange(AMINO_ACID_COUNT):
		aa1 = AminoAcid(aatype(aatype1), 1, acarbon=Point3D.zero())
		aa1.set_axes(Point3D(1.0, 0.0, 0.0), Point3D(0.0, 1.0, 0.0), Point3D(0.0, 0.0, 1.0))
		for aatype2 in xrange(AMINO_ACID_COUNT):
			aa2 = AminoAcid(aatype(aatype2), 0)
			for aatype3 in xrange(AMINO_ACID_COUNT):
				aa3 = AminoAcid(aatype(aatype3), 2)
				#1 is the central amino acid.
				prezones = permissibility.allowed_conformations(aa2, aa1, prior=False)
				for pz in prezones:
					pz.x_axis = pz.x_axis.floor(0.1)
					pz.y_axis = pz.y_axis.floor(0.1)
					pz.z_axis = pz.z_axis.floor(0.1)
					allowed_prezones.add(pz)
				postzones = permissibility.allowed_conformations(aa3, aa1, prior=True)
				for pz in postzones:
					pz.x_axis = pz.x_axis.floor(0.1)
					pz.y_axis = pz.y_axis.floor(0.1)
					pz.z_axis = pz.z_axis.floor(0.1)
					allowed_postzones.add(pz)
				print "zone counts ({}, {}, {}): {} and {}".format(aatype1, aatype2, aatype3, len(allowed_prezones), len(allowed_postzones))
			
	orientations = {}
	aa1 = AminoAcid(0, 1, acarbon=Point3D(0, 0, 0))
	aa1.set_axes(Point3D(1.0, 0.0, 0.0), Point3D(0.0, 1.0, 0.0), Point3D(0.0, 0.0, 1.0))
	for pz1 in prezones:
		aa2 = AminoAcid(0, 0, acarbon=pz1.alpha_zone)
		aa2.set_axes(pz1.x_axis, pz1.y_axis, pz1.z_axis)
		for pz2 in postzones:
			aa3 = AminoAcid(0, 2, acarbon=pz2.alpha_zone)
			aa3.set_axes(pz2.x_axis, pz2.y_axis, pz2.z_axis)
			orientation1 = aa2.tolocal(aa3.acarbon).floor()
			orientation2 = aa3.tolocal(aa2.acarbon).floor()
			retorientation = aa2.aa_position_zone(aa1)
			if orientation1 in orientations:
				if orientation2 in orientations[orientation1]:
					if retorientation not in orientations[orientation1][orientation2]:
						orientations[orientation1][orientation2].append(retorientation)
				else:
					orientations[orientation1][orientation2] = [retorientation]
			else:
				orientations[orientation1] = { orientation2 : [retorientation] }

	print len(orientations), "orientation pairs"
	with open(join(output, "all.txt"), "w") as file:
		for or1, dict in orientations.iteritems():
			for or2, oritems in dict.iteritems():
				file.write("{}, {}, {}, {}, {}, {}; ".format(or1.x, or1.y, or1.z, or2.x, or2.y, or2.z))
				for orient in oritems:
					file.write("{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}; ".format(orient.alpha_zone.x, orient.alpha_zone.y, orient.alpha_zone.z, orient.x_axis.x, orient.x_axis.y, orient.x_axis.z, orient.y_axis.x, orient.y_axis.y, orient.y_axis.z, orient.z_axis.x, orient.z_axis.y, orient.z_axis.z))
				file.write("\n")
	orientations.clear()
	gc.collect()


#MARK: - Hydrophobicity

def format_hydrophobicity_dist(input, output):
	"""This function generates a file for each amino acid type giving the scaled frequency for each "coordination number" as well as the median frequency."""
	dists = [{} for i in xrange(AMINO_ACID_COUNT)]
	if os.path.isdir(input):
		for path in os.listdir(input):
			if ".txt" not in path: continue
			idx = int(path[:path.find(".txt")].split("-")[-1])
			with open(os.path.join(input, path), 'r') as file:
				for line in file:
					line = line.strip()
					if len(line) == 0: continue
					if "," not in line:
						comps = line.split()
					else:
						comps = line.split(",")
					if len(comps) < 2: continue
					coord = int(comps[0])
					dists[idx][coord] = float(comps[1])
	else:
		with open(input, 'r') as file:
			idx = 0
			for line in file:
				line = line.strip()
				if len(line) == 0: continue
				if "," not in line and ":" not in line:
					idx = int(line)
					print idx
				elif ":" in line: continue
				else:
					comps = line.split(",")
					coord = int(comps[0])
					dists[idx][coord] = float(comps[1])

	if not os.path.exists(output): os.mkdir(output)
	#Write them to file and find the median of each one.
	for k, dist in enumerate(dists):
		s = sorted(dist.values())
		if len(s) > 0:
			#median = s[int(len(s) / 2.0)][1]
			median = sum(s) / float(len(s))
		else:
			median = 0.0
		s = sorted(dist.items(), key=lambda x: x[0])
		with open(join(output, "{}.txt".format(k)), 'w') as file:
			for coord, freq in s:
				file.write("{},{}\n".format(coord, freq))
			if median != 0.0:
				file.write(str(median))

def average_coordination_number(input):
	"""Prints the Kyte-Doolittle hydrophobicity scale along with the SPARC-derived average coordination numbers for each amino acid type."""
	kd_hydro = {"ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5, "MET": 1.9, "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8, "TRP": -0.9, "TYR": -1.3, "PRO": -1.6, "HIS": -3.2, "GLU": -3.5, "GLN": -3.5, "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5}
	results = []
	for path in os.listdir(input):
		if ".txt" not in path: continue
		total_coord = 0
		total_occurrences = 0
		with open(os.path.join(input, path), "r") as file:
			for line in file:
				if "," not in line: continue
				comps = line.strip().split(",")
				total_coord += float(comps[0]) * float(comps[1])
				total_occurrences += float(comps[1])
		if aatype(int(path[:path.find('.txt')])) in kd_hydro:
			#print aatype(int(path[:path.find('.txt')]))
			results.append((aatype(int(path[:path.find('.txt')])), kd_hydro[aatype(int(path[:path.find('.txt')]))],(total_coord / total_occurrences)))
	results = sorted(results, key=lambda x: x[1])
	for r in results:
		print r[0] + "\t" + str(r[1]) + "\t" + str(r[2])

#MARK: - Testing Statistical Potential
def sparc_scores(input, dists):
	"""This function iterates through all the files listed in the directory of 'input' and, after reading their structure into memory, evaluates their score using the provided distribution objects. The result is a dictionary containing filepath: score."""
	files = os.listdir(input)
	ret = {}
	for path in files:
		if "pdb" not in path: continue
		peptide = Polypeptide()
		peptide.read(join(input, path))
		s = sum(dist.score(peptide, peptide.aminoacids) for dist in dists)
		print path, s
		ret[path] = s
		del peptide.aminoacids[:]
		peptide = None
		gc.collect()
	return ret

def sparc_scores_file(input, dists, bounds=None, retbounds=False, ignored_aas=None, noeval=False, peptide=None):
	"""This function iterates through the one file at 'input' and, after reading its structure into memory, evaluates their score using the provided distribution objects with all the provided lists of weights. The result is a dictionary containing filepath: score. If you pass a tuple (min, max) for bounds, only those tags (inclusive) will be considered. If you pass True for retbounds, this function will return the min and max of the polypeptide read function in a tuple along with a list of gaps, and the default return value. If you pass True for noeval, the function will only return the bounds and gaps if you so asked, without computing the SPARC scores at all."""
	#files = os.listdir(input)
	ret = None
	if not noeval:
		print os.path.basename(input)
	#if "T" not in input: return ret
	if os.path.basename(input)[0] == ".":
		if retbounds == True:
			return (None, ret)
		else:
			return ret
	if not peptide:
		peptide = Polypeptide()
	gaps, rbounds = peptide.read(input, checkgaps=True, fillgaps=True)
	if bounds is not None:
		bounds = (bounds[0] - rbounds[0], bounds[1] - rbounds[0])
		if bounds[0] < 0 or bounds[-1] >= len(peptide.aminoacids):
			if retbounds == True:
				return (rbounds, ret)
			else:
				return ret
		pre = peptide.aminoacids[:bounds[0]]
		post = peptide.aminoacids[bounds[1] + 1:]
		for aa in pre:
			peptide.hashtable.remove(aa)
		for aa in post:
			peptide.hashtable.remove(aa)
		del peptide.aminoacids[bounds[1] + 1:]
		del peptide.aminoacids[:bounds[0]]
		idx = 0
		for aa in peptide.aminoacids:
			if aa:
				aa.tag = idx
			idx += 1
	if noeval:
		if retbounds == True:
			return (rbounds, gaps)
		else:
			return ret
	for dist in dists:
		dist.weight = 1.0
	considering_aas = peptide.aminoacids
	if ignored_aas:
		considering_aas = [aa for aa in peptide.aminoacids if aa and aa.tag not in ignored_aas]
	ret = [dist.score(peptide, considering_aas) for dist in dists]
	del peptide.aminoacids[:]
	if retbounds == True:
		return (rbounds, gaps, ret)
	else:
		return ret

def _weight_data(input, weights, start=0, basename=None):
	best_combo = 100000000
	best_nm = None
	native = 0.0
	nonnative = []
	scores = []
	with open(input, 'r') as file:
		for line in file:
			if ";" not in line: continue
			line = line.strip()
			nm = line[:line.find(";")]
			if ".pdb" in nm: nm = nm[:nm.find(".pdb")]
			comps = line[line.find(";") + 1:].split(",")
			score = 0.0
			for i, w in enumerate(weights):
				if type(start) is list:
					score += float(comps[start[i]]) * w
				else:
					score += float(comps[i + start]) * w
			if score < best_combo:
				best_combo = score
				best_nm = nm
			if basename and (nm == basename or nm == basename + "_orig"):
				native = score
			else:
				nonnative.append(score)
			scores.append((nm, score))
	scores = sorted(scores, key=lambda x: x[1])
	native_rank = next((i for i in xrange(len(scores)) if basename and (scores[i][0] == basename or scores[i][0] == basename + "_orig")), 0)
	stdev = numpy.std(nonnative)
	mean = sum(nonnative) / float(len(nonnative))
	return (best_nm, float(native - mean) / stdev, native_rank)

def _all_weight_data(input, weights, start=0):
	scores = {}
	with open(input, 'r') as file:
		for line in file:
			if ";" not in line: continue
			line = line.strip()
			nm = line[:line.find(";")]
			if ".pdb" in nm: nm = nm[:nm.find(".pdb")]
			comps = line[line.find(";") + 1:].split(",")
			score = 0.0
			for i, w in enumerate(weights):
				if type(start) is list:
					score += float(comps[start[i]]) * w
				else:
					score += float(comps[i + start]) * w
			scores[nm] = score
	return scores


def best_weights(input, numweights=4, start=0):
	"""This function determines which set of integer weights maximizes the number of correct guesses for the native structure. Correctness is determined by the filename being equivalent to the structure chosen."""
	weights = []
	def genweights(levels, min, max):
		if levels == 1:
			for i in xrange(min, max): yield [i]
		else:
			for w in genweights(levels - 1, min, max):
				for i in xrange(min, max): yield w + [i]
	for weightlist in genweights(numweights, 1, 6):
		weights.append(weightlist)
	#if numweights == 4: weights = [[4, 1, 4, 5]]
	maxcorrect = 0
	maxweights = []
	max_z = -100000
	print os.listdir(input)
	for w in weights:
		correct = 0
		zscores = []
		total_rank = 0
		total_ct = 0
		for path in os.listdir(input):
			if os.path.isdir(os.path.join(input, path)):
				for subpath in os.listdir(os.path.join(input, path)):
					if ".txt" not in subpath: continue
					native_guess, z, rank = _weight_data(os.path.join(input, path, subpath), w, start=start)
					if native_guess == subpath[:subpath.find(".txt")]:
						correct += 1
					zscores.append(z)
					#print subpath[:subpath.find(".txt")], native_guess
				continue
			if ".txt" not in path: continue
			native_guess, z, rank = _weight_data(os.path.join(input, path), w, start=start, basename=path[:path.find(".txt")])
			total_rank += rank
			total_ct += 1
			if native_guess == path[:path.find(".txt")] or native_guess == path[:path.find(".txt")] + "_orig":
				correct += 1
			zscores.append(z)
			#print path[:path.find(".txt")], native_guess
		print str(w) + "\t" + str(correct) + "\t" + str(sum(zscores) / float(len(zscores))) + "\t" + str(total_rank / float(total_ct))
		if correct > maxcorrect:
			maxcorrect = correct
			maxweights = [w]
			for z in zscores:
				print z
			max_z = sum(zscores) / float(len(zscores))
		elif correct == maxcorrect:
			z = sum(zscores) / float(len(zscores))
			if z < max_z:
				for z in zscores:
					print z
				max_z = z
				maxcorrect = correct
				maxweights = [w]
	print "Final: the combos", maxweights, "had a total of", maxcorrect, "correct guesses, with z score", max_z

def best_weights_rmsd(input, rmsd_files, native_paths, numweights=4, start=0, groups=False, tmscore=False):
	"""This method performs exactly the same function as best_weights, but instead of using the Z-score as the rank, it uses the correlation coefficient between the RMSDs of the decoys and the SPARC scores."""
	weights = []
	def genweights(levels, min, max):
		if levels == 1:
			for i in xrange(min, max): yield [i]
		else:
			for w in genweights(levels - 1, min, max):
				for i in xrange(min, max): yield w + [i]
	for weightlist in genweights(numweights, 0, 5):
		weights.append(weightlist)
	maxcorrect = 0
	maxweights = []
	max_rsq = 0.0
	#First compute the RMSD distances for all the decoys
	all_rmsds = {}
	if groups:
		group_rmsds = {}
	if os.path.isdir(rmsd_files):
		rmsd_paths = os.listdir(rmsd_files)
	else:
		rmsd_paths = [rmsd_files]
	for rmsd_file in rmsd_paths:
		if ".txt" not in rmsd_file: continue
		print rmsd_file
		sub_rmsd = {}
		if groups:
			subgroups = {}
		with open(os.path.join(rmsd_files, rmsd_file), "r") as file:
			for line in file:
				if "," not in line or "None" in line: continue
				comps = line.split(",")
				fname = comps[0]
				if ".pdb" in fname: fname = fname[:fname.find(".pdb")]
				sub_rmsd[fname] = float(comps[1])
				if groups:
					nm_comps = fname.split("_")
					if len(nm_comps) == 2:
						if nm_comps[0] in subgroups:
							subgroups[nm_comps[0]][fname] = float(comps[1])
						else:
							subgroups[nm_comps[0]] = { fname: float(comps[1]) }

		all_rmsds[rmsd_file] = sub_rmsd
		if groups:
			group_rmsds[rmsd_file] = subgroups

	print "Stored", len(all_rmsds), "sets of RMSD data"

	for w in weights:
		correct = 0
		correct_model = 0
		total_r2 = 0.0
		num_r2 = 0
		total_best_rmsd = 0.0
		total_score_rank = 0
		total_rmsd_rank = 0
		sparc_scores = []
		rmsd_scores = []
		#if w != [4, 1, 4, 5] and numweights == 4: continue
		if os.path.isdir(input):
			input_paths = os.listdir(input)
		else:
			input_paths = [input]
		for path in input_paths:
			if ".txt" not in path: continue
			if len(input_paths) == 1 and len(all_rmsds) == 1:
				all_rmsds[path] = all_rmsds.values()[0]
			if path not in all_rmsds:
				continue
			scores = _all_weight_data(os.path.join(input, path), w, start=start)
			native_guess = min(scores, key=scores.get)
			if native_guess == path[:path.find(".txt")] or "_orig" in native_guess: #native_guess == path[:path.find(".txt")] + "_orig":
				correct += 1
			#print path[:path.find(".txt")], native_guess
			rmsd_list = all_rmsds[path]
			min_score = 10000000.0
			min_score_rmsd = 0.0
			min_filename = None
			min_rmsd = 0.0
			min_rmsd_score = 0.0
			min_rmsd_filename = None
			for filename, score in scores.iteritems():
				#if filename == path[:path.find(".txt")]: continue
				if filename not in rmsd_list:
					#if filename == path[:path.find(".txt")] or filename == path[:path.find(".txt")] + "_orig":
					#	rmsd_list[filename] = 1.0
					#else:
					continue
				sparc_scores.append(score)
				rmsd_scores.append(rmsd_list[filename])
				if filename != path[:path.find(".txt")] and score < min_score:
					min_filename = filename
					min_score = score
					min_score_rmsd = rmsd_list[filename]
				if filename != path[:path.find(".txt")] and ((tmscore and rmsd_list[filename] > min_rmsd) or (not tmscore and rmsd_list[filename] < min_rmsd)):
					min_rmsd = rmsd_list[filename]
					min_rmsd_score = score
					min_rmsd_filename = filename
			if min_filename == min_rmsd_filename: #min_filename[-1] == "1":
				correct_model += 1
			m, b, r2, p, se = linregress(sparc_scores, rmsd_scores)
			extra = ""
			score_rank = len([x for x in sparc_scores if x < min_rmsd_score])
			if tmscore:
				rmsd_rank = len([x for x in rmsd_scores if x > min_score_rmsd])
			else:
				rmsd_rank = len([x for x in rmsd_scores if x < min_score_rmsd])
			total_score_rank += score_rank
			total_rmsd_rank += rmsd_rank
			if groups:
				group_correct = 0
				subgroups = group_rmsds[path]
				for name, subgroup in subgroups.iteritems():
					min_group = min(subgroup, key=lambda x: (scores[x] if x in scores else 1000000.0))
					if int(min_group.split("_")[1]) < 3: group_correct += 1
				extra = "\t" + str(group_correct) + "\t" + str(len(subgroups))
			if w == [1, 5, 5, 1]:
				best_rmsd = 0
				if tmscore:
					best_rmsd = max(rmsd_scores)
				else:
					best_rmsd = min(rmsd_scores)
				print min_filename + "\t" + str(sum(rmsd_scores) / float(len(rmsd_scores))) + "\t" + str(best_rmsd) + "\t" + str(score_rank) + "\t" + str(min_score_rmsd) + "\t" + str(r2 ** 2) + extra
			total_r2 += r2 ** 2
			total_best_rmsd += min_score_rmsd / min_rmsd
			num_r2 += 1
			del sparc_scores[:], rmsd_scores[:]
		total_r2 /= float(num_r2)
		print "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(w, correct, correct_model, total_r2, total_best_rmsd / float(num_r2), total_score_rank / float(num_r2), total_rmsd_rank / float(num_r2))
		if total_r2 > max_rsq:
			maxweights = w
			maxcorrect = correct
			max_rsq = total_r2

	print "Final: the combos", maxweights, "had a total of", maxcorrect, "correct guesses, with R^2", max_rsq


def z_score(input, weights, structure_files=None):
	"""Prints the z-score of each file in input by computing the score for each line and averaging, subtracting the native score, and dividing by the standard deviation. If structure_files is not None, then this will divide the energies by the number of amino acids in EACH structure before comparing."""
	def file_func(filepath):
		if ".txt" not in filepath: return
		protein_name = os.path.basename(filepath)[:os.path.basename(filepath).find(".txt")]
		nonnative = []
		native = 0
		minscore = 1000000
		print protein_name
		p = Polypeptide()
		with open(filepath, "r") as file:
			for line in file:
				nm = line[:line.find(";")]
				if ".pdb" in nm: nm = nm[:nm.find(".pdb")]
				comps = line[line.find(";") + 1:].split(",")
				score = 0.0
				for i, w in enumerate(weights):
					score += float(comps[i]) * w
				if score != 0.0:
					if structure_files:
						if not os.path.exists(os.path.join(structure_files, protein_name, nm)): break
						p = Polypeptide()
						p.read(os.path.join(structure_files, protein_name, nm))
						score /= len(p.aminoacids)
						del p.aminoacids[:]
						p.aminoacids = []
						p.hashtable = AAHashTable()
					if nm == protein_name or nm == protein_name + "_orig":
						native = score
					else:
						nonnative.append(score)
					if score < minscore:
						minscore = score
					gc.collect()
		stdev = numpy.std(nonnative)
		mean = sum(nonnative) / float(len(nonnative))
		#print protein_name + " & " + str(float(mean - native) / stdev)
		correct = False
		if native == minscore:
			correct = True
		print float(native - mean) / stdev, correct
		return (float(native - mean) / stdev, 1, correct)

	total_zs = 0
	total_zcount = 0
	total_correct = 0
	
	for path in os.listdir(input):
		if os.path.isdir(os.path.join(input, path)):
			for subpath in os.listdir(os.path.join(input, path)):
				ret = file_func(os.path.join(input, path, subpath))
				if ret:
					total_zs += ret[0]
					total_zcount += ret[1]
					total_correct += ret[2]
				gc.collect()
			continue
		ret = file_func(os.path.join(input, path))
		if ret:
			total_zs += ret[0]
			total_zcount += ret[1]
			total_correct += ret[2]
		gc.collect()
	return (total_correct, total_zcount, float(total_zs) / total_zcount)

def tm_sparc_correlation(tm, sparc):
	"""This function reads each file in the tm and sparc directories (should be paired by filename) and finds the correlation coefficient between their values."""
	for path in os.listdir(tm):
		if "txt" not in path: continue
		if not os.path.exists(join(sparc, path)): continue
		values = {}
		with open(join(tm, path), 'r') as file1:
			for line in file1:
				comps = line.strip().split(" ")
				values[comps[0]] = [float(comps[1]), 0.0]
		with open(join(sparc, path), 'r') as file2:
			for line in file2:
				comps = line.strip().split(" ")
				if comps[0] not in values: continue
				values[comps[0]][1] = float(comps[1])

		for struct, scores in values.iteritems():
			print scores[0], scores[1]
		values = values.values()
		slope, intercept, r_value, p_value, std_err = linregress(map(lambda x: x[0], values), map(lambda x: x[1], values))
		print "Correlation:", r_value ** 2

def calc_rmsd(file1, file2, tmpdir):
	origWD = os.getcwd() # remember our original working directory
	os.chdir(tmpdir)
	
	out = None
	try:
		with open(os.devnull, "w") as fnull:
			cmd = "rm \\#*\\#; echo \"3\n3\n\" | gmx confrms -name -f1 \"" + file1 + "\" -f2 \"" + file2 + "\""
			out = subprocess.check_output(cmd, shell=True, stderr=fnull)
	except:
		print "Exception", os.path.basename(file1), os.path.basename(file2)
		os.chdir(origWD) # get back to our original working directory
		return 100000

	os.chdir(origWD) # get back to our original working directory
	if out:
		out = out[out.find("lsq fit = ") + len("lsq fit = "):]
		rms = float(out[:out.find("nm")])
		return rms
	return 100000

def decoys_rmsd(input, native_paths, output, writefile=True):
	"""For every structure in the CASP datasets, calculates the rms distance using GROMACS and puts it in a file in the directory at output (assuming writefile is the default value of False). Also returns a dictionary of filename-RMSD entries."""
	if writefile:
		if not os.path.exists(output): os.mkdir(output)
	
	tmpdir = os.path.join(input, "sp_tmp")
	if os.path.exists(tmpdir):
		shutil.rmtree(tmpdir)
	os.mkdir(tmpdir)
	rmsds = {}
	for path in reversed(os.listdir(input)):
		if not os.path.isdir(os.path.join(input, path)) or path == "sp_tmp": continue
		if output:
			if os.path.exists(os.path.join(output, path)) or os.path.exists(os.path.join(output, path + ".txt")): continue
		if native_paths is not None:
			nativepath = os.path.join(native_paths, path + ".pdb")
		else:
			nativepath = os.path.join(input, path, path + ".pdb")
			assert os.path.exists(nativepath), "No native path for %r!" % path
		print path
		sub_rmsds = {}
		for subpath in os.listdir(os.path.join(input, path)):
			if "T" not in subpath and ".pdb" not in subpath: continue
			#print subpath
			if ".pdb" not in subpath:
				file2 = os.path.join(tmpdir, subpath + ".pdb")
				shutil.copyfile(os.path.join(input, path, subpath), file2)
			else:
				file2 = os.path.join(input, path, subpath)
			if file2 == nativepath: continue
			rmsd = calc_rmsd(nativepath, file2, tmpdir)
			if rmsd > 10000: continue
			sub_rmsds[subpath] = rmsd
			if writefile:
				with open(os.path.join(output, path + ".txt"), "a") as file:
					file.write(subpath + "," + str(rmsd) + "\n")
		rmsds[path] = sub_rmsds
	#shutil.rmtree(tmpdir)
	return rmsds

def link_decoy_rmsd(energy_files, rmsd_files, weights):
	"""energy_files should be a directory of energy files (not nested). rmsd_files should be a directory of rmsd files (not nested)."""
	for path in os.listdir(energy_files):
		if ".txt" not in path: continue
		if not os.path.exists(os.path.join(rmsd_files, path[:path.find(".txt")])):
			print path, "has no match in rmsd."
			continue
		energies = {}
		print "----", path
		with open(os.path.join(energy_files, path), "r") as file:
			for line in file:
				nm = line[:line.find(";")]
				if ".pdb" in nm: nm = nm[:nm.find(".pdb")]
				comps = line[line.find(";") + 1:].split(",")
				score = float(comps[0]) * weights[0] + float(comps[1]) * weights[1] + float(comps[2]) * weights[2]
				energies[nm] = score
		with open(os.path.join(rmsd_files, path[:path.find(".txt")]), "r") as file:
			for line in file:
				nm = line[:line.find(",")]
				if ".pdb" in nm: nm = nm[:nm.find(".pdb")]
				rmsd = float(line[line.find(",") + 1:])
				if nm in energies:
					print str(energies[nm]) + "\t" + str(rmsd)
		gc.collect()

def min_rmsd(input, original_structure, range=None, dists=None, sec_structs=None, separate_scores=False, writeout=False, output=None):
	'''If you pass dists, the SPARC score will be logged next to the rmsd.
		If you pass range, the section of the original structure marked by range will be compared.'''
	tmpdir = os.path.join(os.path.dirname(input), "sp_tmp")
	if os.path.exists(tmpdir): shutil.rmtree(tmpdir)
	os.mkdir(tmpdir)
	min_rmsd = 100000
	min_model = 0
	rmsd_sparc = []
	min_sparc = 0.0
	total_rmsd = 0.0
	num_rmsd = 0
	
	peptide = Polypeptide()
	if range:
		#Write the section of the original structure that corresponds to range
		peptide.read(original_structure)
		original_structure = os.path.join(tmpdir, "orig.pdb")
		with open(original_structure, "w") as file:
			file.write(peptide.pdb(range=range))

	with open(input, "r") as file:
		next_pdb = ""
		modelno = 0
		for line in file:
			next_pdb += line
			if "MODEL" in line:
				modelno = int(line[5:].strip())
			elif "ENDMDL" in line:
				with open(os.path.join(tmpdir, "model.pdb"), "w") as file:
					file.write(next_pdb)
				rmsd = calc_rmsd(original_structure, os.path.join(tmpdir, "model.pdb"), tmpdir)
				if rmsd > 10000: continue
				if dists:
					try:
						peptide.read(os.path.join(tmpdir, "model.pdb"))
					except:
						next_pdb = ""
						del peptide.aminoacids[:]
						peptide.hashtable.clear()
						continue
					if sec_structs:
						if ',' in sec_structs:
							peptide.add_secondary_structures(sec_structs, format='csv')
						else:
							peptide.add_secondary_structures(sec_structs, format='pdb')
					if separate_scores:
						sp_score = ""
						for d in dists:
							sp_score += str(d.score(peptide, peptide.aminoacids)) + "\t"
						sp_score = sp_score[:-1]
						print str(modelno) + "\t" + str(rmsd) + "\t" + sp_score
					else:
						sp_score = sum(d.score(peptide, peptide.aminoacids) for d in dists)
						print str(modelno) + "\t" + str(rmsd) + "\t" + str(sp_score)
					
					rmsd_sparc.append((sp_score, rmsd))
					if rmsd < min_rmsd:
						min_sparc = sp_score
					del peptide.aminoacids[:]
					peptide.hashtable.clear()
				elif writeout:
					print "{}\t{}".format(modelno, rmsd)
				if output:
					with open(output, "a") as outfile:
						outfile.write("{},{}\n".format(modelno, rmsd))
				next_pdb = ""
				if rmsd < min_rmsd:
					min_rmsd = rmsd
					min_model = modelno
				total_rmsd += rmsd
				num_rmsd += 1
	shutil.rmtree(tmpdir)
	if len(rmsd_sparc):
		rank = len([x for x in rmsd_sparc if x[0] <= min_sparc])
		best_items = sorted(rmsd_sparc, key=lambda x: x[0])[:25]
		avg = sum(x[1] for x in best_items) / float(len(best_items))
		maximum = max(x[1] for x in best_items)
		minimum = min(x[1] for x in best_items)
		print "MIN:", min_model, min_rmsd, ". Rank among structures:", rank, minimum, maximum, avg
		return (min_rmsd, min_sparc, sum(x[0] for x in rmsd_sparc) / float(len(rmsd_sparc)), rank, minimum, maximum, avg)
	else:
		if num_rmsd == 0: return 0.0
		print "min:", min_model, min_rmsd, "avg:", total_rmsd / float(num_rmsd)
		return min_rmsd

def sum_squared_difference(input):
	points = []
	initial = 0
	with open(input,"r") as file:
		i = 0
		for line in file:
			if i == 0:
				initial = float(line[:line.find("\t")])
				i = 1
				continue
			line = line.strip()
			points.append(float(line.split()[1]))
	m = (points[-1] - points[0]) / float(len(points))
	b = points[0]
	ssd = 0.0
	for i, p in enumerate(points):
		ssd += (p - m * i - b) ** 2
	return (initial, ssd)

def score_structure_file(input, dists):
	"""Prints the SPARC scores of all structures in the input file."""
	peptide = None
	with open(input, "r") as file:
		next_pdb = []
		modelno = 0
		for line in file:
			next_pdb.append(line)
			if "MODEL" in line:
				modelno = int(line[5:].strip())
			elif "ENDMDL" in line:
				peptide = Polypeptide()
				peptide.read_file(next_pdb)
				newc = peptide.aminoacids[5].toglobal(Point3D(CARBON_BOND_LENGTH, math.pi - math.acos(-1.0 / 2.0) / 2.0, math.acos(-1.0 / 3.0)).tocartesian())
				newn = peptide.aminoacids[5].toglobal(Point3D(NITROGEN_BOND_LENGTH, math.pi + math.acos(-1.0 / 2.0) / 2.0, math.acos(-1.0 / 3.0)).tocartesian())
				print peptide.aminoacids[5].carbon, math.cos(peptide.aminoacids[5].carbon.subtract(peptide.aminoacids[5].acarbon).anglewith(peptide.aminoacids[5].nitrogen.subtract(peptide.aminoacids[5].acarbon))), math.cos(newc.subtract(peptide.aminoacids[5].acarbon).anglewith(newn.subtract(peptide.aminoacids[5].acarbon))), sum(d.score(peptide, peptide.aminoacids) for d in dists)
				del peptide.aminoacids[:]
				peptide.hashtable = None
				peptide = None
				next_pdb = []
	print "\n"


#MARK: Miscellaneous

def read_fasta(input):
	"""Returns a dictionary where each key is the first word after the > in each entry, and the value is the sequence (no newlines)."""
	ret = {}
	with open(input) as file:
		current_sequence = ""
		current_tag = ""
		for line in file:
			if ">" in line:
				ret[current_tag] = current_sequence
				current_tag = line[1: line.find(" ")]
				current_sequence = ""
			else:
				current_sequence += line.strip()
	return ret

#MARK: 2016 new methods for bigger data

def _aggregate_aa(input, output, both, req_aa):
	if os.path.exists(os.path.join(output, "medium", str(req_aa) + ".txt")):
		return
	done_secondary = (req_aa != 0)
	if not done_secondary:
		if not os.path.exists(os.path.join(output, "secondary")): os.mkdir(os.path.join(output, "secondary"))
		secondary_data = {
			secondary_struct_helix + "1" : {}, secondary_struct_helix + "2" : {},
			secondary_struct_helix + "3" : {}, secondary_struct_helix + "4" : {},
			secondary_struct_helix + "5" : {}, secondary_struct_helix + "6" : {},
			secondary_struct_helix + "7" : {}, secondary_struct_helix + "8" : {},
			secondary_struct_helix + "9" : {}, secondary_struct_helix + "10" : {},
			secondary_struct_sheet + "0" : {}, secondary_struct_sheet + "1" : {},
			secondary_struct_sheet + "-1" : {}
		}
	data = {
		"long_range" : [{} for i in range(AMINO_ACID_COUNT)],
		"consec" : [{} for i in range(AMINO_ACID_COUNT)],
		"consec+secondary" : [{} for i in range(AMINO_ACID_COUNT)],
		"short_range" : [{} for i in range(AMINO_ACID_COUNT)],
		"medium" : [0 for i in xrange(100)]
	}
	print "========", req_aa
	def process_line(line, subdata, twopzs=False):
		comps = line.split(";")
		if both:
			coords = comps[0].split(",")
			pz = PositionZone(Point3D(*coords[:3]))
			pz2 = PositionZone(Point3D(*coords[3:]))
			key = (pz, pz2)
		else:
			key = PositionZone(Point3D(*comps[0].split(",")))
		if key in subdata:
			subdata[key] += int(comps[1])
		else:
			subdata[key] = int(comps[1])

	for subdir in os.listdir(input):
		if not os.path.isdir(os.path.join(input, subdir)): continue
		print subdir
		for subsubdir in os.listdir(os.path.join(input, subdir)):
			if "medium" in subsubdir and not os.path.isdir(os.path.join(input, subdir, subsubdir)) and int(subsubdir[7:-4]) == req_aa:
				print "Analyzing medium file", subsubdir, "in", subdir
				with open(os.path.join(input, subdir, subsubdir), "r") as file:
					for line in file:
						comps = line.strip().split()
						if int(comps[0]) < len(data["medium"]):
							data["medium"][int(comps[0])] += int(comps[1])

			if not os.path.isdir(os.path.join(input, subdir, subsubdir)): continue
			if subsubdir == "consec" or subsubdir == "long_range" or subsubdir == "short_range" or subsubdir == "consec+secondary":
				for aacombo in os.listdir(os.path.join(input, subdir, subsubdir)):
					if ".txt" not in aacombo: continue
					aas = [int(x) for x in aacombo.replace(".txt", "").split("-")]
					if aas[0] != req_aa: continue
					with open(os.path.join(input, subdir, subsubdir, aacombo), "r") as file:
						for line in file: process_line(line, data[subsubdir][aas[1]], twopzs=both)
					del file
			elif subsubdir == "secondary" and not done_secondary:
				for sec_struct in os.listdir(os.path.join(input, subdir, subsubdir)):
					if ".txt" not in sec_struct: continue
					struct_type = sec_struct[:-4]
					with open(os.path.join(input, subdir, subsubdir, sec_struct), "r") as file:
						for line in file: process_line(line, secondary_data[struct_type], twopzs=both)
					del file
			if subsubdir == "medium":
				for aafile in os.listdir(os.path.join(input, subdir, subsubdir)):
					if ".txt" not in aafile: continue
					aaid = int(aafile.replace(".txt", ""))
					if aaid != req_aa: continue
					with open(os.path.join(input, subdir, subsubdir, aafile), "r") as file:
						for line in file:
							comps = line.strip().split()
							if len(comps) > 1:
								if int(comps[0]) < len(data["medium"]):
									data["medium"][int(comps[0])] += int(comps[1])
							else:
								if int(line.strip()) < len(data[subsubdir]):
									data[subsubdir][int(line.strip())] += 1
					del file
	gc.collect()

	def write_file(key, twopzs=False):
		has_contents = False
		for j in range(AMINO_ACID_COUNT):
			if len(data[key][j]):
				has_contents = True
				break
		if not has_contents: return
		
		file_path = join(output, key)
		if not os.path.exists(file_path): os.mkdir(file_path)
		i = req_aa
		for j in range(AMINO_ACID_COUNT):
			if not len(data[key][j]) or (i > j and twopzs): continue
			f = open(join(file_path, "%d-%d.txt" % (i, j)), 'w')
			if twopzs:
				for pzs, freq in data[key][j].iteritems():
					f.write(str(int(pzs[0].alpha_zone.x)) + "," + str(int(pzs[0].alpha_zone.y)) + "," + str(int(pzs[0].alpha_zone.z)) + "," + str(int(pzs[1].alpha_zone.x)) + "," + str(int(pzs[1].alpha_zone.y)) + "," + str(int(pzs[1].alpha_zone.z)) + "; " + str(freq) + "\n")
			else:
				for pz, freq in data[key][j].iteritems():
					f.write(str(pz.alpha_zone.x) + ", " + str(pz.alpha_zone.y) + ", " + str(pz.alpha_zone.z) + "; " + str(freq) + "\n")
			f.close()

	write_file("long_range", twopzs=both)
	write_file("consec", twopzs=both)
	write_file("consec+secondary", twopzs=both)
	write_file("short_range", twopzs=both)

	#Medium
	if not both:
		if not os.path.exists(os.path.join(output, "medium")): os.mkdir(os.path.join(output, "medium"))
		medium_path = join(output, "medium", str(req_aa) + ".txt")
		with open(medium_path, "w") as file:
			for i, x in enumerate(data["medium"]):
				file.write(str(i) + ' ' + str(x) + '\n')
		data.clear()
		gc.collect()

	if not done_secondary and len(secondary_data["helix1"]):
		#Secondary
		secondary_path = join(output, "secondary")
		for struct_type in secondary_data:
			if not len(secondary_data[struct_type]): continue
			f = open(join(secondary_path, struct_type + ".txt"), 'w')
			if both:
				for pzs, freq in secondary_data[struct_type].iteritems():
					f.write(str(int(pzs[0].alpha_zone.x)) + "," + str(int(pzs[0].alpha_zone.y)) + "," + str(int(pzs[0].alpha_zone.z)) + "," + str(int(pzs[1].alpha_zone.x)) + "," + str(int(pzs[1].alpha_zone.y)) + "," + str(int(pzs[1].alpha_zone.z)) + "; " + str(freq) + "\n")
			else:
				for pz, freq in secondary_data[struct_type].iteritems():
					f.write(str(pz.alpha_zone.x) + ", " + str(pz.alpha_zone.y) + ", " + str(pz.alpha_zone.z) + "; " + str(freq) + "\n")
			f.close()
		done_secondary = True


def aggregate_networkdata(input, output, both=False):
	"""Aggregates all the data found in the directory specified by input in a folder specified by output. Assumes that the data is structured as follows:
		1
			consec
				0-0.txt
				0-1.txt
				.
				.
				.
			medium
				0.txt
				.
				.
				.
			nonconsec
				0-0.txt
				0-1.txt
				.
				.
				.
		2
		...
		The output is formatted the same way.
		"""
	if not os.path.exists(output): os.mkdir(output)
	pool = multiprocessing.Pool(processes=3, maxtasksperchild=1)
	processor = partial(_aggregate_aa, input, output, both)
	pool.map(processor, xrange(AMINO_ACID_COUNT))
	pool.close()
	pool.join()
	#for req_aa in xrange(0, AMINO_ACID_COUNT):

def _aggregate_filenm(input, output, midpath, filenm):
	def process_line(line, subdata):
		comps = line.split(";")
		key = PositionZone(Point3D(*comps[0].split(",")))
		if key in subdata:
			i = 0
			for x in comps[1].split(","):
				if not len(x.strip()): continue
				subdata[key][i] += int(x)
				i += 1
		else:
			subdata[key] = [int(x) for x in comps[1].split(",") if len(x.strip())]

	print filenm, "==============="
	if ".txt" not in filenm: return
	data = {}
	for subdir in os.listdir(input):
		if not os.path.isdir(os.path.join(input, subdir)): continue
		print subdir
		if midpath:
			desired_path = os.path.join(input, subdir, midpath, filenm)
		else:
			desired_path = os.path.join(input, subdir, filenm)
		if not os.path.exists(desired_path): continue
		with open(desired_path, "r") as file:
			for line in file: process_line(line, data)
		del file
	gc.collect()

	if midpath:
		if not os.path.exists(join(output, midpath)): os.mkdir(join(output, midpath))
		f = open(join(output, midpath, filenm), 'w')
	else:
		f = open(join(output, filenm), 'w')
	for pz, freq in data.iteritems():
		freqstr = ""
		for subfreq in freq:
			freqstr += str(subfreq) + ","
		freqstr = freqstr[:-1]
		f.write(str(pz.alpha_zone.x) + ", " + str(pz.alpha_zone.y) + ", " + str(pz.alpha_zone.z) + "; " + freqstr + "\n")
	f.close()

def aggregate_consolidated_data(input, output):
	if not os.path.exists(output): os.mkdir(output)
	print output
	midpath = None
	filenms = []
	for path in os.listdir(input):
		if os.path.isdir(os.path.join(input, path)):
			filenms = os.listdir(os.path.join(input, path))
			if "default" in filenms:
				midpath = "default"
				filenms = os.listdir(os.path.join(input, path, midpath))
			break
	pool = multiprocessing.Pool(processes=3, maxtasksperchild=1)
	processor = partial(_aggregate_filenm, input, output, midpath)
	pool.map(processor, filenms)
	pool.close()
	pool.join()


def medium_distributions(input, output):
	if not os.path.exists(output): os.mkdir(output)
	for req_aa in xrange(0, AMINO_ACID_COUNT):
		data = [0 for j in xrange(0, 100)]
		print "========", req_aa
		for subdir in os.listdir(input):
			if not os.path.isdir(os.path.join(input, subdir)): continue
			print subdir
			for subsubdir in os.listdir(os.path.join(input, subdir)):
				if subsubdir == "medium":
					for aafile in os.listdir(os.path.join(input, subdir, subsubdir)):
						if ".txt" not in aafile: continue
						aaid = int(aafile.replace(".txt", ""))
						if aaid != req_aa: continue
						with open(os.path.join(input, subdir, subsubdir, aafile), "r") as file:
							for line in file:
								if int(line.strip()) < len(data):
									data[int(line.strip())] += 1
						del file
						gc.collect()

		#Medium
		medium_path = join(output, "medium-" + str(req_aa) + ".txt")
		with open(medium_path, "w") as file:
			for i, x in enumerate(data):
				file.write(str(i) + ' ' + str(x) + '\n')
		del data[:]
		gc.collect()

import pdbstats

def determine_omits(all_pdbs, omitspath):
	"""Returns the total number of PDB files, and the ones which had omits and partial omits as a tuple."""
	with open(all_pdbs, 'r') as file:
		contents = {pdbid.strip() : 0 for pdbid in file}
	num_omits = 0
	partial_omits = 0
	total = 0
	with open(omitspath, 'r') as ofile:
		for line in ofile:
			words = line.strip().split()
			if words[0] == 'Processing':
				if contents[words[1]] == 0: total += 1
				contents[words[1]] = 1
			elif words[0] == '==========================Omit':
				assert words[1].replace("'", "") in contents, "PDBID %r not in contents" % words[1].replace("'", "")
				if contents[words[1].replace("'", "")] < 2: num_omits += 1
				contents[words[1].replace("'", "")] = 2
			elif words[0] == 'Partial':
				assert words[2].replace("'", "") in contents, "PDBID %r not in contents" % words[2].replace("'", "")
				if contents[words[2].replace("'", "")] < 3: partial_omits += 1
				contents[words[2].replace("'", "")] = 3
	for pdbid in contents:
		if contents[pdbid] > 0: continue
		print pdbid
	return (total, num_omits, partial_omits)

def aggregate_secondary_structures(input):
	"""Aggregates the contents of each secondary structure type in each folder into the directory specified in input."""
	data = {}
	for folder in os.listdir(input):
		if not os.path.isdir(os.path.join(input, folder)): continue
		print folder
		filenames = os.listdir(os.path.join(input, folder))
		for struct_file in filenames:
			if ".txt" not in struct_file: continue
			if struct_file not in data:
				data[struct_file] = {}
			with open(os.path.join(input, folder, struct_file), "r") as file:
				for line in file:
					pz = read_pz(line)
					freq = int(line.split(";")[-1])
					if pz in data[struct_file]:
						data[struct_file][pz] += freq
					else:
						data[struct_file][pz] = freq
	for filename in data:
		sorted_points = sorted(data[filename].items(), key=lambda x: -x[1])
		with open(os.path.join(input, filename), "w") as f:
			for pz, freq in sorted_points:
				f.write(str(pz.alpha_zone.x) + ", " + str(pz.alpha_zone.y) + ", " + str(pz.alpha_zone.z) + "; " +
						str(pz.x_axis.x) + ", " + str(pz.x_axis.y) + ", " + str(pz.x_axis.z) + "; " +
						str(pz.y_axis.x) + ", " + str(pz.y_axis.y) + ", " + str(pz.y_axis.z) + "; " +
						str(pz.z_axis.x) + ", " + str(pz.z_axis.y) + ", " + str(pz.z_axis.z) + "; " +
						str(freq) + "\n")

def aggregate_permissible_sequences(input, output, cutoff=0.05):
	"""Aggregates the contents of each permissible sequence folder into the output directory. Recommended to run first with a cutoff of 0 to see how many zones you get, then choose a cutoff accordingly."""
	data = {}
	for folder in os.listdir(input):
		basepath = os.path.join(input, folder)
		if not os.path.isdir(basepath): continue
		print folder
		if os.path.exists(os.path.join(basepath, "permissible_sequences")):
			if os.path.isdir(os.path.join(basepath, "permissible_sequences")):
				basepath = os.path.join(basepath, "permissible_sequences")
		filenames = os.listdir(basepath)
		for struct_file in filenames:
			if ".txt" not in struct_file: continue
			struct_type = struct_file[:-4]
			with open(os.path.join(basepath, struct_file), "r") as file:
				current_zones = None
				for line in file:
					if len(line.strip()) == 0:
						current_zones = None
						continue
					comps = line.split(";")
					if not current_zones:
						comps1 = comps[0].split(",")
						comps2 = comps[1].split(",")
						current_zones = (Point3D(*(float(c) for c in comps1)), Point3D(*(float(c) for c in comps2)))
					else:
						alpha_zone = Point3D(*(float(c) for c in comps[0].split(",")))
						if alpha_zone.magnitude() > 11.0: continue
						freq = int(comps[-1])
						if struct_type not in data:
							data[struct_type] = {}
						if current_zones not in data[struct_type]:
							data[struct_type][current_zones] = {}
						if alpha_zone in data[struct_type][current_zones]:
							data[struct_type][current_zones][alpha_zone] += freq
						else:
							data[struct_type][current_zones][alpha_zone] = freq
	def write_aggregate_freq_data(data_dict, filename):
		min_freq = 0
		if cutoff > 0 and cutoff < 1:
			all_frequencies = [c for coll in data_dict.values() for c in coll.values()]
			all_frequencies = sorted(all_frequencies, key=lambda x: -x)
			min_freq = all_frequencies[int(len(all_frequencies) * cutoff)]
		with open(os.path.join(output, filename), "w") as f:
			for zone1, zone2 in data_dict:
				if min_freq > 0:
					points_to_write = []
					for pz, freq in data_dict[(zone1, zone2)].iteritems():
						if freq > min_freq and next((x for x in data_dict if x[0] == pz or x[1] == pz), None) is not None:
							points_to_write.append((pz, freq))
				else:
					points_to_write = data_dict[(zone1, zone2)].items()
				if len(points_to_write) == 0: continue
				sorted_points = sorted(points_to_write, key=lambda x: -x[1])
				running_freq = 0
				f.write(str(zone1.x) + ", " + str(zone1.y) + ", " + str(zone1.z) + "; " +
						str(zone2.x) + ", " + str(zone2.y) + ", " + str(zone2.z) + "\n")
				for pz, freq in sorted_points:
					f.write(str(pz.x) + ", " + str(pz.y) + ", " + str(pz.z) + "; " + str(freq) + "\n")
					running_freq += freq
					if cutoff > 1 and running_freq > cutoff: break
				f.write("\n")
				del points_to_write, sorted_points
	if not os.path.exists(output): os.mkdir(output)
	if len(data) == 1:
		write_aggregate_freq_data(data.values()[0], "permissible_sequences.txt")
	else:
		for struct_type in data:
			write_aggregate_freq_data(data[struct_type], struct_type + ".txt")

def trim_secondary_structure_pzs(input, fraction=0.9):
	"""Overwrites the secondary structure files by keeping the position zones that contain the top 'fraction' of occurrences (default, top 90%)."""
	for path in os.listdir(input):
		if ".txt" not in path or ("helix" not in path and "sheet" not in path): continue
		print path
		data = []
		total_freq = 0
		with open(os.path.join(input, path), "r") as file:
			for line in file:
				pz = read_pz(line)
				freq = int(line.split(";")[-1])
				data.append((pz, freq))
				total_freq += freq
		cutoff = total_freq * fraction
		with open(os.path.join(input, path), "w") as f:
			cumulative_freq = 0
			for pz, freq in data:
				if cumulative_freq > cutoff: break
				f.write(str(pz.alpha_zone.x) + ", " + str(pz.alpha_zone.y) + ", " + str(pz.alpha_zone.z) + "; " +
					str(pz.x_axis.x) + ", " + str(pz.x_axis.y) + ", " + str(pz.x_axis.z) + "; " +
					str(pz.y_axis.x) + ", " + str(pz.y_axis.y) + ", " + str(pz.y_axis.z) + "; " +
					str(pz.z_axis.x) + ", " + str(pz.z_axis.y) + ", " + str(pz.z_axis.z) + "; " +
					str(freq) + "\n")
				cumulative_freq += freq

def sequence_score_distribution(input, output, weights=[], grain=20.0):
	"""Reads a file of the format
		PDBID,sequence length,scores
		.
		.
		.
		and prints a distribution of the net scores with resolution `grain`, weighted by the weights you provide, separated by sequence length. The output is a file for each segment length listing all the various scores."""
	sparc_scores = {
		3: {},
		4: {},
		5: {},
		6: {},
		7: {},
		8: {},
		9: {},
		10: {}
	}
	with open(input, "r") as file:
		for line in file:
			comps = line.strip().split(",")
			seqlen = int(comps[1])
			tot_score = 0.0
			invalid = False
			for i in xrange(2, len(comps)):
				weight = 1.0
				if i - 2 < len(weights):
					weight = weights[i - 2]
				subscore = float(comps[i])
				if subscore == 0.0:
					invalid = True
					break
				tot_score += subscore * weight
			if seqlen not in sparc_scores or invalid: continue
			with open(os.path.join(output, str(seqlen) + ".txt"), "a") as outfile:
				outfile.write(str(tot_score) + "\n")
			tot_score = math.floor(tot_score / grain) * grain
			if tot_score in sparc_scores[seqlen]:	sparc_scores[seqlen][tot_score] += 1
			else:									sparc_scores[seqlen][tot_score] = 1

	print "Total counts:"
	for i in xrange(3, 11):
		print len(sparc_scores[i])
	for i in xrange(3, 11):
		print "\n", i
		s = sorted(sparc_scores[i].items(), key=lambda x: x[0])
		for score, freq in s:
			if freq > 1:
				print score, freq

#MARK: CHARMM Comparison

def correlate_charmm_sparc(input):
	"""Prints a list of charmm scores and sparc scores separated by commas."""
	data = []
	with open(input, "r") as file:
		idx = 0
		for line in file:
			idx += 1
			if "ERROR" in line or "[" in line: continue
			comps = line.strip().split(",")
			if len(comps):
				data.append((float(comps[0]), 4.0 * float(comps[3]) + float(comps[4]) + 4.0 * float(comps[5]) + 5.0 * float(comps[6]), float(comps[1]), float(comps[2])))
	#data = sorted(data, key=lambda x: x[0])
	for i, d in enumerate(data):
		'''if i > 0:
			comp2 = str(float(d[1] - data[i - 1][1]) / float(d[2] - data[i - 1][2]))
		else: comp2 = ""
		comp1 = d[1] / d[2]
		print str(comp1) + "\t" + comp2'''
		print str(d[0]) + "\t" + str(d[1]) + "\t" + str(d[2]) + "\t" + str(d[3])

#MARK: Reference State

def _compute_ranges(radius, step, appendage=None):
	"""Computes the distance ranges for every position zone."""
	ranges = {}
	for point in Point3D.zero().iteroffsets(radius, step=step):
		if point.x >= 0.0:
			minx = point.x
			maxx = point.x + 1.0
		else:
			minx = point.x + 1.0
			maxx = point.x
		if point.y >= 0.0:
			miny = point.y
			maxy = point.y + 1.0
		else:
			miny = point.y + 1.0
			maxy = point.y
		if point.z >= 0.0:
			minz = point.z
			maxz = point.z + 1.0
		else:
			minz = point.z + 1.0
			maxz = point.z
		mindist = Point3D(minx, miny, minz).magnitude()
		maxdist = Point3D(maxx, maxy, maxz).magnitude()
		if mindist >= radius: continue
		if appendage is not None:
			ranges[point] = [mindist, maxdist, appendage]
		else:
			ranges[point] = (mindist, maxdist)
	return ranges

def generate_distance_constrained_bins(output, radius=10.0, step=1.0):
	"""In the file determined by output, writes one entry for every bucket in a sphere of the provided radius. For each bucket, it writes the number of buckets of side length 'step' that share distances within the range of distances in that bucket."""
	ranges = _compute_ranges(radius, step)
	print "Evaluated all points"
	with open(output, "w") as file:
		currx = 0.0
		for point in Point3D.zero().iteroffsets(radius, step=step):
			if point.x != currx:
				print point.x
				currx = point.x
			# We only write points in Octant I because the other points can be determined from these.
			# For instance, (-2, 4, -7) would have the same number as (1, 4, 6).
			if point not in ranges or point.x < 0.0 or point.y < 0.0 or point.z < 0.0: continue
			baserange = ranges[point]
			count = 0
			for point2, range in ranges.iteritems():
				if point2 != point and (range[0] <= baserange[1] or range[1] <= baserange[1]): #((range[0] >= baserange[0] and range[0] <= baserange[1]) or (range[1] >= baserange[0] and range[1] <= baserange[1])):
					count += 1
			file.write("{}, {}, {}; {}\n".format(point.x, point.y, point.z, count))
	print "Done"

def sparc_distance_constrained_bins(input, output, radius=10.0, step=1.0):
	"""This function performs a similar task as generate_distance_constrained_bins, except that the output is a file for each amino acid combination and for every zone for which there were instances in the dataset. The frequency provided is the number of instances in the opposite direction that could have had the same zone. In other words, for each position zone p_q, the output is the number of occurrences of any zone q_p for which |p_q| ~ |q_p|. The input should be a folder of frequencies."""
	
	if not os.path.exists(output): os.mkdir(output)
	files = os.listdir(input)
	percentage = 0
	alpha_frequencies = [[{} for i in xrange(AMINO_ACID_COUNT)] for j in xrange(AMINO_ACID_COUNT)]
	for n, indfile in enumerate(files):
		if indfile.find(".txt") == -1: continue
		if math.floor(float(n) / len(files) * 10) > percentage:
			percentage = math.floor(float(n) / len(files) * 10)
			print percentage * 10, "percent complete"
		tag1, tag2 = indfile[0:-4].split('-')
		tag1 = int(tag1)
		tag2 = int(tag2)
		with open(join(input, indfile), 'r') as file:
			for line in file:
				if ";" not in line and len(line.strip()) > 0:
					continue
				ptcomps, freq = line.strip().split(";")
				alpha = Point3D(*ptcomps.split(","))
				alpha_frequencies[tag1][tag2][alpha] = [float(freq), 0]
	#Now, given that data, go through each amino acid combo and find its reverses
	print "Done reading"
	ranges = _compute_ranges(radius, step, appendage=0)
	for i in xrange(AMINO_ACID_COUNT):
		for j in xrange(AMINO_ACID_COUNT):
			print i, j
			with open(os.path.join(output, str(i) + "-" + str(j) + ".txt"), "w") as file:
				for alpha, baserange in ranges.iteritems():
					if baserange[2] == 1 or alpha.x < 0.0 or alpha.y < 0.0 or alpha.z < 0.0: continue
					count = 0
					simdist = [alpha]
					for point2, range in ranges.iteritems():
						if ((range[0] >= baserange[0] and range[0] <= baserange[1]) or (range[1] >= baserange[0] and range[1] <= baserange[1])):
							if point2 in alpha_frequencies[j][i]:
								count += alpha_frequencies[j][i][point2][0]
							if point2 != alpha and point2 in alpha_frequencies[i][j] and range == baserange and point2.x >= 0.0 and point2.y >= 0.0 and point2.z >= 0.0:
								simdist.append(point2)
					for point in simdist:
						ranges[point][2] = 1
						file.write("{}, {}, {}; {}\n".format(point.x, point.y, point.z, count))
			for point, range in ranges.iteritems():
				range[2] = 0

def aggregate_possible_interactions(input, output, volume_step=1000.0):
	"""Aggregates the PDB data of the number of interactions each amino acid type pair could have had at consecutive, short-range, and nonconsecutive levels."""
	blocks = os.listdir(input)
	if not os.path.exists(output): os.mkdir(output)
	for i in xrange(AMINO_ACID_COUNT):
		for j in xrange(AMINO_ACID_COUNT):
			if i > j: continue
			print i, j
			data = {}
			for blockpath in blocks:
				if not os.path.isdir(os.path.join(input, blockpath)): continue
				aapath = os.path.join(input, blockpath, "{}-{}.txt".format(i, j))
				if not os.path.exists(aapath): continue
				with open(aapath, "r") as file:
					for line in file:
						comps = line.strip().split(";")
						volume = math.floor(float(comps[0]) / volume_step) * volume_step
						if volume not in data:
							data[volume] = [int(x) for x in comps[1].split(",")]
						else:
							for k, x in enumerate(comps[1].split(",")):
								data[volume][k] += int(x)
			with open(os.path.join(output, "{}-{}.txt".format(i, j)), "w") as file:
				sorted_data = sorted(data.items(), key=lambda x: x[0])
				for v, freqs in sorted_data:
					file.write("{};{},{},{},{},{},{}\n".format(v, freqs[0], freqs[1], freqs[2], freqs[3], freqs[4], freqs[5]))
