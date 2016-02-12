from proteinmath import *
from polypeptide import *
import os
import gc
import multiprocessing
from os.path import join
from functools import partial
from secondary_structure import *

nonconsecutive_mode = 'n'
consecutive_mode = 'c'
hydrophobicity_mode = 'p'
alphacarbon_mode = 'a'
alphacarbon_consec_mode = 'ac'
default_network_mode = 'd'	# Only supported by block_stats_network
sidechain_mode = 'sc'	# Only supported by block_stats_network
secondary_structure_mode = 'ss' # Only supported by block_stats_network

def block_stats(id, mode, input, urls, output):
	"""Calculates a block of PDB statistics based on the mode selected in the mode parameter. Modes are available above."""
	
	if mode == nonconsecutive_mode or mode == consecutive_mode or mode == alphacarbon_mode or mode == alphacarbon_consec_mode:
		data = [[[] for i in range(AMINO_ACID_COUNT)] for j in range(AMINO_ACID_COUNT)]
	elif mode == hydrophobicity_mode:
		data = [[] for i in range(AMINO_ACID_COUNT)]
		
	peptide = Polypeptide()
	for path in urls:
		if path.find("pdb") == -1: continue
		#print threading.currentThread().name, "processing", path, "..."
		print "Processing " + path + " (%d)..." % id
		try:
			peptide.read(join(input, path))
		except (AssertionError, ZeroDivisionError, ValueError):
			print "==========================Omit %r" % path
			continue

		reported = False
		try:
			for i, aa in enumerate(peptide.aminoacids):
				tag1 = aacode(aa.type)
				if mode == nonconsecutive_mode or mode == alphacarbon_mode:
					r = peptide.nearby_aa(peptide.aminoacids[i], 10.0, i, consec=False)
				elif mode == consecutive_mode or mode == alphacarbon_consec_mode:
					r = peptide.nearby_aa(peptide.aminoacids[i], 10.0, i, consec=True)
				elif mode == hydrophobicity_mode:
					r = peptide.nearby_aa(aa, 10.0, i)
					if tag1 < 22:
						data[tag1].append(len(r))
					else:
						print "Partial omit %r (unknown aa)." % path
						reported = True
						break
					continue
				for aa2 in r:
					tag2 = aacode(aa2.type)
					if tag1 < 22 and tag2 < 22:
						if mode == alphacarbon_mode or mode == alphacarbon_consec_mode:
							data[tag1][tag2].append(aa.tolocal(aa2.acarbon))
						else:
							data[tag1][tag2].append(aa.aa_position_zone(aa2))
					elif not reported:
						print "Partial omit %r (unknown aa)." % path
						reported = True
						break
		except (AssertionError, ZeroDivisionError, ValueError):
			if not reported:
				print "Partial omit %r (exception)." % path

		del peptide.aminoacids[:]
		peptide.hashtable = None
		gc.collect()
	del peptide

	if mode == nonconsecutive_mode or mode == consecutive_mode or mode == alphacarbon_mode or mode == alphacarbon_consec_mode:
		if output is not None and not os.path.exists(output): os.mkdir(output)
		for i in range(AMINO_ACID_COUNT):
			for j in range(AMINO_ACID_COUNT):
				if output is not None:
					f = open(join(output, "%d-%d.txt" % (i, j)), 'w')
				for datum in data[i][j]:
					if output is not None:
						f.write(str(datum).translate(None, "()") + "\n")
				if output is not None:
					f.close()
		del data
	elif mode == hydrophobicity_mode:
		if output is not None and not os.path.exists(output): os.mkdir(output)
		for i in range(AMINO_ACID_COUNT):
			if output is not None:
				f = open(join(output, "%d.txt" % i), 'w')
			for datum in data[i]:
				if output is not None:
					f.write(str(datum) + '\n')
			if output is not None:
				f.close()
		del data


def pool_block_stats(input, output, mode, (urls, i)=None):
	block_stats(i, mode, input, urls, join(output, str(i)))
	print "Done with batch", i

#12/8/14 8:44 pm: calculating the first 10 blocks (nullified)
#12/10/14 4:41 pm: calculated the first 3 blocks with improved memory management, starting blocks 3-22 inclusive
#12/10/14 9pm: it works! Now, blocks 23-end
#12/11/14 6:48 am: starting all blocks, consecutive

def pool_initializer():
	print "Starting", multiprocessing.current_process().name

def calculate_pdb_stats(input, output, mode=nonconsecutive_mode):
	
	if output is not None and not os.path.exists(output): os.mkdir(output)
	
	contents = os.listdir(input)
	block_size = 100
	i = 0
	blocks = int(math.ceil(len(contents) / block_size))
	print blocks, "blocks"
	
	# Single thread program used to debug using memory_profiler
	'''block_stats(0, input, [contents[1]], None)
	block_stats(1, input, [contents[2]], None)
	block_stats(2, input, [contents[3]], None)
	return'''
	pool = multiprocessing.Pool(processes=3, initializer=pool_initializer, maxtasksperchild=1)
	partial_block_stats = partial(pool_block_stats, input, output, mode)
	zipped = [(contents[block_size * k : min(block_size * (k + 1), len(contents))], k) for k in xrange(i, i + blocks)]
	#print zipped
	pool.map(partial_block_stats, zipped)
	pool.close()
	pool.join()
	print "done"

def calculate_hydrophobicity_stats(input, output):
	
	if not os.path.exists(output): os.mkdir(output)
	
	contents = os.listdir(input)
	block_size = 100
	i = 0
	blocks = int(math.ceil(len(contents) / block_size))
	print blocks, "blocks"
	
	# Single thread program used to debug using memory_profiler
	'''block_stats(0, input, [contents[1]], None)
		block_stats(1, input, [contents[2]], None)
		block_stats(2, input, [contents[3]], None)
		return'''
	pool = multiprocessing.Pool(processes=3, initializer=pool_initializer, maxtasksperchild=1)
	partial_block_stats = partial(pool_block_stats, input, output, hydrophobicity_mode)
	zipped = [(contents[block_size * k : min(block_size * (k + 1), len(contents))], k) for k in xrange(i, i + blocks)]
	#print zipped
	pool.map(partial_block_stats, zipped)
	pool.close()
	pool.join()
	print "done"

def calculate_pdb_alphacarbons(input, output, mode=alphacarbon_mode):
	
	if not os.path.exists(output): os.mkdir(output)
	
	contents = os.listdir(input)[:1000]
	block_size = 100
	i = 0
	blocks = int(math.ceil(len(contents) / block_size))
	print blocks, "blocks"
	
	pool = multiprocessing.Pool(processes=3, initializer=pool_initializer, maxtasksperchild=1)
	partial_block_stats = partial(pool_block_stats, input, output, mode)
	zipped = [(contents[block_size * k : min(block_size * (k + 1), len(contents))], k) for k in xrange(i, i + blocks)]
	#print zipped
	pool.map(partial_block_stats, zipped)
	pool.close()
	pool.join()
	print "done"

def calculate_amino_acid_counts(input):
	contents = os.listdir(input)
	data = [0 for i in xrange(AMINO_ACID_COUNT)]
	
	peptide = Polypeptide()
	for path in contents:
		if path.find("pdb") == -1: continue
		print "Processing %r..." % path
		try:
			peptide.read(join(input, path))
		except (AssertionError, ZeroDivisionError, ValueError):
			print "==========================Omit %r" % path
			continue
				
		for aa in peptide.aminoacids:
			tag = aacode(aa.type)
			if tag < 22: data[tag] += 1

		del peptide.aminoacids[:]
		peptide.hashtable = None
		gc.collect()
	del peptide

	for x in data:
		print x
	print "Sum:", sum(data)

# MARK: Fall 2015, new code for a larger sample size

import urllib2, time, numpy

def block_stats_network(id, pdbids, output, mode='d'):
	"""Calculates a block of PDB statistics, and puts all the data into a directory specified by output. For mode, you may pass default_network_mode, sidechain_mode, or secondary_structure_mode."""
	
	if mode == default_network_mode:
		data = {
			"nonconsec" : [[{} for i in range(AMINO_ACID_COUNT)] for j in range(AMINO_ACID_COUNT)],
			"short-range" : [[{} for i in range(AMINO_ACID_COUNT)] for j in range(AMINO_ACID_COUNT)],
			"consec" : [[{} for i in range(AMINO_ACID_COUNT)] for j in range(AMINO_ACID_COUNT)],
			"medium" : [[] for i in range(AMINO_ACID_COUNT)],
			"secondary" : {
				secondary_struct_helix + "1" : {}, secondary_struct_helix + "2" : {},
				secondary_struct_helix + "3" : {}, secondary_struct_helix + "4" : {},
				secondary_struct_helix + "5" : {}, secondary_struct_helix + "6" : {},
				secondary_struct_helix + "7" : {}, secondary_struct_helix + "8" : {},
				secondary_struct_helix + "9" : {}, secondary_struct_helix + "10" : {},
				secondary_struct_sheet + "0" : {}, secondary_struct_sheet + "1" : {},
				secondary_struct_sheet + "-1" : {}
			}
		}
	elif mode == sidechain_mode:
		data = [{} for i in range(AMINO_ACID_COUNT)]
	elif mode == secondary_structure_mode:
		data = {
			secondary_struct_helix + "1" : {}, secondary_struct_helix + "2" : {},
			secondary_struct_helix + "3" : {}, secondary_struct_helix + "4" : {},
			secondary_struct_helix + "5" : {}, secondary_struct_helix + "6" : {},
			secondary_struct_helix + "7" : {}, secondary_struct_helix + "8" : {},
			secondary_struct_helix + "9" : {}, secondary_struct_helix + "10" : {},
			secondary_struct_sheet + "0" : {}, secondary_struct_sheet + "1" : {},
			secondary_struct_sheet + "-1" : {}
		}
	
	def add_to_data(aa, aa2, tag1, tag2, key):
		if tag1 < 22 and tag2 < 22:
			pz = aa.aa_position_zone(aa2).alpha_zone
			if pz in data[key][tag1][tag2]:
				data[key][tag1][tag2][pz] += 1
			else:
				data[key][tag1][tag2][pz] = 1
			return True
		else:
			return False


	peptide = Polypeptide()
	for pdbid in pdbids:
		pdbid = pdbid.strip()
		print "Processing " + pdbid + " (%d)..." % id
		try:
			response = urllib2.urlopen('http://www.rcsb.org/pdb/files/' + pdbid + '.pdb')
			if mode == default_network_mode:
				peptide.read_file(response, secondary_structure=True, fillgaps=True)
			elif mode == sidechain_mode:
				peptide.read_file(response, otheratoms=True)
			elif mode == secondary_structure_mode:
				peptide.read_file(response, secondary_structure=True, fillgaps=True)
		except:
			print "==========================Omit %r" % pdbid
			continue
		response.close()
		del response

		reported = False
		try:
			if mode == secondary_structure_mode:
				for sec_struct in peptide.secondary_structures:
					for strand in sec_struct.strands:
						prev_aa = peptide.aminoacids[strand.start]
						for aa in peptide.aminoacids[strand.start + 1 : strand.end + 1]:
							if not prev_aa or not aa:
								prev_aa = aa
								continue
							pz = prev_aa.aa_position_zone(aa)
							if sec_struct.type + str(strand.identifiers[0]) not in data:
								continue
							if pz in data[sec_struct.type + str(strand.identifiers[0])]:
								data[sec_struct.type + str(strand.identifiers[0])][pz] += 1
							else:
								data[sec_struct.type + str(strand.identifiers[0])][pz] = 1
							prev_aa = aa
			else:
				for i, aa in enumerate(peptide.aminoacids):
					if aa == None: continue
					tag1 = aacode(aa.type)
					if mode == default_network_mode:
						'''#Nonconsec
						r = peptide.nearby_aa(peptide.aminoacids[i], 10.0, i, consec=False)
						for aa2 in r:
							tag2 = aacode(aa2.type)
							key = "nonconsec"
							if math.fabs(aa2.tag - aa.tag) <= 5:
								key = "short-range"
							if not add_to_data(aa, aa2, tag1, tag2, key):
								print "Partial omit %r (unknown aa)." % pdbid
								reported = True
								break'''
							
						#Consec
						r = peptide.nearby_aa(peptide.aminoacids[i], 10.0, i, consec=True)
						for aa2 in r:
							tag2 = aacode(aa2.type)
							key = "consec"
							#sec_struct = peptide.secondary_structure_aa(i)
							if False: #sec_struct:
								key = "secondary"
								pz = aa.aa_position_zone(aa2).alpha_zone
								sec_name = sec_struct[0].type + str(sec_struct[1].identifiers[0])
								if sec_name in data[key]:
									if pz in data[key][sec_name]:
										data[key][sec_name][pz] += 1
									else:
										data[key][sec_name][pz] = 1
							else:
								if not add_to_data(aa, aa2, tag1, tag2, key):
									print "Partial omit %r (unknown aa)." % pdbid
									reported = True
									break

						'''#Medium
						r = peptide.nearby_aa(aa, 10.0, i)
						if tag1 < 22:
							data["medium"][tag1].append(len(r))
						elif not reported:
							print "Partial omit %r (unknown aa)." % pdbid
							reported = True'''
				
					elif mode == sidechain_mode:
						if tag1 >= 22 and not reported:
							print "Partial omit %r (unknown aa)." % pdbid
							reported = True
							continue
						for atomname, location in aa.otheratoms.iteritems():
							if atomname in data[tag1]:
								data[tag1][atomname].append(location)
							else:
								data[tag1][atomname] = [location]

		except (AssertionError, ZeroDivisionError, ValueError):
			if not reported:
				print "Partial omit %r (exception)." % pdbid

		del peptide.aminoacids[:]
		peptide.hashtable = None
		gc.collect()
	del peptide

	if output is None: return
	if not os.path.exists(output): os.mkdir(output)

	if mode == default_network_mode:
		def write_data(key):
			write_path = join(output, key)
			if not os.path.exists(write_path): os.mkdir(write_path)
			for i in range(AMINO_ACID_COUNT):
				for j in range(AMINO_ACID_COUNT):
					f = open(join(write_path, "%d-%d.txt" % (i, j)), 'w')
					for pz, freq in data[key][i][j].iteritems():
						f.write(str(pz.x) + ", " + str(pz.y) + ", " + str(pz.z) + "; " + str(freq) + "\n")
					f.close()

		#write_data("nonconsec")
		write_data("consec")
		#write_data("short-range")
		
		'''#Secondary
		write_path = join(output, "secondary")
		if not os.path.exists(write_path): os.mkdir(write_path)
		for sec_struct_type in data["secondary"]:
			f = open(join(write_path, sec_struct_type + ".txt"), 'w')
			for pz, freq in data["secondary"][sec_struct_type].iteritems():
				f.write(str(pz.x) + ", " + str(pz.y) + ", " + str(pz.z) + "; " + str(freq) + "\n")
			f.close()

		#Medium
		medium_path = join(output, "medium")
		if not os.path.exists(medium_path): os.mkdir(medium_path)
		for i in range(AMINO_ACID_COUNT):
			f = open(join(medium_path, "%d.txt" % i), 'w')
			for datum in data["medium"][i]:
				f.write(str(datum) + '\n')
			f.close()'''
	elif mode == sidechain_mode:
		for i in range(AMINO_ACID_COUNT):
			f = open(join(output, "%d.txt" % i), 'w')
			for atomname, locations in data[i].iteritems():
				totalpoint = Point3D.zero()
				for location in locations:
					totalpoint = totalpoint.add(location)
				totalpoint = totalpoint.multiply(1.0 / len(locations))
				distance = 0.0
				for location in locations:
					distance += totalpoint.distanceto(location)
				distance /= float(len(locations))
				f.write(atomname + ',' + str(totalpoint.x) + ',' + str(totalpoint.y) + ',' + str(totalpoint.z) + ',' + str(distance) + '\n')
			f.close()
	elif mode == secondary_structure_mode:
		for struct_type in data:
			out_path = join(output, struct_type + ".txt")
			with open(out_path, "w") as f:
				for pz, freq in data[struct_type].iteritems():
					f.write(str(pz.alpha_zone.x) + ", " + str(pz.alpha_zone.y) + ", " + str(pz.alpha_zone.z) + "; " +
							str(pz.x_axis.x) + ", " + str(pz.x_axis.y) + ", " + str(pz.x_axis.z) + "; " +
							str(pz.y_axis.x) + ", " + str(pz.y_axis.y) + ", " + str(pz.y_axis.z) + "; " +
							str(pz.z_axis.x) + ", " + str(pz.z_axis.y) + ", " + str(pz.z_axis.z) + "; " +
							str(freq) + "\n")

def pool_block_stats_network(input, output, mode, (urls, i)=None):
	if output:
		if os.path.exists(join(output, str(i))): return
		block_stats_network(i, urls, join(output, str(i)), mode=mode)
	else:
		block_stats_network(i, urls, None, mode=mode)
	print "Done with batch", i
	gc.collect()
	time.sleep(60)


def calculate_pdb_stats_network(input, output):
	"""Calculates pdb stats by retrieving files with the PDB IDs listed in the newline-separated input file."""
	
	if output is not None and not os.path.exists(output): os.mkdir(output)
	
	with open(input, 'r') as file:
		contents = file.readlines()
	#contents = numpy.random.choice(contents, 500)
	block_size = 50
	i = 0
	blocks = int(math.ceil(len(contents) / block_size))
	print blocks, "blocks"
	
	pool = multiprocessing.Pool(processes=3, initializer=pool_initializer, maxtasksperchild=1)
	partial_block_stats = partial(pool_block_stats_network, input, output, default_network_mode)
	zipped = [(contents[block_size * k : min(block_size * (k + 1), len(contents))], k) for k in xrange(i, blocks)]
	#print zipped
	pool.map(partial_block_stats, zipped)
	pool.close()
	pool.join()
	print "done"

#MARK: Structure pairs

def generate_structurepairs(pdbsample, output):
	"""Generates PDB files according with each pair of amino acid types from the pdb sample. These PDB files contain all sidechain atoms from that structure."""
	
	if not os.path.exists(output):
		os.mkdir(output)

	with open(pdbsample, 'r') as file:
		contents = file.readlines()
	
	allcombos = [(i, j) for l in [[(i, j) for j in xrange(i, AMINO_ACID_COUNT)] for i in xrange(AMINO_ACID_COUNT)] for i, j in l if i != 3 and i != 8 and j != 3 and j != 8]
	
	
	for pdbid in contents:
		pdbid = pdbid.strip()
		print "Processing " + pdbid, ",", len(allcombos), "combinations left..."
		try:
			response = urllib2.urlopen('http://www.rcsb.org/pdb/files/' + pdbid + '.pdb')
			peptide = Polypeptide()
			peptide.read_file(response, otheratoms=True)
			response.close()
		except:
			print "==========================Omit %r" % pdbid
			continue
		del response

		for i, aa in enumerate(peptide.aminoacids):
			for aa2 in peptide.aminoacids[i + 1:]:
				type1 = aacode(aa.type)
				type2 = aacode(aa2.type)
				if type1 < type2: combo = (type1, type2)
				else:			  combo = (type2, type1)
				if combo in allcombos:
					print combo
					newpeptide = Polypeptide(preaas=[aa, aa2])
					newpeptide.center()
					with open(os.path.join(output, str(type1) + "-" + str(type2) + ".pdb"), "w") as file:
						file.write(newpeptide.pdb())
					allcombos.remove(combo)
					del newpeptide
		del peptide
		gc.collect()
		if len(allcombos) == 0:
			break
	print "Done"

#MARK: Secondary Structure Analysis

def analyze_secondary_structure(input, output, sample=-1):
	"""Retrieves files with the PDB IDs listed in the newline-separated input file, then determines the relative orientations going forward through each secondary structure in the A chain. The output is a file for each secondary structure type (10 types of helices, and initial, parallel, and anti-parallel sheets)."""

	if output is not None and not os.path.exists(output): os.mkdir(output)
	
	with open(input, 'r') as file:
		contents = file.readlines()
	
	if sample > 0:
		contents = random.sample(contents, sample)

	block_size = 50
	i = 0
	blocks = int(math.ceil(len(contents) / block_size))
	print blocks, "blocks"

	pool = multiprocessing.Pool(processes=3, initializer=pool_initializer, maxtasksperchild=1)
	partial_block_stats = partial(pool_block_stats_network, input, output, secondary_structure_mode)
	zipped = [(contents[block_size * k : min(block_size * (k + 1), len(contents))], k) for k in xrange(i, blocks)]

	pool.map(partial_block_stats, zipped)
	pool.close()
	pool.join()
	print "done"
