from proteinmath import *
from polypeptide import *
import os
import gc
import multiprocessing
from os.path import join
from functools import partial
from secondary_structure import *
from sparc_distribution import *
from charmm import *

nonconsecutive_mode = 'n'
consecutive_mode = 'c'
hydrophobicity_mode = 'p'
alphacarbon_mode = 'a'
alphacarbon_consec_mode = 'ac'
default_network_mode = 'd'	# Only supported by block_stats_network
sidechain_mode = 'sc'	# Only supported by block_stats_network
secondary_structure_mode = 'ss' # Only supported by block_stats_network
permissible_sequences_mode = 'ps' # Only supported by block_stats_network
permissible_sequences_ss_mode = 'pss'
score_sequences_mode = 'scs' # Only supported by block_stats_network
both_orientations_mode = 'bo'
num_interacting_mode = 'ni'
random_coil_mode = 'rc'

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
	pass #print "Starting", multiprocessing.current_process().name

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

import requests, time, numpy, urllib2, sys

def block_stats_network(id, pdbids, output, mode='d', **kwargs):
	"""Calculates a block of PDB statistics, and puts all the data into a directory specified by output. For mode, you may pass default_network_mode, sidechain_mode, or secondary_structure_mode."""
	
	if mode == default_network_mode:
		data = {
			"nonconsec" : [[{} for i in range(AMINO_ACID_COUNT)] for j in range(AMINO_ACID_COUNT)],
			"short-range" : [[{} for i in range(AMINO_ACID_COUNT)] for j in range(AMINO_ACID_COUNT)],
			"consec" : [[{} for i in range(AMINO_ACID_COUNT)] for j in range(AMINO_ACID_COUNT)],
			"consec+secondary" : [[{} for i in range(AMINO_ACID_COUNT)] for j in range(AMINO_ACID_COUNT)],
			"medium" : [[] for i in range(AMINO_ACID_COUNT)],
			"secondary" : {
				secondary_struct_helix + "1" : {}, secondary_struct_helix + "2" : {},
				secondary_struct_helix + "3" : {}, secondary_struct_helix + "4" : {},
				secondary_struct_helix + "5" : {}, secondary_struct_helix + "6" : {},
				secondary_struct_helix + "7" : {}, secondary_struct_helix + "8" : {},
				secondary_struct_helix + "9" : {}, secondary_struct_helix + "10" : {},
				secondary_struct_sheet + "0" : {}, secondary_struct_sheet + "1" : {},
				secondary_struct_sheet + "-1" : {}
			},
			"permissible" : {
				secondary_struct_helix + "1" : {}, secondary_struct_helix + "2" : {},
				secondary_struct_helix + "3" : {}, secondary_struct_helix + "4" : {},
				secondary_struct_helix + "5" : {}, secondary_struct_helix + "6" : {},
				secondary_struct_helix + "7" : {}, secondary_struct_helix + "8" : {},
				secondary_struct_helix + "9" : {}, secondary_struct_helix + "10" : {},
				secondary_struct_sheet + "0" : {}, secondary_struct_sheet + "1" : {},
				secondary_struct_sheet + "-1" : {}, "default": {}, "all": {}
			}
		}
	elif mode == num_interacting_mode:
		data = [[{} for i in range(AMINO_ACID_COUNT)] for j in range(AMINO_ACID_COUNT)]
	elif mode == both_orientations_mode:
		data = {
			"nonconsec" : [[{} for i in range(AMINO_ACID_COUNT)] for j in range(AMINO_ACID_COUNT)],
			"short-range" : [[{} for i in range(AMINO_ACID_COUNT)] for j in range(AMINO_ACID_COUNT)],
			"consec" : [[{} for i in range(AMINO_ACID_COUNT)] for j in range(AMINO_ACID_COUNT)],
			"consec+secondary" : [[{} for i in range(AMINO_ACID_COUNT)] for j in range(AMINO_ACID_COUNT)],
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
	elif mode == permissible_sequences_mode:
		#It will contain tuples (zone1, zone2) for an amino acid pair as keys, and dictionaries of possible zones for the pair's neighbors as values.
		data = {}
	elif mode == permissible_sequences_ss_mode:
		data = {
			secondary_struct_helix + "1" : {}, secondary_struct_helix + "2" : {},
			secondary_struct_helix + "3" : {}, secondary_struct_helix + "4" : {},
			secondary_struct_helix + "5" : {}, secondary_struct_helix + "6" : {},
			secondary_struct_helix + "7" : {}, secondary_struct_helix + "8" : {},
			secondary_struct_helix + "9" : {}, secondary_struct_helix + "10" : {},
			secondary_struct_sheet + "0" : {}, secondary_struct_sheet + "1" : {},
			secondary_struct_sheet + "-1" : {}, "default": {}, "all": {}
		}
	elif mode == random_coil_mode:
		data = {
			secondary_struct_helix + "1" : {}, secondary_struct_helix + "2" : {},
			secondary_struct_helix + "3" : {}, secondary_struct_helix + "4" : {},
			secondary_struct_helix + "5" : {}, secondary_struct_helix + "6" : {},
			secondary_struct_helix + "7" : {}, secondary_struct_helix + "8" : {},
			secondary_struct_helix + "9" : {}, secondary_struct_helix + "10" : {},
			secondary_struct_sheet + "0" : {}, secondary_struct_sheet + "1" : {},
			secondary_struct_sheet + "-1" : {}, "default": {}, "all": {}
		}
	
	def add_to_data(aa, aa2, tag1, tag2, key, both=False):
		if tag1 < 22 and tag2 < 22:
			pz = aa.aa_position_zone(aa2).alpha_zone
			if both:
				pz2 = aa2.aa_position_zone(aa).alpha_zone
				if (pz, pz2) in data[key][tag1][tag2]:
					data[key][tag1][tag2][(pz, pz2)] += 1
				else:
					data[key][tag1][tag2][(pz, pz2)] = 1
			else:
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
			#response = urllib2.urlopen('http://www.rcsb.org/pdb/files/' + pdbid + '.pdb')
			r = requests.get('http://www.rcsb.org/pdb/files/' + pdbid + '.pdb', headers={'user-agent': 'Mozilla/5.0'})
			response = r.text.split('\n')
			if len(response) <= 1:
				assert False, "Got one-liner for PDB {}: {}".format(pdbid, r.text)
			if mode == default_network_mode or mode == permissible_sequences_mode or mode == permissible_sequences_ss_mode or mode == both_orientations_mode:
				peptide.read_file(response, secondary_structure=True, fillgaps=True)
			elif mode == sidechain_mode:
				peptide.read_file(response, otheratoms=True)
			elif mode == secondary_structure_mode:
				peptide.read_file(response, secondary_structure=True, fillgaps=True, fillends=True)
			elif mode == random_coil_mode:
				peptide.read_file(response, secondary_structure=True)
			else:
				peptide.read_file(response, secondary_structure=True, fillgaps=True)
			del response
		except Exception as e:
			print "==========================Omit %r - %r" % (pdbid, e)
			continue
		#response.close()
		#del response

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
							pz = prev_aa.tolocal(aa.acarbon).floor()
							pz2 = aa.tolocal(prev_aa.acarbon).floor()
							if sec_struct.type + str(strand.identifiers[0]) not in data:
								continue
							if pz in data[sec_struct.type + str(strand.identifiers[0])]:
								data[sec_struct.type + str(strand.identifiers[0])][pz] += 1
							else:
								data[sec_struct.type + str(strand.identifiers[0])][pz] = 1
							if pz2 in data[sec_struct.type + str(strand.identifiers[0])]:
								data[sec_struct.type + str(strand.identifiers[0])][pz2] += 1
							else:
								data[sec_struct.type + str(strand.identifiers[0])][pz2] = 1
							prev_aa = aa
			elif mode == num_interacting_mode:
				#For each amino acid type, the dictionary will have the floored volume as its key, and a list of [observed interactions consec, obs short-range, obs long-range, possible interactions consec, possible short-range, possible long-range] as the value.
				volume = math.floor(peptide.volume())
				for aa in peptide.aminoacids:
					if aa == None: continue
					tag1 = aacode(aa.type)
					#Find possible interactions as well as actually contacting residues
					for aa2 in peptide.aminoacids:
						if not aa2: continue
						tag2 = aacode(aa2.type)
						if tag1 > tag2 or (tag1 == tag2 and aa.tag > aa2.tag): continue
						try:
							if volume not in data[tag1][tag2]:
								data[tag1][tag2][volume] = [0, 0, 0, 0, 0, 0]
							idx = 0
							if math.fabs(aa.tag - aa2.tag) == 1:
								idx = 0
							elif math.fabs(aa.tag - aa2.tag) <= 5:
								idx = 1
							else:
								idx = 2
							
							if aa.acarbon.distanceto(aa2.acarbon) <= 10.0:
								data[tag1][tag2][volume][idx] += 1
							data[tag1][tag2][volume][idx + 3] += 1
						except:
							print "Partial omit -", pdbid
			elif mode == random_coil_mode:
				iterations = 100
				if "iterations" in kwargs: iterations = kwargs["iterations"]
				assert "permissions" in kwargs, "No permissions object provided to block stats. Use the keyword argument 'permissions'."
				permissions = kwargs["permissions"]
				sec_struct_permissions = None
				if "sec_struct_permissions" in kwargs: sec_struct_permissions = kwargs["sec_struct_permissions"]
				for n in xrange(iterations):
					#if n == iterations / 2: print n, "of", iterations, "iterations"
					worked = peptide.randomcoil(permissions=permissions, struct_permissions=sec_struct_permissions, repeat_cutoff=10)
					if not worked: continue
					for i, aa in enumerate(peptide.aminoacids):
						r = peptide.nearby_aa(aa, 10.0, i)
						for j, aa2 in enumerate(r):
							pz = aa.aa_position_zone(aa2).alpha_zone
							separation = min(int(math.fabs(j - i)), 6) - 1
							sec_struct = peptide.secondary_structure_aa(i)
							sec_struct_2 = peptide.secondary_structure_aa(j)
							if sec_struct and sec_struct_2 and sec_struct[1].start == sec_struct_2[1].start:
								sec_name = sec_struct[0].type + str(sec_struct[1].identifiers[0])
								if pz in data[sec_name]:
									data[sec_name][pz][separation] += 1
								else:
									data[sec_name][pz] = [0 for i in xrange(1, 7)]
									data[sec_name][pz][separation] = 1
							else:
								if pz in data["default"]:
									data["default"][pz][separation] += 1
								else:
									data["default"][pz] = [0 for i in xrange(1, 7)]
									data["default"][pz][separation] = 1
							if pz in data["all"]:
								data["all"][pz][separation] += 1
							else:
								data["all"][pz] = [0 for i in xrange(1, 7)]
								data["all"][pz][separation] = 1
			else:
				for i, aa in enumerate(peptide.aminoacids):
					if aa == None: continue
					tag1 = aacode(aa.type)
					if mode == permissible_sequences_mode or mode == permissible_sequences_ss_mode or mode == default_network_mode:
						if i == 0 or i == len(peptide.aminoacids) - 1: continue
						
						#Figure out the orientation pair between the last consecutive residue and this one
						prev_aa = peptide.aminoacids[i - 1]
						if not prev_aa: continue
						zone1 = prev_aa.aa_position_zone(aa).alpha_zone
						zone2 = aa.aa_position_zone(prev_aa).alpha_zone
					
						sec_struct = peptide.secondary_structure_aa(i)
						do_next = True
						do_prev = True
						if mode == default_network_mode:
							cover_data = data["permissible"]
						else:
							cover_data = data
						if mode == permissible_sequences_ss_mode or mode == default_network_mode:
							if sec_struct:
								if sec_struct[1].start > i - 1: continue
								if sec_struct[1].start > i - 2: do_prev = False
								if sec_struct[1].end < i + 1: do_next = False
								sec_name = sec_struct[0].type + str(sec_struct[1].identifiers[0])
								if sec_name not in cover_data: continue
								data_lists = [cover_data[sec_name]]
							else:
								data_lists = [cover_data["default"]]
							data_lists.append(cover_data["all"])
						else:
							data_list = [cover_data]
						if peptide.aminoacids[i + 1] and do_next:
							#Figure out the orientations between the next consecutive residue and this one
							next_zone = aa.aa_position_zone(peptide.aminoacids[i + 1]).alpha_zone
							for data_list in data_lists:
								if (zone1, zone2) not in data_list:
									data_list[(zone1, zone2)] = {}
								if next_zone in data_list[(zone1, zone2)]:
									data_list[(zone1, zone2)][next_zone] += 1
								else:
									data_list[(zone1, zone2)][next_zone] = 1
				
						if peptide.aminoacids[i - 2] and do_prev:
							#Figure out the orientations between the last two consecutive residues
							prev_zone = prev_aa.aa_position_zone(peptide.aminoacids[i - 2]).alpha_zone
							for data_list in data_lists:
								if (zone1, zone2) not in data_list:
									data_list[(zone1, zone2)] = {}
								if prev_zone in data_list[(zone1, zone2)]:
									data_list[(zone1, zone2)][prev_zone] += 1
								else:
									data_list[(zone1, zone2)][prev_zone] = 1
				
					if mode == default_network_mode:
						#Nonconsec
						r = peptide.nearby_aa(peptide.aminoacids[i], 10.0, i, consec=False)
						for aa2 in r:
							tag2 = aacode(aa2.type)
							key = "nonconsec"
							if math.fabs(aa2.tag - aa.tag) <= 5:
								key = "short-range"
							if not add_to_data(aa, aa2, tag1, tag2, key):
								print "Partial omit %r (unknown aa)." % pdbid
								reported = True
								break
							
						#Consec
						r = peptide.nearby_aa(peptide.aminoacids[i], 10.0, i, consec=True)
						for aa2 in r:
							tag2 = aacode(aa2.type)
							key = "consec"
							sec_struct = peptide.secondary_structure_aa(i)
							if sec_struct:
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
							if not add_to_data(aa, aa2, tag1, tag2, "consec+secondary"):
								print "Partial omit %r (unknown aa)." % pdbid
								reported = True
								break

						#Medium
						r = peptide.nearby_aa(aa, 10.0, i)
						if tag1 < 22:
							data["medium"][tag1].append(len(r))
						elif not reported:
							print "Partial omit %r (unknown aa)." % pdbid
							reported = True

					if mode == both_orientations_mode:
						#Nonconsec
						r = peptide.nearby_aa(peptide.aminoacids[i], 10.0, i, consec=False)
						for aa2 in r:
							tag2 = aacode(aa2.type)
							if tag1 > tag2 or (tag1 == tag2 and aa.tag > aa2.tag): continue
							key = "nonconsec"
							if math.fabs(aa2.tag - aa.tag) <= 5:
								key = "short-range"
							if not add_to_data(aa, aa2, tag1, tag2, key, both=True):
								print "Partial omit %r (unknown aa)." % pdbid
								reported = True
								break
							
						#Consec
						r = peptide.nearby_aa(peptide.aminoacids[i], 10.0, i, consec=True)
						for aa2 in r:
							tag2 = aacode(aa2.type)
							if tag1 > tag2 or (tag1 == tag2 and aa.tag > aa2.tag): continue
							key = "consec"
							sec_struct = peptide.secondary_structure_aa(i)
							if sec_struct:
								key = "secondary"
								pz = aa.aa_position_zone(aa2).alpha_zone
								pz2 = aa2.aa_position_zone(aa).alpha_zone
								sec_name = sec_struct[0].type + str(sec_struct[1].identifiers[0])
								if sec_name in data[key]:
									if (pz, pz2) in data[key][sec_name]:
										data[key][sec_name][(pz, pz2)] += 1
									else:
										data[key][sec_name][(pz, pz2)] = 1
							else:
								if not add_to_data(aa, aa2, tag1, tag2, key, both=True):
									print "Partial omit %r (unknown aa)." % pdbid
									reported = True
									break
							if not add_to_data(aa, aa2, tag1, tag2, "consec+secondary", both=True):
								print "Partial omit %r (unknown aa)." % pdbid
								reported = True
								break

					if mode == sidechain_mode:
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

	def write_data(key, twopzs=False):
		write_path = join(output, key)
		if not os.path.exists(write_path): os.mkdir(write_path)
		for i in range(AMINO_ACID_COUNT):
			for j in range(AMINO_ACID_COUNT):
				f = open(join(write_path, "%d-%d.txt" % (i, j)), 'w')
				if twopzs:
					for pzs, freq in data[key][i][j].iteritems():
						f.write(str(pzs[0].x) + ", " + str(pzs[0].y) + ", " + str(pzs[0].z) + ", " + str(pzs[1].x) + ", " + str(pzs[1].y) + ", " + str(pzs[1].z) + "; " + str(freq) + "\n")
				else:
					for pz, freq in data[key][i][j].iteritems():
						f.write(str(pz.x) + ", " + str(pz.y) + ", " + str(pz.z) + "; " + str(freq) + "\n")
				f.close()

	if mode == default_network_mode:
		write_data("nonconsec")
		write_data("consec")
		write_data("short-range")
		write_data("consec+secondary")
		
		#Secondary
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
			f.close()

		#Permissible sequences
		permissibles_path = os.path.join(output, "permissible_sequences")
		if not os.path.exists(permissibles_path): os.mkdir(permissibles_path)
		for sec_struct_type in data["permissible"]:
			with open(os.path.join(permissibles_path, sec_struct_type + ".txt"), "w") as file:
				for zone1, zone2 in data["permissible"][sec_struct_type]:
					file.write("{}, {}, {}; {}, {}, {}\n".format(zone1.x, zone1.y, zone1.z, zone2.x, zone2.y, zone2.z))
					for pz, freq in data["permissible"][sec_struct_type][(zone1, zone2)].iteritems():
						file.write(str(pz.x) + ", " + str(pz.y) + ", " + str(pz.z) + "; " + str(freq) + "\n")
					file.write("\n")

	elif mode == random_coil_mode:
		for sec_struct_type in data:
			f = open(join(output, sec_struct_type + ".txt"), "w")
			for pz, freqs in data[sec_struct_type].iteritems():
				freqstr = ""
				for freq in freqs: freqstr += str(freq) + ","
				freqstr = freqstr[:-1]
				f.write("{:.0f},{:.0f},{:.0f};{}\n".format(pz.x, pz.y, pz.z, freqstr))
			f.close()
	elif mode == both_orientations_mode:
		write_data("nonconsec", twopzs=True)
		write_data("consec", twopzs=True)
		write_data("short-range", twopzs=True)
		write_data("consec+secondary", twopzs=True)
		
		#Secondary
		write_path = join(output, "secondary")
		if not os.path.exists(write_path): os.mkdir(write_path)
		for sec_struct_type in data["secondary"]:
			f = open(join(write_path, sec_struct_type + ".txt"), 'w')
			for pzs, freq in data["secondary"][sec_struct_type].iteritems():
				f.write(str(pzs[0].x) + ", " + str(pzs[0].y) + ", " + str(pzs[0].z) + ", " + str(pzs[1].x) + ", " + str(pzs[1].y) + ", " + str(pzs[1].z) + "; " + str(freq) + "\n")
			f.close()
	elif mode == num_interacting_mode:
		for i in range(AMINO_ACID_COUNT):
			for j in range(AMINO_ACID_COUNT):
				f = open(join(output, "%d-%d.txt" % (i, j)), 'w')
				for volume, counts in data[i][j].iteritems():
					f.write("{};{},{},{},{},{},{}\n".format(volume, counts[0], counts[1], counts[2], counts[3], counts[4], counts[5]))
				f.close()

	elif mode == permissible_sequences_mode:
		write_path = join(output, "permissible_sequences.txt")
		with open(write_path, "w") as file:
			for zone1, zone2 in data:
				file.write("{}, {}, {}; {}, {}, {}\n".format(zone1.x, zone1.y, zone1.z, zone2.x, zone2.y, zone2.z))
				for pz, freq in data[(zone1, zone2)].iteritems():
					file.write(str(pz.x) + ", " + str(pz.y) + ", " + str(pz.z) + "; " + str(freq) + "\n")
				file.write("\n")
	elif mode == permissible_sequences_ss_mode:
		for sec_struct_type in data:
			with open(os.path.join(output, sec_struct_type + ".txt"), "w") as file:
				for zone1, zone2 in data[sec_struct_type]:
					file.write("{}, {}, {}; {}, {}, {}\n".format(zone1.x, zone1.y, zone1.z, zone2.x, zone2.y, zone2.z))
					for pz, freq in data[sec_struct_type][(zone1, zone2)].iteritems():
						file.write(str(pz.x) + ", " + str(pz.y) + ", " + str(pz.z) + "; " + str(freq) + "\n")
					file.write("\n")
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
					f.write(str(pz.x) + ", " + str(pz.y) + ", " + str(pz.z) + "; " + str(freq) + "\n")

def pool_block_stats_network(input, output, mode, (urls, i)=None):
	if output:
		if os.path.exists(join(output, str(i))): return
		block_stats_network(i, urls, join(output, str(i)), mode=mode)
	else:
		block_stats_network(i, urls, None, mode=mode)
	print "Done with batch", i
	gc.collect()
	time.sleep(60)


def calculate_pdb_stats_network(input, output, mode=default_network_mode):
	"""Calculates pdb stats by retrieving files with the PDB IDs listed in the newline-separated input file."""
	
	if output is not None and not os.path.exists(output): os.mkdir(output)
	
	with open(input, 'r') as file:
		contents = file.readlines()
	#contents = numpy.random.choice(contents, 500)
	block_size = 50
	i = 0
	blocks = int(math.ceil(len(contents) / block_size))
	if mode == default_network_mode:
		time = 3.0 * len(contents)
		time_unit = "seconds"
		if time > 60:
			time /= 60.0
			time_unit = "minutes"
			if time > 60:
				time /= 60.0
				time_unit = "hours"
		print "NOTE: The calculation process will be performed in 50-structure blocks. Generally this takes roughly 3 seconds per structure, so with {} blocks your job will likely take about {:.2f} {}.".format(blocks, time, time_unit)
	
	pool = multiprocessing.Pool(processes=3, initializer=pool_initializer, maxtasksperchild=1)
	partial_block_stats = partial(pool_block_stats_network, input, output, mode)
	zipped = [(contents[block_size * k : min(block_size * (k + 1), len(contents))], k) for k in xrange(i, blocks)]
	#print zipped
	pool.map(partial_block_stats, zipped)
	pool.close()
	pool.join()
	print "done"

def num_possible_interacting(pdbsample, output):
	"""Writes the number of possible interactions each amino acid pair could have in the PDB sample."""
	calculate_pdb_stats_network(pdbsample, output, mode=num_interacting_mode)

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

def sparc_score_sequences(pdb_path, output, distributions, stop_pt=None, charmm_scores=False):
	
	tmpdir = os.path.join(os.path.dirname(pdb_path), "sp_tmp")
	if charmm_scores:
		if os.path.exists(tmpdir):
			shutil.rmtree(tmpdir)
		os.mkdir(tmpdir)

	peptide = Polypeptide()
	pdbids = []
	with open(pdb_path, "r") as file:
		found_stopping_pt = (stop_pt == None or len(stop_pt) == 0)
		for line in file:
			if len(line.strip()) == 0: continue
			if not found_stopping_pt and line.strip() == stop_pt:
				found_stopping_pt = True
				continue
			if not found_stopping_pt: continue
			pdbids.append(line.strip())

	for pdbid in pdbids:
		pdbid = pdbid.strip()
		if random.randint(0, 100) < 80: continue
		print "Processing " + pdbid + "..."
		try:
			response = urllib2.urlopen('http://www.rcsb.org/pdb/files/' + pdbid + '.pdb')
			peptide.read_file(response, secondary_structure=True, fillgaps=True, otheratoms=charmm_scores)
		except:
			print "==========================Omit %r" % pdbid
			continue
		response.close()
		del response

		for i in xrange(len(peptide.aminoacids)):
			if random.randint(0, 100) < 90: continue
			length = random.randint(3, 10)
			length = min(len(peptide.aminoacids), i + length) - i
			try:
				print "Evaluating segment", i, "-", i + length
				scores = [d.score(peptide, peptide.aminoacids[i:i + length], isolate=True) for d in distributions]
				if charmm_scores == True:
					npeptide = Polypeptide()
					npeptide.add_aas(peptide.aminoacids[i:i + length])
					cscores = charmm_score(npeptide, os.path.join(tmpdir, "sequence.pdb"), mdpfile="../../default_nonminimized.mdp")
					for c in cscores:
						scores.append(c)
				with open(output, "a") as file:
					file.write("{},{},".format(pdbid, length))
					scorestr = ""
					for s in scores: scorestr += str(s) + ","
					file.write(scorestr[:-1] + "\n")
			except:
				print "Partial omit -", pdbid

		del peptide.aminoacids[:]
		peptide.hashtable = None
		gc.collect()
	del peptide
	if charmm_scores:
		shutil.rmtree(tmpdir)

def pool_pdb_reference_state(input, output, mode, iterations, permissions, sec_struct_permissions, (urls, i)=None):
	if output:
		if os.path.exists(join(output, str(i))): return
		block_stats_network(i, urls, join(output, str(i)), mode=mode, iterations=iterations, permissions=permissions, sec_struct_permissions=sec_struct_permissions)
	else:
		block_stats_network(i, urls, None, mode=mode, iterations=iterations, permissions=permissions, sec_struct_permissions=sec_struct_permissions)
	print "Done with batch", i
	gc.collect()
	time.sleep(60)



def pdb_reference_state(input, output, permissions, sec_struct_permissions, iterations=1):
	"""Calculates the SPARC reference state by retrieving files with the PDB IDs listed in the newline-separated input file. Each structure undergoes a default of 20 random coils, and the orientations found in each structure are recorded."""
	
	if output is not None and not os.path.exists(output): os.mkdir(output)
	
	with open(input, 'r') as file:
		contents = file.readlines()
	#contents = numpy.random.choice(contents, 500)
	block_size = 5
	i = 0
	blocks = int(math.ceil(len(contents) / block_size))
	
	pool = multiprocessing.Pool(processes=3, initializer=pool_initializer, maxtasksperchild=1)
	partial_block_stats = partial(pool_pdb_reference_state, input, output, random_coil_mode, iterations, permissions, sec_struct_permissions)
	zipped = [(contents[block_size * k : min(block_size * (k + 1), len(contents))], k) for k in xrange(i, 1)] #blocks
	#print zipped
	map(partial_block_stats, zipped)
	pool.close()
	pool.join()
	print "done"

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

def generate_permissible_sequences(input, output, sample=-1):
	"""Retrieves files with the PDB IDs listed in the newline-separated input file, then determines the sets of permissible orientations between triplets of amino acids. The format of the return files will be 
		zone1; zone2
		pz; freq
		pz; freq
		pz; freq
		.
		.
		.
		
		zone1; zone2
		pz; freq
		...
		For example, given amino acids ABCD, zone1 and zone2 refer to B and C's orientation w.r.t. each other. pz could refer to the orientation between C and D or B and A."""

	if output is not None and not os.path.exists(output): os.mkdir(output)
	
	with open(input, 'r') as file:
		contents = file.readlines()
	
	if sample > 0:
		contents = random.sample(contents, sample)

	block_size = 50
	i = 0
	blocks = int(math.ceil(len(contents) / float(block_size)))

	pool = multiprocessing.Pool(processes=3, initializer=pool_initializer, maxtasksperchild=1)
	partial_block_stats = partial(pool_block_stats_network, input, output, permissible_sequences_ss_mode)
	zipped = [(contents[block_size * k : min(block_size * (k + 1), len(contents))], k) for k in xrange(i, blocks)]
	pool.map(partial_block_stats, zipped)
	pool.close()
	pool.join()
	print "done"

