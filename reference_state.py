"""This module represents a singleton for working with amino acid reference states."""

import os
from aminoacids import *
import loading_indicator

# This dictionary stores all the position zones in the first octant along with the number of
# other position zones which contain points at the same distances as each one.
radial_matches = {}

possible_interactions = [[{} for i in xrange(AMINO_ACID_COUNT)] for j in xrange(AMINO_ACID_COUNT)]

num_buckets = 5111

def load_reference_state(helper_path):
	"""This function performs a one-time initialization for the reference state. Currently, this version uses a text file that provides the number of zones with the same distance range as each position zone in the first octant."""
	global radial_matches
	radial_matches = {}
	loading_indicator.add_loading_data(1)
	with open(helper_path, 'r') as file:
		for line in file:
			comps = line.split(";")
			point = Point3D(*(float(x) for x in comps[0].split(",")))
			matches = int(comps[1])
			radial_matches[point] = matches
	loading_indicator.update_progress(1)

def load_possible_interactions(interactions_path):
	"""Performs a one-time initialization for the reference contact frequencies of various amino acid types."""
	global possible_interactions
	paths = os.listdir(interactions_path)
	loading_indicator.add_loading_data(len(paths))
	for path in paths:
		loading_indicator.update_progress(1)
		if ".txt" not in path or path[0] == ".": continue
		tag1, tag2 = path[:-4].split('-')
		tag1 = int(tag1)
		tag2 = int(tag2)
		with open(os.path.join(interactions_path, path), "r") as file:
			for line in file:
				comps = line.split(";")
				if len(comps) < 2: continue
				possible_interactions[tag1][tag2][float(comps[0])] = [int(x) for x in comps[1].split(",")]

def is_initialized():
	"""Returns True if the radial_matches dictionary has been populated, indicating that the reference_state is ready to be computed, and False if not."""
	global radial_matches, possible_interactions
	return len(radial_matches) > 0 and len(possible_interactions[0][0]) > 0

def position_ref(pz):
	"""Calculates the reference state for two amino acids given the location of one alpha carbon in the other's local coordinate system. Currently, that reference state refers to the expression
			1
		--------
		N * N(d)
		where N is the total number of buckets, d is the distance between the two amino acids, and N(d) is the number of buckets that contain points with distance d from the origin."""
	global radial_matches, num_buckets
	if len(radial_matches) == 0:
		print "Asked for a reference state without providing the helper files. Call load_reference_state first!"
		return 0
	pz = pz.floor()
	if pz.x < 0.0: pz.x = -pz.x - 1.0
	if pz.y < 0.0: pz.y = -pz.y - 1.0
	if pz.z < 0.0: pz.z = -pz.z - 1.0
	if pz not in radial_matches: return 1.0
	nd = radial_matches[pz]
	return 1.0 / (num_buckets * nd)

def contact_probability(type1, type2, volume, interaction_type):
	"""Returns the probability of two amino acids of the provided types being in contact, given that they are both found within a protein of the provided volume.
		Interaction types: 0 - consecutive, 1 - short-range, 2 - long-range"""
	global possible_interactions
	volume = math.floor(volume / 1000.0) * 1000.0
	if volume not in possible_interactions[type1][type2]: return 1.0
	contacts = possible_interactions[type1][type2][volume]
	if contacts[interaction_type + 3] == 0: return 1.0
	return contacts[interaction_type] / float(contacts[interaction_type + 3])