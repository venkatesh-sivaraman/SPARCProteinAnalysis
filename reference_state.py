"""This module represents a singleton for working with amino acid reference states."""

import os
from aminoacids import *

# This dictionary stores all the position zones in the first octant along with the number of
# other position zones which contain points at the same distances as each one.
radial_matches = {}

num_buckets = 8000

def load_reference_state(helper_path):
	"""This function performs a one-time initialization for the reference state. Currently, this version uses a text file that provides the number of zones with the same distance range as each position zone in the first octant."""
	global radial_matches
	print "Reading reference state data"
	radial_matches = {}
	with open(helper_path, 'r') as file:
		for line in file:
			comps = line.split(";")
			point = Point3D(*(float(x) for x in comps[0].split(",")))
			matches = int(comps[1])
			radial_matches[point] = matches
	print "Done reading reference state", len(radial_matches)

def is_initialized():
	"""Returns True if the radial_matches dictionary has been populated, indicating that the reference_state is ready to be computed, and False if not."""
	global radial_matches
	return len(radial_matches) > 0

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
	nd = radial_matches[pz]
	return 1.0 / (num_buckets * nd)