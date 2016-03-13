from aminoacids import *
import os
from os.path import join
import random
from numpy import random as nprand
from secondary_structure import *

def _compute_bounds(p, mag):
	''' 1. x = p.x
		2. x = p.x + 1
		3. y = p.y
		4. y = p.y + 1
		5. z = p.z
		6. z = p.z + 1
		
		x^2 + y^2 + z^2 = a^2
		1/3: z^2 = a^2 - p.x^2 - p.y^2
		2/3: z^2 = a^2 - (p.x + 1)^2 - p.y^2
		1/4: z^2 = a^2 - p.x^2 - (p.y + 1)^2
		2/4: z^2 = a^2 - (p.x + 1)^2 - (p.y + 1)^2
		1/5: y^2 = a^2 - p.x^2 - p.z^2
		2/5: y^2 = a^2 - (p.x + 1)^2 - p.z^2
		1/6: y^2 = a^2 - p.x^2 - (p.z + 1)^2
		2/6: y^2 = a^2 - (p.x + 1)^2 - (p.z + 1)^2
		3/5: x^2 = a^2 - p.y^2 - p.z^2
		4/5: x^2 = a^2 - (p.y + 1)^2 - p.z^2
		3/6: x^2 = a^2 - p.y^2 - (p.z + 1)^2
		4/6: x^2 = a^2 - (p.y + 1)^2 - (p.z + 1)^2'''
	print p, mag
	xs = [mag ** 2 - p.y ** 2 - p.z ** 2,
		  mag ** 2 - (p.y + 1.0) ** 2 - p.z ** 2,
		  mag ** 2 - p.y ** 2 - (p.z + 1.0) ** 2,
		  mag ** 2 - (p.y + 1.0) ** 2 - (p.z + 1.0) ** 2]
	if p.x < 0.0: c = -1.0
	else: c = 1.0
	xs = [min(max(c * math.sqrt(x), p.x), p.x + 1.0) for x in xs if x >= 0.0]
	idx = 0
	while idx < len(xs):
		if xs[idx] not in xs[:idx]: idx += 1
		else: del xs[idx]
	ys = [mag ** 2 - p.x ** 2 - p.z ** 2,
		  mag ** 2 - (p.x + 1.0) ** 2 - p.z ** 2,
		  mag ** 2 - p.x ** 2 - (p.z + 1.0) ** 2,
		  mag ** 2 - (p.x + 1.0) ** 2 - (p.z + 1.0) ** 2]
	print ys
	if p.y < 0.0: c = -1.0
	else: c = 1.0
	ys = [min(max(c * math.sqrt(y), p.y), p.y + 1.0) for y in ys if y >= 0.0]
	idx = 0
	while idx < len(ys):
		if ys[idx] not in ys[:idx]: idx += 1
		else: del ys[idx]
	idx = 0
	while idx < len(xs):
		if next((y for y in ys if mag ** 2 - xs[idx] ** 2 - y ** 2 >= 0.0), None) is None: del xs[idx]
		else: idx += 1
	idx = 0
	while idx < len(ys):
		if next((x for x in xs if mag ** 2 - ys[idx] ** 2 - x ** 2 >= 0.0), None) is None: del ys[idx]
		else: idx += 1
	print xs, ys
	return (xs, ys)


def xybounds(zone, mag):
	"""This helper function computes the minimum and maximum x/y values that have vectors with magnitude mag in the cube defined by zone. Returns (minx, maxx, miny, maxy)."""
	xs, ys = _compute_bounds(zone, mag)
	minx = min(xs)
	maxx = max(xs)
	miny = min(ys)
	maxy = max(ys)
	return (minx, maxx, miny, maxy)

class PermissionsManager(object):
	"""PermissionsManager objects control which conformations are allowed and restricted. They should not only be able to verify the permissions for a specific conformation, but also provide on demand a set of allowed conformations for a simulation.

		This class is abstract and must be subclassed for use in a simulation."""

	def allowed_conformations(self, *args, **kwargs):
		"""This function should return a list of allowed conformations given the information in args and kwargs. The format of input and output depends on the subclass and usage."""
		raise NotImplementedError("Subclasses of PermissionsManager must implement allowed_conformations().")

	def is_valid(self, *args, **kwargs):
		"""This function should return a Boolean value representing whether or not the given conformation is allowed. Parameters depend on the subclass and usage."""
		raise NotImplementedError("Subclasses of PermissionsManager must implement is_valid().")

class AAPermissionsManager(PermissionsManager):
	"""AAPermissionsManager is a concrete subclass of PermissionsManager that provides means for interpreting protein conformations and discriminating between allowed and restricted consecutive amino acid orientations."""
	
	def __init__(self, source, triplets_source):
		"""For source, pass in a directory of files (one for each combination of amino acid types) specifying the allowed alpha carbon zones for that type of consecutive interaction. triplets_source should be a text file specifying the permissible orientations given another pair of orientations."""
		super(AAPermissionsManager, self).__init__()
		self.allowed_zones = [[[] for i in xrange(AMINO_ACID_COUNT)] for j in xrange(AMINO_ACID_COUNT)]
		self.permissible_sequences = {}
		self.load_permissions_data(source)
		self.load_sequences_data(triplets_source)

	def load_permissions_data(self, source):
		if not os.path.exists(source):
			return
		paths = os.listdir(source)
		for path in paths:
			if ".txt" not in path: continue
			tag1, tag2 = path[:-4].split("-")
			tag1 = int(tag1)
			tag2 = int(tag2)
			del self.allowed_zones[tag1][tag2][:]
			with open(join(source, path), 'r') as file:
				for line in file:
					pt = Point3D(*line.strip().split(","))
					self.allowed_zones[tag1][tag2].append(pt)

	def load_sequences_data(self, source):
		assert ".txt" in source, "Permissible sequences data not in the form of a .txt file: {}".format(source)
		self.permissible_sequences = {}
		with open(source, "r") as file:
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
					if current_zones not in self.permissible_sequences:
						self.permissible_sequences[current_zones] = []
					self.permissible_sequences[current_zones].append(alpha_zone)

	def allowed_conformations(self, aminoacid, reference_aa, prior=True, opposite_aa=None):
		"""This overridden function accepts an amino acid whose allowed position zones are requested, and a reference_aa to which the amino acid's position should be favorable. Set prior to True if reference_aa comes before aminoacid, and False if not. Return values are position zones in the global coordinate system.
			If you pass the residue on the opposite side of reference_aa from aminoacid to opposite_aa, this method searches in the list of permissible sequences for the orientations defined between reference_aa and opposite_aa."""
		type1 = aacode(aminoacid.type)
		type2 = aacode(reference_aa.type)
		ret = []
		candidate_zones = []
		opposite_permissible = True
		if opposite_aa is not None:
			if prior == True:
				zone1 = opposite_aa.aa_position_zone(reference_aa).alpha_zone
				zone2 = reference_aa.aa_position_zone(opposite_aa).alpha_zone
			else:
				zone1 = reference_aa.aa_position_zone(opposite_aa).alpha_zone
				zone2 = opposite_aa.aa_position_zone(reference_aa).alpha_zone
			if (zone1, zone2) in self.permissible_sequences:
				candidate_zones = self.permissible_sequences[(zone1, zone2)]
			else:
				opposite_permissible = False
		if len(candidate_zones) == 0:
			if prior:
				candidate_zones = [zones[0] for zones in self.permissible_sequences]
			else:
				candidate_zones = [zones[1] for zones in self.permissible_sequences]
	
		for alpha_zone in candidate_zones:
			if (alpha_zone.y <= 0.0 and prior == True) or (alpha_zone.y > 0.0 and prior == False):
				continue
			retro_candidates = []
			if opposite_aa is not None and opposite_permissible:
				#Search for pairs within the top-level keys of permissible_sequences that the new amino acid could use
				for zonea, zoneb in self.permissible_sequences:
					if prior == True and zonea == alpha_zone:
						if zone2 in self.permissible_sequences[(zonea, zoneb)]:
							retro_candidates.append(zoneb)
					elif prior == False and zoneb == alpha_zone:
						if zone1 in self.permissible_sequences[(zonea, zoneb)]:
							retro_candidates.append(zonea)
			if len(retro_candidates) == 0:
				if prior:
					retro_candidates = [zones[1] for zones in self.permissible_sequences if zones[0] == alpha_zone]
				else:
					retro_candidates = [zones[0] for zones in self.permissible_sequences if zones[1] == alpha_zone]

			alpha_zone = alpha_zone.add(Point3D(0.5, 0.5, 0.5))
			alocation = reference_aa.toglobal(alpha_zone)
			alpha_zone = reference_aa.toglobal(alpha_zone).subtract(reference_aa.acarbon).multiply(-1.0)
			#Imagine that we are changing the location of the reference_aa, but in fact we are just rotating the aminoacid to see at which axes it will have allowed interactions with the reference amino acid.
			for retro_zone in retro_candidates:
				if (retro_zone.y <= 0.0 and prior == False) or (retro_zone.y > 0.0 and prior == True):
					continue
				axes = aminoacid.axes_for_zone(alpha_zone, retro_zone)
				if not axes:
					continue
				hypo = aminoacid.hypothetical(PositionZone(alocation, axes[0], axes[1], axes[2]), True)
				if (prior == True and not self.is_valid(hypo, reference_aa)) or (prior == False and not self.is_valid(hypo, reference_aa, prior)):
					continue
				ret.append(PositionZone(alocation, axes[0], axes[1], axes[2]))
		return ret

	def is_valid(self, aminoacid, reference_aa, prior=True):
		"""This overridden function accepts an amino acid whose location's permissions are requested, and a reference_aa to which the amino acid's position should be compared. Set prior to True if reference_aa comes before aminoacid, and False if not. Return values are position zones in the global coordinate system."""
		type1 = aacode(aminoacid.type)
		type2 = aacode(reference_aa.type)
		if (prior == True and reference_aa.carbon.distanceto(aminoacid.nitrogen) > 3.0) or (prior == False and reference_aa.nitrogen.distanceto(aminoacid.carbon) > 3.0):
			return False
		if prior:
			loc = reference_aa.tolocal(aminoacid.acarbon).floor()
			loc2 = aminoacid.tolocal(reference_aa.acarbon).floor()
		else:
			loc = aminoacid.tolocal(reference_aa.acarbon).floor()
			loc2 = reference_aa.tolocal(aminoacid.acarbon).floor()

		allowed = [zones[0] for zones in self.permissible_sequences]
		if loc in allowed: #self.allowed_zones[type2][type1]:
			allowed = [zones[1] for zones in self.permissible_sequences if zones[0] == loc]
			if loc2 in allowed:
				return True
			else:
				return False
		else:
			return False

class AASecondaryStructurePermissionsManager(PermissionsManager):
	"""AASecondaryStructurePermissionsManager is a concrete subclass of PermissionsManager that manages the allowed orientations for helices and sheets."""
	
	def __init__(self, source):
		"""For source, pass in a directory of files (one for each secondary structure type) specifying the allowed FULL position zones for that type of structure."""
		super(AASecondaryStructurePermissionsManager, self).__init__()
		self.permissible_sequences = {
			secondary_struct_helix + "1" : {}, secondary_struct_helix + "2" : {},
			secondary_struct_helix + "3" : {}, secondary_struct_helix + "4" : {},
			secondary_struct_helix + "5" : {}, secondary_struct_helix + "6" : {},
			secondary_struct_helix + "7" : {}, secondary_struct_helix + "8" : {},
			secondary_struct_helix + "9" : {}, secondary_struct_helix + "10" : {},
			secondary_struct_sheet + "0" : {}, secondary_struct_sheet + "1" : {},
			secondary_struct_sheet + "-1" : {}
		}
		self.load_permissions_data(source)

	def load_permissions_data(self, source):
		paths = os.listdir(source)
		for path in paths:
			if ".txt" not in path: continue
			structure_type = path[:-4]
			if structure_type not in self.permissible_sequences: continue
			del self.permissible_sequences[structure_type]
			self.permissible_sequences[structure_type] = {}
			with open(join(source, path), 'r') as file:
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
						if current_zones not in self.permissible_sequences[structure_type]:
							self.permissible_sequences[structure_type][current_zones] = []
						self.permissible_sequences[structure_type][current_zones].append(alpha_zone)

	def allowed_conformations(self, aminoacid, reference_aa, structure_type, subtype, prior=True, opposite_aa=None, size=10):
		"""This overridden function accepts an amino acid whose allowed position zones are requested, and a reference_aa (prior amino acid or not depending on the value of prior) to which the amino acid's position should be favorable. The structure_type and subtype can be obtained from the secondary structure and strand objects in the protein. Return values are position zones in the global coordinate system."""
		info_key = structure_type + str(subtype)
		if info_key not in self.permissible_sequences:
			print "Secondary structure type", info_key, "not supported."
			return []
		data = self.permissible_sequences[info_key]
		if len(data) == 0:
			print "No data available for secondary structure type", info_key
			return []
		ret = []

		candidate_zones = []
		opposite_permissible = True
		if opposite_aa is not None:
			if prior == True:
				zone1 = opposite_aa.aa_position_zone(reference_aa).alpha_zone
				zone2 = reference_aa.aa_position_zone(opposite_aa).alpha_zone
			else:
				zone1 = reference_aa.aa_position_zone(opposite_aa).alpha_zone
				zone2 = opposite_aa.aa_position_zone(reference_aa).alpha_zone
			if (zone1, zone2) in data:
				candidate_zones = data[(zone1, zone2)]
			else:
				opposite_permissible = False
		if len(candidate_zones) == 0:
			if prior:
				candidate_zones = [zones[0] for zones in data]
			else:
				candidate_zones = [zones[1] for zones in data]

		for alpha_zone in candidate_zones:
			if (alpha_zone.y <= 0.0 and prior == True) or (alpha_zone.y > 0.0 and prior == False):
				continue
			
			retro_candidates = []
			if opposite_aa is not None and opposite_permissible:
				#Search for pairs within the top-level keys of permissible_sequences that the new amino acid could use
				for zonea, zoneb in data:
					if prior == True and zonea == alpha_zone:
						if zone2 in data[(zonea, zoneb)]:
							retro_candidates.append(zoneb)
					elif prior == False and zoneb == alpha_zone:
						if zone1 in data[(zonea, zoneb)]:
							retro_candidates.append(zonea)
			if len(retro_candidates) == 0:
				if prior:
					retro_candidates = [zones[1] for zones in data if zones[0] == alpha_zone]
				else:
					retro_candidates = [zones[0] for zones in data if zones[1] == alpha_zone]

			
			alpha_zone = alpha_zone.add(Point3D(0.5, 0.5, 0.5))
			alocation = reference_aa.toglobal(alpha_zone)
			alpha_zone = reference_aa.toglobal(alpha_zone).subtract(reference_aa.acarbon).multiply(-1.0)
			#mag = alpha_zone.magnitude()
			#Imagine that we are changing the location of the reference_aa, but in fact we are just rotating the aminoacid to see at which axes it will have allowed interactions with the reference amino acid.
			for retro_zone in retro_candidates:
				if (retro_zone.y <= 0.0 and prior == False) or (retro_zone.y > 0.0 and prior == True):
					continue
				axes = aminoacid.axes_for_zone(alpha_zone, retro_zone)
				if not axes: continue
				hypo = aminoacid.hypothetical(PositionZone(alocation, axes[0], axes[1], axes[2]), True)
				if (prior == True and not self.is_valid(hypo, reference_aa, structure_type, subtype)) or (prior == False and not self.is_valid(hypo, reference_aa, prior, structure_type, subtype)):
					continue
				ret.append(PositionZone(alocation, axes[0], axes[1], axes[2]))
		return ret
			

	def is_valid(self, aminoacid, reference_aa, structure_type, subtype, prior=True):
		"""This overridden function accepts an amino acid whose location's permissions are requested, and a reference_aa to which the amino acid's position should be compared. Set prior to True if reference_aa comes before aminoacid, and False if not. Return values are position zones in the global coordinate system."""
		info_key = structure_type + str(subtype)
		if info_key not in self.permissible_sequences:
			print "Secondary structure type", info_key, "not supported."
			return []
		data = self.permissible_sequences[info_key]
		if len(data) == 0:
			print "No data available for secondary structure type", info_key
			return []

		if (prior == True and reference_aa.carbon.distanceto(aminoacid.nitrogen) > 3.0) or (prior == False and reference_aa.nitrogen.distanceto(aminoacid.carbon) > 3.0):
			return False
		if prior:
			loc = reference_aa.tolocal(aminoacid.acarbon).floor()
			loc2 = aminoacid.tolocal(reference_aa.acarbon).floor()
		else:
			loc = aminoacid.tolocal(reference_aa.acarbon).floor()
			loc2 = reference_aa.tolocal(aminoacid.acarbon).floor()

		allowed = [zones[0] for zones in data]
		if loc in allowed: #self.allowed_zones[type2][type1]:
			allowed = [zones[1] for zones in data if zones[0] == loc]
			if loc2 in allowed:
				return True
			else:
				return False
		else:
			return False
