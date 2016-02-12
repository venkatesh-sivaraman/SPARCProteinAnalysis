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
	'''mins = [math.sqrt(mag ** 2 - zone.z ** 2), -math.sqrt(mag ** 2 - zone.z ** 2)]
	maxs = [math.sqrt(mag ** 2 - (zone.z + 1.0) ** 2), -math.sqrt(mag ** 2 - (zone.z + 1.0) ** 2)]
	if mins[0] > zone.x and mins[0] < zone.x + 1.0:
		minx = mins[0]
	elif mins[1] > zone.x and mins[1] < zone.x + 1.0:
		minx = mins[1]
	else:
		assert False, "Neither min is within the zone (%r): %r" % (zone, mins)
	if mins[0] > zone.y and mins[0] < zone.y + 1.0:
		miny = mins[0]
	elif mins[1] > zone.y and mins[1] < zone.y + 1.0:
		miny = mins[1]
	else:
		assert False, "Neither min is within the zone (%r): %r" % (zone, mins)
	if maxs[0] > zone.x and maxs[0] < zone.x + 1.0:
		maxx = maxs[0]
	elif maxs[1] > zone.x and maxs[1] < zone.x + 1.0:
		maxx = maxs[1]
	else:
		assert False, "Neither max is within the zone (%r): %r" % (zone, maxs)
	if maxs[0] > zone.y and maxs[0] < zone.y + 1.0:
		maxy = maxs[0]
	elif maxs[1] > zone.y and maxs[1] < zone.y + 1.0:
		maxy = maxs[1]
	else:
		assert False, "Neither max is within the zone (%r): %r" % (zone, maxs)
	if minx > maxx:
		temp = minx
		minx = maxx
		maxx = temp
	if miny > maxy:
		temp = miny
		miny = maxy
		maxy = temp'''
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
	
	def __init__(self, source):
		"""For source, pass in a directory of files (one for each combination of amino acid types) specifying the allowed alpha carbon zones for that type of consecutive interaction."""
		super(AAPermissionsManager, self).__init__()
		self.allowed_zones = [[[] for i in xrange(AMINO_ACID_COUNT)] for j in xrange(AMINO_ACID_COUNT)]
		self.load_permissions_data(source)

	def load_permissions_data(self, source):
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

	def allowed_conformations(self, aminoacid, reference_aa, prior=True):
		"""This overridden function accepts an amino acid whose allowed position zones are requested, and a reference_aa to which the amino acid's position should be favorable. Set prior to True if reference_aa comes before aminoacid, and False if not. Return values are position zones in the global coordinate system."""
		type1 = aacode(aminoacid.type)
		type2 = aacode(reference_aa.type)
		ret = []
		for alpha_zone in self.allowed_zones[type2][type1]:
			if (alpha_zone.y <= 0.0 and prior == True) or (alpha_zone.y > 0.0 and prior == False):
				continue
			alpha_zone = alpha_zone.add(Point3D(0.5, 0.5, 0.5))
			alocation = reference_aa.toglobal(alpha_zone)
			alpha_zone = reference_aa.toglobal(alpha_zone).subtract(reference_aa.acarbon).multiply(-1.0)
			#mag = alpha_zone.magnitude()
			#Imagine that we are changing the location of the reference_aa, but in fact we are just rotating the aminoacid to see at which axes it will have allowed interactions with the reference amino acid.
			for retro_zone in self.allowed_zones[type1][type2]:
				if (retro_zone.y <= 0.0 and prior == False) or (retro_zone.y > 0.0 and prior == True):
					continue
				axes = aminoacid.axes_for_zone(alpha_zone, retro_zone)
				if not axes: continue
				hypo = aminoacid.hypothetical(PositionZone(alocation, axes[0], axes[1], axes[2]), True)
				if (prior == True and reference_aa.carbon.distanceto(hypo.nitrogen) > 3.0) or (prior == False and reference_aa.nitrogen.distanceto(hypo.carbon) > 3.0):
					continue
				if prior is True:
					assert self.is_valid(hypo, reference_aa), "Yielding invalid orientation for post: {} and {} ({})".format(hypo, reference_aa, reference_aa.tolocal(alocation))
				else:
					assert self.is_valid(hypo, reference_aa, prior), "Yielding invalid orientation for prior: {} and {} ({})".format(hypo, reference_aa, reference_aa.tolocal(alocation))
				ret.append(PositionZone(alocation, axes[0], axes[1], axes[2]))
		return ret

	def is_valid(self, aminoacid, reference_aa, prior=True):
		"""This overridden function accepts an amino acid whose location's permissions are requested, and a reference_aa to which the amino acid's position should be compared. Set prior to True if reference_aa comes before aminoacid, and False if not. Return values are position zones in the global coordinate system."""
		type1 = aacode(aminoacid.type)
		type2 = aacode(reference_aa.type)
		loc = reference_aa.tolocal(aminoacid.acarbon).floor()
		if (prior == True and reference_aa.carbon.distanceto(aminoacid.nitrogen) > 3.0) or (prior == False and reference_aa.nitrogen.distanceto(aminoacid.carbon) > 3.0):
			return False
		if loc in self.allowed_zones[type2][type1]:
			loc2 = aminoacid.tolocal(reference_aa.acarbon).floor()
			if loc2 in self.allowed_zones[type1][type2]:
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
		self.allowed_zones = {
			secondary_struct_helix + "1" : [], secondary_struct_helix + "2" : [],
			secondary_struct_helix + "3" : [], secondary_struct_helix + "4" : [],
			secondary_struct_helix + "5" : [], secondary_struct_helix + "6" : [],
			secondary_struct_helix + "7" : [], secondary_struct_helix + "8" : [],
			secondary_struct_helix + "9" : [], secondary_struct_helix + "10" : [],
			secondary_struct_sheet + "0" : [], secondary_struct_sheet + "1" : [],
			secondary_struct_sheet + "-1" : []
		}
		self.load_permissions_data(source)

	def load_permissions_data(self, source):
		paths = os.listdir(source)
		for path in paths:
			if ".txt" not in path: continue
			structure_type = path[:-4]
			del self.allowed_zones[structure_type]
			self.allowed_zones[structure_type] = []
			with open(join(source, path), 'r') as file:
				for line in file:
					pt = Point3D(*line.strip().split(","))
					self.allowed_zones[structure_type].append(pt)

	def allowed_conformations(self, aminoacid, reference_aa, structure_type, subtype, prior=True, size=10):
		"""This overridden function accepts an amino acid whose allowed position zones are requested, and a reference_aa (prior amino acid or not depending on the value of prior) to which the amino acid's position should be favorable. The structure_type and subtype can be obtained from the secondary structure and strand objects in the protein. Return values are position zones in the global coordinate system."""
		info_key = structure_type + str(subtype)
		if info_key not in self.allowed_zones:
			print "Secondary structure type", info_key, "not supported."
			return []
		data = self.allowed_zones[info_key]
		if len(data) == 0:
			print "No data available for secondary structure type", info_key
			return []

		ret = []
		for alpha_zone in data:
			if (alpha_zone.y <= 0.0 and prior == True) or (alpha_zone.y > 0.0 and prior == False):
				continue
			alpha_zone = alpha_zone.add(Point3D(0.5, 0.5, 0.5))
			alocation = reference_aa.toglobal(alpha_zone)
			alpha_zone = reference_aa.toglobal(alpha_zone).subtract(reference_aa.acarbon).multiply(-1.0)
			#mag = alpha_zone.magnitude()
			#Imagine that we are changing the location of the reference_aa, but in fact we are just rotating the aminoacid to see at which axes it will have allowed interactions with the reference amino acid.
			for retro_zone in self.allowed_zones[info_key]:
				if (retro_zone.y <= 0.0 and prior == False) or (retro_zone.y > 0.0 and prior == True):
					continue
				axes = aminoacid.axes_for_zone(alpha_zone, retro_zone)
				if not axes: continue
				hypo = aminoacid.hypothetical(PositionZone(alocation, axes[0], axes[1], axes[2]), True)
				if (prior == True and reference_aa.carbon.distanceto(hypo.nitrogen) > 3.0) or (prior == False and reference_aa.nitrogen.distanceto(hypo.carbon) > 3.0):
					continue
				if prior is True:
					assert self.is_valid(hypo, reference_aa, structure_type, subtype), "Yielding invalid orientation for post: {} and {} ({})".format(hypo, reference_aa, reference_aa.tolocal(alocation))
				else:
					assert self.is_valid(hypo, reference_aa, structure_type, subtype, prior), "Yielding invalid orientation for prior: {} and {} ({})".format(hypo, reference_aa, reference_aa.tolocal(alocation))
				ret.append(PositionZone(alocation, axes[0], axes[1], axes[2]))
		return ret
			

	def is_valid(self, aminoacid, reference_aa, structure_type, subtype, prior=True):
		"""This overridden function accepts an amino acid whose location's permissions are requested, and a reference_aa to which the amino acid's position should be compared. Set prior to True if reference_aa comes before aminoacid, and False if not. Return values are position zones in the global coordinate system."""
		info_key = structure_type + str(subtype)
		if info_key not in self.allowed_zones:
			print "Secondary structure type", info_key, "not supported."
			return []
		data = self.allowed_zones[info_key]
		if len(data) == 0:
			print "No data available for secondary structure type", info_key
			return []

		loc = reference_aa.tolocal(aminoacid.acarbon).floor()
		if (prior == True and reference_aa.carbon.distanceto(aminoacid.nitrogen) > 3.0) or (prior == False and reference_aa.nitrogen.distanceto(aminoacid.carbon) > 3.0):
			return False
		if loc in self.allowed_zones[info_key]:
			loc2 = aminoacid.tolocal(reference_aa.acarbon).floor()
			if loc2 in self.allowed_zones[info_key]:
				return True
			else:
				return False
		else:
			return False
