from proteinmath import *

AMINO_ACID_COUNT = 22

amino_acid_alanine = "ALA"
amino_acid_arginine = "ARG"
amino_acid_asparagine = "ASN"
amino_acid_aspartic_acid = "ASP"
amino_acid_asparagine_or_aspartic_acid = "ASX"
amino_acid_cysteine = "CYS"
amino_acid_glutamine = "GLN"
amino_acid_glutamic_acid = "GLU"
amino_acid_glutamine_or_glutamic_acid = "GLX"
amino_acid_glycine = "GLY"
amino_acid_histidine = "HIS"
amino_acid_isoleucine = "ILE"
amino_acid_leucine = "LEU"
amino_acid_lysine = "LYS"
amino_acid_methionine = "MET"
amino_acid_phenylalanine = "PHE"
amino_acid_proline = "PRO"
amino_acid_serine = "SER"
amino_acid_threonine = "THR"
amino_acid_tryptophan = "TRP"
amino_acid_tyrosine = "TYR"
amino_acid_valine = "VAL"

def aacode(name):
	if (name is amino_acid_alanine):
		return 0
	elif (name is amino_acid_arginine):
		return 1
	elif (name is amino_acid_asparagine):
		return 2
	elif (name is amino_acid_asparagine_or_aspartic_acid):
		return 3
	elif (name is amino_acid_aspartic_acid):
		return 4
	elif (name is amino_acid_cysteine):
		return 5
	elif (name is amino_acid_glutamic_acid):
		return 6
	elif (name is amino_acid_glutamine):
		return 7
	elif (name is amino_acid_glutamine_or_glutamic_acid):
		return 8
	elif (name is amino_acid_glycine):
		return 9
	elif (name is amino_acid_histidine):
		return 10
	elif (name is amino_acid_isoleucine):
		return 11
	elif (name is amino_acid_leucine):
		return 12
	elif (name is amino_acid_lysine):
		return 13
	elif (name is amino_acid_methionine):
		return 14
	elif (name is amino_acid_phenylalanine):
		return 15
	elif (name is amino_acid_proline):
		return 16
	elif (name is amino_acid_serine):
		return 17
	elif (name is amino_acid_threonine):
		return 18
	elif (name is amino_acid_tryptophan):
		return 19
	elif (name is amino_acid_tyrosine):
		return 20
	elif (name is amino_acid_valine):
		return 21
	else:
		return 22

def aatypec(name):
	if (name is "A"):
		return amino_acid_alanine
	elif (name is "R"):
		return amino_acid_arginine
	elif (name is "N"):
		return amino_acid_asparagine
	elif (name is "B"):
		return amino_acid_asparagine_or_aspartic_acid
	elif (name is "D"):
		return amino_acid_aspartic_acid
	elif (name is "C"):
		return amino_acid_cysteine
	elif (name is "E"):
		return amino_acid_glutamic_acid
	elif (name is "Q"):
		return amino_acid_glutamine
	elif (name is "Z"):
		return amino_acid_glutamine_or_glutamic_acid
	elif (name is "G"):
		return amino_acid_glycine
	elif (name is "H"):
		return amino_acid_histidine
	elif (name is "I"):
		return amino_acid_isoleucine
	elif (name is "L"):
		return amino_acid_leucine
	elif (name is "K"):
		return amino_acid_lysine
	elif (name is "M"):
		return amino_acid_methionine
	elif (name is "F"):
		return amino_acid_phenylalanine
	elif (name is "P"):
		return amino_acid_proline
	elif (name is "S"):
		return amino_acid_serine
	elif (name is "T"):
		return amino_acid_threonine
	elif (name is "W"):
		return amino_acid_tryptophan
	elif (name is "Y"):
		return amino_acid_tyrosine
	elif (name is "V"):
		return amino_acid_valine
	else:
		return ""

class AminoAcid(object):
	'Represents a single amino acid.'

	def __init__(self, type):
		self.type = type
		self.tag = 0
		self.acarbon = Point3D.zero()
		self.bcarbon = Point3D.zero()
	
	def __init__(self, type, tag, acarbon, bcarbon):
		self.type = type
		self.tag = tag
		self.acarbon = acarbon
		self.bcarbon = bcarbon
	
	def __str__(self):
		return "Amino acid %s: %s" % (self.type, self.acarbon.__str__())

	@property
	def type(self):
		return self._type

	@type.setter
	def type(self, value):
		self._type = value
		radii = (2.0745,
				 5.6371,
				 2.6128,
				 0.0000,
				 2.3381,
				 2.1038,
				 3.6782,
				 3.7643,
				 0.0000,
				 1.2853,
				 4.2895,
				 2.9068,
				 3.1256,
				 4.7853,
				 3.6317,
				 4.2265,
				 2.2682,
				 2.0983,
				 2.0542,
				 5.3609,
				 5.8212,
				 2.0446,
				 999999)
		self.radius = radii[aacode(value)]

	def aa_position_zone(self, aa):
		alpha_zone = 0
		beta_zone = 0
		alpha_loc = aa.acarbon
		beta_loc = aa.bcarbon
		print self.acarbon, self.bcarbon, alpha_loc, beta_loc
	
		#This method for finding the locations of the other amino acid's carbons in self's coordinate space uses vectors to define a new coordinate system. The in_coordinate_system() function expresses each carbon's location in the coordinate system defined by i, j, and k. These vectors should be orthogonal unit vectors, with i corresponding to the alpha-beta carbon bond in self. Penultimate has diagrams on how this works.
	
		i = self.bcarbon.subtract(self.acarbon).normalize()
		spherical = i.tospherical()
		j = Point3D(1.0, math.pi / 2.0, math.pi / 2.0).add(Point3D(1.0, 0.0, math.pi / 2.0).subtract(spherical))
		j = j.tocartesian()

		#Cross product of i and j gives perpendicular vector axis
		k = Point3D(i.y * j.z - i.z * j.y, i.z * j.x - i.x * j.z, i.x * j.y - i.y * j.x).normalize()
		
		myBetaCarbon = self.bcarbon.subtract(self.acarbon).in_coordinate_system(i, j, k)
		alpha_loc = alpha_loc.subtract(self.acarbon).in_coordinate_system(i, j, k)
		beta_loc = beta_loc.subtract(self.acarbon).in_coordinate_system(i, j, k)
	
		print alpha_loc, beta_loc
		alpha_zone = zone_point(alpha_loc)
		beta_zone = zone_point(beta_loc)
	
		if alpha_zone < 0 or beta_zone < 0:
			print "Abnormal zone calculation"
		else:
			return PositionZone(alpha_zone, beta_zone)
		return None

	def position_zone_from_other(self, aa, zp):
		'This method takes a zone pair that is defined with respect to the other amino acid, aa, and transforms it into a zone pair that is defined with respect to self. Of course, the zone pair must be within 10 angstroms of self in order to be valid.'

		alpha_loc_aa = Point3D((zp.alpha_zone % (POSITION_ZONE_EXTENSION * 2)) - POSITION_ZONE_EXTENSION + 0.5,
							   math.floor((zp.alpha_zone % ((POSITION_ZONE_EXTENSION * 2) * (POSITION_ZONE_EXTENSION * 2))) / (POSITION_ZONE_EXTENSION * 2.0)) - POSITION_ZONE_EXTENSION + 0.5,
							   math.floor(zp.alpha_zone / ((POSITION_ZONE_EXTENSION * 2.0) * (POSITION_ZONE_EXTENSION * 2.0))) - POSITION_ZONE_EXTENSION + 0.5)
		beta_loc_aa = Point3D((zp.beta_zone % (POSITION_ZONE_EXTENSION * 2)) - POSITION_ZONE_EXTENSION + 0.5,
							  math.floor((zp.beta_zone % ((POSITION_ZONE_EXTENSION * 2) * (POSITION_ZONE_EXTENSION * 2))) / (POSITION_ZONE_EXTENSION * 2.0)) - POSITION_ZONE_EXTENSION + 0.5,
							  math.floor(zp.beta_zone / ((POSITION_ZONE_EXTENSION * 2.0) * (POSITION_ZONE_EXTENSION * 2.0))) - POSITION_ZONE_EXTENSION + 0.5)

		#First, find out the location of the unit vector parallel to the positive x-axis with aa's alpha carbon at the center. Remember that in the zp coordinate system, aa's beta carbon is at (x, 0, 0).
		bc_location = aa.bcarbon.subtract(aa.acarbon).normalize().tospherical()
		print "bc: ", bc_location
		#get the location of the z-axis in the zp coordinate system to use bc_location.phi properly.
		zaxis = Point3D(math.cos(bc_location.z), 0.0, math.sin(bc_location.z))
		#the angle from the z-axis to the beta carbon in the zp coordinate system, which is phi, should equal the angle from the 
		i = Point3D(1.0, -bc_location.y, bc_location.z).tocartesian()
		print "i:", i
		'''xaxis = Point3D(1, 0, 0)
		angle1 = math.acos(dotproduct(i, xaxis) / i.magnitude())
		angle2 = math.acos(dotproduct(xaxis, aa.bcarbon.subtract(aa.acarbon)) / aa.bcarbon.subtract(aa.acarbon).magnitude())
		assert math.fabs(angle1 - angle2) <= 0.001, "Incorrect math in position_zone_from_other: %r, %r" % (angle1, angle2)'''
		
		'''i = aa.bcarbon.subtract(aa.acarbon).normalize()
		j = Point3D(1.0, math.pi / 2.0, math.pi / 2.0).add(Point3D(1.0, 0.0, math.pi / 2.0).subtract(i.tospherical())).tocartesian()
		print "i: {i} j: {j}".format(i=i,j=j)'''



bucket_count = 8000
location_min = -500.0
location_max = 500.0
box_dimension = 5.0

class AAHashTable(object):
	'Represents a hash table of amino acids that can be retrieved based on proximity.'
	
	def hash_function(self, point):
		cubeWidth = (location_max - location_min) / (10.0)
		cubePoint = Point3D(math.fmod(point.x - location_min, cubeWidth), math.fmod(point.y - location_min, cubeWidth), math.fmod(point.z - location_min, cubeWidth))
		cubeWidth /= box_dimension
		hash = int(math.floor(cubePoint.z / box_dimension) * cubeWidth * cubeWidth + math.floor(cubePoint.y / box_dimension) * cubeWidth + math.floor(cubePoint.x / box_dimension))
		if hash > 8000:
			print point, ", ", cubePoint, ", ", hash
		return hash
	
	def __init__(self):
		self.buckets = [[] for x in range(bucket_count)]

	def add(self, aa):
		bucket = self.buckets[self.hash_function(aa.acarbon)]
		if (aa not in bucket):
			bucket.append(aa)
		bucket = self.buckets[self.hash_function(aa.bcarbon)]
		if (aa not in bucket):
			bucket.append(aa)

	def remove(self, aa):
		bucket = self.buckets[self.hash_function(aa.acarbon)]
		if (aa in bucket): bucket.remove(aa)
		bucket = self.buckets[self.hash_function(aa.bcarbon)]
		if (aa in bucket): bucket.remove(aa)

	def contains(self, aa):
		bucket = self.buckets[self.hash_function(aa.acarbon)]
		bucket2 = self.buckets[self.hash_function(aa.bcarbon)]
		return aa in bucket or aa in bucket2

	def nearby_aa(self, aa, distance):
		bucket = self.buckets[self.hash_function(aa.acarbon)]
		n = int(math.ceil(distance / box_dimension))
		retVal = []
		for z in range(-n, n + 1):
			for y in range(-n, n + 1):
				for x in range(-n, n + 1):
					searchBucket = self.hash_function(Point3D(aa.acarbon.x + x * box_dimension,
															  aa.acarbon.y + y * box_dimension,
															  aa.acarbon.z + z * box_dimension))
					for candidate in self.buckets[searchBucket]:
						if candidate.tag == aa.tag:	continue;
						firstDistance = aa.acarbon.distanceto(candidate.acarbon)
						if firstDistance <= distance and firstDistance > 0.001:
							retVal.append(candidate)
						else:
							firstDistance = aa.acarbon.distanceto(candidate.bcarbon)
							if firstDistance <= distance and firstDistance > 0.001:
								retVal.append(candidate)
		return retVal

