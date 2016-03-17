from proteinmath import *
import resource

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

CARBON_BOND_LENGTH = 1.53
NITROGEN_BOND_LENGTH = 1.47
BOLTZMANN_CONSTANT = 1.3806488E-23
TEMPERATURE = 298

running_lcs_calculations = 0
running_gcs_calculations = 0

aa_masses = {
	"ILE" : 131.1736,
	"LEU" : 131.1736,
	"LYS" : 146.1882,
	"MET" : 149.2124,
	"PHE" : 165.1900,
	"THR" : 119.1197,
	"TRP" :	204.2262,
	"VAL" : 117.1469,
	"ARG" : 174.2017,
	"HIS" : 155.1552,
	"ALA" : 89.0935,
	"ASN" : 132.1184,
	"ASP" : 133.1032,
	"CYS" : 121.1590,
	"GLU" : 147.1299,
	"GLN" : 146.1451,
	"GLY" : 75.0669,
	"PRO" : 115.1310,
	"SER" : 105.0930,
	"TYR" : 181.1894
}

def reset_stats():
	global running_lcs_calculations, running_gcs_calculations
	#print "Local:", running_lcs_calculations, "Global:", running_gcs_calculations
	running_gcs_calculations = 0
	running_lcs_calculations = 0

def aacode(name):
	name = name.upper()
	if (name == amino_acid_alanine):
		return 0
	elif (name == amino_acid_arginine):
		return 1
	elif (name == amino_acid_asparagine):
		return 2
	elif (name == amino_acid_asparagine_or_aspartic_acid):
		return 3
	elif (name == amino_acid_aspartic_acid):
		return 4
	elif (name == amino_acid_cysteine):
		return 5
	elif (name == amino_acid_glutamic_acid):
		return 6
	elif (name == amino_acid_glutamine):
		return 7
	elif (name == amino_acid_glutamine_or_glutamic_acid):
		return 8
	elif (name == amino_acid_glycine):
		return 9
	elif (name == amino_acid_histidine):
		return 10
	elif (name == amino_acid_isoleucine):
		return 11
	elif (name == amino_acid_leucine):
		return 12
	elif (name == amino_acid_lysine):
		return 13
	elif (name == amino_acid_methionine):
		return 14
	elif (name == amino_acid_phenylalanine):
		return 15
	elif (name == amino_acid_proline):
		return 16
	elif (name == amino_acid_serine):
		return 17
	elif (name == amino_acid_threonine):
		return 18
	elif (name == amino_acid_tryptophan):
		return 19
	elif (name == amino_acid_tyrosine):
		return 20
	elif (name == amino_acid_valine):
		return 21
	else:
		return 22

def aatype(code):
	if code == 0:
		return amino_acid_alanine
	elif code == 1:
		return amino_acid_arginine
	elif code == 2:
		return amino_acid_asparagine
	elif code == 3:
		return amino_acid_asparagine_or_aspartic_acid
	elif code == 4:
		return amino_acid_aspartic_acid
	elif code == 5:
		return amino_acid_cysteine
	elif code == 6:
		return amino_acid_glutamic_acid
	elif code == 7:
		return amino_acid_glutamine
	elif code == 8:
		return amino_acid_glutamine_or_glutamic_acid
	elif code == 9:
		return amino_acid_glycine
	elif code == 10:
		return amino_acid_histidine
	elif code == 11:
		return amino_acid_isoleucine
	elif code == 12:
		return amino_acid_leucine
	elif code == 13:
		return amino_acid_lysine
	elif code == 14:
		return amino_acid_methionine
	elif code == 15:
		return amino_acid_phenylalanine
	elif code == 16:
		return amino_acid_proline
	elif code == 17:
		return amino_acid_serine
	elif code == 18:
		return amino_acid_threonine
	elif code == 19:
		return amino_acid_tryptophan
	elif code == 20:
		return amino_acid_tyrosine
	elif code == 21:
		return amino_acid_valine
	else:
		return "UNK"

def aatypec(name):
	if (name == "A"):
		return amino_acid_alanine
	elif (name == "R"):
		return amino_acid_arginine
	elif (name == "N"):
		return amino_acid_asparagine
	elif (name == "B"):
		return amino_acid_asparagine_or_aspartic_acid
	elif (name == "D"):
		return amino_acid_aspartic_acid
	elif (name == "C"):
		return amino_acid_cysteine
	elif (name == "E"):
		return amino_acid_glutamic_acid
	elif (name == "Q"):
		return amino_acid_glutamine
	elif (name == "Z"):
		return amino_acid_glutamine_or_glutamic_acid
	elif (name == "G"):
		return amino_acid_glycine
	elif (name == "H"):
		return amino_acid_histidine
	elif (name == "I"):
		return amino_acid_isoleucine
	elif (name == "L"):
		return amino_acid_leucine
	elif (name == "K"):
		return amino_acid_lysine
	elif (name == "M"):
		return amino_acid_methionine
	elif (name == "F"):
		return amino_acid_phenylalanine
	elif (name == "P"):
		return amino_acid_proline
	elif (name == "S"):
		return amino_acid_serine
	elif (name == "T"):
		return amino_acid_threonine
	elif (name == "W"):
		return amino_acid_tryptophan
	elif (name == "Y"):
		return amino_acid_tyrosine
	elif (name == "V"):
		return amino_acid_valine
	else:
		return ""


class AminoAcid(object):
	'''Represents a single amino acid. Important: amino acid tags between protein objects should not be close together.'''

	def __init__(self, type, tag=0, acarbon=Point3D.zero(), nitrogen=Point3D.zero(), carbon=Point3D.zero(), sidechain=None):
		self.mass = 0.0
		self.type = type
		self.tag = tag
		self.acarbon = acarbon
		self.sidechain = sidechain
		self.nitrogen = nitrogen
		self.carbon = carbon
		self.i = Point3D.zero()
		self.j = Point3D.zero()
		self.k = Point3D.zero()
		self.localscore = 0
		self.otheratoms = {}
		self.cluster = (self.tag, self.tag + 1)
		self.clusterscore = 0
		self.has_break = False
		self._observers = {}
		self._tmp_atoms = None
	
	def withtag(self, newtag):
		"""Creates a new amino acid with the specified tag number. Note: The new amino acid does NOT automatically gain the observers that were assigned to this object."""
		aa = AminoAcid(self.type, newtag, self.acarbon, self.nitrogen, self.carbon, self.sidechain)
		aa.i = self.i
		aa.j = self.j
		aa.k = self.k
		aa.otheratoms = self.otheratoms
		aa.localscore = self.localscore
		aa.cluster = self.cluster
		aa.clusterscore = self.clusterscore
		aa.has_break = self.has_break
		return aa
	
	def set_axes(self, i, j, k, normalized=False, calculate_atoms=True):
		old = (self.i, self.j, self.k)
		if normalized:
			self.i = i
			self.j = j
			self.k = k
		else:
			self.i = i.normalize()
			self.j = j.normalize()
			self.k = k.normalize()
		self.localscore = 0.0
		if calculate_atoms is True:
			self.nitrogen = self.toglobal(Point3D(NITROGEN_BOND_LENGTH, math.pi + math.acos(-1.0 / 2.0) / 2.0, math.acos(-1.0 / 3.0)).tocartesian())
			self.carbon = self.toglobal(Point3D(CARBON_BOND_LENGTH, math.pi - math.acos(-1.0 / 2.0) / 2.0, math.acos(-1.0 / 3.0)).tocartesian())
		if "axes" in self._observers:
			for callback in self._observers["axes"]:
				callback(self, "axes", old)
	
	def __str__(self):
		ret = "%s (%d): %s" % (self.type, self.tag, self.acarbon.__str__())
		if self.has_break: ret += " <end-of-chain>"
		return ret

	def __repr__(self):
		ret = "%s (%d): %s" % (self.type, self.tag, self.acarbon.__str__())
		if self.has_break: ret += " <end-of-chain>"
		return ret
	
	@classmethod
	def nlocation(cls, pz):
		globali = Point3D(1.0, 0.0, 0.0).in_coordinate_system(pz.x_axis, pz.y_axis, pz.z_axis)
		globalj = Point3D(0.0, 1.0, 0.0).in_coordinate_system(pz.x_axis, pz.y_axis, pz.z_axis)
		globalk = Point3D(0.0, 0.0, 1.0).in_coordinate_system(pz.x_axis, pz.y_axis, pz.z_axis)
		return Point3D(NITROGEN_BOND_LENGTH, math.pi + math.acos(-1.0 / 2.0) / 2.0, math.acos(-1.0 / 3.0)).tocartesian().in_coordinate_system(globali, globalj, globalk).add(pz.alpha_zone)

	@classmethod
	def clocation(cls, pz):
		globali = Point3D(1.0, 0.0, 0.0).in_coordinate_system(pz.x_axis, pz.y_axis, pz.z_axis)
		globalj = Point3D(0.0, 1.0, 0.0).in_coordinate_system(pz.x_axis, pz.y_axis, pz.z_axis)
		globalk = Point3D(0.0, 0.0, 1.0).in_coordinate_system(pz.x_axis, pz.y_axis, pz.z_axis)
		return Point3D(CARBON_BOND_LENGTH, math.pi - math.acos(-1.0 / 2.0) / 2.0, math.acos(-1.0 / 3.0)).tocartesian().in_coordinate_system(globali, globalj, globalk).add(pz.alpha_zone)
	
	@property
	def acarbon(self):
		return self._acarbon
	
	@acarbon.setter
	def acarbon(self, value):
		has = False
		if hasattr(self, "_acarbon"):
			old = self._acarbon
			has = True
		self.localscore = 0.0
		self._acarbon = value
		if has is True:
			self.nitrogen = self.nitrogen.subtract(old).add(self._acarbon)
			self.carbon = self.carbon.subtract(old).add(self._acarbon)
		if hasattr(self, "_tmp_atoms") and self._acarbon != Point3D.zero():
			if self._tmp_atoms:
				for key, value in self._tmp_atoms.iteritems():
					self.otheratoms[key] = self.tolocal(value)
				self._tmp_atoms = None
		if has is True and "acarbon" in self._observers:
			for callback in self._observers["acarbon"]:
				callback(self, "acarbon", old)

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
		if self.type in aa_masses:
			self.mass = aa_masses[self.type]
		else:
			self.mass = 0.0

	def add_observer(self, key, callback):
		"""Pass in a string key ("acarbon" and "axes" at present) and a function object to be notified when the value changes.
			Callback should accept three arguments: an AminoAcid object, key, and old value. (If key="axes", old value is a tuple of axes.)"""
		if key in self._observers:
			self._observers[key].append(callback)
		else:
			self._observers[key] = [callback]

	def remove_observer(self, key, callback):
		if key in self._observers:
			self._observers[key].remove(callback)
		else:
			print "Could not remove observer because {} is not being observed by anyone.".format(key)

	def hypothetical(self, pz, calculate_atoms=True):
		aa = AminoAcid(self.type, self.tag, pz.alpha_zone, self.nitrogen, self.carbon, self.sidechain)
		aa.set_axes(pz.x_axis, pz.y_axis, pz.z_axis, calculate_atoms=calculate_atoms)
		return aa

	def _calc_x(self, j, k, c):
		"""Use the cross products of j and k to determine two possible vectors perpendicular to their plane, and use the one that maximizes the angle with c."""
		crossproduct1 = Point3D(j.y * k.z - j.z * k.y, j.z * k.x - j.x * k.z, j.x * k.y - j.y * k.x)
		crossproduct2 = Point3D(k.y * j.z - k.z * j.y, k.z * j.x - k.x * j.z, k.x * j.y - k.y * j.x)
		return (crossproduct1 if (crossproduct1.anglewith(c) > crossproduct2.anglewith(c)) else crossproduct2).normalize()
	
	def compute_coordinate_system_vectors(self):
		"""Compute the vectors i, j, and k that produce a right-handed coordinate system such that the tetrahedral structure of the amino acid is represented in a constant manner: the N and C atoms are located in the -x, -z region, while the sidechain points upward and the hydrogen atom points forward along the xz-plane."""

		n = self.nitrogen.subtract(self.acarbon).normalize()
		c = self.carbon.subtract(self.acarbon).normalize()
		axes = self.coordinate_axes(n, c)
		if len(axes) > 0:
			self.set_axes(*axes, calculate_atoms=False)
	
	def coordinate_axes(self, n, c):
		#First, find the y-axis by finding the candidate that minimizes the angle with the c-vector.
		candidate1 = Point3D(c.x - n.x, c.y - n.y, c.z - n.z)
		candidate2 = Point3D(n.x - c.x, n.y - c.y, n.z - c.z)
		vectors = [ Point3D.zero(), Point3D.zero(), Point3D.zero() ]
		vectors[1] = (candidate1 if (candidate1.anglewith(c) < candidate2.anglewith(c)) else candidate2).normalize()
		
		#Next, find the z-axis using the quadratic equation described in Notability, "Science Fair Calculations."
		alpha = 3 * (c.y * n.z - n.y * c.z)
		beta = c.z - n.z
		phi = 3 * (c.x * n.z - n.x * c.z)
		lam = n.y - c.y
		mu = 3 * (c.y * n.x - c.x * n.y)
		if alpha == 0: vectors[1] = Point3D.zero(); return []
		ak = alpha ** 2 + phi ** 2 + mu ** 2
		bk = -2 * (phi * beta + lam * mu)
		ck = beta ** 2 + lam ** 2 - alpha ** 2
		discriminant = max(0.0, bk ** 2 - 4 * ak * ck)
		can1 = (-bk + math.sqrt(discriminant)) / (2 * ak)
		can2 = (-bk - math.sqrt(discriminant)) / (2 * ak)
		candidate1 = Point3D(can1, (beta - phi * can1) / alpha, (lam - mu * can1) / alpha)
		candidate2 = Point3D(can2, (beta - phi * can2) / alpha, (lam - mu * can2) / alpha)
		i1 = self._calc_x(vectors[1], candidate1, c)
		i2 = self._calc_x(vectors[1], candidate2, c)
		if math.fabs(vectors[1].in_coordinate_system(i1, vectors[1], candidate1).y - 1.0) <= 0.001:
			vectors[0] = i1
			vectors[2] = candidate1
		elif math.fabs(vectors[1].in_coordinate_system(i2, vectors[1], candidate2).y - 1.0) <= 0.001:
			vectors[0] = i2
			vectors[2] = candidate2
		else:
			return []
			#assert False, "no valid y-axis found: %.6f, %.6f" % (math.fabs(vectors[1].in_coordinate_system(i1, vectors[1], candidate1).y - 1.0), math.fabs(vectors[1].in_coordinate_system(i2, vectors[1], candidate2).y - 1.0))
		return vectors

	def aa_position_zone(self, aa):
		"""This method for finding the locations of the other amino acid's carbons in self's coordinate space uses vectors to define a new coordinate system. The in_coordinate_system() function expresses each carbon's location in the coordinate system defined by i, j, and k. Penultimate has diagrams on how this works."""
		if self.i == Point3D.zero() or self.j == Point3D.zero() or self.k == Point3D.zero():
			self.compute_coordinate_system_vectors()
		assert self.i != Point3D.zero() and self.j != Point3D.zero() and self.k != Point3D.zero(), "amino acid with no LCS: %r" % self
		alpha_loc = aa.acarbon.subtract(self.acarbon).in_coordinate_system(self.i, self.j, self.k).floor()
		axisoffset = aa.acarbon.subtract(self.acarbon)
		i_loc = aa.i.in_coordinate_system(self.i, self.j, self.k).floor(0.1)
		j_loc = aa.j.in_coordinate_system(self.i, self.j, self.k).floor(0.1)
		k_loc = aa.k.in_coordinate_system(self.i, self.j, self.k).floor(0.1)
		return PositionZone(alpha_loc, i_loc, j_loc, k_loc)

	def tolocal(self, point):
		"""Converts a point from the global coordinate system (in which self.acarbon is defined) to the local coordinate system defined by the i, j, and k vectors."""
		if self.i == Point3D.zero() or self.j == Point3D.zero() or self.k == Point3D.zero():
			self.compute_coordinate_system_vectors()
		pt = point.subtract(self.acarbon)
		global running_lcs_calculations
		running_lcs_calculations += 1
		return pt.in_coordinate_system(self.i, self.j, self.k)

	def localpz(self, pz):
		intermediate = PositionZone(pz.alpha_zone, pz.x_axis.add(pz.alpha_zone), pz.y_axis.add(pz.alpha_zone), pz.z_axis.add(pz.alpha_zone))
		local_a = self.tolocal(intermediate.alpha_zone)
		return PositionZone(local_a,
							self.tolocal(intermediate.x_axis).subtract(local_a),
							self.tolocal(intermediate.y_axis).subtract(local_a),
							self.tolocal(intermediate.z_axis).subtract(local_a))

	def toglobal(self, point):
		"""Converts a point from the local coordinate system to the global coordinate system."""
		if self.i == Point3D.zero() or self.j == Point3D.zero() or self.k == Point3D.zero():
			self.compute_coordinate_system_vectors()
		globali = Point3D(1.0, 0.0, 0.0).in_coordinate_system(self.i, self.j, self.k)
		globalj = Point3D(0.0, 1.0, 0.0).in_coordinate_system(self.i, self.j, self.k)
		globalk = Point3D(0.0, 0.0, 1.0).in_coordinate_system(self.i, self.j, self.k)
		global running_gcs_calculations
		running_gcs_calculations += 1
		return point.in_coordinate_system(globali, globalj, globalk).add(self.acarbon)

	def globalpz(self, pz):
		intermediate = PositionZone(pz.alpha_zone, pz.x_axis.add(pz.alpha_zone), pz.y_axis.add(pz.alpha_zone), pz.z_axis.add(pz.alpha_zone))
		global_a = self.toglobal(intermediate.alpha_zone)
		return PositionZone(global_a,
							self.toglobal(intermediate.x_axis).subtract(global_a),
							self.toglobal(intermediate.y_axis).subtract(global_a),
							self.toglobal(intermediate.z_axis).subtract(global_a))

	def add_other_atom(self, atomname, location, convert=True):
		"""Adds the atom of type `atomname` to the amino acid. `location` is a global position which is converted into a local coordinate in the otheratoms dictionary by default; if it is already local, then pass in False for the `convert` kwarg. These atoms are NOT updated when the alpha carbon location changes."""
		if self.acarbon == Point3D.zero():
			if self._tmp_atoms: self._tmp_atoms[atomname] = location
			else: self._tmp_atoms = { atomname : location }
			return
		else:
			if self._tmp_atoms:
				for key, value in self._tmp_atoms.iteritems():
					self.otheratoms[key] = self.tolocal(value)
				self._tmp_atoms = None
		if convert:
			self.otheratoms[atomname] = self.tolocal(location)
		else:
			self.otheratoms[atomname] = location

	def axes_for_zone(self, alpha_zone, retro_zone, timeout=False):
		"""This method computes a set of coordinate axes such that the global vector from this amino acid to some aa2 is alpha_zone, and the LCS position zone of aa2 is retro_zone. It essentially reverse engineers the process of choosing coordinate axes to create a certain relative orientation. Returns a tuple (i, j, k). If timeout is True, then the algorithm will choose a default value instead of a random vector within the position zone."""
		#Only select retro zones inside which there is a vector with magnitude equal to the magnitude of alpha_zone.
		mag = alpha_zone.magnitude()
		minim = retro_zone
		maxim = Point3D(retro_zone.x + 1, retro_zone.y + 1, retro_zone.z + 1)
		if (minim.magnitude() > mag and maxim.magnitude() > mag) or (minim.magnitude() < mag and maxim.magnitude() < mag):
			return None
		#Find the coordinate axes that yield the vector of equal magnitude given retro zone.
		a = Point3D.zero()
		mid = retro_zone.add(Point3D(0.5, 0.5, 0.5))
		z = (mag ** 2 + mid.magnitude() ** 2 - 0.5) / (2 * mid.magnitude())
		if mag ** 2 - z ** 2 < 0.0: z = (mag ** 2 + mid.magnitude() ** 2 - 0.75) / (2 * mid.magnitude())
		x = math.sqrt(mag ** 2 - z ** 2)
		k_axis = mid.normalize()
		i_axis = Point3D(random.uniform(-1.0, 1.0), random.uniform(-1.0, 1.0), 0.0)
		i_axis.z = -(k_axis.x * i_axis.x + k_axis.y * i_axis.y) / k_axis.z
		i_axis = i_axis.normalize()
		j_axis = crossproduct(k_axis, i_axis).normalize()
		globali = Point3D(1.0, 0.0, 0.0).in_coordinate_system(i_axis, j_axis, k_axis)
		globalj = Point3D(0.0, 1.0, 0.0).in_coordinate_system(i_axis, j_axis, k_axis)
		globalk = Point3D(0.0, 0.0, 1.0).in_coordinate_system(i_axis, j_axis, k_axis)
		zpt = mid.multiply(z / mid.magnitude())
		#print zpt.add(i_axis.multiply(x)).distanceto(mid), zpt.add(i_axis.multiply(x)).magnitude()
		#Loosely using the sphere-sphere plane of intersection formula from Notability, "Science Fair Calculations."
		#print retro_zone
		idx = 0
		while (a == Point3D.zero() or a.floor() != retro_zone) and (not timeout or idx < 100):
			r = random.uniform(0.0, x)
			theta = random.uniform(0.0, 2 * math.pi)
			a = Point3D(r * math.cos(theta), r * math.sin(theta), 0.0)
			a = a.in_coordinate_system(globali, globalj, globalk).add(zpt)
			#Project this point onto the circular intersection region defined by center and normal mid.
			a = a.multiply(mag / a.magnitude())
			#print a.floor(), retro_zone
			idx += 1
		if idx == 100 and timeout:
			return None
		axes = alpha_zone.coordinate_system_for_transform(a)
		return axes

	def pz_representation(self):
		return PositionZone(self.acarbon, self.i, self.j, self.k)

	#MARK: Save/Restore
	def save(self, id):
		"""This function saves the current location and orientation of the amino acid to the stack."""
		if not hasattr(self, "_savestack"):
			self._savestack = {}
		self._savestack[id] = self.pz_representation()

	def restore(self, id):
		"""This function restores the location and orientation stored at id. Returns the location and orientation BEFORE the restoration."""
		if not hasattr(self, "_savestack"): return
		assert len(self._savestack) > 0, "Nothing to restore"
		curr = self.pz_representation()
		if id in self._savestack:
			loc = self._savestack[id]
			self.acarbon = loc.alpha_zone
			self.set_axes(loc.x_axis, loc.y_axis, loc.z_axis)
		return curr

	def discard_save(self, id):
		"""This function pops the location for id without restoring it."""
		del self._savestack[id]

	def clear_save(self):
		"""Removes all save caches."""
		self._savestack.clear()

#MARK: Helpers

def _rotate_vectors_pz(ptvectors, rot, length_ratio=1.0, translation=Point3D.zero()):
	"""Returns an array of position zones for the given pt vectors. They should be grouped in 4s for the position zone conversion: acarbon, x, y, z axes. The length ratio offers a means to stretch the bonds before """
	newpoints = np.dot([x.nparray() for x in ptvectors], rot.T).tolist()
	ret = []
	for i in xrange(0, len(newpoints), 4):
		ac = Point3D(newpoints[i][0] * length_ratio + translation.x,
					 newpoints[i][1] * length_ratio + translation.y,
					 newpoints[i][2] * length_ratio + translation.z)
		ret.append(PositionZone(ac,
								Point3D(newpoints[i + 1][0] * length_ratio + translation.x - ac.x,
										newpoints[i + 1][1] * length_ratio + translation.y - ac.y,
										newpoints[i + 1][2] * length_ratio + translation.z - ac.z).normalize(),
								Point3D(newpoints[i + 2][0] * length_ratio + translation.x - ac.x,
										newpoints[i + 2][1] * length_ratio + translation.y - ac.y,
										newpoints[i + 2][2] * length_ratio + translation.z - ac.z).normalize(),
								Point3D(newpoints[i + 3][0] * length_ratio + translation.x - ac.x,
										newpoints[i + 3][1] * length_ratio + translation.y - ac.y,
										newpoints[i + 3][2] * length_ratio + translation.z - ac.z).normalize()))
	return ret

def rotate_segment(aminoacids, newstart, newend):
	"""This helper function returns a list of alpha carbon locations for the amino acids in the list aminoacids when the first and last residues are moved to newstart and newend, respectively."""
	start = aminoacids[0].acarbon
	end = aminoacids[-1].acarbon
	orig_vec = end.subtract(start).normalize()
	new_vec = newend.subtract(newstart).normalize()
	if orig_vec == new_vec:
		#Translation
		offset = newend.subtract(end)
		return map(lambda x: PositionZone(x.acarbon.add(offset), x.i, x.j, x.k), aminoacids)
	else:
		rot = orig_vec.rotation_matrix(new_vec)
		ptvectors = []
		for x in aminoacids:
			ac = x.acarbon.subtract(start)
			ptvectors.append(ac)
			ptvectors.append(ac.add(x.i))
			ptvectors.append(ac.add(x.j))
			ptvectors.append(ac.add(x.k))
		old_length = end.subtract(start).magnitude()
		new_length = newend.subtract(newstart).magnitude()
		return _rotate_vectors_pz(ptvectors, rot, new_length / old_length, newstart)

def rotate_segment_anchor(aminoacids, anchor, beginning=True):
	"""anchor is a point representing the center of rotation. Pass in a PositionZone object to anchor, and a Boolean to beginning to specify whether the anchor comes at the beginning or the end of the segment. The provided position zone should indicate the new alpha carbon and axes of the amino acid closest to the anchor holding the segment to the chain - so, it would include the acarbon and axes provided by a PermissionsManager."""
	if beginning is True:
		anchor_aa = aminoacids[0]
	else:
		anchor_aa = aminoacids[-1]
	if isinstance(anchor_aa, PositionZone):
		aa = AminoAcid(amino_acid_alanine, -1, anchor_aa.alpha_zone)
		aa.set_axes(anchor_aa.x_axis, anchor_aa.y_axis, anchor_aa.z_axis)
		anchor_aa = aa
	#Find all segment points in the local coordinate system of the original anchor_aa.
	lcs_points = []
	for aa in aminoacids:
		if isinstance(aa, AminoAcid):
			lcs_points.append(anchor_aa.tolocal(aa.acarbon))
			lcs_points.append(anchor_aa.tolocal(aa.acarbon.add(aa.i)))
			lcs_points.append(anchor_aa.tolocal(aa.acarbon.add(aa.j)))
			lcs_points.append(anchor_aa.tolocal(aa.acarbon.add(aa.k)))
		elif isinstance(aa, PositionZone):
			lcs_points.append(anchor_aa.tolocal(aa.alpha_zone))
			lcs_points.append(anchor_aa.tolocal(aa.alpha_zone.add(aa.x_axis)))
			lcs_points.append(anchor_aa.tolocal(aa.alpha_zone.add(aa.y_axis)))
			lcs_points.append(anchor_aa.tolocal(aa.alpha_zone.add(aa.z_axis)))
	#Now create a hypothetical amino acid with the new anchor's orientation, and convert all those points into the global coordinate system from the hypothetical.
	hypothetical = anchor_aa.hypothetical(anchor, False)
	lcs_points = [hypothetical.toglobal(pt) for pt in lcs_points]
	ret = []
	for i in xrange(0, len(lcs_points), 4):
		ac = lcs_points[i]
		ret.append(PositionZone(ac, lcs_points[i + 1].subtract(ac), lcs_points[i + 2].subtract(ac), lcs_points[i + 3].subtract(ac)))
	if beginning is True:
		assert ret[0] == anchor, "Rotated amino acid is not the same: {} and {}.".format(ret[0], anchor)
	else:
		assert ret[-1] == anchor, "Rotated amino acid is not the same: {} and {}.".format(ret[-1], anchor)
	return ret

class AAAnchor(AminoAcid):
	"""AAAnchor objects represent pulling forces on a folding simulation. These can be used to ensure that a chain structure is maintained during simulation."""
	def __init__(self, type, tag=0, acarbon=Point3D.zero(), nitrogen=Point3D.zero(), carbon=Point3D.zero(), sidechain=None, weight=1.0, hook=0):
		"""Pass 0 for hook to apply to the first amino acid, or -1 to apply to the last amino acid. Any other number is treated as the index of the amino acid."""
		super(AAAnchor, self).__init__(type, tag, acarbon, nitrogen, carbon, sidechain)
		self.weight = weight
		self.hook = hook

	@classmethod
	def make(cls, aa, weight=1.0, hook=0):
		a2 = AAAnchor(aa.type, aa.tag, aa.acarbon, aa.nitrogen, aa.carbon, aa.sidechain, weight, hook)
		a2.set_axes(aa.i, aa.j, aa.k)
		return a2
		
bucket_count = 8000
location_min = -500.0
location_max = 500.0
box_dimension = 5.0
cubeWidth = 100.0

class AAHashTable(object):
	'Represents a hash table of amino acids that can be retrieved based on proximity.'
	
	def hash_function(self, point):
		modfn = math.fmod
		cubePoint = Point3D(modfn(point.x - location_min, cubeWidth), modfn(point.y - location_min, cubeWidth), modfn(point.z - location_min, cubeWidth))
		w = int(cubeWidth / box_dimension)
		hash = int(cubePoint.z / box_dimension) * w * w + int(cubePoint.y / box_dimension) * w + int(cubePoint.x / box_dimension)
		if hash > 8000:
			print point, ", ", cubePoint, ", ", hash
		return hash
	
	def __init__(self):
		self.buckets = [[] for x in xrange(bucket_count)]
	
	def __del__(self):
		for b in self.buckets:
			for aa in b:
				aa.remove_observer("acarbon", self.update_aa)
		del self.buckets[:]

	def clear(self):
		for b in self.buckets:
			for aa in b:
				aa.remove_observer("acarbon", self.update_aa)
			del b[:]
	
	def update_aa(self, aa, key, old):
		if key != "acarbon": return
		bucket = self.buckets[self.hash_function(old)]
		if (aa in bucket):
			bucket.remove(aa)
		else:
			print "Tried to observe an amino acid that's not in this hash table."
		bucket = self.buckets[self.hash_function(aa.acarbon)]
		if (aa not in bucket):
			bucket.append(aa)
	
	def add(self, aa):
		if not aa: return
		bucket = self.buckets[self.hash_function(aa.acarbon)]
		if (aa not in bucket):
			aa.add_observer("acarbon", self.update_aa)
			bucket.append(aa)

	def remove(self, aa):
		if not aa: return
		bucket = self.buckets[self.hash_function(aa.acarbon)]
		if (aa in bucket):
			aa.remove_observer("acarbon", self.update_aa)
			bucket.remove(aa)

	def contains(self, aa):
		if not aa: return
		bucket = self.buckets[self.hash_function(aa.acarbon)]
		return aa in bucket

	def nearby_aa(self, aa, distance, consec=2):
		"""You can also pass in a point for aa to check for simple proximity, with no consecutive/nonconsecutive considerations."""
		if not aa: return []
		check_tags = False
		if isinstance(aa, AminoAcid):
			check_tags = True
			pt = aa.acarbon
		elif isinstance(aa, Point3D): pt = aa
		else: return []
		bucket = self.buckets[self.hash_function(pt)]
		n = int(math.ceil(distance / box_dimension))
		retVal = []
		for z in xrange(-n, n + 1):
			for y in xrange(-n, n + 1):
				for x in xrange(-n, n + 1):
					searchBucket = self.hash_function(Point3D(pt.x + x * box_dimension,
															  pt.y + y * box_dimension,
															  pt.z + z * box_dimension))
					for candidate in self.buckets[searchBucket]:
						if check_tags and candidate.tag == aa.tag:	continue;
						firstDistance = aa.acarbon.distanceto(candidate.acarbon)
						if firstDistance <= distance and (not check_tags or consec == 2 or ((consec == True or math.fabs(candidate.tag - aa.tag) > 1) and (consec == False or math.fabs(candidate.tag - aa.tag) == 1))):
							retVal.append(candidate)
		return retVal

class PointHashTable(object):
	'Represents a hash table of points as keys to values that can be retrieved based on proximity.'
	
	def hash_function(self, point):
		cubePoint = Point3D(math.fmod(point.x - location_min, cubeWidth), math.fmod(point.y - location_min, cubeWidth), math.fmod(point.z - location_min, cubeWidth))
		w = cubeWidth / box_dimension
		hash = int(math.floor(cubePoint.z / box_dimension) * w * w + math.floor(cubePoint.y / box_dimension) * w + math.floor(cubePoint.x / box_dimension))
		if hash > 8000:
			print point, ", ", cubePoint, ", ", hash
		return hash
	
	def __init__(self, dict={}):
		self.buckets = [{} for x in range(bucket_count)]
		for pt, value in dict.iteritems():
			self.add(pt, value)
	
	def add(self, pt, value):
		bucket = self.buckets[self.hash_function(pt)]
		if (pt not in bucket):
			bucket[pt] = value

	def remove(self, pt):
		bucket = self.buckets[self.hash_function(pt)]
		if (pt in bucket): bucket.pop(pt)

	def contains(self, pt):
		bucket = self.buckets[self.hash_function(pt)]
		return pt in bucket

	def nearby_points(self, pt, distance):
		bucket = self.buckets[self.hash_function(pt)]
		n = int(math.ceil(distance / box_dimension))
		retVal = {}
		for z in range(-n, n + 1):
			for y in range(-n, n + 1):
				for x in range(-n, n + 1):
					searchBucket = self.hash_function(Point3D(pt.x + x * box_dimension,
															  pt.y + y * box_dimension,
															  pt.z + z * box_dimension))
					for candidate, value in self.buckets[searchBucket].iteritems():
						firstDistance = pt.distanceto(candidate)
						if firstDistance <= distance:
							retVal[candidate] = value
		return retVal