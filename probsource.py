"""The probsource module provides a class ProbabilitySource whose instances return cumulative distribution functions for use by the folding module. This module also includes functions for interconverting frequency, score, and probability.
	
	The provided subclass of ProbabilitySource, DistributionProbabilitySource, acts as a coordinator of DistributionManager objects. Initialize a DistributionProbabilitySource object with instances of the provided subclasses of DistributionManager, then pass it to folding_iteration() to apply frequency weighting."""

from molecular_systems import *
from functools import partial
import numpy as np
from numpy import matlib
from randomcoil import random_axes
import random
import datetime
import distributions
from memory_profiler import profile
import resource

#MARK: Helper functions

def apply_weight(probabilities, weight_function=lambda x: 1.0, passprob=False):
	"""This function applies a weight described by weight_function to every element in probabilities. weight_function is passed the conformation and asked for a weight to multiply the probability by. Pass in True for passprob to pass the probability to the weight function as a keyword argument 'prob'."""
	for x in probabilities:
		if passprob == True:
			x[1] *= weight_function(x[0], prob=x[1])
		else:
			x[1] *= weight_function(x[0])
		assert x[1] >= 0.0, "apply_weight has detected a negative weight for conformation %s. Please ensure that no negative weights are returned by your weight_function." % str(x[0])

def to_cdf(probabilities):
	"""This function generates a cumulative distribution function by setting the ith probability to the sum of all probabilities from 0 to i. probabilities should be sorted before calling this function."""
	sum = 0.0
	for x in probabilities:
		sum += x[1]
		x[1] = sum
	return probabilities

def score_to_probability(score):
	"""This function converts a score (which is negative for favorable results and positive for non-favorable ones) to a probability (which is always >= 0). This probability is not scaled to 0-1."""
	return math.exp(-score / 100.0)

def sample_cdf(probabilities):
	"""This helper function takes a list of probabilities in tuple form (conformation, probability number), representing a cumulative distribution function (cdf). It chooses a random one and returns its conformation."""
	rand = random.uniform(0.0, probabilities[-1][1])
	selected_conformation = next((x for x in probabilities if x[1] >= rand), None)
	assert selected_conformation is not None, "No valid conformation found with probability number greater than %.3f" % rand
	return selected_conformation[0]

#MARK: - ProbabilitySource

class ProbabilitySource(object):
	def __init__(self, protein, system=None):
		"""The default ProbabilitySource object has a rudimentary set of weights. Subclass __init__ to use your own weights. These weights should be applied last in the distribution-generating process (see probabilities)."""
		self.protein = protein
		self.system = system
		self.steric_cutoff = 4.0
		self.steric_consec_diff = 5

		def distance_weight(conformation, aminoacids=[], proximity=1):
			distance = sum(conformation[i].alpha_zone.distanceto(aminoacids[i].acarbon) for i in xrange(len(aminoacids))) / len(aminoacids)
			if proximity == 0.0: proximity = 0.01
			return max(0.0, -distance / proximity + 1.5)
		
		def bond_weight(conformation, point=None, point2=None):
			points = [x.alpha_zone for x in conformation]
			if point is not None: points.insert(0, point)
			if point2 is not None: points.append(point2)
			distance = 2
			for i in xrange(len(points) - 1):
				d = points[i].distanceto(points[i + 1])
				distance = min(-0.3333 * (d - 1.0) * (d - 3.5), distance)
			if distance == 2: distance = 1
			return max(0.0, distance)
		
		def steric_weight(conformation, aminoacids):
			for i, pz in enumerate(conformation):
				if self.system.check_steric_clash(aminoacids[i].hypothetical(pz), self.protein, 2.0):
					return 0.0
			return 1.0
		self.weights = { "distance" : distance_weight, "bond" : bond_weight, "steric" : steric_weight }
	
	def random_axis_rotation(self, aminoacids, pzs):
		"""If aminoacids is of length one, performs a rotation around the y-axis of the position zones. If it has more than one amino acid, the entire segment is rotated around the y-axis of the first amino acid. The return value is a new list of position zones."""
		if len(aminoacids) == 0:
			return
		theta = 0.0 # random.uniform(-math.pi / 3.0, math.pi / 3.0)
		#if random.randint(0, 2) == 0: theta += math.pi
		#The rotated axis vectors should be along the xz-plane in the amino acid's local coordinate system.
		aa0 = aminoacids[0].hypothetical(pzs[0])
		i = aa0.toglobal(Point3D(math.cos(theta), 0.0, math.sin(theta))).subtract(aa0.acarbon)
		k = aa0.toglobal(Point3D(-math.sin(theta), 0.0, math.cos(theta))).subtract(aa0.acarbon)
		if len(aminoacids) == 1:
			return [PositionZone(aa0.acarbon, i, aa0.j, k)]
		else:
			first_pz = PositionZone(aa0.acarbon, i, aa0.j, k)
			return rotate_segment_anchor(pzs, first_pz)
	
	def iter_ensemble_alpha(self, aminoacids, proximity=1, step=1):
		"""This method is an iterator over a series of possible destinations for the aminoacids (alpha carbons only). The return value is an array of position zones. Subclasses may want to use this method to assist in implementing probabilities()."""
		if len(aminoacids) == 1:
			residue = aminoacids[0]
			for alphapt in residue.acarbon.iter_randomoffsets(proximity, count=10):
				if alphapt.distanceto(residue.acarbon) > proximity: continue
				if self.system.check_steric_clash(residue.hypothetical(PositionZone(alphapt, residue.i, residue.j, residue.k)), self.protein, self.steric_cutoff, consec=False, mindiff=self.steric_consec_diff):
					continue
				pz_list = [PositionZone(alphapt, residue.i, residue.j, residue.k)]
				yield pz_list
				if self.mode == psource_erratic_mode and random.randint(0, 2) == 0:
					pz_list = self.random_axis_rotation(aminoacids, pz_list)
					yield pz_list
		else:
			#Conformations will be defined by the locations of the start and end points.
			start = aminoacids[0].acarbon
			end = aminoacids[-1].acarbon
			length = start.distanceto(end)
			for startpt in start.iter_randomoffsets(proximity, count=10):
				for endpt in end.iter_randomoffsets(proximity, count=10):
					#print start, end, startpt, endpt
					if start.distanceto(startpt) > proximity: continue
					if math.fabs(startpt.distanceto(endpt) - start.distanceto(end)) >= 0.5: continue
					endpt = endpt.subtract(startpt)
					endpt = endpt.multiply(length / endpt.magnitude()).add(startpt)
					pzs = rotate_segment(aminoacids, startpt, endpt)
					valid = True
					for i, pz in enumerate(pzs):
						if self.system.check_steric_clash(aminoacids[i].hypothetical(pz), self.protein, self.steric_cutoff, consec=False, mindiff=self.steric_consec_diff):
							valid = False
							break
					if valid != True: continue
					yield pzs
					if self.mode == psource_erratic_mode and random.randint(0, 2) == 0:
						pzs = self.random_axis_rotation(aminoacids, pzs)
						yield pzs
						
	def iter_ensemble_axis(self, aminoacid, proximity=0.1, step=0.1):
		"""This method is an iterator over a series of possible orientations for the aminoacid. The return value is an array of position zones. Subclasses may want to use this method to assist in implementing angleprobabilities()."""
		#for alphapt in aminoacid.acarbon.iteroffsets(1, 1):
		for x in np.arange(max(aminoacid.i.x - proximity, -1.0), min(aminoacid.i.x + proximity, 1.0) + step, step):
			for y in np.arange(max(aminoacid.i.y - proximity, -1.0), min(aminoacid.i.y + proximity, 1.0) + step, step):
				if math.sqrt(x ** 2 + y ** 2) > 1.0: continue
				#x = random.uniform(max(x - step / 2.0, -1.0), min(x + step / 2.0, 1.0))
				#y = random.uniform(max(y - step / 2.0, -1.0), min(y + step / 2.0, 1.0))
				z1 = math.sqrt(max(1.0 - x ** 2 - y ** 2, 0.0))
				z2 = -math.sqrt(max(1.0 - x ** 2 - y ** 2, 0.0))
				if abs(z1 - aminoacid.i.z) < abs(z2 - aminoacid.i.z):
					z = z1
				else:
					z = z2
				if Point3D(x, y, z).distanceto(aminoacid.i) < 0.01:
					yield [PositionZone(aminoacid.acarbon, aminoacid.i, aminoacid.j, aminoacid.k)]
					continue
				v = crossproduct(aminoacid.i, Point3D(x, y, z))
				s = v.magnitude()
				c = dotproduct(aminoacid.i, Point3D(x, y, z))
				skewsymmetric_cp = np.matrix( ((0, -v.z, v.y),
											   (v.z, 0, -v.x),
											   (-v.y, v.x, 0)) )
				rot = matlib.eye(3) + skewsymmetric_cp + skewsymmetric_cp * skewsymmetric_cp * ((1 - c) / (s ** 2))
				ptvectors = [aminoacid.i, aminoacid.j, aminoacid.k]
				newpoints = np.dot([pt.nparray() for pt in ptvectors], rot.T).tolist()
				assert len(newpoints) == 3, "Mismatched points: %r" % newpoints
				yield [PositionZone(aminoacid.acarbon,
									Point3D(newpoints[0][0], newpoints[0][1], newpoints[0][2]).normalize(),
									Point3D(newpoints[1][0], newpoints[1][1], newpoints[1][2]).normalize(),
									Point3D(newpoints[2][0], newpoints[2][1], newpoints[2][2]).normalize())]


	def probabilities(self, aminoacids, anchors=[], **kwargs):
		"""This function is called by folding simulators to determine a cumulative distribution function for the possible new locations of a segment of the protein. This default implementation returns position zones around the amino acids with a uniform probability. You should subclass ProbabilitySource to provide a better implementation."""
		proximity = 1
		probabilities = []
		if len(aminoacids) == 1:
			residue = aminoacids[0]
			#Set up the initial probabilities using the iterator helper function.
			for conformation in self.iter_ensemble_alpha(aminoacids, proximity):
				probabilities.append([conformation, 1.0])
			#Apply default distance and bond weights
			apply_weight(probabilities, partial(self.weights["distance"], aminoacids=[residue], proximity=proximity))
			if residue.tag > 0:
				apply_weight(probabilities, partial(self.weights["bond"], point=self.protein.aminoacids[residue.tag - 1].acarbon))
			if residue.tag < len(self.protein.aminoacids) - 1:
				apply_weight(probabilities, partial(self.weights["bond"], point=self.protein.aminoacids[residue.tag + 1].acarbon))
			apply_weight(probabilities, partial(self.weights["steric"], aminoacids=[residue]))
		else:
			#Set up the initial probabilities using the iterator helper function.
			for conformation in self.iter_ensemble_alpha(aminoacids, proximity):
				probabilities.append([conformation, 1.0])
			#Apply default distance and bond weights
			apply_weight(probabilities, partial(self.weights["distance"], aminoacids=aminoacids, proximity=proximity * 1.5))
			if aminoacids[0].tag > 0 and aminoacids[-1].tag < len(self.protein.aminoacids) - 1:
				apply_weight(probabilities, partial(self.weights["bond"],
													point=self.protein.aminoacids[aminoacids[0].tag - 1].acarbon,
													point2=self.protein.aminoacids[aminoacids[-1].tag + 1].acarbon))
			elif aminoacids[0].tag > 0:
				apply_weight(probabilities, partial(self.weights["bond"],
													point=self.protein.aminoacids[aminoacids[0].tag - 1].acarbon))
			elif aminoacids[-1].tag < len(self.protein.aminoacids) - 1:
				apply_weight(probabilities, partial(self.weights["bond"],
													point2=self.protein.aminoacids[aminoacids[-1].tag + 1].acarbon))
			else:
				apply_weight(probabilities, self.weights["bond"])
			apply_weight(probabilities, partial(self.weights["steric"], aminoacids=aminoacids))
		probabilities = [x for x in probabilities if x[1] > 0.0]
		if len(probabilities) == 0:
			probabilities = [[[PositionZone(x.acarbon, x.i, x.j, x.k) for x in aminoacids], 1.0]]
		return to_cdf(sorted(probabilities, key=lambda x: x[1]))

	def angleprobabilities(self, aminoacid, **kwargs):
		"""This function is called by folding simulators to determine a cumulative distribution function for the possible new orientations of a single amino acid in the protein (in its initial location)."""
		return [[[PositionZone(aminoacid.acarbon, aminoacid.i, aminoacid.j, aminoacid.k)], 1.0]]

	def randomization_margin(self, data, aminoacids):
		"""This function is called by folding simulators to determine how much the method should randomly mutate a conformation within its vicinity. It should pass a number that reflects how viable the conformation is."""
		return 0.1
	
	def choose_segment(self, segment_length):
		"""This function is called by folding simulators to choose a segment of the protein to be mutated. The returned segment should be of the length specified in segment_length. This default implementation returns an unweighted random segment.
			IMPORTANT: The return value for this function is not guaranteed to have the same length as segment_length."""
		random_start = random.randrange(len(self.protein.aminoacids))
		segment = self.protein.aminoacids[random_start : random_start + segment_length]
		return segment

	def cdf_from_hypotheticals(self, cases):
		"""This convenience method for folding simulators creates a cdf by evaluating the score function for each array of hypothetical amino acids in 'cases'. IMPORTANT: This method's returned probabilities have amino acid hypotheticals as members, not PositionZones."""
		probabilities = []
		for conformation in cases:
			probabilities.append([conformation, 1.0])
		#Apply default distance and bond weights
		apply_weight(probabilities, partial(self.weights["distance"], aminoacids=[residue], proximity=proximity))
		if residue.tag > 0:
			apply_weight(probabilities, partial(self.weights["bond"], point=self.protein.aminoacids[residue.tag - 1].acarbon))
		if residue.tag < len(self.protein.aminoacids) - 1:
			apply_weight(probabilities, partial(self.weights["bond"], point=self.protein.aminoacids[residue.tag + 1].acarbon))
		apply_weight(probabilities, partial(self.weights["steric"], aminoacids=[residue]))
		probabilities = [x for x in probabilities if x[1] > 0.0]
		return to_cdf(sorted(probabilities, key=lambda x: x[1]))

#MARK: - AAProbabilitySource

#These define the modes for AAProbabilitySource operation.
psource_erratic_mode = 1
psource_gentle_mode = 2

class AAProbabilitySource(ProbabilitySource):
	"""AAProbabilitySource provides the interface for coordinating DistributionManager objects to return probabilities for each possible conformation visited by a folding simulator."""
	
	def __init__(self, protein, distributions, permissions=None, sec_struct_permissions=None, system=None):
		super(AAProbabilitySource,self).__init__(protein, system)
		self.distributions = distributions
		self.permissions = permissions
		self.sec_struct_permissions = sec_struct_permissions
		def distance_weight(conformation, aminoacids=[], proximity=0.01):
			a = 100.0
			b = math.atan(a * proximity) + math.pi / 2.0
			c = 0.5 * math.pi / b
			distance = (conformation[0].alpha_zone.distanceto(aminoacids[0].acarbon) + conformation[-1].alpha_zone.distanceto(aminoacids[-1].acarbon)) / 2.0
			return -math.atan(a * (distance - 2 * proximity)) / b + c

		def bond_weight(conformation, point=None, point2=None):
			points = [x.alpha_zone for x in conformation]
			if point is not None: points.insert(0, point)
			if point2 is not None: points.append(point2)
			distance = 2
			for i in xrange(len(points) - 1):
				d = points[i].distanceto(points[i + 1])
				distance = min(-4 * (d - 1.5) * (d - 2.5), distance)
			if distance == 2: distance = 1
			return max(0.0, distance)

		def angle_distance_weight(conformation, aminoacid=None, proximity=0.15):
			distance = conformation[0].x_axis.distanceto(aminoacid.i)
			return max(0.0, -distance / proximity + 1.0)
		
		def comparison_weight(conformation, currentprob, prob=0.0):
			if currentprob == 0: return 1.0
			else:
				#print conformation, (prob / currentprob) ** 2
				return max(0.0, (prob / currentprob) ** 3)
	
		def anchor_weight(conformation, anchors=[], negweight=1.0):
			"""Pass in a float value for negweight to use in addition to the regular anchor weight if the distance constraint is not satisfied."""
			weight = 1.0
			for anchor in anchors:
				if anchor.hook == 0:
					d = AminoAcid.nlocation(conformation[anchor.hook]).distanceto(anchor.carbon)
				elif anchor.hook == -1:
					d = AminoAcid.clocation(conformation[anchor.hook]).distanceto(anchor.nitrogen)
				else:
					d = conformation[anchor.hook].alpha_zone.distanceto(anchor.acarbon)

				if anchor.hook == 0 and d > 1.5:
					target = anchor.carbon
					d = conformation[0].y_axis.multiply(-1.0).add(conformation[0].alpha_zone).distanceto(target) - conformation[0].alpha_zone.distanceto(target)
					if d > -0.6: weight /= anchor.weight * (d + 0.6)
				elif anchor.hook == -1 and d > 1.5:
					target = anchor.nitrogen
					d = conformation[-1].y_axis.add(conformation[-1].alpha_zone).distanceto(target) - conformation[-1].alpha_zone.distanceto(target)
					if d > -0.6: weight /= anchor.weight * (d + 0.6)
				elif anchor.hook != 0 and anchor.hook != -1 and (d < 2.0 or d > 3.5):
					if anchor.weight >= 1e9: weight = 0
					else:
						weight /= anchor.weight * d
						weight *= negweight / d
			return max(0.0, weight)
		
		self.weights = { "distance": distance_weight, "bond" : bond_weight, "angle_distance" : angle_distance_weight, "comparison" : comparison_weight, "anchor" : anchor_weight }
		self.mode = psource_erratic_mode
		self.erratic_proximity = 4.0
	
	def iter_ensemble_pivot(self, aminoacids, anchor, prior=True, proximity=1.0, maxsamples=-1, connected=False):
		"""This new method (2/6/15) iterates the possible conformations created by pivoting a segment around the anchor (should be an AminoAcid). This function does not take into account the connectivity with the next unmutated segment. However, the connected parameter specifies whether the anchor is connected to the protein on the other side. This helps to produce realistic angles during a series of pivots."""
		assert self.permissions is not None, "Cannot perform a pivot without a PermissionsManager."
		if prior is True:
			aa = aminoacids[0]
		else:
			aa = aminoacids[-1]
		sec_struct = self.protein.secondary_structure_aa(aa.tag)
		opposite = None
		if prior and anchor.tag > 0 and connected:
			opposite = self.protein.aminoacids[anchor.tag - 1]
		elif not prior and anchor.tag < len(self.protein.aminoacids) - 1 and connected:
			opposite = self.protein.aminoacids[anchor.tag + 1]
		if sec_struct and self.sec_struct_permissions and ((prior and sec_struct[1].start < aa.tag) or (not prior and sec_struct[1].end > aa.tag)):
			allowed_conformations = self.sec_struct_permissions.allowed_conformations(aa, anchor, sec_struct[0].type, sec_struct[1].identifiers[0], prior=prior, opposite_aa=opposite)
		else:
			allowed_conformations = self.permissions.allowed_conformations(aa, anchor, prior, opposite_aa=opposite)
		
		random.shuffle(allowed_conformations)
		numsamples = 0
		for pz in allowed_conformations:
			if len(aminoacids) > 1:
				conformation = rotate_segment_anchor(aminoacids, pz, prior)
			else:
				conformation = [pz]
			if prior is True:
				if len(aminoacids) == 1 and conformation[0].x_axis.distanceto(aminoacids[0].i) > proximity: continue
				elif conformation[-1].alpha_zone.distanceto(aminoacids[-1].acarbon.subtract(aminoacids[0].acarbon).add(conformation[0].alpha_zone)) > proximity * 3.0: continue
			else:
				if len(aminoacids) == 1 and conformation[-1].x_axis.distanceto(aminoacids[-1].i) > proximity: continue
				elif conformation[0].alpha_zone.distanceto(aminoacids[0].acarbon.subtract(aminoacids[-1].acarbon).add(conformation[-1].alpha_zone)) > proximity * 3.0: continue

			if not self.aa_connected(aa.hypothetical(pz), anchor):
				continue
			valid = True
			for i, pz in enumerate(conformation):
				if self.system.check_steric_clash(aminoacids[i].hypothetical(pz), self.protein, self.steric_cutoff, consec=False, mindiff=self.steric_consec_diff, excluded=[a.tag for a in aminoacids]):
					valid = False
					break
			if valid is not True: continue
			yield conformation
			numsamples += 1
			if maxsamples > 0 and numsamples > maxsamples:
				return

	def probabilities(self, aminoacids, anchors=[], **kwargs):
		"""This overridden function generates a list of possible conformations and their probabilities, as does its superclass's implementation. However, this implementation consults its list of DistributionManagers, polling their score() method for each amino acid in the segment. These scores are summed for each configuration, and the resulting score is transformed into a probability.
			NOTES: Pass in an integer for keyword argument "primanchor" to signify which anchor is used as a pivot point for mutation. Pivoting will not be used unless you specify a value for "primanchor". Also, you should provide a Boolean for keyword argument "prior" to signify whether or not primanchor is before the segment. If not provided, "prior" is assumed True. Pass the keyword "connected" in a cascade of modifications to show whether the anchor(s) are connected to their outer segments."""
		probabilities = []
		current_score = sum(d.score(self.protein, aminoacids, system=self.system) for d in self.distributions)
		if self.mode == psource_erratic_mode:
			proximity = self._randomization_margin(current_score / len(aminoacids)) * self.erratic_proximity
			step = proximity
		else:
			proximity = 0.1
			step = 0.1

		if "primanchor" in kwargs:
			if "prior" in kwargs: prior = kwargs["prior"]
			else: prior = True
			current_score = 0.0
			for d in self.distributions:
				if d.type != distributions.frequency_consec_disttype:
					current_score += d.score(self.protein, aminoacids, system=self.system)
				else:
					current_score += d.score(self.protein, aminoacids, system=self.system, prior=prior)
			connected = False
			if "connected" in kwargs: connected = kwargs["connected"]
			for conformation in self.iter_ensemble_pivot(aminoacids, anchors[kwargs["primanchor"]], prior, proximity * 2.0, connected=connected):
				hypotheticals = [aminoacids[i].hypothetical(conformation[i]) for i in xrange(len(aminoacids))]
				score = 0.0
				invalid = False
				for d in self.distributions:
					if d.type != distributions.frequency_consec_disttype:
						subscore = d.score(self.protein, hypotheticals, system=self.system)
					else:
						subscore = d.score(self.protein, hypotheticals, prior=prior, system=self.system)
					if d.type == distributions.frequency_nonconsec_disttype and d.short_range == True:
						if subscore > 5.0:
							invalid = True
							break
					score += subscore
				if invalid: continue
				if self.mode != psource_gentle_mode or score < current_score:
					probabilities.append([conformation, score_to_probability(score / len(hypotheticals))])
				del hypotheticals[:]
			if len(probabilities) == 0:
				for conformation in self.iter_ensemble_pivot(aminoacids, anchors[kwargs["primanchor"]], prior, 1000.0, connected=connected):
					hypotheticals = [aminoacids[i].hypothetical(conformation[i]) for i in xrange(len(aminoacids))]
					score = 0.0
					invalid = False
					for d in self.distributions:
						if d.type != distributions.frequency_consec_disttype:
							subscore = d.score(self.protein, hypotheticals, system=self.system)
						else:
							subscore = d.score(self.protein, hypotheticals, prior=prior, system=self.system)
						if d.type == distributions.frequency_nonconsec_disttype and d.short_range == True:
							if subscore > 5.0:
								invalid = True
								break
						score += subscore
					if invalid: continue
					if self.mode != psource_gentle_mode or score < current_score:
						probabilities.append([conformation, score_to_probability(score / len(hypotheticals))])
					del hypotheticals[:]
		else:
			current_score = sum(d.score(self.protein, aminoacids, system=self.system) for d in self.distributions if d.type != distributions.frequency_consec_disttype)
			for conformation in self.iter_ensemble_alpha(aminoacids, proximity, step=step):
				hypotheticals = [aminoacids[i].hypothetical(conformation[i]) for i in xrange(len(aminoacids))]
				score = sum(d.score(self.protein, hypotheticals, system=self.system) for d in self.distributions if d.type != distributions.frequency_consec_disttype)
				if self.mode != psource_gentle_mode or score < current_score:
					probabilities.append([conformation, score_to_probability(score / len(hypotheticals))])
				del hypotheticals[:]
			orig = [[PositionZone(aa.acarbon, aa.i, aa.j, aa.k) for aa in aminoacids], score_to_probability(current_score / len(aminoacids))]
			if orig not in probabilities:
				probabilities.append(orig)
		apply_weight(probabilities, partial(self.weights["distance"], aminoacids=aminoacids, proximity=self._randomization_margin(current_score / len(aminoacids)) * 4.0))
		apply_weight(probabilities, partial(self.weights["comparison"], currentprob=score_to_probability(current_score / len(aminoacids))), passprob=True)
		apply_weight(probabilities, partial(self.weights["anchor"], anchors=anchors, negweight=proximity))

		probabilities = [x for x in probabilities if x[1] > 0.0]
		if len(probabilities) == 0 and "primanchor" not in kwargs:
			probabilities.append([[PositionZone(aa.acarbon, aa.i, aa.j, aa.k) for aa in aminoacids], score_to_probability(current_score / len(aminoacids))])
		#We're only returning the top 25 probabilities! Need to add a control over this for the user
		return to_cdf(sorted(probabilities, key=lambda x: x[1])[max(len(probabilities) - 25, 0):])
			
	def _iter_permissible_randomcoils(self, aminoacids, start, end, chain=[], yieldct=0):
		"""This helper function is a recursive generator for permissible self-avoiding walks."""
		if len(chain) == len(aminoacids) - 1:
			#Connect the last amino acid in a permissible orientation to both chain[-1] and end.
			if start.tag < end.tag:
				second_last = None
				if len(chain) > 0:
					last = aminoacids[len(chain) - 1].hypothetical(chain[-1])
					if len(chain) > 1:
						second_last = aminoacids[len(chain) - 2].hypothetical(chain[-2])
				else:
					last = start
					if start.tag > 0:
						second_last = self.protein.aminoacids[start.tag - 1]
				sec_struct = self.protein.secondary_structure_aa(aminoacids[-1].tag)
				if sec_struct and self.sec_struct_permissions and sec_struct[1].start < aminoacids[-1].tag:
					confs = self.sec_struct_permissions.allowed_conformations(aminoacids[-1], last, sec_struct[0].type, sec_struct[1].identifiers[0], opposite_aa=second_last)
					confs2 = self.sec_struct_permissions.allowed_conformations(aminoacids[-1], end, sec_struct[0].type, sec_struct[1].identifiers[0], prior=False)
				else:
					confs = self.permissions.allowed_conformations(aminoacids[-1], last, opposite_aa=second_last)
					confs2 = self.permissions.allowed_conformations(aminoacids[-1], end, prior=False)
				random.shuffle(confs)
				random.shuffle(confs2)
				for pz in confs:
					for pz2 in confs2:
						if pz2.alpha_zone.distanceto(pz.alpha_zone) > 0.5: continue
						pz = PositionZone(pz.alpha_zone.add(pz2.alpha_zone).multiply(0.5),
										  pz.x_axis, pz.y_axis, pz.z_axis)
						hyp = aminoacids[-1].hypothetical(pz)
						if self.aa_connected(hyp, end) and self.aa_connected(hyp, last):
							if self.system.check_steric_clash(hyp, self.protein, self.steric_cutoff, consec=False, mindiff=self.steric_consec_diff, excluded=[aa.tag for aa in aminoacids]):
								continue
							yieldct += 1
							yield chain + [pz]
							if (yieldct >= 10 and not sec_struct) or yieldct == 20: return
			else:
				for saw in self._iter_permissible_randomcoils(aminoacids, end, start, chain, yieldct):
					yieldct += 1
					yield saw
		else:
			#Try to add one more amino acid in a permissible orientation.
			if start.tag < end.tag:
				second_last = None
				if len(chain) > 0:
					last = aminoacids[len(chain) - 1].hypothetical(chain[-1])
					if len(chain) > 1:
						second_last = aminoacids[len(chain) - 2].hypothetical(chain[-2])
				else:
					last = start
					if start.tag > 0:
						second_last = self.protein.aminoacids[start.tag - 1]
				sec_struct = self.protein.secondary_structure_aa(aminoacids[len(chain)].tag)
				if sec_struct and self.sec_struct_permissions and sec_struct[1].start < aminoacids[len(chain)].tag:
					confs = self.sec_struct_permissions.allowed_conformations(aminoacids[len(chain)], last, sec_struct[0].type, sec_struct[1].identifiers[0], opposite_aa=second_last)
				else:
					confs = self.permissions.allowed_conformations(aminoacids[len(chain)], last, opposite_aa=second_last)
				#I'm using a constant number here to save on performance. Only 10 of the given conformations will be chosen, which still gives a total of 10^(n - 1) iterations before the last residue.
				if len(confs) > 30:
					confs = random.sample(confs, 30)
				random.shuffle(confs)
				iters = 0
				for pz in confs:
					hyp = aminoacids[len(chain)].hypothetical(pz, True)
					if not self.connectivity_possible(hyp, end, len(aminoacids) - len(chain) - 1): continue
					if self.system.check_steric_clash(hyp, self.protein, self.steric_cutoff, consec=False, mindiff=self.steric_consec_diff, excluded=[aa.tag for aa in aminoacids]):
						continue
					iters += 1
					for saw in self._iter_permissible_randomcoils(aminoacids, start, end, chain + [pz], yieldct):
						yieldct += 1
						yield saw
						if (yieldct >= 10 and not sec_struct) or yieldct == 20: return
				#print len(chain), "did", iters
			else:
				for saw in self._iter_permissible_randomcoils(aminoacids, end, start, chain, yieldct):
					yieldct += 1
					yield saw

	def randomcoil_probabilities(self, aminoacids, start, end):
		"""Pass in an array of amino acids, a start anchor amino acid (not included in aminoacids), and an end anchor amino acid. Each entry in the returned probability distribution will be a permissible SAW from start to end. If no permissible SAW is found that connects start and end, this method will return no probabilities, in which case you should use probabilities() instead."""
		probabilities = []
		current_score = sum(d.score(self.protein, aminoacids, system=self.system) for d in self.distributions)

		for conformation in self._iter_permissible_randomcoils(aminoacids, start, end):
			hypotheticals = [aminoacids[i].hypothetical(conformation[i]) for i in xrange(len(aminoacids))]
			probabilities.append([conformation, score_to_probability(sum(d.score(self.protein, hypotheticals, system=self.system) for d in self.distributions))])
			if not self.aa_connected(hypotheticals[0], start):
				print "AA not connected to start random coil", hypotheticals[0], start
			if not self.aa_connected(hypotheticals[-1], end):
				print "AA not connected to end random coil", hypotheticals[-1], end
			del hypotheticals[:]


		apply_weight(probabilities, partial(self.weights["comparison"], currentprob=score_to_probability(current_score)), passprob=True)
		probabilities = [x for x in probabilities if x[1] > 0.0]
		if len(probabilities) == 0:
			return []
		return to_cdf(sorted(probabilities, key=lambda x: x[1]))

	
	def angleprobabilities(self, aminoacid, **kwargs):
		"""This overridden function generates a certain number of random axis configurations for aminoacid, then computes the probability of each one occurring using the DistributionManager score() method."""
		probabilities = []
		current_score = sum(d.score(self.protein, [aminoacid], system=self.system) for d in self.distributions)
		proximity = max(self._randomization_margin(current_score), 0.001)
		step = proximity
		for conformation in self.iter_ensemble_axis(aminoacid, proximity=proximity, step=step):
			hypothetical = aminoacid.hypothetical(conformation[0])
			probabilities.append([conformation, score_to_probability(sum(d.score(self.protein, [hypothetical], system=self.system) for d in self.distributions))])
		
		probabilities.append([[PositionZone(aminoacid.acarbon, aminoacid.i, aminoacid.j, aminoacid.k)], score_to_probability(current_score)])
		apply_weight(probabilities, partial(self.weights["distance"], aminoacids=[aminoacid]))
		apply_weight(probabilities, partial(self.weights["angle_distance"], aminoacid=aminoacid))
		apply_weight(probabilities, partial(self.weights["comparison"], currentprob=score_to_probability(current_score)), passprob=True)

		probabilities = [x for x in probabilities if x[1] > 0.0]
		if len(probabilities) == 0:
			probabilities.append([[PositionZone(aminoacid.acarbon, aminoacid.i, aminoacid.j, aminoacid.k)], score_to_probability(current_score)])
		return to_cdf(sorted(probabilities, key=lambda x: x[1]))

	def connectivity_possible(self, start, end, numaa):
		"""Pass in AminoAcid objects for start and end, and an int for numaa. This function will determine whether it is possible to build a permissible SAW from start to end using numaa amino acids. Note that a True return value does not guarantee the existence of a permissible SAW, but a False value does guarantee the non-existence of one.
			Updated 2/22/15: This method returns FALSE if the segment is at either end of the chain. It does so to notify folding simulators that they can use normal pivoting procedures to determine the location of this segment."""

		#Updated 2/22/15 to use a new method. We will check if the distance between start and end is at an acceptable level relative to the distance between their outer neighbors.
		if start.tag == 0 or end.tag == len(self.protein.aminoacids) - 1:
			return False
		else:
			l = self.protein.aminoacids[start.tag - 1].acarbon.distanceto(self.protein.aminoacids[end.tag + 1].acarbon)
			min = l - 6.0 * math.sqrt(3.0)
			max = 7.0 * numaa #2 * math.sqrt(l ** 2 / 4.0 + 27.0)
			d = start.nitrogen.distanceto(end.carbon)
			if d >= min and d <= max:
				return True
			else:
				return False

	def aa_connected(self, aa1, aa2):
		"""This function detects whether the amino acids at the beginning and end of segment are connected to their neighbors outside the segment."""
		if aa1.tag < aa2.tag:
			pre = aa1
			post = aa2
		else:
			pre = aa2
			post = aa1
		ss_pre = self.protein.secondary_structure_aa(pre.tag)
		ss_post = self.protein.secondary_structure_aa(post.tag)
		if ss_pre and ss_post and ss_pre[1].start == ss_post[1].start:
			return self.sec_struct_permissions.is_valid(post, pre, ss_pre[0].type, ss_pre[1].identifiers[0], prior=True)
		else:
			return self.permissions.is_valid(post, pre, prior=True)

	def is_connected(self, segment):
		"""This function detects whether the amino acids at the beginning and end of segment are connected to their neighbors outside the segment."""
		if segment[0].tag == 0 and segment[-1].tag == len(self.protein.aminoacids) - 1:
			return True
		elif segment[0].tag == 0:
			return self.aa_connected(segment[-1], self.protein.aminoacids[segment[-1].tag + 1])
		elif segment[-1].tag == len(self.protein.aminoacids) - 1:
			return self.aa_connected(segment[0], self.protein.aminoacids[segment[0].tag - 1])
		else:
			return self.aa_connected(segment[0], self.protein.aminoacids[segment[0].tag - 1]) and self.aa_connected(segment[-1], self.protein.aminoacids[segment[-1].tag + 1])

	def choose_segment(self, segment_length):
		#First we have to cluster out the chain into groups that have similarly stable consecutive scores.
		clusters = [[0, 0, 0, 100000, -100000]] #Start, end, total/avg score, min score, max score
		cluster_range = [0.0, 0.0]
		min_bound = 0.9
		max_bound = 1.1
		for i, aa in enumerate(self.protein.aminoacids):
			if aa.localscore == 0.0:
				aa.localscore = sum(d.score(self.protein, [aa], system=self.system) for d in self.distributions)
		
			#Find a secondary structure if this amino acid is in one
			current_struct = None
			struct_info = self.protein.secondary_structure_aa(i)
			if struct_info:
				current_struct = struct_info[1]
			
			def matches_score(ls, clus_r):
				return ls >= clus_r[0] and ls <= clus_r[1] and ls <= -125.0
			
			if clusters[-1][1] == 0:
				clusters[-1][1] += 1
				clusters[-1][2] += aa.localscore
				clusters[-1][3] = aa.localscore
				clusters[-1][4] = aa.localscore
			elif current_struct and random.randint(0, 100) < 60:
				if current_struct.start == i and not matches_score(aa.localscore, cluster_range):
					clusters[-1][2] /= float(clusters[-1][1] - clusters[-1][0])
					clusters.append([i, i + 1, aa.localscore, aa.localscore, aa.localscore])
				else:
					clusters[-1][1] += 1
					clusters[-1][2] += aa.localscore
				if aa.localscore < clusters[-1][3]: clusters[-1][3] = aa.localscore
				if aa.localscore > clusters[-1][4]: clusters[-1][4] = aa.localscore
			elif matches_score(aa.localscore, cluster_range):
				clusters[-1][1] += 1
				clusters[-1][2] += aa.localscore
				if aa.localscore < clusters[-1][3]: clusters[-1][3] = aa.localscore
				if aa.localscore > clusters[-1][4]: clusters[-1][4] = aa.localscore
			else:
				clusters[-1][2] /= float(clusters[-1][1] - clusters[-1][0])
				clusters.append([i, i + 1, aa.localscore, aa.localscore, aa.localscore])
			cluster_range[0] = min(clusters[-1][3] * min_bound, clusters[-1][4] * max_bound)
			cluster_range[1] = max(clusters[-1][3] * min_bound, clusters[-1][3] * max_bound)
			
		clusters[-1][2] /= float(clusters[-1][1] - clusters[-1][0])
		#print "---", len(clusters), "clusters for", len(self.protein.aminoacids), "amino acids"
		for clust in clusters:
			tupcluster = (clust[0], clust[1])
			for i in xrange(*tupcluster):
				self.protein.aminoacids[i].cluster = tupcluster
				self.protein.aminoacids[i].clusterscore = clust[2]
		self.clusters = clusters
		weights = [score_to_probability(-1 * sum(d.score(self.protein, self.protein.aminoacids[i:i + segment_length], system=self.system) for d in self.distributions)) for i in xrange(len(self.protein.aminoacids) - 1)]
		s = sum(weights)
		weights = [w / s for w in weights]
		random_start = np.random.choice(range(len(self.protein.aminoacids) - 1), p=weights)
		selected_cluster = self.protein.aminoacids[random_start].cluster
		if selected_cluster[1] - selected_cluster[0] == 1:
			segment = self.protein.aminoacids[random_start : random_start + segment_length]
		else:
			segment = self.protein.aminoacids[selected_cluster[0] : selected_cluster[1]]
		return segment
	
	def _randomization_margin(self, score):
		"""This function yields 0 for a score of -25, and 0.1 for a score of 5.
			Updated 2/13/15: this function yields 0 for a score of -180, and 0.4 for a score of -50."""
		if self.mode == psource_gentle_mode:
			return 10.8 ** (score / 40.0)
		else:
			return max(0.0, 1.0 / 325.0 * (score + 180.0))
	
	def randomization_margin(self, data, aminoacids):
		"""Pass in an array of position zones for data, or pass in None for data and an array of amino acids."""
		if data is not None:
			hypotheticals = [aminoacids[i].hypothetical(data[i]) for i in xrange(len(aminoacids))]
		else:
			hypotheticals = aminoacids
		score = sum(d.score(self.protein, hypotheticals, system=self.system) for d in self.distributions)
		return self._randomization_margin(score)

	def bond_breach(self, aa1, aa2):
		return aa1.acarbon.distanceto(aa2.acarbon) - self.steric_cutoff

	def cdf_from_hypotheticals(self, cases):
		probabilities = []
		for hypotheticals in cases:
			probabilities.append([hypotheticals, score_to_probability(sum(d.score(self.protein, hypotheticals, system=self.system) for d in self.distributions))])
		probabilities = [x for x in probabilities if x[1] > 0.0]
		return to_cdf(sorted(probabilities, key=lambda x: x[1]))

#MARK: - AAConstructiveProbabilitySource

class AAConstructiveProbabilitySource(AAProbabilitySource):
	"""This class copies all the functionality from AAProbabilitySource except that it allows you to run simulations with two big clusters of amino acids, mutating the orientations between them."""

	def __init__(self, protein, segment1, segment2, distributions, permissions=None, sec_struct_permissions=None, system=None):
		"""For segment1 and segment2, pass in tuples (min, max). For instance, (0, 5) represents amino acids 0-5 for 5 amino acids total."""
		super(AAConstructiveProbabilitySource,self).__init__(protein, distributions, permissions, sec_struct_permissions, system)
		#These are tuples defining the range of each cluster.
		self.cluster1 = segment1
		self.cluster2 = segment2
		#These are lists of tuples: (list of position zones, SPARC score).
		self.c1_conformations = []
		self.c2_conformations = []
	
	def load_cluster_conformations(self, cluster_num, path, n=25):
		'''This new function loads the conformations for one of the segments of the simulation. If the SPARC scores are stored in the comment lines of the multi-model PDB file you pass in at path, they will be read in with the format "SPARC xxxx". Otherwise, the SPARC scores will be computed within this method.'''

		curr_score = 0.0
		peptide = Polypeptide()
		model_count = n
		confs = [[0, 10000000, 0] for i in xrange(model_count)]
		backups = [[0, 10000000, 0] for i in xrange(model_count * 3)]

		for modelno in peptide.iter_models(path):
			curr_score = sum(d.score(peptide, peptide.aminoacids) for d in self.distributions)
			add_to_backup = True
			for k in xrange(model_count):
				if curr_score < confs[k][1] * 1.1:
					for m in reversed(xrange(max(k, 1), model_count)):
						confs[m] = confs[m - 1]
					confs[k] = [modelno, curr_score, score_to_probability(curr_score)]
					add_to_backup = False
					break
				elif curr_score < confs[k][1]:
					break
			if add_to_backup:
				for k in xrange(model_count):
					if curr_score < backups[k][1]:
						for m in reversed(xrange(max(k, 1), model_count)):
							backups[m] = backups[m - 1]
						backups[k] = [modelno, curr_score, score_to_probability(curr_score)]
						break

			#confs.append([modelno, curr_score, score_to_probability(curr_score)])
			if modelno % 50 == 0:
				print "Loaded", modelno, "conformations..."

		for i, c in enumerate(confs):
			if c[1] > 10000.0:
				confs[i] = backups[0]
				del backups[0]
		i = 0
		while i < len(confs):
			if confs[i][1] > 10000.0:
				del confs[i]
			else:
				i += 1
		
		confs2 = sorted(confs, key=lambda x: x[1])
		#Go back and read through the file, saving the conformations that we noted above
		last_read = 0
		modelno = 0
		peptide = None
		with open(path, 'r') as file:
			lines = file.readlines()
			cache = None
			for i, l in enumerate(lines):
				if "MODEL" in l:
					modelno = int(l[5:].strip())
				elif "ENDMDL" in l:
					idx = next((k for k, x in enumerate(confs2) if x[0] == modelno), -1)
					if idx >= 0:
						if not peptide:
							peptide = Polypeptide()
						try:
							peptide.read_file(lines[last_read : i + 1], cache_aas=cache)
						except AssertionError as e:
							print "Exception reading structure!"
							last_read = i + 1
							del peptide.aminoacids[:]
							peptide.hashtable = None
							del confs2[idx]
							continue
						confs2[idx][0] = [aa.pz_representation() for aa in peptide.aminoacids]
						cache = peptide.aminoacids
					last_read = i + 1
			del lines, cache
		del confs
		print "Loaded", len(confs2), "conformations for segment", cluster_num
		for i, c in enumerate(confs2):
			assert len(c[0]), "Invalid conf {} {}".format(i, c)
		if cluster_num == 1:
			self.c1_conformations = confs2
		elif cluster_num == 2:
			self.c2_conformations = confs2
		else:
			print "AAConstructiveProbabilitySource doesn't support segment number", cluster_num, ". Sorry!"

	def generate_structure_from_segments(self, sequence, c1_struct, c2_struct):
		'''This method combines the structures in c1_struct and c2_struct with a random connection orientation.'''
		aminoacids = []
		hashtable = AAHashTable()
		
		assert len(c1_struct) == self.cluster1[1] - self.cluster1[0], "Mismatched lengths: {} position zones for cluster {}".format(len(c1_struct), self.cluster1)
		for pz in c1_struct:
			aminoacid = AminoAcid(aatypec(sequence[len(aminoacids)]), tag=len(aminoacids), acarbon=pz.alpha_zone)
			aminoacid.set_axes(pz.x_axis, pz.y_axis, pz.z_axis)
			hashtable.add(aminoacid)
			aminoacids.append(aminoacid)

		#Then get permissible orientations for the last residue in c1 and the first residue in c2
		pivot = AminoAcid(aatypec(sequence[len(aminoacids)]), tag=len(aminoacids))
		second_last = None
		if len(aminoacids) > 1: second_last = aminoacids[-2]
		candidates = self.permissions.allowed_conformations(pivot, aminoacids[-1], opposite_aa=second_last)
		if len(candidates) == 0:
			candidates = self.permissions.allowed_conformations(pivot, aminoacids[-1])
		assert len(candidates) > 0, "No permissible candidates for pivot"
		while len(candidates) > 0 and (pivot.acarbon == Point3D.zero() or len(hashtable.nearby_aa(pivot, self.steric_cutoff, consec=False)) > 0):
			zone = random.choice(candidates)
			pivot.acarbon = zone.alpha_zone
			pivot.set_axes(zone.x_axis, zone.y_axis, zone.z_axis)
			candidates.remove(zone)
		hashtable.add(pivot)
		aminoacids.append(pivot)

		assert len(c2_struct) == self.cluster2[1] - self.cluster2[0], "Mismatched lengths: {} position zones for cluster {}".format(len(c2_struct), self.cluster2)
		pivot_hypo = pivot.hypothetical(c2_struct[0])
		for	pz in c2_struct[1:]:
			aminoacid = AminoAcid(aatypec(sequence[len(aminoacids)]), tag=len(aminoacids))
			pivoted_pz = pivot.globalpz(pivot_hypo.localpz(pz))
			aminoacid.acarbon = pivoted_pz.alpha_zone
			aminoacid.set_axes(pivoted_pz.x_axis, pivoted_pz.y_axis, pivoted_pz.z_axis)
			hashtable.add(aminoacid)
			aminoacids.append(aminoacid)

		return (aminoacids, hashtable)
	
	def generate_initial_structure(self, sequence):
		'''This method creates an initial structure by combining the two segments in a random way. Provide the amino acid sequence of the entire segment pair so that we know what amino acids to add.'''
		#First choose a random structure in c1_conformations.
		weights = [c[2] for c in self.c1_conformations]
		total = float(sum(weights))
		weights = [w / total for w in weights]
		c1_struct = self.c1_conformations[np.random.choice(len(self.c1_conformations), p=weights)][0]

		#Choose a random structure in c2_conformations that does not cause a steric clash.
		weights = [c[2] for c in self.c2_conformations]
		total = float(sum(weights))
		weights = [w / total for w in weights]
		c2_struct = self.c2_conformations[np.random.choice(len(self.c2_conformations), p=weights)][0]

		test_idx = 0
		aas, hashtable = self.generate_structure_from_segments(sequence, c1_struct, c2_struct)
		
		#Check for steric violation
		while next((aa for aa in aas if len(hashtable.nearby_aa(aa, self.steric_cutoff, consec=False))), None):
			if test_idx == 50:
				self.generate_initial_structure(sequence)
				return
			for i in xrange(self.cluster2[0] + 1, self.cluster2[1]):
				hashtable.remove(aminoacids[i])
			del aminoacids[self.cluster2[0] + 1 : self.cluster2[1]]
			c2_struct = self.c2_conformations[np.random.choice(len(self.c2_conformations), p=weights)][0]
			aas, hashtable = self.generate_structure_from_segments(sequence, c1_struct, c2_struct)
			test_idx += 1
		self.protein.add_aas(aas)


	def choose_segment(self, segment_length):
		'''This subclassed method ignores the parameter segment_length, instead choosing one of the pivot amino acids to return.'''
		cluster1_aa = self.protein.aminoacids[self.cluster1[0] : self.cluster1[1]]
		cluster2_aa = self.protein.aminoacids[self.cluster2[0] : self.cluster2[1]]
		cluster1_score = sum(d.score(self.protein, cluster1_aa, system=self.system) for d in self.distributions)
		cluster2_score = sum(d.score(self.protein, cluster2_aa, system=self.system) for d in self.distributions)
		for aa in cluster1_aa:
			if aa.localscore == 0.0:
				aa.localscore = sum(d.score(self.protein, [aa], system=self.system) for d in self.distributions)
			aa.cluster = self.cluster1
			aa.clusterscore = cluster1_score
		for aa in cluster2_aa:
			if aa.localscore == 0.0:
				aa.localscore = sum(d.score(self.protein, [aa], system=self.system) for d in self.distributions)
			aa.cluster = self.cluster2
			aa.clusterscore = cluster2_score

		weights = [score_to_probability(-cluster1_score), score_to_probability(-cluster2_score)]
		s = sum(weights)
		weights = [w / s for w in weights]
		random_start = np.random.choice([self.cluster1[1] - 1, self.cluster2[0]], p=weights)
		selected_cluster = self.protein.aminoacids[random_start].cluster
		return self.protein.aminoacids[selected_cluster[0] : selected_cluster[1]]

	def probabilities(self, aminoacids, anchors=[], **kwargs):
		"""This overridden function generates a list of possible conformations and their probabilities weighted by the distribution scores, as does its superclass's implementation. However, this implementation considers variants on the segment structures based on the possible conformations provided earlier.
			NOTES: Pass in an integer for keyword argument "primanchor" to signify which anchor is used as a pivot point for mutation. Pivoting will not be used unless you specify a value for "primanchor". Also, you should provide a Boolean for keyword argument "prior" to signify whether or not primanchor is before the segment. If not provided, "prior" is assumed True.
			If you pass the keyword argument 'numconfs', it will be used to determine a number of possible conformations from the previously-generated structures."""
		probabilities = []
		current_score = sum(d.score(self.protein, aminoacids, system=self.system) for d in self.distributions)
		if self.mode == psource_erratic_mode:
			proximity = self._randomization_margin(current_score / len(aminoacids)) * 4.0
			step = proximity
		else:
			proximity = 0.1
			step = 0.1
		assert "primanchor" in kwargs, "Need a primary anchor to calculate probabilities"
		if "prior" in kwargs: prior = kwargs["prior"]
		else: prior = True
		current_score = 0.0
		for d in self.distributions:
			if d.type != distributions.frequency_consec_disttype:
				current_score += d.score(self.protein, aminoacids, system=self.system)
			else:
				current_score += d.score(self.protein, aminoacids, system=self.system, prior=prior)

		#Choose some predetermined conformations to try as well
		if "numconfs" in kwargs: numconfs = kwargs["numconfs"]
		else:					 numconfs = 25
		if numconfs > 0:
			if prior:
				weights = [c[2] for c in self.c2_conformations]
				total = float(sum(weights))
				weights = [w / total for w in weights]
				seg_confs = [self.c2_conformations[c][0] for c in np.random.choice(len(self.c2_conformations), numconfs, p=weights)]
			else:
				weights = [c[2] for c in self.c1_conformations]
				total = float(sum(weights))
				weights = [w / total for w in weights]
				seg_confs = [self.c1_conformations[c][0] for c in np.random.choice(len(self.c1_conformations), numconfs, p=weights)]
		else: seg_confs = []

		def _iter_confs(seg_conf=None):
			starters = aminoacids
			if seg_conf:
				second_last = None
				if prior:
					pivot = starters[0]
					if starters[0].tag > 0:
						second_last = self.protein.aminoacids[starters[0].tag - 1]
				else:
					pivot = starters[-1]
					if starters[-1].tag < len(self.protein.aminoacids) - 1:
						second_last = self.protein.aminoacids[starters[-1].tag + 1]
				candidates = self.permissions.allowed_conformations(pivot, anchors[kwargs["primanchor"]], prior, opposite_aa=second_last)
				if not len(candidates):
					print "No permissible candidates for pivot"
					return
				zone = random.choice(candidates)
				pivot.acarbon = zone.alpha_zone
				pivot.set_axes(zone.x_axis, zone.y_axis, zone.z_axis)
				if prior:
					original_pivot = pivot.hypothetical(seg_conf[0])
				else:
					original_pivot = pivot.hypothetical(seg_conf[-1])

				for	i, pz in enumerate(seg_conf):
					hypo = starters[i].hypothetical(pivot.globalpz(original_pivot.localpz(pz)))
					starters[i] = hypo
		
			maxsamples = -1
			if seg_confs: maxsamples = 5
			for conformation in self.iter_ensemble_pivot(starters, anchors[kwargs["primanchor"]], prior, proximity * 2.0, maxsamples=maxsamples):
				hypotheticals = [starters[i].hypothetical(conformation[i]) for i in xrange(len(starters))]
				score = 0.0
				for d in self.distributions:
					if d.type != distributions.frequency_consec_disttype:
						score += d.score(self.protein, hypotheticals, system=self.system)
					else:
						score += d.score(self.protein, hypotheticals, system=self.system, prior=prior)
				if self.mode != psource_gentle_mode or score < current_score:
					probabilities.append([conformation, score_to_probability(score / len(hypotheticals))])
				del hypotheticals[:]

		_iter_confs()
		for sc in seg_confs:
			_iter_confs(sc)
		if len(probabilities) == 0:
			for conformation in self.iter_ensemble_pivot(aminoacids, anchors[kwargs["primanchor"]], prior, 1000.0):
				hypotheticals = [aminoacids[i].hypothetical(conformation[i]) for i in xrange(len(aminoacids))]
				score = 0.0
				for d in self.distributions:
					if d.type != distributions.frequency_consec_disttype:
						score += d.score(self.protein, hypotheticals, system=self.system)
					else:
						score += d.score(self.protein, hypotheticals, system=self.system, prior=prior)
				if self.mode != psource_gentle_mode or score < current_score:
					probabilities.append([conformation, score_to_probability(score / len(hypotheticals))])
				del hypotheticals[:]

		apply_weight(probabilities, partial(self.weights["distance"], aminoacids=aminoacids, proximity=self._randomization_margin(current_score / len(aminoacids)) * 4.0))
		apply_weight(probabilities, partial(self.weights["comparison"], currentprob=score_to_probability(current_score / len(aminoacids))), passprob=True)
		apply_weight(probabilities, partial(self.weights["anchor"], anchors=anchors, negweight=proximity))

		probabilities = [x for x in probabilities if x[1] > 0.0]
		if len(probabilities) == 0 and "primanchor" not in kwargs:
			probabilities.append([[PositionZone(aa.acarbon, aa.i, aa.j, aa.k) for aa in aminoacids], score_to_probability(current_score / len(aminoacids))])
		return to_cdf(sorted(probabilities, key=lambda x: x[1]))
