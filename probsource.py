"""The probsource module provides a class ProbabilitySource whose instances return cumulative distribution functions for use by the folding module. This module also includes functions for interconverting frequency, score, and probability.
	
	The provided subclass of ProbabilitySource, DistributionProbabilitySource, acts as a coordinator of DistributionManager objects. Initialize a DistributionProbabilitySource object with instances of the provided subclasses of DistributionManager, then pass it to folding_iteration() to apply frequency weighting."""

from aminoacids import *
from functools import partial
import numpy as np
from numpy import matlib
from randomcoil import random_axes
import random
import datetime
import distributions

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

#MARK: - ProbabilitySource

class ProbabilitySource(object):
	def __init__(self, protein):
		"""The default ProbabilitySource object has a rudimentary set of weights. Subclass __init__ to use your own weights. These weights should be applied last in the distribution-generating process (see probabilities)."""
		self.protein = protein
		
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
				if len(self.protein.hashtable.nearby_aa(aminoacids[i].hypothetical(pz), 2.0)) > 0:
					return 0.0
			return 1.0
		self.weights = { "distance" : distance_weight, "bond" : bond_weight, "steric" : steric_weight }
		self.steric_cutoff = 4.0
		self.steric_consec_diff = 5
	
	def iter_ensemble_alpha(self, aminoacids, proximity=1, step=1):
		"""This method is an iterator over a series of possible destinations for the aminoacids (alpha carbons only). The return value is an array of position zones. Subclasses may want to use this method to assist in implementing probabilities()."""
		if len(aminoacids) == 1:
			residue = aminoacids[0]
			'''if step <= 0.001:
				yield [PositionZone(alpha=residue.acarbon, x=residue.i, y=residue.j, z = residue.k)]
				return'''
			for alphapt in residue.acarbon.iter_randomoffsets(proximity, count=10):
				if alphapt.distanceto(residue.acarbon) > proximity: continue
				if len(self.protein.nearby_aa(residue.hypothetical(PositionZone(alphapt, residue.i, residue.j, residue.k)), self.steric_cutoff, consec=False, mindiff=self.steric_consec_diff)) > 0:
					continue
				yield [PositionZone(alpha=alphapt, x=residue.i, y=residue.j, z = residue.k)]
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
						if len(self.protein.nearby_aa(aminoacids[i].hypothetical(pz), self.steric_cutoff, consec=False, mindiff=self.steric_consec_diff)) > 0:
							valid = False
							break
					if valid != True: continue
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

	def random_vicinity(self, point, tag, distance=0.1):
		"""This helper function returns a random point within a cube of side distance * 2 centered at point."""
		point = Point3D(random.uniform(point.x - distance, point.x + distance),
						random.uniform(point.y - distance, point.y + distance),
						random.uniform(point.z - distance, point.z + distance))
		return point
	
	def random_vicinity_axes(self, conformation, distance=0.01):
		if math.fabs(distance) <= 0.0001: return (conformation.x_axis, conformation.y_axis, conformation.z_axis)
		i = self.random_vicinity(conformation.x_axis, distance).normalize()
		j = Point3D(random.uniform(conformation.y_axis.x - distance, conformation.y_axis.x + distance), random.uniform(conformation.y_axis.y - distance, conformation.y_axis.y + distance), 0.0)
		j.z = -(i.x * j.x + i.y * j.y) / i.z
		k = crossproduct(i, j)
		return (i, j.normalize(), k.normalize())

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
	def __init__(self, protein, distributions, permissions=None):
		super(AAProbabilitySource,self).__init__(protein)
		self.distributions = distributions
		self.permissions = permissions
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
		self.erratic_proximity = 0.6
	
	def iter_ensemble_pivot(self, aminoacids, anchor, prior=True, proximity=1.0):
		"""This new method (2/6/15) iterates the possible conformations created by pivoting a segment around the anchor (should be an AminoAcid). This function does not take into account the connectivity with the next unmutated segment."""
		assert self.permissions is not None, "Cannot perform a pivot without a PermissionsManager."
		if prior is True:	aa = aminoacids[0]
		else:				aa = aminoacids[-1]
		allowed_conformations = self.permissions.allowed_conformations(aa, anchor, prior)
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

			valid = True
			for i, pz in enumerate(conformation):
				if len(self.protein.nearby_aa(aminoacids[i].hypothetical(pz), self.steric_cutoff, consec=False, mindiff=self.steric_consec_diff, excluded=[aa.tag for aa in aminoacids])) > 0:
					valid = False
					break
			if valid is not True: continue
			yield conformation

	def probabilities(self, aminoacids, anchors=[], **kwargs):
		"""This overridden function generates a list of possible conformations and their probabilities, as does its superclass's implementation. However, this implementation consults its list of DistributionManagers, polling their score() method for each amino acid in the segment. These scores are summed for each configuration, and the resulting score is transformed into a probability.
			NOTES: Pass in an integer for keyword argument "primanchor" to signify which anchor is used as a pivot point for mutation. Pivoting will not be used unless you specify a value for "primanchor". Also, you should provide a Boolean for keyword argument "prior" to signify whether or not primanchor is before the segment. If not provided, "prior" is assumed True."""
		probabilities = []
		current_score = sum(d.score(self.protein, aminoacids) for d in self.distributions)
		if self.mode == psource_erratic_mode:
			proximity = self._randomization_margin(current_score / len(aminoacids)) * 4.0 #sum(d.score(self.protein, self.protein.aminoacids) for d in self.distributions) / len(self.protein.aminoacids) #self.erratic_proximity
			step = proximity #max(proximity / 2.0, min(0.1, proximity))
		else:
			proximity = 0.1
			step = 0.1
		#print proximity, "Now:", repr(aminoacids), score_to_probability(current_score)
		#print len(anchors), aminoacids
		if "primanchor" in kwargs:
			if "prior" in kwargs: prior = kwargs["prior"]
			else: prior = True
			current_score = 0.0
			for d in self.distributions:
				if d.type != distributions.frequency_consec_disttype:
					current_score += d.score(self.protein, aminoacids)
				else:
					current_score += d.score(self.protein, aminoacids, prior=prior)
			for conformation in self.iter_ensemble_pivot(aminoacids, anchors[kwargs["primanchor"]], prior, proximity * 2.0):
				hypotheticals = [aminoacids[i].hypothetical(conformation[i]) for i in xrange(len(aminoacids))]
				score = 0.0
				for d in self.distributions:
					if d.type != distributions.frequency_consec_disttype:
						score += d.score(self.protein, hypotheticals)
					else:
						score += d.score(self.protein, hypotheticals, prior=prior)
				if self.mode != psource_gentle_mode or score < current_score:
					probabilities.append([conformation, score_to_probability(score / len(hypotheticals))])
				del hypotheticals[:]
			if len(probabilities) == 0:
				for conformation in self.iter_ensemble_pivot(aminoacids, anchors[kwargs["primanchor"]], prior, 1000.0):
					hypotheticals = [aminoacids[i].hypothetical(conformation[i]) for i in xrange(len(aminoacids))]
					score = 0.0
					for d in self.distributions:
						if d.type != distributions.frequency_consec_disttype:
							score += d.score(self.protein, hypotheticals)
						else:
							score += d.score(self.protein, hypotheticals, prior=prior)
					if self.mode != psource_gentle_mode or score < current_score:
						probabilities.append([conformation, score_to_probability(score / len(hypotheticals))])
					del hypotheticals[:]
		else:
			current_score = sum(d.score(self.protein, aminoacids) for d in self.distributions if d.type != distributions.frequency_consec_disttype)
			for conformation in self.iter_ensemble_alpha(aminoacids, proximity, step=step):
				hypotheticals = [aminoacids[i].hypothetical(conformation[i]) for i in xrange(len(aminoacids))]
				score = sum(d.score(self.protein, hypotheticals) for d in self.distributions if d.type != distributions.frequency_consec_disttype)
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
		return to_cdf(sorted(probabilities, key=lambda x: x[1]))
			
	def _iter_permissible_randomcoils(self, aminoacids, start, end, chain=[], yieldct=0):
		"""This helper function is a recursive generator for permissible self-avoiding walks."""
		if len(chain) == len(aminoacids) - 1:
			#Connect the last amino acid in a permissible orientation to both chain[-1] and end.
			if start.tag < end.tag:
				if len(chain) > 0:
					last = aminoacids[len(chain) - 1].hypothetical(chain[-1])
				else:
					last = start
				confs = self.permissions.allowed_conformations(aminoacids[-1], last)
				random.shuffle(confs)
				confs2 = self.permissions.allowed_conformations(aminoacids[-1], end, prior=False)
				random.shuffle(confs2)
				for pz in confs:
					for pz2 in confs2:
						if pz2.alpha_zone.distanceto(pz.alpha_zone) > 0.5: continue
						pz = PositionZone(pz.alpha_zone.add(pz2.alpha_zone).multiply(0.5),
										  pz.x_axis, pz.y_axis, pz.z_axis)
						hyp = aminoacids[-1].hypothetical(pz)
						if self.permissions.is_valid(hyp, end, prior=False) and self.permissions.is_valid(hyp, last):
							if len(self.protein.nearby_aa(hyp, self.steric_cutoff, consec=False, mindiff=self.steric_consec_diff, excluded=[aa.tag for aa in aminoacids])) > 0:
								continue
							yieldct += 1
							yield chain + [pz]
							if yieldct == 10: return
			else:
				for saw in self._iter_permissible_randomcoils(aminoacids, end, start, chain, yieldct):
					yieldct += 1
					yield saw
		else:
			#Try to add one more amino acid in a permissible orientation.
			if start.tag < end.tag:
				if len(chain) > 0:
					last = aminoacids[len(chain) - 1].hypothetical(chain[-1])
				else:
					last = start
				confs = self.permissions.allowed_conformations(aminoacids[len(chain)], last)
				#I'm using a constant number here to save on performance. Only 10 of the given conformations will be chosen, which still gives a total of 10^(n - 1) iterations before the last residue.
				if len(confs) > 30:
					confs = random.sample(confs, 30)
				random.shuffle(confs)
				iters = 0
				for pz in confs:
					hyp = aminoacids[len(chain)].hypothetical(pz, True)
					if not self.connectivity_possible(hyp, end, len(aminoacids) - len(chain) - 1): continue
					if len(self.protein.nearby_aa(hyp, self.steric_cutoff, consec=False, mindiff=self.steric_consec_diff, excluded=[aa.tag for aa in aminoacids])) > 0:
						continue
					iters += 1
					for saw in self._iter_permissible_randomcoils(aminoacids, start, end, chain + [pz], yieldct):
						yieldct += 1
						yield saw
						if yieldct == 10: return
				#print len(chain), "did", iters
			else:
				for saw in self._iter_permissible_randomcoils(aminoacids, end, start, chain, yieldct):
					yieldct += 1
					yield saw

	def randomcoil_probabilities(self, aminoacids, start, end):
		"""Pass in an array of amino acids, a start anchor amino acid (not included in aminoacids), and an end anchor amino acid. Each entry in the returned probability distribution will be a permissible SAW from start to end. If no permissible SAW is found that connects start and end, this method will return no probabilities, in which case you should use probabilities() instead."""
		probabilities = []
		current_score = sum(d.score(self.protein, aminoacids) for d in self.distributions)

		for conformation in self._iter_permissible_randomcoils(aminoacids, start, end):
			hypotheticals = [aminoacids[i].hypothetical(conformation[i]) for i in xrange(len(aminoacids))]
			probabilities.append([conformation, score_to_probability(sum(d.score(self.protein, hypotheticals) for d in self.distributions))])
			del hypotheticals[:]

		apply_weight(probabilities, partial(self.weights["comparison"], currentprob=score_to_probability(current_score)), passprob=True)
		probabilities = [x for x in probabilities if x[1] > 0.0]
		if len(probabilities) == 0:
			return []
		return to_cdf(sorted(probabilities, key=lambda x: x[1]))

	
	def angleprobabilities(self, aminoacid, **kwargs):
		"""This overridden function generates a certain number of random axis configurations for aminoacid, then computes the probability of each one occurring using the DistributionManager score() method."""
		probabilities = []
		current_score = sum(d.score(self.protein, [aminoacid]) for d in self.distributions)
		proximity = max(self._randomization_margin(current_score), 0.001)
		step = proximity
		for conformation in self.iter_ensemble_axis(aminoacid, proximity=proximity, step=step):
			hypothetical = aminoacid.hypothetical(conformation[0])
			probabilities.append([conformation, score_to_probability(sum(d.score(self.protein, [hypothetical]) for d in self.distributions))])
		
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

	def is_connected(self, segment):
		"""This function detects whether the amino acids at the beginning and end of segment are connected to their neighbors outside the segment."""
		if segment[0].tag == 0 and segment[-1].tag == len(self.protein.aminoacids) - 1:
			return True
		elif segment[0].tag == 0:
			if self.permissions.is_valid(self.protein.aminoacids[segment[-1].tag + 1], segment[-1]):
				return True
			else:
				return False
		elif segment[-1].tag == len(self.protein.aminoacids) - 1:
			if self.permissions.is_valid(segment[0], self.protein.aminoacids[segment[0].tag - 1]):
				return True
			else:
				return False
		else:
			pre = self.protein.aminoacids[segment[0].tag - 1]
			post = self.protein.aminoacids[segment[-1].tag + 1]
			#print pre, segment[0], "...", segment[-1], post
			if self.permissions.is_valid(segment[0], pre, prior=True) and self.permissions.is_valid(segment[-1], post, prior=False):
				return True
			else:
				return False

	def choose_segment(self, segment_length):
		#First we have to cluster out the chain into groups that have similarly stable consecutive scores.
		clusters = [[0, 0, 0]]
		ranges = []
		cluster_range = [0.0, 0.0]
		min_bound = 0.9
		max_bound = 1.1
		for i, aa in enumerate(self.protein.aminoacids):
			#assert aa.localscore != 0.0, "The amino acid {} was mutated and the local score not updated before calling choose_segment.".format(aa)
			if aa.localscore == 0.0:
				aa.localscore = sum(d.score(self.protein, [aa]) for d in self.distributions)
			if clusters[-1][1] == 0:
				clusters[-1][1] += 1
				clusters[-1][2] += aa.localscore
				cluster_range[0] = min(aa.localscore * min_bound, aa.localscore * max_bound)
				cluster_range[1] = max(aa.localscore * min_bound, aa.localscore * max_bound)
				ranges.append([cluster_range[0], cluster_range[1]])
			elif aa.localscore >= cluster_range[0] and aa.localscore <= cluster_range[1] and aa.localscore <= -125.0:
				clusters[-1][1] += 1
				clusters[-1][2] += aa.localscore
			else:
				clusters[-1][2] /= float(clusters[-1][1] - clusters[-1][0])
				clusters.append([i, i + 1, aa.localscore])
				cluster_range[0] = min(aa.localscore * min_bound, aa.localscore * max_bound)
				cluster_range[1] = max(aa.localscore * min_bound, aa.localscore * max_bound)
				ranges.append([cluster_range[0], cluster_range[1]])
		print "---", len(clusters), "clusters for", len(self.protein.aminoacids), "amino acids"
		for clust in clusters:
			tupcluster = (clust[0], clust[1])
			for i in xrange(*tupcluster):
				self.protein.aminoacids[i].cluster = tupcluster
				self.protein.aminoacids[i].clusterscore = clust[2]
		self.clusters = clusters
		weights = [score_to_probability(-1 * sum(d.score(self.protein, self.protein.aminoacids[i:i + segment_length]) for d in self.distributions)) for i in xrange(len(self.protein.aminoacids) - 1)]
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
		score = sum(d.score(self.protein, hypotheticals) for d in self.distributions)
		return self._randomization_margin(score)

	def bond_breach(self, aa1, aa2):
		return aa1.acarbon.distanceto(aa2.acarbon) - self.steric_cutoff

	def cdf_from_hypotheticals(self, cases):
		probabilities = []
		for hypotheticals in cases:
			probabilities.append([hypotheticals, score_to_probability(sum(d.score(self.protein, hypotheticals) for d in self.distributions))])
		probabilities = [x for x in probabilities if x[1] > 0.0]
		return to_cdf(sorted(probabilities, key=lambda x: x[1]))
