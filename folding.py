from distributions import *
from molecular_systems import *
from probsource import *
import random
from main import load_dists
import os, sys
import numpy
from multiprocessing import Process, Queue
import datetime, time

#MARK: Helpers

def mutate_aa_orientation(protein, psource, residue):
	#Find angle probabilities.
	probabilities = psource.angleprobabilities(residue)
	
	#Sample the cumulative distribution function and execute the change.
	selected_conformation = sample_cdf(probabilities)
	if psource.distributions[0].score(protein, [residue.hypothetical(selected_conformation[0])]) < 0:
		print "Favorable conformation."
	residue.set_axes(random_vicinity_axes(selected_conformation[0], distance=psource.randomization_margin(selected_conformation, [residue])), normalized=True)
	
	#Now radiate outward from the selected amino acid and adjust the neighbors' positions.
	delta = 1000000.0
	preresidue = postresidue = residue
	while delta > 1.0 and (preresidue.tag > 0 or postresidue.tag < len(protein.aminoacids) - 1):
		if preresidue.tag > 0:
			preresidue = protein.aminoacids[preresidue.tag - 1]
			
			#Find angle probabilities
			probabilities = psource.angleprobabilities(preresidue)
			
			#Sample the cumulative distribution function and execute the change.
			selected_conformation = sample_cdf(probabilities)
			preresidue.set_axes(random_vicinity_axes(selected_conformation[0], distance=psource.randomization_margin(selected_conformation, [preresidue])), normalized=True)
		#preresidue.set_axes(selected_conformation[0].x_axis, selected_conformation[0].y_axis, selected_conformation[0].z_axis)
		
		if postresidue.tag < len(protein.aminoacids) - 1:
			postresidue = protein.aminoacids[postresidue.tag + 1]
			
			#Find angle probabilities
			probabilities = psource.angleprobabilities(postresidue)
			
			#Sample the cumulative distribution function and execute the change.
			selected_conformation = sample_cdf(probabilities)
			postresidue.set_axes(random_vicinity_axes(selected_conformation[0], distance=psource.randomization_margin(selected_conformation, [postresidue])), normalized=True)
			#postresidue.set_axes(selected_conformation[0].x_axis, selected_conformation[0].y_axis, selected_conformation[0].z_axis)

def mutate_aa_pos_orientation(protein, psource, residue):
	#Get the probability list for the residue. It should be sorted by probability.
	probabilities = psource.probabilities([residue])
	
	#Sample the cumulative distribution function and execute the change.
	selected_conformation = sample_cdf(probabilities)
	residue.acarbon = selected_conformation[0].alpha_zone.random_vicinity(distance=psource.randomization_margin(selected_conformation, [residue]))
	
	#Find angle probabilities
	probabilities = psource.angleprobabilities(residue)
	
	#Sample the cumulative distribution function and execute the change.
	selected_conformation = sample_cdf(probabilities)
	residue.set_axes(random_vicinity_axes(selected_conformation[0], distance=psource.randomization_margin(selected_conformation, [residue])), normalized=True)
	#residue.set_axes(selected_conformation[0].x_axis, selected_conformation[0].y_axis, selected_conformation[0].z_axis)
	
	#Now radiate outward from the selected amino acid and adjust the neighbors' positions.
	delta = 1000000.0
	preresidue = postresidue = residue
	while delta > 1.0 and (preresidue.tag > 0 or postresidue.tag < len(protein.aminoacids) - 1):
		if preresidue.tag > 0:
			preresidue = protein.aminoacids[preresidue.tag - 1]
			#Get the probability list for the residue.
			anchors = [AAAnchor.make(protein.aminoacids[preresidue.tag + 1])]
			if preresidue.tag > 0:
				anchors.append(AAAnchor.make(protein.aminoacids[preresidue.tag - 1], weight=0.5))
			probabilities = psource.probabilities([preresidue], anchors=anchors)
			
			#Sample the cumulative distribution function and execute the change.
			selected_conformation = sample_cdf(probabilities)
			#for anchor in anchors:
			#	d = min(x.alpha_zone.distanceto(anchor.acarbon) for x in selected_conformation)
			#	print "Preselected:", d
			new_loc = selected_conformation[0].alpha_zone.random_vicinity(distance=psource.randomization_margin(selected_conformation, [preresidue]))
			delta = min(delta, preresidue.acarbon.distanceto(new_loc))
			preresidue.acarbon = new_loc
			
			#Find angle probabilities
			probabilities = psource.angleprobabilities(preresidue)
			
			#Sample the cumulative distribution function and execute the change.
			selected_conformation = sample_cdf(probabilities)
			preresidue.set_axes(random_vicinity_axes(selected_conformation[0], distance=psource.randomization_margin(selected_conformation, [preresidue])), normalized=True)
			#preresidue.set_axes(selected_conformation[0].x_axis, selected_conformation[0].y_axis, selected_conformation[0].z_axis)
		
		if postresidue.tag < len(protein.aminoacids) - 1:
			postresidue = protein.aminoacids[postresidue.tag + 1]
			#Get the probability list for the residue.
			anchors = [AAAnchor.make(protein.aminoacids[postresidue.tag - 1])]
			if postresidue.tag + 1 < len(protein.aminoacids):
				anchors.append(AAAnchor.make(protein.aminoacids[postresidue.tag + 1], weight=0.5))
			probabilities = psource.probabilities([postresidue], anchors=anchors)
			
			#Sample the cumulative distribution function and execute the change.
			selected_conformation = sample_cdf(probabilities)
			#for anchor in anchors:
			#	d = min(x.alpha_zone.distanceto(anchor.acarbon) for x in selected_conformation)
			#	print "Postselected:", d
			new_loc = selected_conformation[0].alpha_zone.random_vicinity(distance=psource.randomization_margin(selected_conformation, [postresidue]))
			delta = min(delta, postresidue.acarbon.distanceto(new_loc))
			postresidue.acarbon = new_loc
			
			#Find angle probabilities
			probabilities = psource.angleprobabilities(postresidue)
			
			#Sample the cumulative distribution function and execute the change.
			selected_conformation = sample_cdf(probabilities)
			postresidue.set_axes(random_vicinity_axes(selected_conformation[0], distance=psource.randomization_margin(selected_conformation, [postresidue])), normalized=True)
			#postresidue.set_axes(selected_conformation[0].x_axis, selected_conformation[0].y_axis, selected_conformation[0].z_axis)

def _apply_conformation_recursive(protein, psource, segment_length, segment, selected_conformation, restore=False, wave=0, file=None, cutoff_wave=5.0):
	"""wave: 0=no wave, 1=wave left, 2=wave right, 3=wave both.
		Returns: if no probabilities are found for the wave, returns None. If the conformation has been saved, returns a tuple (min, max) where min is the lowest tag mutated and max is the lowest tag not mutated (e.g., if amino acids 1-3 were mutated, min=1 and max=4."""
	begin_offset = end_offset = Point3D.zero()
	
	for i in xrange(len(segment)):
		#if restore == True: residue.save()
		residue = protein.aminoacids[i + segment[0].tag]
		old = residue.acarbon
		residue.acarbon = selected_conformation[i].alpha_zone
		if i == 0:
			begin_offset = residue.acarbon.subtract(old)
		elif i == len(segment) - 1:
			end_offset = residue.acarbon.subtract(old)
		residue.set_axes(selected_conformation[i].x_axis,
						 selected_conformation[i].y_axis,
						 selected_conformation[i].z_axis)

	#print protein.xyz(escaped=False, highlight=range(segment[0].tag, segment[-1].tag + 1))
	mutation_list = [segment[0].tag, segment[-1].tag + 1]
	if wave == 0: return mutation_list

	#Now radiate outward from the selected segment and adjust the neighbor segments' positions.
	segment_length = 1

	if segment[0].tag > 0 and not psource.is_connected(segment) and wave != 2:
		cluster = protein.aminoacids[segment[0].tag - 1].cluster
		if cluster[1] - cluster[0] <= segment_length or random.uniform(-160.0, 0.0) <= protein.aminoacids[segment[0].tag - 1].clusterscore:
			presegment = protein.aminoacids[max(segment[0].tag - segment_length, 0) : segment[0].tag]
		else:
			presegment = protein.aminoacids[cluster[0] : segment[0].tag]
		anchors = [AAAnchor.make(protein.aminoacids[presegment[-1].tag + 1], weight=4, hook=-1)]
		if presegment[0].tag > 0:
			anchors.append(AAAnchor.make(protein.aminoacids[presegment[0].tag - 1], weight=8, hook=0))
		#Translate all the amino acids to put them closer to the recently changed segment.
		begin_offset = protein.aminoacids[presegment[-1].tag + 1].acarbon.subtract(presegment[-1].acarbon)
		begin_offset = begin_offset.multiply((begin_offset.magnitude() - random.uniform(2.5, 3.5)) / begin_offset.magnitude())
		for aa in presegment:
			if restore == True: aa.save()
			aa.acarbon = aa.acarbon.add(begin_offset)
		#Get the probability list for the segment. It should be sorted by probability.
		probabilities = []
		if len(anchors) >= 2 and len(presegment) == 1:
			if psource.connectivity_possible(anchors[0], anchors[1], 1):
				probabilities = psource.randomcoil_probabilities(presegment, anchors[1], anchors[0])
		if len(probabilities) == 0:
			probabilities = psource.probabilities(presegment, anchors=anchors, primanchor=0, prior=False, connected=(wave == 1))
		if len(probabilities) == 0:
			if restore == True:
				for aa in presegment:
					aa.restore()
			return None

		#Sample the cumulative distribution function and execute the change.
		application_ret = None
		selected_conformation = None
		while application_ret is None and len(probabilities) > 0:
			selected_conformation = sample_cdf(probabilities)
			application_ret = _apply_conformation_recursive(protein, psource, segment_length, presegment, selected_conformation, restore, wave=1, file=file)
			if application_ret is None:
				entry = next(i for i in probabilities if i[0] == selected_conformation)
				probabilities.remove(entry)
		if application_ret is None:
			if restore == True:
				for aa in presegment:
					aa.restore()
			return None
		else:
			mutation_list[0] = application_ret[0]
		assert psource.permissions.is_valid(presegment[-1], anchors[0], prior=False), "Invalid orientation for presegment: {} and {} ({})".format(presegment, anchors, anchors[0].tolocal(presegment[0].acarbon))
		#print protein.xyz(escaped=False, highlight=range(presegment[0].tag, presegment[-1].tag + 1))
	
	if segment[-1].tag < len(protein.aminoacids) - 1 and not psource.is_connected(segment) and wave != 1:
		cluster = protein.aminoacids[segment[-1].tag + 1].cluster
		if cluster[1] - cluster[0] <= segment_length or random.uniform(-160.0, 0.0) <= protein.aminoacids[segment[-1].tag + 1].clusterscore:
			postsegment = protein.aminoacids[segment[-1].tag + 1 : min(segment[-1].tag + 1 + segment_length, len(protein.aminoacids))]
		else:
			postsegment = protein.aminoacids[segment[-1].tag + 1 : cluster[1]]
		anchors = [AAAnchor.make(protein.aminoacids[postsegment[0].tag - 1], weight=4, hook=0)]
		if postsegment[-1].tag + 1 < len(protein.aminoacids):
			anchors.append(AAAnchor.make(protein.aminoacids[postsegment[-1].tag + 1], weight=8, hook=-1))
		#Translate all the amino acids to put them closer to the recently changed segment.
		end_offset = protein.aminoacids[postsegment[0].tag - 1].acarbon.subtract(postsegment[0].acarbon)
		end_offset = end_offset.multiply((end_offset.magnitude() - random.uniform(2.5, 3.5)) / end_offset.magnitude())
		for aa in postsegment:
			if restore == True: aa.save()
			aa.acarbon = aa.acarbon.add(end_offset)
		#Get the probability list for the segment. It should be sorted by probability.
		probabilities = []
		if len(anchors) >= 2 and len(postsegment) == 1:
			if psource.connectivity_possible(anchors[0], anchors[1], 1):
				probabilities = psource.randomcoil_probabilities(postsegment, anchors[0], anchors[1])
		if len(probabilities) == 0:
			probabilities = psource.probabilities(postsegment, anchors=anchors, primanchor=0, prior=True, connected=(wave == 2))
		if len(probabilities) == 0:
			if restore == True:
				for aa in postsegment:
					aa.restore()
				for aa in protein.aminoacids[mutation_list[0] : segment[0].tag]:
					aa.restore()
			return None

		#Sample the cumulative distribution function and execute the change.
		application_ret = None
		selected_conformation = None
		while application_ret is None and len(probabilities) > 0:
			selected_conformation = sample_cdf(probabilities)
			application_ret = _apply_conformation_recursive(protein, psource, segment_length, postsegment, selected_conformation, restore, wave=2, file=file)
			if application_ret is None:
				entry = next(i for i in probabilities if i[0] == selected_conformation)
				probabilities.remove(entry)
		if application_ret is None:
			if restore == True:
				for aa in postsegment:
					aa.restore()
				for aa in protein.aminoacids[mutation_list[0] : segment[0].tag]:
					aa.restore()
			return None
		else:
			mutation_list[1] = application_ret[1]
		assert psource.permissions.is_valid(postsegment[0], anchors[0]), "Invalid orientation for postsegment: {} and {} ({})".format(postsegment, anchors, anchors[0].tolocal(postsegment[0].acarbon))
	
	if not psource.is_connected(protein.aminoacids[mutation_list[0] : mutation_list[-1]]):
		return None

	if file is not None:
		highlight = []
		for clust in psource.clusters:
			if clust[1] - clust[0] > 1: highlight += range(clust[0], clust[1])
		file.write(protein.xyz(escaped=False, highlight=highlight)) #, highlight=range(mutation_list[0], mutation_list[1])))
	return mutation_list

def apply_conformation(protein, psource, segment_length, segment, selected_conformation, restore=False, file=None):
	"""If no valid permissible conformation is found using the selected_conformation, this function returns None.
		If restore=False, this function returns True and simply adjusts the location of the amino acids (including a wave of adjustments).
		If restore=True, this function returns an array of hypothetical amino acids representing the final locations of the changed residues (including the wave of adjustments). All the amino acids will be restored back to their original positions at the end."""
	if restore == True:
		for aa in segment:
			aa.save()
	application_ret = _apply_conformation_recursive(protein, psource, segment_length, segment, selected_conformation, restore, wave=3, file=file)
	if application_ret is None:
		if restore == True:
			for aa in segment:
				aa.restore()
		return None
	else:
		if restore == True:
			hypotheticals = []
			for aa in protein.aminoacids[application_ret[0] : application_ret[1]]:
				hypotheticals.append(aa.hypothetical(aa.restore(), True))
			return hypotheticals
		else:
			return True

#MARK: - Iterations

def folding_iteration(system, psources, segment_length=1, file=None):
	"""The backbone of the protein folding simulator. The gist of the algorithm is to choose a random subset of the amino acids in 'protein' (which is a Polypeptide) and move them randomly. The function returns no value, just updates the protein.
		
		folding_iteration uses the ProbabilitySource object's probabilities() and angleprobabilities() methods to determine where to move the segment. psource should accept a list of relevant amino acids (the residues for which probabilities are required) and return a list of lists (form: [[conformation, probability], ...]) representing a cumulative distribution function and sorted by probability. (The conformation should be specified as a list of position zones, one for each amino acid.)
		
		The segment_length parameter can be used to manipulate the scale/granularity of the mutation. If 0 is passed, only the orientation of a single amino acid will be changed. For 1, the position and orientation will be changed. For any higher integer, a segment of that length will be pivoted/rotated/translated.
		
		Note that the final transformation is random and not guaranteed to match one of the probabilities exactly (no lattice is involved). Psource provides a randomization_margin function that uses score to determine how much randomization to introduce.
		
		Specify an open file object to write intermediate output to it."""

	#Choose a psource at random to use in this iteration
	psource = random.choice(psources)
	protein = psource.protein
	
	#Save a copy of the chain just in case there's a steric violation.
	for aa in protein.aminoacids:
		aa.save()

	if segment_length == 0:
		mutate_aa_orientation(protein, psource, random.choice(protein.aminoacids))
		reset_stats()
	else:
		#Choose a random segment out of the protein.
		segment = psource.choose_segment(segment_length)
		segment_length = len(segment)
		#print psource.is_connected(segment)

		#Get the probability list for the segment. It should be sorted by probability.
		probabilities = psource.probabilities(segment)
		
		#Sample the cumulative distribution function once and execute the change.
		selected_conformation = None
		application_ret = None
		while application_ret is None and len(probabilities) > 0:
			selected_conformation = sample_cdf(probabilities)
			application_ret = apply_conformation(protein, psource, segment_length, segment, selected_conformation, file=file)
			if application_ret is None:
				entry = next(i for i in probabilities if i[0] == selected_conformation)
				probabilities.remove(entry)
		#assert application_ret is not None, "No valid permissible conformation found from this location."
		reset_stats()
		if application_ret is None:
			return

	violated = False
	for aa in protein.aminoacids:
		if system.check_steric_clash(aa, protein, steric_cutoff=psource.steric_cutoff, consec=False, mindiff=psource.steric_consec_diff):
			violated = True
			break

	if violated is True:
		print "Steric violation."
		for aa in protein.aminoacids:
			aa.restore()
	else:
		aa.clear_save()

def constructive_folding_iteration(protein, psources, system=None, file=None):
	"""The backbone of the protein folding simulator, but for segmented, constructive simulations. The gist of the algorithm is to choose a random subset of the amino acids in 'protein' (which is a Polypeptide) and move them randomly. The function returns no value, just updates the protein.
		
		folding_iteration uses the ProbabilitySource object's probabilities() method to determine where to move the segment. psource should accept a list of relevant amino acids (the residues for which probabilities are required) and return a list of lists (form: [[conformation, probability], ...]) representing a cumulative distribution function and sorted by probability. (The conformation should be specified as a list of position zones, one for each amino acid.)
		
		Note that the final transformation is random and not guaranteed to match one of the probabilities exactly (no lattice is involved). Psource provides a randomization_margin function that uses score to determine how much randomization to introduce.
		
		Specify an open file object to write intermediate output to it."""
	
	#Save a copy of the chain just in case there's a steric violation.
	for aa in protein.aminoacids:
		aa.save()
	
	#Choose a random segment out of the protein.
	segment = psource.choose_segment(0)
	#print psource.is_connected(segment)

	#Get the probability list for the segment. It should be sorted by probability.
	#For the anchor, use the amino acid immediately before or after the segment.
	prior = True
	if segment[0].tag == 0:
		anchor = AAAnchor.make(protein.aminoacids[segment[-1].tag + 1], weight=4, hook=-1)
		prior = False
	elif segment[-1].tag == len(protein.aminoacids) - 1:
		anchor = AAAnchor.make(protein.aminoacids[segment[0].tag - 1], weight=4, hook=0)
	else:
		assert False, "Folding doesn't understand this segment: %r" % segment
	probabilities = psource.probabilities(segment, anchors=[anchor], primanchor=0, prior=prior, numconfs=0)

	#Sample the cumulative distribution function once and execute the change.
	selected_conformation = None
	application_ret = None
	while application_ret is None and len(probabilities) > 0:
		selected_conformation = sample_cdf(probabilities)
		application_ret = apply_conformation(protein, psource, 1, segment, selected_conformation, file=file)
		if application_ret is None:
			entry = next(i for i in probabilities if i[0] == selected_conformation)
			probabilities.remove(entry)
	#assert application_ret is not None, "No valid permissible conformation found from this location."
	reset_stats()
	if application_ret is None:
		for aa in protein.aminoacids:
			aa.restore()
		return

	violated = False
	for aa in protein.aminoacids:
		if (system and system.check_steric_clash(aa, protein, steric_cutoff=psource.steric_cutoff, consec=False, mindiff=psource.steric_consec_diff)) or (not system and len(protein.nearby_aa(aa, psource.steric_cutoff, consec=False, mindiff=psource.steric_consec_diff)) > 0):
			violated = True
			break
	if violated is True:
		print "Steric violation."
		for aa in protein.aminoacids:
			aa.restore()
	else:
		aa.clear_save()

#MARK: - Simulators

def segment_length(avg_score):
	"""This function needs a home."""
	#v = max(5 - math.exp(4.0 / 225.0 * avg_score + 1.6), 1) This decreases as score becomes less stable
	v = max(math.exp(4.0 / 225.0 * avg_score + 2), 1) #This increases as score becomes less stable
	weights = [(-1.0 / 9.0 * (x - v) ** 2 + 4.0) for x in xrange(1, 6)]
	'''v = max(5.0 / math.cosh(0.114622 * (avg_score + 60.0)), 1.0)
		weights = [(-((((x - v) / 6.0) ** 2) ** (1.0 / 7.0)) + 1.0) for x in xrange(1, 6)]'''
	s = sum(weights)
	weights = [w / s for w in weights]
	return numpy.random.choice(range(1, 4))#, p=weights)

def test_segment_combo(q, dists, seg_prob, seq1, seq2, conf1, conf2):
	aas, hashtable = seg_prob.generate_structure_from_segments(seq1 + seq2, conf1, conf2)
	test_no = 0
	while next((aa for aa in aas if len(hashtable.nearby_aa(aa, seg_prob.steric_cutoff, consec=False))), None):
		if test_no == 25:
			test_no = 100
			break
		aas, hashtable = seg_prob.generate_structure_from_segments(seq1 + seq2, conf1, conf2)
		test_no += 1
	if test_no > 25:
		q.put(0)
	protein = Polypeptide(aas)
	q.put(sum(dist.score(protein, protein.aminoacids) for dist in dists))

def segment_fold(dists, seq, range1, range2, infiles, output, sec_structs=None, outname="simulation_test.pdb", cluster_confs=25, sims=75, candidates=20):
	cluster_confs = int(cluster_confs)
	sims = int(sims)
	candidates = int(candidates)
	
	permissions = AAPermissionsManager(os.path.join(sparc_dir, "permissions"), os.path.join(sparc_dir, "permissible_sequences", "all.txt"))
	peptide = Polypeptide()
	system = MolecularSystem([peptide])
	seq1 = seq[range1[0] - 1 : range1[1]]
	seq2 = seq[range2[0] - 1 : range2[1]]
	seg_prob = AAConstructiveProbabilitySource(peptide, (0, len(seq1)), (len(seq1), len(seq1) + len(seq2)), dists, permissions, system=system)
	for i, inf in enumerate(infiles):
		seg_prob.load_cluster_conformations(i + 1, inf, n=cluster_confs)
	
	#First, test all possible combos of the segments with a random linking orientation
	i = 0
	j = 0
	scores = []
	print "Preliminary conformation testing..."
	
	queue = Queue()
	for i in xrange(len(seg_prob.c1_conformations)):
		print "Testing c1", i
		for j in xrange(len(seg_prob.c2_conformations)):
			p = Process(target=test_segment_combo, args=(queue, dists, seg_prob, seq1, seq2, seg_prob.c1_conformations[i][0], seg_prob.c2_conformations[j][0]))
			p.start()
			p.join() # this blocks until the process terminates
			result = queue.get()
			if result != 0:
				scores.append([i, j, result])

	scores = sorted(scores, key=lambda x: x[2])
	scores = scores[:min(len(scores), candidates)]
	print "The range of", len(scores), "scores is", scores[0][2], "to", scores[-1][2]
	gc.collect()

	file = open(os.path.join(output, outname), 'w')
	scoresfile = open(os.path.join(output, outname[:-4] + "-scores.txt"), 'w')
	pdb_model_idx = 1
	prob = AAProbabilitySource(peptide, dists, permissions) #sec_struct_permissions
	prob.mode = psource_gentle_mode
	model_count = 10
	best_models = [[] for i in xrange(model_count)]
	best_scores = [1000 for i in xrange(model_count)]

	a = datetime.datetime.now()
	for i, j, score in scores:
		print "Testing combination {}-{} ({})...".format(i, j, score)
		aas, hashtable = seg_prob.generate_structure_from_segments(seq1 + seq2, seg_prob.c1_conformations[i][0], seg_prob.c2_conformations[j][0])
		peptide.add_aas(aas)
		if sec_structs:
			if ',' in sec_structs:
				peptide.add_secondary_structures(sec_structs, format='csv', range=(range1[0], range2[1]))
			else:
				peptide.add_secondary_structures(sec_structs, format='pdb', range=(range1[0], range2[1]))
		print sec_structs
		system.center()

		scores = []
		t_scores = ""
		curscore = 0.0
		for dist in dists:
			sc = dist.score(peptide, peptide.aminoacids)
			t_scores += str(sc) + " (" + dist.identifier + ") "
			curscore += sc
		scores.append(curscore)
		scoresfile.write(str(pdb_model_idx) + " " + t_scores + "\n")
		file.write(system.pdb(modelno=pdb_model_idx))
		pdb_model_idx += 1
		for n in xrange(sims):
			#seglen = segment_length(scores[-1] / len(peptide.aminoacids))
			constructive_folding_iteration(peptide, seg_prob, system=system)
			system.center()
			curscore = 0.0
			t_scores = ""
			for dist in dists:
				sc = dist.score(peptide, peptide.aminoacids, system=system)
				t_scores += str(sc) + " (" + dist.identifier + ") "
				curscore += sc
			scores.append(curscore)
			scoresfile.write(str(pdb_model_idx) + " " + t_scores + "\n")
			file.write(system.pdb(modelno=pdb_model_idx))
			pdb_model_idx += 1
			'''#Save the conformation if it is the best so far.
			for k in xrange(model_count):
				if scores[-1] < best_scores[k] * 1.2:
					for m in reversed(xrange(max(k, 1), model_count)):
						best_scores[m] = best_scores[m - 1]
						best_models[m] = best_models[m - 1]
					best_scores[k] = scores[-1]
					best_models[k] = [PositionZone(aa.acarbon, aa.i, aa.j, aa.k) for aa in peptide.aminoacids]
					break
				elif scores[-1] <= best_scores[k]:
					break'''

		del peptide.aminoacids[:]
		peptide.hashtable = None
		gc.collect()

	'''new_best_models = []
	new_best_scores = [1000 for n in xrange(model_count)]
	for k, model in enumerate(best_models):
		print "Refining model {}: {}".format(k, model)
		for i, aa in enumerate(peptide.aminoacids):
			aa.acarbon = model[i].alpha_zone
			aa.set_axes(model[i].x_axis, model[i].y_axis, model[i].z_axis)
		#Keep track of how long the model has had this score.
		running_count = 0
		last_score = 0.0
		new_best_models.append([])
		print dists, peptide.aminoacids
		currscore = sum(dist.score(peptide, peptide.aminoacids) for dist in dists)
		print "First has score", currscore
		while running_count < 50:
			seglen = segment_length(scores[-1] / len(peptide.aminoacids))
			folding_iteration(peptide, prob, seglen)
			peptide.center()
			newscore = sum(dist.score(peptide, peptide.aminoacids) for dist in dists)
			if last_score == newscore:
				running_count += 1
			else:
				running_count = 0
				last_score = newscore
			if newscore < new_best_scores[k]:
				print "New best score"
				new_best_scores[k] = newscore
				new_best_models[k] = [PositionZone(aa.acarbon, aa.i, aa.j, aa.k) for aa in peptide.aminoacids]
			elif newscore > new_best_scores[k] + 20.0:
				print "Too unstable"
				break
	for model in new_best_models:
		for i, aa in enumerate(peptide.aminoacids):
			aa.acarbon = model[i].alpha_zone
			aa.set_axes(model[i].x_axis, model[i].y_axis, model[i].z_axis)
		file.write(system.pdb(modelno=pdb_model_idx))
		pdb_model_idx += 1
	print "Best model scores:"
	for i in xrange(model_count):
		print "{} ({})".format(i, new_best_scores[i])'''


	'''gentle_cutoff = 0
	pdb_model_idx = 2
	a = datetime.datetime.now()
	for i in xrange(500):
		constructive_folding_iteration(peptide, prob)
		#peptide.center()
		curscore = 0.0
		for aa in peptide.aminoacids:
			aa.localscore = sum(dist.score(peptide, [aa]) for dist in dists)
			curscore += aa.localscore
		scores.append(curscore)
		t_scores = ""
		for dist in dists: t_scores += str(dist.score(peptide, peptide.aminoacids)) + " (" + dist.identifier + ") "
		print i, t_scores
		#print i, scores[-1] / len(peptide.aminoacids)
		file.write(system.pdb(modelno=pdb_model_idx))
		pdb_model_idx += 1'''
		
	b = datetime.datetime.now()
	print "Finished segment fold in {0:.4f} sec.".format((b - a).total_seconds())
	del peptide.aminoacids[:]
	peptide = None
	del scores
	del system
	gc.collect()

	file.close()
	scoresfile.close()


def simulate_fold(sparc_dir, dists, seq, range, output, outname="simulation.pdb", sec_structs=None, model_count=5, n=500):
	"""For seq, pass the sequence of the entire protein. Only the range of amino acids specified in the tuple 'range' (inclusive) will be simulated. Pass the entire protein worth of secondary structures to sec_structs if it is available."""
	#"GRYRRCIPGMFRAYCYMD" (2LWT - GRY...MD, 2MDB - KWC...CR, 1QLQ - RPDF...GGA, insulin MALW...YCN)
	#"KWCFRVCYRGICYRRCR"
	#"TTCCPSIVARSNFNVCRLPGTPSEALICATYTGCIIIPGATCPGDYAN"
	#"RPDFCLEPPYAGACRARIIRYFYNAKAGLCQTFVYGGCRAKRNNFKSAEDCLRTCGGA"
	model_count = int(model_count)
	n = int(n)
	
	permissions = AAPermissionsManager(os.path.join(sparc_dir, "permissions"), os.path.join(sparc_dir, "permissible_sequences", "all.txt"))
	sec_struct_permissions = AASecondaryStructurePermissionsManager(os.path.join(sparc_dir, "permissible_sequences"))
	peptide = Polypeptide()
	#peptide.read("/Users/venkatesh-sivaraman/Downloads/1QLQ.pdb")
	if sec_structs:
		if ',' in sec_structs:
			peptide.add_secondary_structures(sec_structs, format='csv', range=range)
		else:
			peptide.add_secondary_structures(sec_structs, format='pdb', range=range)
	peptide.randomcoil(seq[range[0] - 1 : range[1]], permissions=permissions, struct_permissions=sec_struct_permissions)
	print seq[range[0] - 1 : range[1]], peptide.secondary_structures
	system = MolecularSystem([peptide])
	system.center()
	file = open(os.path.join(output, outname), 'w')
	scoresfile = open(os.path.join(output, outname[:-4] + "-scores.txt"), 'w')
	#file.write(peptide.xyz(escaped=False))
	scores = []
	t_scores = ""
	curscore = 0.0
	for dist in dists:
		sc = dist.score(peptide, peptide.aminoacids, system=system)
		t_scores += str(sc) + " (" + dist.identifier + ") "
		curscore += sc
	scores.append(curscore)
	scoresfile.write("-1 " + t_scores + "\n")
	file.write(system.pdb(modelno=1))
	
	prob = AAProbabilitySource(peptide, dists, permissions, sec_struct_permissions, system=system)
	gentle_cutoff = 0
	best_models = [[] for i in xrange(model_count)]
	best_scores = [1000 for i in xrange(model_count)]
	pdb_model_idx = 2
	a = datetime.datetime.now()
	time_wasted = 0.0
	for i in xrange(n):
		print "Iteration", i
		seglen = segment_length(scores[-1] / len(peptide.aminoacids))
		folding_iteration(system, [prob], seglen)
		peptide.center()
		for aa in peptide.aminoacids:
			aa.localscore = sum(dist.score(peptide, [aa], system=system) for dist in dists)
		t_scores = ""
		curscore = 0.0
		for dist in dists:
			sc = dist.score(peptide, peptide.aminoacids, system=system)
			t_scores += str(sc) + " (" + dist.identifier + ") "
			curscore += sc
		scores.append(curscore)
		#print i, t_scores
		scoresfile.write(str(i) + " " + t_scores + "\n")
		#print i, scores[-1] / len(peptide.aminoacids), dists[-1].score(peptide, peptide.aminoacids)
		#if scores[-1] / len(peptide.aminoacids) <= -100.0:
			#file.write(peptide.xyz(escaped=False))
		if curscore / len(peptide.aminoacids) <= 0.0:
			file.write(system.pdb(modelno=pdb_model_idx))
			pdb_model_idx += 1
		#print peptide.xyz(escaped=False)
		#Save the conformation if it is the best so far.
		for k in xrange(model_count):
			if scores[-1] < best_scores[k] * 1.2:
				for m in reversed(xrange(max(k, 1), model_count)):
					best_scores[m] = best_scores[m - 1]
					best_models[m] = best_models[m - 1]
				best_scores[k] = scores[-1]
				best_models[k] = [PositionZone(aa.acarbon, aa.i, aa.j, aa.k) for aa in peptide.aminoacids]
				break
			elif scores[-1] <= best_scores[k]:
				break
		'''if i == 24:
			#Save the best of the first ten scores as a comparison for later iterations.
			gentle_cutoff = min(scores)
		elif i > 24:
			if scores[-1] < gentle_cutoff and prob.mode == psource_erratic_mode:
				prob.mode = psource_gentle_mode
				print "Switching to gentle mode"
			elif scores[-1] >= gentle_cutoff and prob.mode == psource_gentle_mode:
				prob.mode = psource_erratic_mode
				print "Switching to erratic mode"'''

		#time_wasted += 0.5
		#time.sleep(0.5)

	print "Iterating over the best scores:", best_scores
	prob.mode = psource_gentle_mode
	new_best_models = []
	new_best_scores = [1000 for n in xrange(model_count)]
	for k, model in enumerate(best_models):
		if len(model) == 0: continue
		print "Refining model {}".format(k)
		for i, aa in enumerate(peptide.aminoacids):
			aa.acarbon = model[i].alpha_zone
			aa.set_axes(model[i].x_axis, model[i].y_axis, model[i].z_axis)
		#Keep track of how long the model has had this score.
		running_count = 0
		last_score = 0.0
		new_best_models.append([])
		currscore = sum(dist.score(peptide, peptide.aminoacids, system=system) for dist in dists)
		print "First has score", currscore
		while running_count < 50:
			seglen = segment_length(scores[-1] / len(peptide.aminoacids))
			folding_iteration(system, [prob], seglen)
			peptide.center()
			newscore = sum(dist.score(peptide, peptide.aminoacids, system=system) for dist in dists)
			if last_score == newscore:
				running_count += 1
			else:
				running_count = 0
				last_score = newscore
			print "{} ({})".format(newscore, running_count)
			if newscore < new_best_scores[k]:
				print "New best score"
				new_best_scores[k] = newscore
				new_best_models[k] = [PositionZone(aa.acarbon, aa.i, aa.j, aa.k) for aa in peptide.aminoacids]
			elif newscore > new_best_scores[k] + 20.0:
				print "Too unstable"
				break
	print "\nFinal:"
	for model in new_best_models:
		for i, aa in enumerate(peptide.aminoacids):
			aa.acarbon = model[i].alpha_zone
			aa.set_axes(model[i].x_axis, model[i].y_axis, model[i].z_axis)
		file.write(system.pdb(modelno=pdb_model_idx))
		pdb_model_idx += 1
	print "\n\n"
	print "Best model scores:"
	for i in xrange(model_count):
		print "{} ({})".format(i, new_best_scores[i])
	del scores
	del peptide.aminoacids[:]
	peptide = None
	b = datetime.datetime.now()
	print "Finished 750 iterations in {0:.4f} sec.".format((b - a).total_seconds() - time_wasted)
	scoresfile.write("\n" + str((b - a).total_seconds() - time_wasted))
	file.close()
	scoresfile.close()
	gc.collect()

def run_simulation(directives, output):
	"""This method takes the simulation directives found at the input path and performs the simulations."""
	seq = None
	sec_structs = None
	runs = []
	with open(directives, "r") as file:
		lines = file.readlines()
		seq = lines[0].strip()
		del lines[0]
		processing_runs = False
		for line in lines:
			if line[0] == "#": continue
			if len(line.strip()) == 0:
				processing_runs = True
				continue
			if processing_runs:
				comps = line.strip().split(";")
				runs.append([[y.strip() for y in x.split(",")] for x in comps])
			else:
				if not sec_structs: sec_structs = ""
				sec_structs += line
	sec_structs = sec_structs.strip()
	print sec_structs

	sparc_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "potential")
	weights = { "consec": 3.0, "secondary": 3.0, "short-range": 2.0, "nonconsec": 2.0, "medium": 0.0 }
	distributions = load_dists(sparc_dir, secondary=True, weights=weights)

	for run in runs:
		assert len(run) >= 2, "Run directive is invalid: {}".format(run)
		if len(run[1]) > 1:
			# run[1] must be the input paths, and run[2] must be the output path name
			extra_args = { "sec_structs": sec_structs }
			if len(run) > 3:
				# Get extra parameters
				for arg in run[3:]:
					kv = arg.split("=")
					extra_args[kv[0]] = kv[1]
			range1 = [int(x) for x in run[0][0].split("-")]
			range2 = [int(x) for x in run[0][1].split("-")]
			segment_fold(sparc_dir, distributions, seq, range1, range2, output, outname=run[2][0], **extra_args)
		else:
			# run[1] must be the output path name
			extra_args = { "sec_structs": sec_structs }
			if len(run) > 2:
				# Get extra parameters
				for arg in run[2:]:
					kv = arg.split("=")
					extra_args[kv[0]] = kv[1]
			range = [int(x) for x in run[0][0].split("-")]
			simulate_fold(sparc_dir, distributions, seq, range, output, outname=run[1][0], **extra_args)

if __name__ == '__main__':
	args = sys.argv[1:]
	input = None
	output = None
	i = 0
	while i < len(args):
		if args[i].lower() == "-i":
			assert len(args) > i + 1, "Not enough arguments"
			input = args[i + 1]
			i += 2
		elif args[i].lower() == "-o":
			assert len(args) > i + 1, "Not enough arguments"
			output = args[i + 1]
			i += 2
		else:
			assert False, "Unexpected command-line argument {}".format(args[i])
	run_simulation(input, output)