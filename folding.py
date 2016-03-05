from distributions import *
from polypeptide import *
from probsource import *
import random

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
			print "Going back"
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
			print "Going back"
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
			print "Restoring"
			for aa in segment:
				aa.restore()
		return None
	else:
		if restore == True:
			hypotheticals = []
			print "Doing something with hypotheticals"
			for aa in protein.aminoacids[application_ret[0] : application_ret[1]]:
				hypotheticals.append(aa.hypothetical(aa.restore(), True))
			return hypotheticals
		else:
			return True

def folding_iteration(protein, psource, segment_length=1, file=None):
	"""The backbone of the protein folding simulator. The gist of the algorithm is to choose a random subset of the amino acids in 'protein' (which is a Polypeptide) and move them randomly. The function returns no value, just updates the protein.
		
		folding_iteration uses the ProbabilitySource object's probabilities() and angleprobabilities() methods to determine where to move the segment. psource should accept a list of relevant amino acids (the residues for which probabilities are required) and return a list of lists (form: [[conformation, probability], ...]) representing a cumulative distribution function and sorted by probability. (The conformation should be specified as a list of position zones, one for each amino acid.)
		
		The segment_length parameter can be used to manipulate the scale/granularity of the mutation. If 0 is passed, only the orientation of a single amino acid will be changed. For 1, the position and orientation will be changed. For any higher integer, a segment of that length will be pivoted/rotated/translated.
		
		Note that the final transformation is random and not guaranteed to match one of the probabilities exactly (no lattice is involved). Psource provides a randomization_margin function that uses score to determine how much randomization to introduce.
		
		Specify an open file object to write intermediate output to it."""
	
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
			for aa in protein.aminoacids:
				aa.restore()
			return

	violated = False
	for aa in protein.aminoacids:
		if len(protein.nearby_aa(aa, psource.steric_cutoff, consec=False, mindiff=psource.steric_consec_diff)) > 0:
			violated = True
			break
	if violated is True:
		print "Steric violation."
		for aa in protein.aminoacids:
			aa.restore()
	else:
		aa.discard_save()

def constructive_folding_iteration(protein, psource, file=None):
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
		if len(protein.nearby_aa(aa, psource.steric_cutoff, consec=False, mindiff=psource.steric_consec_diff)) > 0:
			violated = True
			break
	if violated is True:
		print "Steric violation."
		for aa in protein.aminoacids:
			aa.restore()
	else:
		aa.discard_save()