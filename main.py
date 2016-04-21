from pdbstats import *
from pdbanalysis import *
from generic_distributions import *
from central_distributions import *
from permissions import *
from randomcoil import *
from molecular_systems import *
import os, sys
import numpy
import scipy
from probsource import *
import cProfile
from multiprocessing import Process, Queue
from memory_profiler import profile
import datetime, time
from tmscore import *
from charmm import *
from sparc_distribution import *
from maxwell_scores import *
import loading_indicator

def test_probsource(prob, peptide):
	for conformation in prob._iter_permissible_randomcoils(peptide.aminoacids[4:8], peptide.aminoacids[3], peptide.aminoacids[8]):
		print "Found valid conformation."
		'''old = []
		for i, residue in enumerate(peptide.aminoacids[4:7]):
			old.append(PositionZone(residue.acarbon, residue.i, residue.j, residue.k))
			residue.acarbon = conformation[i].alpha_zone
			residue.set_axes(conformation[i].x_axis,
							 conformation[i].y_axis,
							 conformation[i].z_axis)
		print peptide.xyz(escaped=False, highlight=range(4, 7))
		for i, aa in enumerate(peptide.aminoacids[4:7]):
			aa.acarbon = old[i].alpha_zone
			aa.set_axes(old[i].x_axis, old[i].y_axis, old[i].z_axis)'''

def load_dists(basepath, weights={}, concurrent=True, secondary=True):
	print "Loading SPARC from {}...".format(basepath)
	if not reference_state.is_initialized():
		if os.path.exists(os.path.join(basepath, "reference_states.txt")):
			reference_state.load_reference_state(os.path.join(basepath, "reference_states.txt"))
		if os.path.exists(os.path.join(basepath, "possible_interactions")):
			reference_state.load_possible_interactions(os.path.join(basepath, "possible_interactions"))
	nonconsec = os.path.join(basepath, "long_range")
	if secondary or not os.path.exists(os.path.join(basepath, "consec+secondary")):
		consec = os.path.join(basepath, "consec") #+secondary
	else:
		consec = os.path.join(basepath, "consec+secondary")
	medium = os.path.join(basepath, "medium")
	short_range = os.path.join(basepath, "short_range")
	secondary_path = os.path.join(basepath, "secondary")
	DistClass = SPARCBasicDistributionManager #SPARCBothOrientationDistributionManager
	if concurrent == False:
		dist1 = MediumDistributionManager(medium)
		dist2 = DistClass(consec, True, blocks_sec_struct=secondary)
		if os.path.exists(short_range):
			dist3 = DistClass(nonconsec, False, short_range=False)
			dist5 = DistClass(short_range, False, short_range=True)
		else:
			dist3 = DistClass(nonconsec, False)
			dist5 = None
		if os.path.exists(secondary_path) and secondary:
			dist4 = SPARCSecondaryDistributionManager(secondary_path)
		else:
			dist4 = None
		dists = [dist1, dist2, dist3, dist4, dist5]
		for d in dists:
			if d and d.identifier in weights:
				d.weight = weights[d.identifier]

		#Always return in the order (consec, secondary, short-range, nonconsec, medium)
		loading_indicator.clear_loading_data()
		print "Finished loading."
		if secondary:
			return [x for x in [dist2, dist4, dist5, dist3, dist1] if x]
		else:
			return [x for x in [dist2, dist5, dist3, dist1] if x]
				
	processes = []
	queue = multiprocessing.Queue()
	def generate_distmanager(cls, q, *args):
		dist = cls(*args)
		if dist.identifier in weights:
			dist.weight = weights[dist.identifier]
		q.put(dist)
	p1 = multiprocessing.Process(target=generate_distmanager, args=(SPARCBasicDistributionManager, queue, nonconsec, False, False, False))
	processes.append(p1)
	p1.start()
	p2 = multiprocessing.Process(target=generate_distmanager, args=(SPARCBasicDistributionManager, queue, consec, True, True))
	processes.append(p2)
	p2.start()
	p3 = multiprocessing.Process(target=generate_distmanager, args=(MediumDistributionManager, queue, medium))
	processes.append(p3)
	p3.start()
	p4 = multiprocessing.Process(target=generate_distmanager, args=(SPARCSecondaryDistributionManager, queue, secondary_path))
	processes.append(p4)
	p4.start()
	p5 = multiprocessing.Process(target=generate_distmanager, args=(SPARCBasicDistributionManager, queue, short_range, False, False, True))
	processes.append(p5)
	p5.start()

	dists = []
	dists.append(queue.get())
	dists.append(queue.get())
	dists.append(queue.get())
	dists.append(queue.get())
	dists.append(queue.get())
	p1.join()
	p2.join()
	p3.join()
	p4.join()
	p5.join()

	if secondary:
		distributions = ["", "", "", "", ""]
		for dist in dists:
			if dist.identifier == "consec": distributions[0] = dist
			elif dist.identifier == "secondary": distributions[1] = dist
			elif dist.identifier == "short_range": distributions[2] = dist
			elif dist.identifier == "long_range": distributions[3] = dist
			elif "medium" in dist.identifier: distributions[4] = dist
	else:
		distributions = ["", "", "", ""]
		for dist in dists:
			if dist.identifier == "consec+secondary": distributions[0] = dist
			elif dist.identifier == "short_range": distributions[1] = dist
			elif dist.identifier == "long_range": distributions[2] = dist
			elif "medium" in dist.identifier: distributions[3] = dist
	print "Finished loading."
	loading_indicator.clear_loading_data()
	return distributions

def load_central_dist(basepath, secondary=True):
	central_manager = SPARCCentralDistributionManager(os.path.join(basepath, "default"), os.path.join(basepath, "random_coil_ref"))
	medium = os.path.join(basepath, "medium")
	medium_dist = MediumDistributionManager(medium)
	if secondary:
		return [SPARCCentralDistributionPuppet(central_manager, sparc_consecutive_mode), SPARCCentralDistributionPuppet(central_manager, sparc_secondary_mode), SPARCCentralDistributionPuppet(central_manager, sparc_short_range_mode), SPARCCentralDistributionPuppet(central_manager, sparc_long_range_mode), medium_dist]
	else:
		return [SPARCCentralDistributionPuppet(central_manager, sparc_consec_secondary_mode), SPARCCentralDistributionPuppet(central_manager, sparc_short_range_mode), SPARCCentralDistributionPuppet(central_manager, sparc_long_range_mode), medium_dist]

def extract_dist_weights(dists):
	ret = {}
	for d in dists:
		ret[d.identifier] = d.weight
	return ret

def apply_dist_weights(dists, w):
	for d in dists:
		if d.identifier in w:
			d.weight = w[d.identifier]

#{ "consec": 4.0, "secondary": 4.0, "short_range": 1.0, "long_range": 4.0, "medium": 5.0 }
def func(weights={ "consec": 3.0, "secondary": 3.0, "short_range": 2.0, "long_range": 2.0, "medium": 3.0 }, base="refined-bpti/"):
	#Weights used to be 2, 4, 8
	
	dists = load_dists(weights=weights) #load_dists(weights={frequency_nonconsec_disttype: 9.0, frequency_consec_disttype: 4.0, medium_disttype: 3.0})
	sec_struct_weights = { "consec": 3.0, "secondary": 3.0, "short_range": 1.0, "long_range": 1.0, "medium": 0.0 }
	#cProfile.runctx('simulate_fold(dists, seq="RPDFCLE", outname="segments/seg1.pdb")', {'dists': dists, 'simulate_fold': simulate_fold}, {})
	#Insulin - GIVEQCCTSICSLYQLENYCN, FVNQHLCGSHLVEALYLVCGERGFFYTPKT
	'''apply_dist_weights(dists, sec_struct_weights)
	simulate_fold(dists, seq="GIVEQCC", outname=base + "seg1.pdb", sec_structs="helix,1,1,7")
	apply_dist_weights(dists, weights)
	simulate_fold(dists, seq="TSIC", outname=base + "seg2.pdb")
	apply_dist_weights(dists, sec_struct_weights)
	simulate_fold(dists, seq="SLYQLEN", outname=base + "seg3.pdb", sec_structs="helix,1,1,7")
	apply_dist_weights(dists, weights)
	simulate_fold(dists, seq="YCN", outname=base + "seg4.pdb")
	simulate_fold(dists, seq="FVNQHLC", outname=base + "seg5.pdb", sec_structs="helix,1,1,7")
	apply_dist_weights(dists, sec_struct_weights)
	simulate_fold(dists, seq="GSHLVEAL", outname=base + "seg6.pdb", sec_structs="helix,1,1,8")
	simulate_fold(dists, seq="YLVCG", outname=base + "seg7.pdb", sec_structs="helix,1,1,5")
	simulate_fold(dists, seq="ERG", outname=base + "seg8.pdb", sec_structs="helix,5,1,3")
	apply_dist_weights(dists, weights)
	simulate_fold(dists, seq="FFYTPKT", outname=base + "seg9.pdb")'''
	'''segment_fold(dists, seq1="GIVEQCC", seq2="TSIC", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg1.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg2.pdb"], outname=base + "seg12.pdb", sec_structs="helix,1,1,7")
	segment_fold(dists, seq1="SLYQLEN", seq2="YCN", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg3.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg4.pdb"], outname=base + "seg34.pdb", sec_structs="helix,1,1,7")
	segment_fold(dists, seq1="FVNQHLC", seq2="GSHLVEAL", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg5.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg6.pdb"], outname=base + "seg56.pdb", sec_structs="helix,1,1,15")
	segment_fold(dists, seq1="YLVCG", seq2="ERG", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg7.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg8.pdb"], outname=base + "seg78.pdb", sec_structs="helix,1,1,5\nhelix,5,6,8")
	segment_fold(dists, seq1="GIVEQCCTSIC", seq2="SLYQLENYCN", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg12.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg34.pdb"], outname=base + "seg1234.pdb", sec_structs="helix,1,1,7\nhelix,1,12,18")
	segment_fold(dists, seq1="YLVCGERG", seq2="FFYTPKT", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg78.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg9.pdb"], outname=base + "seg789.pdb", sec_structs="helix,1,1,5\nhelix,5,6,8")'''
	#segment_fold(dists, seq1="FVNQHLCGSHLVEAL", seq2="YLVCGERGFFYTPKT", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg56.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg789.pdb"], outname=base + "seg56789.pdb", sec_structs="helix,1,1,15\nhelix,1,16,20\nhelix,5,21,23")

	#simulate_fold(dists, outname=base + "seg_all.pdb", model_count=20, n=2500, sec_structs="HELIX    1   1 PRO A    2  GLU A    7  5                                   6\nHELIX    2   2 SER A   47  GLY A   56  1                                  10\nSHEET    1   A 2 ILE A  18  ASN A  24  0\nSHEET    2   A 2 LEU A  29  TYR A  35 -1  N  TYR A  35   O  ILE A  18")
	
	print "Moving to segment 1"
	#apply_dist_weights(dists, sec_struct_weights)
	#simulate_fold(dists, seq="RPDFCLE", outname=base + "seg1.pdb", sec_structs="helix,5,2,7")
	'''print "Moving to segment 2"
	apply_dist_weights(dists, weights)
	simulate_fold(dists, seq="PPYAG", outname=base + "seg2.pdb")
	print "Moving to segment 3"
	simulate_fold(dists, seq="ACRAR", outname=base + "seg3.pdb")'''
	print "Moving to segment 4"
	apply_dist_weights(dists, sec_struct_weights)
	simulate_fold(dists, seq="IIRYFYN", outname=base + "seg4.pdb", sec_structs="sheet,0,1,7", n=1000)
	'''print "Moving to segment 5"
	apply_dist_weights(dists, weights)
	simulate_fold(dists, seq="AKAG", outname=base + "seg5.pdb")'''
	print "Moving to segment 6"
	apply_dist_weights(dists, sec_struct_weights)
	simulate_fold(dists, seq="LCQTFVY", outname=base + "seg6.pdb", sec_structs="sheet,0,1,7", n=1000)
	'''print "Moving to segment 7"
	apply_dist_weights(dists, weights)
	simulate_fold(dists, seq="GGCRA", outname=base + "seg7.pdb")
	print "Moving to segment 8"
	simulate_fold(dists, seq="KRNNFK", outname=base + "seg8.pdb")'''
	print "Moving to segment 9"
	apply_dist_weights(dists, sec_struct_weights)
	simulate_fold(dists, seq="SAEDC", outname=base + "seg9.pdb", sec_structs="helix,1,1,5")
	print "Moving to segment 10"
	simulate_fold(dists, seq="LRTCGGA", outname=base + "seg10.pdb", sec_structs="helix,1,1,5")

def test_folding_parameters(dists):
	permissions = AAPermissionsManager("/Users/venkatesh-sivaraman/Desktop/sciencefair/allowed-zones")
	peptide = Polypeptide()
	#"GRYRRCIPGMFRAYCYMD" (2LWT - GRY...MD, 2MDB - KWC...CR)
	#"KWCFRVCYRGICYRRCR"
	
	prob = AAProbabilitySource(peptide, dists, permissions)
	for cutoff in numpy.arange(3.0, 5.5, 1.0):
		for proximity in numpy.arange(0.4, 1.1, 0.2):
			prob.steric_cutoff = cutoff
			prob.erratic_proximity = proximity
			print "Starting to analyze {} and {}".format(cutoff, proximity)
			peptide.randomcoil("TTCCPSIVARSNFNVCRLPGTPSEALICATYTGCIIIPGATCPGDYAN", permissions)
			center = Point3D.zero()
			for aa in peptide.aminoacids:
				center = center.add(aa.acarbon)
			center = center.multiply(1.0 / len(peptide.aminoacids))
			for aa in peptide.aminoacids:
				aa.acarbon = aa.acarbon.subtract(center)
			scores = [sum(dist.score(peptide, peptide.aminoacids) for dist in dists)]
			print scores[-1], "Avg:", scores[-1] / len(peptide.aminoacids)
			for i in xrange(100):
				seglen = segment_length(scores[-1] / len(peptide.aminoacids))
				folding.folding_iteration(peptide, prob, seglen)
				center = Point3D.zero()
				for aa in peptide.aminoacids:
					center = center.add(aa.acarbon)
				center = center.multiply(1.0 / len(peptide.aminoacids))
				for aa in peptide.aminoacids:
					aa.acarbon = aa.acarbon.subtract(center)
				scores.append(sum(dist.score(peptide, peptide.aminoacids) for dist in dists))
				print i, scores[-1] / len(peptide.aminoacids)
			print "Average score: {}".format(float(sum(scores)) / float(len(scores)))
			'''for score in xrange(int(math.floor(min(scores))), int(math.ceil(max(scores)))):
				print score, sum((s >= score and s < score + 1.0) for s in scores)'''
			del scores[:]

def randomcoiltest():
	peptide = Polypeptide()
	peptide.randomcoil("ARDRFGANMILILGGA")
	#print "movie.addFrame([ChemDoodle.readXYZ(\'" + peptide.xyz() + "\')],[]);"
	print peptide.xyz(escaped=False)
	
	prob = ProbabilitySource(peptide)
	for i in xrange(400):
		folding.folding_iteration(peptide, prob, random.randint(1,4))
		#print "movie.addFrame([ChemDoodle.readXYZ(\'" + peptide.xyz() + "\')],[]);"
		#print i
		print peptide.xyz(escaped=False)

def analysis():
	root = "/Volumes/External Hard Drive/School archives/Science Fair/2015/sciencefair/filtered-frequencies-consec" #"/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC/nonconsec"
	total_zones = 0
	total_count = 0
	for path in os.listdir(root):
		if path.find(".txt") == -1: continue
		print path
		zones = quantiles(os.path.join(root, path))
		if path == "all.txt":
			#total_zones += zones * total_count
			pass
		else:
			total_zones += zones
			total_count += 1
	print "Average:", float(total_zones) / float(total_count)

def list_alphacarbons():
	root = "/Volumes/External Hard Drive/Science Fair 2014-15/non-redundant-pdb-data-2"
	output = "/Volumes/External Hard Drive/Science Fair 2014-15/test-nonconsec"
	calculate_pdb_alphacarbons(root, output, mode=alphacarbon_mode)

def test_directedrandomwalk():
	startpt = Point3D(0.0, 0.0, 0.0)
	endpt = Point3D(6.0, 6.0, 6.0)
	for i in xrange(100):
		zones = directed_randomwalk(startpt, endpt, 3, 3.0)
		print "%d\nNo comment\n" % (len(self.aminoacids) * 4)
		zones.insert(0, PositionZone(startpt))
		zones.append(PositionZone(endpt))
		for zone in zones:
			print "C\t%.4f\t%.4f\t%.4f\n" % (zone.alpha_zone.x, zone.alpha_zone.y, zone.alpha_zone.z)
		print "\n"

def illustrate_permissible():
	protein = Polypeptide()
	aa1 = AminoAcid("ALA", 0, acarbon=Point3D(0.0, 0.0, 0.0))
	aa1.set_axes(Point3D(1.0, 0.0, 0.0), Point3D(0.0, 1.0, 0.0), Point3D(0.0, 0.0, 1.0))
	aa2 = AminoAcid("ALA", 1)
	protein.aminoacids = [aa1, aa2]
	permissions = AASecondaryStructurePermissionsManager("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC 3/permissible_sequences")
	for pz in permissions.allowed_conformations(aa2, aa1, "helix", 1):
		aa2.acarbon = pz.alpha_zone
		aa2.set_axes(pz.x_axis, pz.y_axis, pz.z_axis)
		print protein.xyz(escaped=False)

def supplement_natives(input):
	nativepath = "/Users/venkatesh-sivaraman/Downloads/casp11.targets_unsplitted.release11242014"
	distributions = load_dists(concurrent=False)
	dists = ["", "", ""]
	for dist in distributions:
		if dist.type == frequency_nonconsec_disttype: dists[0] = dist
		elif dist.type == frequency_consec_disttype: dists[1] = dist
		elif dist.type == medium_disttype: dists[2] = dist

	for path in os.listdir(input):
		if path == "T0781.txt" or path == "T0786.txt" or path == "T0791.txt" or path == "T0801.txt" or path == "T0800.txt": continue
		if ".txt" not in path: continue
		with open(os.path.join(input, path), "r") as file:
			contents = file.read()
			nativename = path[:path.find(".txt")] + "_orig.pdb"
			nativefile = path[:path.find(".txt")] + ".pdb"
			if nativename not in contents and nativefile not in contents:
				print path
				time.sleep(5)
				if os.path.exists(join(nativepath, nativefile)):
					bounds, scores = sparc_scores_file(join(nativepath, nativefile), dists, retbounds=True)
					if scores is not None:
						with open(os.path.join(input, path), "a") as file:
							file.write("{}; {}, {}, {}\n".format(nativename, scores[0], scores[1], scores[2]))
				else:
					print join(nativepath, nativefile), "does not exist."
		del contents
		gc.collect()

def analyze_relative_orientations(path):
	peptide = Polypeptide()
	peptide.read(path, otheratoms=True)
	for i, aa in enumerate(peptide.aminoacids):
		if len(peptide.aminoacids) > i + 1:
			zone = aa.aa_position_zone(peptide.aminoacids[i + 1]).alpha_zone
			retro_zone = peptide.aminoacids[i + 1].aa_position_zone(aa).alpha_zone
			print "{}\t{}\t{}".format(i, zone, retro_zone)

def aminoacid_type_variation():
	aa1 = AminoAcid(amino_acid_alanine, 1)
	aa1.set_axes(Point3D(1, 0, 0), Point3D(0, 1, 0), Point3D(0, 0, 1))
	aa2 = AminoAcid(amino_acid_alanine, 2)
	aa2.set_axes(Point3D(1, 0, 0), Point3D(0, 1, 0), Point3D(0, 0, 1))
	basepath = "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC 3"
	if not reference_state.is_initialized() and os.path.exists(os.path.join(basepath, "reference_states.txt")):
		reference_state.load_reference_state(os.path.join(basepath, "reference_states.txt"))
	nonconsec = os.path.join(basepath, "short-range")
	dist = SPARCBasicDistributionManager(nonconsec, False, short_range=True)
	for i in xrange(AMINO_ACID_COUNT):
		aa2.type = aatype(i)
		for point in Point3D.zero().iteroffsets(10.0):
			if point.x != 0.0: continue
			aa2.acarbon = point
			freq = dist.alpha_frequency(aacode(aa1.type), aacode(aa2.type), aa1.tolocal(aa2.acarbon))
			print "{}\t{}\t{}".format(point.y, point.z, freq)
		print "\n"

if __name__ == '__main__':
	#func()
	#generate_distance_constrained_bins("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/bins_test.txt")
	#print z_score("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/tasser", w, structure_files="/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoys/tasser-decoys")
	#print z_score("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/casp", w, structure_files="/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoys/casp-decoys")
	#print determine_omits("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Nonredundant/all_pdb_ids.txt", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Nonredundant/omits.txt")
	#best_weights_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/bpti-analysis/sparc_scores", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/bpti-analysis/rmsds", None, numweights=4, start=[0, 2, 3, 4])
	#load_central_dist("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/consolidated-sparc/SPARC 4")
	#best_weights_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/corrected", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/TM-scores/casp", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/casp-natives", numweights=4, tmscore=True)
	#best_weights_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/corrected", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/TM-scores/casp", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/casp-natives", numweights=4, tmscore=True)
	#best_weights("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/rw/casp-rw", numweights=1)
	#best_weights("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/rw/tasser-rw", numweights=1)
	#best_weights("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/goap/casp-goap", numweights=1)
	#best_weights("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/goap/tasser-goap", numweights=1)
	#best_weights_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/ref-tests/average", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/TM-scores/casp", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/casp-natives", numweights=4)
	#Average: Final: the combos [1, 1, 5, 1] had a total of 1 correct guesses, with R^2 0.486152814813
	#Interaction Median: Final: the combos [1, 5, 5, 1] had a total of 1 correct guesses, with R^2 0.335183666746
	#best_weights_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/tasser-correct-orient", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/TM-scores/tasser", None, numweights=4, tmscore=True) #"/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/casp-natives"
	#best_weights_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/bpti-long-refine-scores/sparc", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/bpti-long-refine-scores/tm", None, numweights=5)
	#analyze_relative_orientations("/Users/venkatesh-sivaraman/Downloads/1QLQ.pdb")
	#decoys_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/tasser-decoys", None, "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/rmsd") #/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoys/casp-decoys,
	#min1 = min_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/bpti-laptop/seg12.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments/native.pdb", range=[1, 12], writeout=True)
	#min2 = min_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/bpti-laptop/seg56.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments/native.pdb", range=[36, 46])
	#print min1, min2
	#print "We did main"
	'''lowest_scores = [ [10000, None], [10000, None], [10000, None], [10000, None], [10000, None] ]
	basepath = "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segment_weight_test"
	for fname in os.listdir(basepath):
		if "pdb" not in fname or "seg" not in fname: continue
		idx = int(fname[fname.find("seg") + 3 : fname.find(".")])
		weights = [int(x) for x in fname[:3]]
		print idx, weights
		min = min_rmsd(os.path.join(basepath, fname), os.path.join("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments", "native_seg" + str(idx) + ".pdb"))
		if min < lowest_scores[idx - 1][0]:
			lowest_scores[idx - 1][0] = min
			lowest_scores[idx - 1][1] = weights
			print "New best"
	print lowest_scores'''
	'''sparc_dir = "/Users/venkatesh-sivaraman/Documents/Xcode Projects/PythonProteins/potential"
	dists_noref = load_dists(sparc_dir, concurrent=True, secondary=True)
	for d in dists_noref: d.refstate = False
	dists_yesref = load_dists(sparc_dir, concurrent=True, secondary=True)
	for d in dists_yesref: d.refstate = True
	distributions = dists_noref + dists_yesref
	min_rmsd("/Users/venkatesh-sivaraman/Desktop/bpti/seg4.pdb", "/Users/venkatesh-sivaraman/Downloads/1QLQ.pdb", range=[18, 35], dists=distributions, separate_scores=True)'''
	'''mins = []
	for i in range(1, 9):
		mins.append(min_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/refined-bpti/seg" + str(i) + ".pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/refined-bpti/native_seg" + str(i) + ".pdb"))
	print mins'''
	#link_decoy_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/casp", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/rmsd-casp-backbone", [9.0, 4.0, 3.0])
	#protein_protein_energies("/Users/venkatesh-sivaraman/Downloads/1DUM.pdb", load_dists(secondary=False), "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/magainin_test_yz_gromos.txt")
	'''basepath = "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments"
	for decoyset in os.listdir(basepath):
		if "segments" in basepath and ("native" in decoyset or not os.path.exists(os.path.join(basepath, "native_" + decoyset))): continue
		if decoyset[0] == ".": continue
		if os.path.exists(os.path.join("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/TM-scores/simul", decoyset + ".txt")): continue
		calculate_tm_scores(os.path.join(basepath, decoyset), os.path.join("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/TM-scores/simul", decoyset + ".txt"), natives=os.path.join(basepath, "native_" + decoyset)) #, natives="/Users/venkatesh-sivaraman/Downloads/casp11.targets_unsplitted.release11242014")'''
	#batch_compare_charmm_sparc("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/casp-decoys", dists, "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/casp_charmm.txt")
	#score_structure_file("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments-test/seg12-3.pdb", distributions)
	#sparc_score_sequences("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Nonredundant/all_pdb_ids.txt", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/sequence_scores_charmm.txt", dists, charmm_scores=False)
	'''for x in xrange(1, 4):
		for y in xrange(1, 4):
			for z in xrange(1, 4):
				print "Testing weights", x, y, z
				func(weights={ "consec": x, "secondary": x, "short-range": y, "nonconsec": 0.0, "medium": z }, base="segment_weight_test_sec/" + str(x) + str(y) + str(z))'''