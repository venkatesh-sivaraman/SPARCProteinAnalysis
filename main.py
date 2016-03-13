from pdbstats import *
from pdbanalysis import *
from generic_distributions import *
from permissions import *
from randomcoil import *
import os, sys
import folding
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

def load_dists(weights={}, basepath="/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC 3", concurrent=True, secondary=True):
	print "Loading SPARC..."
	if not reference_state.is_initialized() and os.path.exists(os.path.join(basepath, "reference_states.txt")):
		reference_state.load_reference_state(os.path.join(basepath, "reference_states.txt"))
	nonconsec = os.path.join(basepath, "nonconsec")
	if secondary or not os.path.exists(os.path.join(basepath, "consec+secondary")):
		consec = os.path.join(basepath, "consec") #+secondary
	else:
		consec = os.path.join(basepath, "consec+secondary")
	medium = os.path.join(basepath, "medium")
	short_range = os.path.join(basepath, "short-range")
	secondary_path = os.path.join(basepath, "secondary")
	if concurrent == False:
		dist1 = MediumDistributionManager(medium)
		dist2 = SPARCBasicDistributionManager(consec, True, blocks_sec_struct=secondary)
		if os.path.exists(short_range):
			dist3 = SPARCBasicDistributionManager(nonconsec, False, short_range=False)
			dist5 = SPARCBasicDistributionManager(short_range, False, short_range=True)
		else:
			dist3 = SPARCBasicDistributionManager(nonconsec, False)
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
			elif dist.identifier == "short-range": distributions[2] = dist
			elif dist.identifier == "nonconsec": distributions[3] = dist
			elif "medium" in dist.identifier: distributions[4] = dist
	else:
		distributions = ["", "", "", ""]
		for dist in dists:
			if dist.identifier == "consec+secondary": distributions[0] = dist
			elif dist.identifier == "short-range": distributions[1] = dist
			elif dist.identifier == "nonconsec": distributions[2] = dist
			elif "medium" in dist.identifier: distributions[3] = dist
	print "Finished loading."
	loading_indicator.clear_loading_data()
	return distributions

def process_decoys((input, weights, path)):
	if not os.path.isdir(join(input, path)):
		process_decoys.q.put("")
		return
	dists = load_dists(weights={frequency_nonconsec_disttype: weights[0], frequency_consec_disttype: weights[1], medium_disttype: weights[2]}, concurrent=False)
	raw = sparc_scores(join(input, path), dists)
	scores = sorted(raw.items(), key=lambda x: x[1])
	ret = "{}, {}, {}:".format(*weights) + scores[0][0] + "\n"
	process_decoys.q.put(ret)

def process_decoys_file((input, output, nativepath)):
	if not os.path.isdir(input): return
	if os.path.exists(output): return
	if os.path.basename(input) == "doc": return
	protein_name = os.path.basename(input)
	print protein_name
	dists_old = load_dists(concurrent=False, secondary=False, basepath="/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC")
	for d in dists_old: d.refstate = False
	dists_noref = load_dists(concurrent=False, secondary=False, basepath="/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC 3")
	for d in dists_noref: d.refstate = False
	dists_yesref = load_dists(concurrent=False, secondary=False)
	for d in dists_yesref: d.refstate = True
	distributions = dists_old + dists_noref + dists_yesref #load_dists(concurrent=False, secondary=False)
	paths = os.listdir(input)
	allpaths = [os.path.join(input, path) for path in paths]

	if "-" in protein_name:
		path = protein_name[:protein_name.find("-")] + ".pdb"
	else:
		path = protein_name + ".pdb"
	if nativepath: allpaths.append(os.path.join(nativepath, path))
	'''max_min = -10000000	#Maximum min bound
	min_max = 1000000	#Minimum max bound
	gaps = []
	for boundpath in allpaths:
		if os.path.basename(boundpath).count("_") > 1: continue
		bounds, newgaps = sparc_scores_file(boundpath, distributions, retbounds=True, noeval=True)
		if not bounds: continue
		if bounds[0] > max_min: max_min = bounds[0]
		if bounds[1] < min_max: min_max = bounds[1]
		for gap in newgaps:
			if gap not in gaps: gaps.append(gap)
	bounds = (max_min, min_max)'''

	nativescores = None
	bounds = None
	if nativepath:
		if os.path.exists(join(nativepath, path)):
			bounds, gaps, scores = sparc_scores_file(join(nativepath, path), distributions, retbounds=True) #, ignored_aas=gaps
			print "native scores:", scores
			nativescores = scores
			if output and scores is not None:
				scorestr = ""
				for s in scores: scorestr += str(s) + ","
				scorestr = scorestr[:-1]
				with open(output, "w") as file:
					file.write("{}; {}\n".format(protein_name + "_orig.pdb", scorestr))
		else:
			print join(nativepath, path), "does not exist."
	else:
		scores = None
	for path in paths:
		if path == "list" or path == "rmsds": continue
		try:
			scores = sparc_scores_file(join(input, path), distributions, bounds=bounds) #, ignored_aas=gaps
		except:
			print path, "exception"
			continue
		if output and scores is not None:
			scorestr = ""
			for s in scores: scorestr += str(s) + ","
			scorestr = scorestr[:-1]
			with open(output, "a") as file:
				file.write("{}; {}\n".format(path, scorestr))
		del scores
		gc.collect()
	ret = ""
	print "Done"
	del paths
	del distributions
	gc.collect()


def test_sparc(input, output):
	files = os.listdir(input)
	print len(files), "files"
	if "casp11_seqs.txt" in files:
		files.remove("casp11_seqs.txt")
	'''for file in files:
		if not os.path.isdir(join(input, file)): continue
		if os.path.exists(join(output, file + ".txt")): continue
		if file == "doc": continue
		process_decoys_file((join(input, file), join(output, file + ".txt"), None))
		gc.collect()'''
	pool = multiprocessing.Pool(processes=2, maxtasksperchild=1)
	zipped = [(join(input, file), join(output, file + ".txt"), "/Users/venkatesh-sivaraman/Downloads/casp11.targets_unsplitted.release11242014") for file in reversed(files)]
	#print zipped
	pool.map(process_decoys_file, zipped)
	pool.close()
	pool.join()
	print "done"

def apply_dist_weights(dists, w):
	for d in dists:
		if d.identifier in w:
			d.weight = w[d.identifier]

#{ "consec": 4.0, "secondary": 4.0, "short-range": 1.0, "nonconsec": 4.0, "medium": 5.0 }
def func(weights={ "consec": 3.0, "secondary": 3.0, "short-range": 2.0, "nonconsec": 2.0, "medium": 3.0 }, base="insulin/"):
	#Weights used to be 2, 4, 8
	
	dists = load_dists(weights=weights) #load_dists(weights={frequency_nonconsec_disttype: 9.0, frequency_consec_disttype: 4.0, medium_disttype: 3.0})
	sec_struct_weights = { "consec": 3.0, "secondary": 3.0, "short-range": 1.0, "nonconsec": 1.0, "medium": 0.0 }
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
	segment_fold(dists, seq1="GIVEQCC", seq2="TSIC", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg1.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg2.pdb"], outname=base + "seg12.pdb", sec_structs="helix,1,1,7")
	segment_fold(dists, seq1="SLYQLEN", seq2="YCN", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg3.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg4.pdb"], outname=base + "seg34.pdb", sec_structs="helix,1,1,7")
	segment_fold(dists, seq1="FVNQHLC", seq2="GSHLVEAL", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg5.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg6.pdb"], outname=base + "seg56.pdb", sec_structs="helix,1,1,15")
	segment_fold(dists, seq1="YLVCG", seq2="ERG", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg7.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg8.pdb"], outname=base + "seg78.pdb", sec_structs="helix,1,1,5\nhelix,5,6,8")
	segment_fold(dists, seq1="GIVEQCCTSIC", seq2="SLYQLENYCN", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg12.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg34.pdb"], outname=base + "seg1234.pdb", sec_structs="helix,1,1,7\nhelix,1,12,18")
	segment_fold(dists, seq1="YLVCGERG", seq2="FFYTPKT", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg78.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg9.pdb"], outname=base + "seg789.pdb", sec_structs="helix,1,1,5\nhelix,5,6,8")
	#segment_fold(dists, seq1="FVNQHLCGSHLVEAL", seq2="YLVCGERGFFYTPKT", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg56.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/insulin/seg789.pdb"], outname=base + "seg56789.pdb", sec_structs="helix,1,1,15\nhelix,1,16,20\nhelix,5,21,23")

	#simulate_fold(dists, outname=base + "seg_all.pdb", model_count=20, n=2500, sec_structs="HELIX    1   1 PRO A    2  GLU A    7  5                                   6\nHELIX    2   2 SER A   47  GLY A   56  1                                  10\nSHEET    1   A 2 ILE A  18  ASN A  24  0\nSHEET    2   A 2 LEU A  29  TYR A  35 -1  N  TYR A  35   O  ILE A  18")
	'''
	print "Moving to segment 1"
	apply_dist_weights(dists, sec_struct_weights)
	simulate_fold(dists, seq="RPDFCLE", outname=base + "seg1.pdb", sec_structs="helix,5,2,7")
	print "Moving to segment 2"
	apply_dist_weights(dists, weights)
	simulate_fold(dists, seq="PPYAG", outname=base + "seg2.pdb")
	print "Moving to segment 3"
	simulate_fold(dists, seq="ACRAR", outname=base + "seg3.pdb")
	print "Moving to segment 4"
	apply_dist_weights(dists, sec_struct_weights)
	simulate_fold(dists, seq="IIRYFYN", outname=base + "seg4.pdb", sec_structs="sheet,0,1,7")
	print "Moving to segment 5"
	apply_dist_weights(dists, weights)
	simulate_fold(dists, seq="AKAG", outname=base + "seg5.pdb")
	print "Moving to segment 6"
	apply_dist_weights(dists, sec_struct_weights)
	simulate_fold(dists, seq="LCQTFVY", outname=base + "seg6.pdb", sec_structs="sheet,0,1,7")
	print "Moving to segment 7"
	apply_dist_weights(dists, weights)
	simulate_fold(dists, seq="GGCRA", outname=base + "seg7.pdb")
	print "Moving to segment 8"
	simulate_fold(dists, seq="KRNNFK", outname=base + "seg8.pdb")
	print "Moving to segment 9"
	apply_dist_weights(dists, sec_struct_weights)
	simulate_fold(dists, seq="SAEDC", outname=base + "seg9.pdb", sec_structs="helix,1,1,5")
	print "Moving to segment 10"
	simulate_fold(dists, seq="LRTCGGA", outname=base + "seg10.pdb", sec_structs="helix,1,1,5")
	#test_folding_parameters(dists)
	apply_dist_weights(dists, weights)'''
	'''segment_fold(dists, seq1="RPDFCLE", seq2="PPYAG", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments-test/seg1.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments-test/seg2.pdb"], outname=base + "seg12.pdb")
	segment_fold(dists, seq1="ACRAR", seq2="IIRYFYN", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments-test/seg3.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments-test/seg4.pdb"], outname=base + "seg34.pdb")
	segment_fold(dists, seq1="AKAG", seq2="LCQTFVY", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments-test/seg5.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments-test/seg6.pdb"], outname=base + "seg56.pdb")
	segment_fold(dists, seq1="GGCRA", seq2="KRNNFK", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments-test/seg7.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments-test/seg8.pdb"], outname=base + "seg78.pdb")
	segment_fold(dists, seq1="SAEDC", seq2="LRTCGGA", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments-test/seg9.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments-test/seg10.pdb"], outname=base + "seg910.pdb")
	segment_fold(dists, seq1="RPDFCLEPPYAG", seq2="ACRARIIRYFYN", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/new-segments/seg12.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/new-segments/seg34.pdb"], outname=base + "seg1234.pdb", sec_structs="helix,5,2,7\nsheet,0,18,24")
	segment_fold(dists, seq1="AKAGLCQTFVY", seq2="GGCRAKRNNFK", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/new-segments/seg56.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/new-segments/seg78.pdb"], outname=base + "seg5678.pdb", sec_structs="sheet,0,5,11")
	segment_fold(dists, seq1="AKAGLCQTFVYGGCRAKRNNFK", seq2="SAEDCLRTCGGA", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/new-segments/seg5678.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/new-segments/seg910.pdb"], outname=base + "seg5678910.pdb", sec_structs="sheet,0,5,11\nhelix,1,23,32")
	segment_fold(dists, seq1="RPDFCLEPPYAGACRARIIRYFYN", seq2="AKAGLCQTFVYGGCRAKRNNFKSAEDCLRTCGGA", infiles=["/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/new-segments/seg1234.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/new-segments/seg5678910.pdb"], outname=base + "seg12345678910.pdb", sec_structs="helix,5,2,7\nsheet,0,18,24\nsheet,0,29,35\nhelix,1,47,56")'''
	'''simulate_fold(dists, sec_structs="HELIX    1 AA1 PRO A  121  ILE A  125  5                                   5 \n\
		HELIX    2 AA2 ASP A  156  ASP A  161  5                                   6 \n\
		HELIX    3 AA3 ILE A  162  ARG A  166  5                                   5 \n\
		HELIX    4 AA4 LEU A  198  ALA A  204  1                                   7 \n\
		HELIX    5 AA5 THR A  209  HIS A  230  1                                  22 \n\
		HELIX    6 AA6 LYS A  238  PRO A  240  5                                   3 \n\
		HELIX    7 AA7 ASP A  263  GLY A  272  1                                  10 \n\
		HELIX    8 AA8 THR A  273  MET A  277  5                                   5 \n\
		HELIX    9 AA9 ALA A  278  ARG A  283  1                                   6 \n\
		HELIX   10 AB1 GLU A  289  GLY A  306  1                                  18 \n\
		HELIX   11 AB2 ASP A  314  SER A  324  1                                  11 \n\
		HELIX   12 AB3 PRO A  336  TRP A  347  1                                  12 \n\
		HELIX   13 AB4 LYS A  350  ARG A  354  5                                   5 \n\
		HELIX   14 AB5 SER A  356  SER A  373  1                                  18 \n\
		HELIX   15 AB6 PRO A  375  ILE A  397  1                                  23 \n\
		SHEET    1 AA1 5 LEU A 126  GLY A 134  0 \n\
		SHEET    2 AA1 5 GLY A 137  PHE A 144 -1  O  LEU A 141   N  GLN A 129 \n\
		SHEET    3 AA1 5 GLU A 147  LYS A 153 -1  O  LYS A 153   N  ALA A 138 \n\
		SHEET    4 AA1 5 CYS A 187  GLU A 191 -1  O  ILE A 188   N  LYS A 152 \n\
		SHEET    5 AA1 5 PHE A 176  CYS A 180 -1  N  LYS A 177   O  LEU A 189 \n\
		SHEET    1 AA2 3 GLY A 196  GLN A 197  0 \n\
		SHEET    2 AA2 3 MET A 242  ILE A 244 -1  O  ILE A 244   N  GLY A 196 \n\
		SHEET    3 AA2 3 VAL A 250  ILE A 252 -1  O  LYS A 251   N  LEU A 243           ")'''

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

def segment_fold(dists, seq1="", seq2="", infiles=[], sec_structs=None, outname="simulation_test.pdb", cluster_confs=25, sims=75, candidates=20):
	permissions = AAPermissionsManager("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC 3/permissions", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC 3/permissible_sequences/all.txt")
	peptide = Polypeptide()
	seg_prob = AAConstructiveProbabilitySource(peptide, (0, len(seq1)), (len(seq1), len(seq1) + len(seq2)), dists, permissions)
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

	file = open(os.path.join("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations", outname), 'w')
	scoresfile = open(os.path.join("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations", outname[:-4] + "-scores.txt"), 'w')
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
				peptide.add_secondary_structures(sec_structs, format='csv')
			else:
				peptide.add_secondary_structures(sec_structs, format='pdb')
		peptide.center()

		scores = []
		t_scores = ""
		curscore = 0.0
		for dist in dists:
			sc = dist.score(peptide, peptide.aminoacids)
			t_scores += str(sc) + " (" + dist.identifier + ") "
			curscore += sc
		scores.append(curscore)
		scoresfile.write(str(pdb_model_idx) + " " + t_scores + "\n")
		file.write(peptide.pdb(modelno=pdb_model_idx))
		pdb_model_idx += 1
		for n in xrange(sims):
			#seglen = segment_length(scores[-1] / len(peptide.aminoacids))
			folding.constructive_folding_iteration(peptide, seg_prob)
			peptide.center()
			curscore = 0.0
			t_scores = ""
			for dist in dists:
				sc = dist.score(peptide, peptide.aminoacids)
				t_scores += str(sc) + " (" + dist.identifier + ") "
				curscore += sc
			scores.append(curscore)
			scoresfile.write(str(pdb_model_idx) + " " + t_scores + "\n")
			file.write(peptide.pdb(modelno=pdb_model_idx))
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
			folding.folding_iteration(peptide, prob, seglen)
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
		file.write(peptide.pdb(modelno=pdb_model_idx))
		pdb_model_idx += 1
	print "Best model scores:"
	for i in xrange(model_count):
		print "{} ({})".format(i, new_best_scores[i])'''


	'''gentle_cutoff = 0
	pdb_model_idx = 2
	a = datetime.datetime.now()
	for i in xrange(500):
		folding.constructive_folding_iteration(peptide, prob)
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
		file.write(peptide.pdb(modelno=pdb_model_idx))
		pdb_model_idx += 1'''
		
	b = datetime.datetime.now()
	print "Finished segment fold in {0:.4f} sec.".format((b - a).total_seconds())
	del peptide.aminoacids[:]
	peptide = None
	del scores
	gc.collect()

	file.close()
	scoresfile.close()


def simulate_fold(dists, seq="RPDFCLEPPYAGACRARIIRYFYNAKAGLCQTFVYGGCRAKRNNFKSAEDCLRTCGGA", outname="simulation_final_5.pdb", sec_structs=None, model_count=5, n=500):
	#"GRYRRCIPGMFRAYCYMD" (2LWT - GRY...MD, 2MDB - KWC...CR, 1QLQ - RPDF...GGA, insulin MALW...YCN)
	#"KWCFRVCYRGICYRRCR"
	#"TTCCPSIVARSNFNVCRLPGTPSEALICATYTGCIIIPGATCPGDYAN"
	#"RPDFCLEPPYAGACRARIIRYFYNAKAGLCQTFVYGGCRAKRNNFKSAEDCLRTCGGA"
	#"HELIX    1   1 PRO A    2  GLU A    7  5                                   6\nHELIX    2   2 SER A   47  GLY A   56  1                                  10\nSHEET    1   A 2 ILE A  18  ASN A  24  0\nSHEET    2   A 2 LEU A  29  TYR A  35 -1  N  TYR A  35   O  ILE A  18"
	
	#if os.path.exists(os.path.join("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations", outname)):
	#	return

	permissions = AAPermissionsManager("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC 3/permissions", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC 3/permissible_sequences/default.txt")
	sec_struct_permissions = AASecondaryStructurePermissionsManager("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC 3/permissible_sequences")
	peptide = Polypeptide()
	#peptide.read("/Users/venkatesh-sivaraman/Downloads/1QLQ.pdb")
	if sec_structs:
		if ',' in sec_structs:
			peptide.add_secondary_structures(sec_structs, format='csv')
		else:
			peptide.add_secondary_structures(sec_structs, format='pdb')
	peptide.randomcoil(seq, permissions=permissions, struct_permissions=sec_struct_permissions)
	print peptide.secondary_structures
	peptide.center()
	file = open(os.path.join("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations", outname), 'w')
	scoresfile = open(os.path.join("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations", outname[:-4] + "-scores.txt"), 'w')
	#file.write(peptide.xyz(escaped=False))
	scores = []
	t_scores = ""
	curscore = 0.0
	for dist in dists:
		sc = dist.score(peptide, peptide.aminoacids)
		t_scores += str(sc) + " (" + dist.identifier + ") "
		curscore += sc
	scores.append(curscore)
	scoresfile.write("-1 " + t_scores + "\n")
	file.write(peptide.pdb(modelno=1))
	
	prob = AAProbabilitySource(peptide, dists, permissions) #sec_struct_permissions
	gentle_cutoff = 0
	best_models = [[] for i in xrange(model_count)]
	best_scores = [1000 for i in xrange(model_count)]
	pdb_model_idx = 2
	a = datetime.datetime.now()
	time_wasted = 0.0
	for i in xrange(n):
		seglen = segment_length(scores[-1] / len(peptide.aminoacids))
		folding.folding_iteration(peptide, prob, seglen)
		peptide.center()
		for aa in peptide.aminoacids:
			aa.localscore = sum(dist.score(peptide, [aa]) for dist in dists)
		t_scores = ""
		curscore = 0.0
		for dist in dists:
			sc = dist.score(peptide, peptide.aminoacids)
			t_scores += str(sc) + " (" + dist.identifier + ") "
			curscore += sc
		scores.append(curscore)
		#print i, t_scores
		scoresfile.write(str(i) + " " + t_scores + "\n")
		#print i, scores[-1] / len(peptide.aminoacids), dists[-1].score(peptide, peptide.aminoacids)
		#if scores[-1] / len(peptide.aminoacids) <= -100.0:
			#file.write(peptide.xyz(escaped=False))
		if curscore / len(peptide.aminoacids) <= 0.0:
			file.write(peptide.pdb(modelno=pdb_model_idx))
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
		currscore = sum(dist.score(peptide, peptide.aminoacids) for dist in dists)
		print "First has score", currscore
		while running_count < 50:
			seglen = segment_length(scores[-1] / len(peptide.aminoacids))
			folding.folding_iteration(peptide, prob, seglen)
			peptide.center()
			newscore = sum(dist.score(peptide, peptide.aminoacids) for dist in dists)
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
		file.write(peptide.pdb(modelno=pdb_model_idx))
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
	permissions = AAPermissionsManager("/Users/venkatesh-sivaraman/Desktop/sciencefair/allowed-zones")
	for pz in permissions.allowed_conformations(aa2, aa1):
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
	#generate_distance_constrained_bins("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC 3/reference_states.txt")
	#sparc_distance_constrained_bins("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC 3/consec+secondary", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC 3/consec+secondary-ref")
	#test_sparc("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/tasser-decoys", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/comparison"),
	#supplement_natives("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/casp")
	#tm_sparc_correlation("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/TM-scores", "/Volumes/External Hard Drive/Science Fair 2014-15/Decoy Output/casp_hard")
	#permissions = AAPermissionsManager("/Users/venkatesh-sivaraman/Desktop/sciencefair/allowed-zones")
	#func()
	#w = [9,4,3]
	'''print z_score("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/4state_reduced", w)
	print z_score("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/fisa_casp3", w)
	print z_score("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/lattice_ssfit", w)
	print z_score("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/lmds", w)
	print z_score("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/hg_structal", w)
	print z_score("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/ig_structal", w)
	print z_score("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/vhp_mcmd", w)'''
	#print z_score("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/tasser", w, structure_files="/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoys/tasser-decoys")
	#print z_score("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/casp", w, structure_files="/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoys/casp-decoys")
	#print determine_omits("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Nonredundant/all_pdb_ids.txt", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Nonredundant/omits.txt")
	'''best_weights("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/comparison", numweights=4, start=[3, 4, 5, 10])
	best_weights("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/comparison", numweights=4, start=7)
	best_weights("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/comparison", numweights=4, start=3)
	best_weights("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/comparison", numweights=3)'''
	#best_weights("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/rw/casp-rw", numweights=1)
	#best_weights("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/rw/tasser-rw", numweights=1)
	#best_weights("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/goap/casp-goap", numweights=1)
	#best_weights("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/goap/tasser-goap", numweights=1, start=1)
	#best_weights_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/tasser-decoys-new", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/TM-scores/tasser", None) #, "/Users/venkatesh-sivaraman/Downloads/casp11.targets_unsplitted.release11242014")
	#decoys_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/tasser-decoys", None, "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/rmsd") #/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoys/casp-decoys,
	#min1 = min_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments-test/seg12345678910.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments/native.pdb")
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
	'''mins = []
	for i in [12, 34, 56, 78, 910, 1234, 5678, 5678910]: #range(1, 11) +
		mins.append(min_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments-test/seg" + str(i) + ".pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments/native_seg" + str(i) + ".pdb"))
	print mins'''
	#min2 = min_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments-test/seg12345678910_back.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments/native.pdb")
	#link_decoy_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/casp", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/rmsd-casp-backbone", [9.0, 4.0, 3.0])
	#compare_charmm_sparc_sp("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/structure_pairs", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/comparison_results.txt") #, pairwise=False, idealized=False, minimize=True)
	#protein_protein_energies("/Users/venkatesh-sivaraman/Downloads/1DUM.pdb", load_dists(secondary=False), "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/magainin_test_yz_gromos.txt")
	#analysis()
	#average_coordination_number("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC/medium")
	#batch_compare_charmm_sparc("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/tasser-decoys", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/solvated_comparison.txt")
	#trim_secondary_structure_pzs("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/secondary_structures", fraction=0.5)
	#sequence_score_distribution("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/sequence_scores.txt", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Sequence Scores", weights=[3.0, 3.0, 2.0, 0.0, 0.0])
	#for i in xrange(3, 11):
	#i = 9
	#fit_maxwell("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Sequence Scores/" + str(i) + ".txt")
	'''basepath = "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments"
	for decoyset in os.listdir(basepath):
		if "segments" in basepath and ("native" in decoyset or not os.path.exists(os.path.join(basepath, "native_" + decoyset))): continue
		if decoyset[0] == ".": continue
		if os.path.exists(os.path.join("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/TM-scores/simul", decoyset + ".txt")): continue
		calculate_tm_scores(os.path.join(basepath, decoyset), os.path.join("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/TM-scores/simul", decoyset + ".txt"), natives=os.path.join(basepath, "native_" + decoyset)) #, natives="/Users/venkatesh-sivaraman/Downloads/casp11.targets_unsplitted.release11242014")'''