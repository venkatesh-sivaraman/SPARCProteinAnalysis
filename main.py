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
import datetime, time
from tmscore import *
from charmm import *
from sparc_distribution import *

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
	return numpy.random.choice(range(1, 6), p=weights)

def load_dists(weights={frequency_nonconsec_disttype: 1.0, frequency_consec_disttype: 1.0, medium_disttype: 1.0}, concurrent=True):
	nonconsec = "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC/nonconsec"
	consec = "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC/consec"
	medium = "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC/medium"
	if concurrent == False:
		dist1 = MediumDistributionManager(medium)
		dist1.weight = weights[medium_disttype]
		dist2 = SPARCBasicDistributionManager(consec, True)
		dist2.weight = weights[frequency_consec_disttype]
		dist3 = SPARCBasicDistributionManager(nonconsec, False)
		dist3.weight = weights[frequency_nonconsec_disttype]
		return [dist1, dist2, dist3]
	processes = []
	queue = multiprocessing.Queue()
	def generate_distmanager(cls, q, *args):
		dist = cls(*args)
		dist.weight = weights[dist.type]
		q.put(dist)
	p1 = multiprocessing.Process(target=generate_distmanager, args=(SPARCBasicDistributionManager, queue, nonconsec, False))
	processes.append(p1)
	p1.start()
	p2 = multiprocessing.Process(target=generate_distmanager, args=(SPARCBasicDistributionManager, queue, consec, True))
	processes.append(p2)
	p2.start()
	p3 = multiprocessing.Process(target=generate_distmanager, args=(MediumDistributionManager, queue, medium))
	processes.append(p3)
	p3.start()
	dists = []
	dists.append(queue.get())
	dists.append(queue.get())
	dists.append(queue.get())
	p1.join()
	p2.join()
	p3.join()
	return dists

def func():
	#Weights used to be 2, 4, 8
	dists = load_dists(weights={frequency_nonconsec_disttype: 9.0, frequency_consec_disttype: 4.0, medium_disttype: 3.0})
	#cProfile.runctx('simulate_fold(dists)', {'dists': dists, 'simulate_fold': simulate_fold}, {})
	simulate_fold(dists)
	#test_folding_parameters(dists)

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
	distributions = load_dists(concurrent=False)
	dists = ["", "", ""]
	for dist in distributions:
		if dist.type == frequency_nonconsec_disttype: dists[0] = dist
		elif dist.type == frequency_consec_disttype: dists[1] = dist
		elif dist.type == medium_disttype: dists[2] = dist
	paths = os.listdir(input)

	if "-" in protein_name:
		path = protein_name[:protein_name.find("-")] + ".pdb"
	else:
		path = protein_name + ".pdb"
	nativescores = None
	if nativepath:
		if os.path.exists(join(nativepath, path)):
			bounds, scores = sparc_scores_file(join(nativepath, path), dists, retbounds=True)
			nativescores = scores
			if output and scores is not None:
				with open(output, "w") as file:
					file.write("{}; {}, {}, {}\n".format(protein_name + "_orig.pdb", scores[0], scores[1], scores[2]))
		else:
			print join(nativepath, path), "does not exist."
	else:
		bounds = None
		scores = None
	for path in paths:
		try:
			scores = sparc_scores_file(join(input, path), dists)
		except:
			print path, "exception"
			continue
		if output and scores is not None:
			with open(output, "a") as file:
				file.write("{}; {}, {}, {}\n".format(path, scores[0], scores[1], scores[2]))
		del scores
		gc.collect()
	ret = ""
	print "Done"
	del dists
	del paths


def decoy_initializer():
	print "Starting"

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
	pool = multiprocessing.Pool(processes=2, initializer=decoy_initializer, maxtasksperchild=1)
	zipped = [(join(input, file), join(output, file + ".txt"), "/Users/venkatesh-sivaraman/Downloads/casp11.targets_unsplitted.release11242014") for file in files]
	#print zipped
	pool.map(process_decoys_file, zipped)
	pool.close()
	pool.join()
	print "done"


def simulate_fold(dists):
	permissions = AAPermissionsManager("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC/permissions")
	sec_struct_permissions = AASecondaryStructurePermissionsManager("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC/secondary")
	peptide = Polypeptide()
	#peptide.read("/Users/venkatesh-sivaraman/Downloads/1QLQ.pdb")
	peptide.add_secondary_structures("HELIX    1   1 PRO A    2  GLU A    7  5                                   6\nHELIX    2   2 SER A   47  GLY A   56  1                                  10\nSHEET    1   A 2 ILE A  18  ASN A  24  0\nSHEET    2   A 2 LEU A  29  TYR A  35 -1  N  TYR A  35   O  ILE A  18", format='pdb')
	peptide.randomcoil("RPDFCLEPPYAGACRARIIRYFYNAKAGLCQTFVYGGCRAKRNNFKSAEDCLRTCGGA", permissions=permissions, struct_permissions=sec_struct_permissions) #"GRYRRCIPGMFRAYCYMD" (2LWT - GRY...MD, 2MDB - KWC...CR, 1QLQ - RPDF...GGA)
	print peptide.secondary_structures
	#"KWCFRVCYRGICYRRCR"
	#"TTCCPSIVARSNFNVCRLPGTPSEALICATYTGCIIIPGATCPGDYAN"
	peptide.center()
	file = open("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/simulation_test.pdb", 'w')
	#file.write(peptide.xyz(escaped=False))
	scores = [sum(dist.score(peptide, peptide.aminoacids) for dist in dists)]
	print scores[-1], "Avg:", scores[-1] / len(peptide.aminoacids)
	file.write(peptide.pdb(modelno=1))
	
	prob = AAProbabilitySource(peptide, dists, permissions)
	gentle_cutoff = 0
	model_count = 5
	best_models = [[] for i in xrange(model_count)]
	best_scores = [1000 for i in xrange(model_count)]
	pdb_model_idx = 2
	a = datetime.datetime.now()
	for i in xrange(2500):
		seglen = segment_length(scores[-1] / len(peptide.aminoacids))
		folding.folding_iteration(peptide, prob, seglen)
		center = Point3D.zero()
		for aa in peptide.aminoacids:
			center = center.add(aa.acarbon)
		center = center.multiply(1.0 / len(peptide.aminoacids))
		for aa in peptide.aminoacids:
			aa.acarbon = aa.acarbon.subtract(center)
		curscore = 0.0
		for aa in peptide.aminoacids:
			aa.localscore = sum(dist.score(peptide, [aa]) for dist in dists)
			curscore += aa.localscore
		scores.append(curscore)
		print i, scores[-1] / len(peptide.aminoacids)
		if scores[-1] / len(peptide.aminoacids) <= -100.0:
			#file.write(peptide.xyz(escaped=False))
			file.write(peptide.pdb(modelno=pdb_model_idx))
			pdb_model_idx += 1
			#print peptide.xyz(escaped=False)
		#Save the conformation if it is the best so far.
		for k in xrange(model_count):
			if scores[-1] == best_scores[k]:
				break
			if scores[-1] < best_scores[k]:
				for m in reversed(xrange(max(k, 1), model_count)):
					best_scores[m] = best_scores[m - 1]
					best_models[m] = best_models[m - 1]
				best_scores[k] = scores[-1]
				best_models[k] = [PositionZone(aa.acarbon, aa.i, aa.j, aa.k) for aa in peptide.aminoacids]
				break
		if i == 9:
			#Save the best of the first ten scores as a comparison for later iterations.
			gentle_cutoff = min(scores)
		elif i > 9:
			if scores[-1] < gentle_cutoff and prob.mode == psource_erratic_mode:
				prob.mode = psource_gentle_mode
				print "Switching to gentle mode"
			elif scores[-1] >= gentle_cutoff and prob.mode == psource_gentle_mode:
				prob.mode = psource_erratic_mode
				print "Switching to erratic mode"
	b = datetime.datetime.now()
	print "Finished 1000 iterations in {0:.4f} sec.".format((b - a).total_seconds())
	del peptide.aminoacids[:]
	peptide = None
	del scores
	gc.collect()

	print "Iterating over the best scores:", best_scores
	prob.mode = psource_gentle_mode
	new_best_models = []
	new_best_scores = [1000 for n in xrange(model_count)]
	for k, model in enumerate(best_models):
		print "Refining model {}".format(k)
		for i, aa in enumerate(peptide.aminoacids):
			aa.acarbon = model[i].alpha_zone
			aa.set_axes(model[i].x_axis, model[i].y_axis, model[i].z_axis)
		#Keep track of how long the model has had this score.
		running_count = 0
		last_score = 0.0
		new_best_models.append([])
		while running_count < 100:
			seglen = segment_length(scores[-1] / len(peptide.aminoacids))
			folding.folding_iteration(peptide, prob, seglen)
			center = Point3D.zero()
			for aa in peptide.aminoacids:
				center = center.add(aa.acarbon)
			center = center.multiply(1.0 / len(peptide.aminoacids))
			for aa in peptide.aminoacids:
				aa.acarbon = aa.acarbon.subtract(center)
			newscore = sum(dist.score(peptide, peptide.aminoacids) for dist in dists)
			if last_score == newscore:
				running_count += 1
			else:
				running_count = 0
				last_score = newscore
			print "{} ({})".format(newscore, running_count)
			if newscore < new_best_scores[k]:
				new_best_scores[k] = newscore
				new_best_models[k] = [PositionZone(aa.acarbon, aa.i, aa.j, aa.k) for aa in peptide.aminoacids]
			if newscore > best_scores[k] + 10.0:
				break
	print "\nFinal:"
	for i in xrange(5):
		file.write("MODEL        %d\nENDMDL\n".format(pdb_model_idx))
		pdb_model_idx += 1
	for model in new_best_models:
		for i, aa in enumerate(peptide.aminoacids):
			aa.acarbon = model[i].alpha_zone
			aa.set_axes(model[i].x_axis, model[i].y_axis, model[i].z_axis)
		file.write(peptide.pdb(modelno=pdb_model_idx))
		pdb_model_idx += 1
	file.close()
	print "\n\n"
	print "Best model scores:"
	for i in xrange(model_count):
		print "{} ({})".format(i, new_best_scores[i])

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


if __name__ == '__main__':
	#test_sparc("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoys/casp-decoys", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/casp"),
	#supplement_natives("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/casp")
	#tm_sparc_correlation("/Volumes/External Hard Drive/Science Fair 2014-15/decoy-tm", "/Volumes/External Hard Drive/Science Fair 2014-15/decoy-scores")
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
	#best_weights("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output")
	#decoys_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoys/casp-decoys", "/Users/venkatesh-sivaraman/Downloads/casp11.targets_unsplitted.release11242014", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/rmsd-casp")
	#min_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/simulation_long.pdb", "/Users/venkatesh-sivaraman/Downloads/1QLQ.pdb")
	#link_decoy_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Decoy Output/casp", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/rmsd-casp-backbone", [9.0, 4.0, 3.0])
	#compare_charmm_sparc_sp("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/structure_pairs", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/comparison_results.txt") #, pairwise=False, idealized=False, minimize=True)
	#protein_protein_energies("/Users/venkatesh-sivaraman/Downloads/1DUM.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/magainin_test_yz.txt")
	#analysis()
	#average_coordination_number("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC/medium")
	#batch_compare_charmm_sparc("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/tasser-decoys", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/solvated_comparison.txt")
	#trim_secondary_structure_pzs("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/secondary_structures", fraction=0.5)