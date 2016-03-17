from main import *
import gc

if __name__ == '__main__':
	'''root = "/Volumes/External Hard Drive/Science Fair 2014-15/tasser-decoys-2/"
	output = "/Volumes/External Hard Drive/Science Fair 2014-15/decoy-tm/"
	for struct in os.listdir(root):
		if not os.path.isdir(join(root, struct)): continue
		print struct
		calculate_tm_scores(join(root, struct), join(output, struct + ".txt"))'''
	#analyze_secondary_structure("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Nonredundant/all_pdb_ids.txt", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/secondary_structures")
	#dists = load_dists(concurrent=True, secondary=False, weights={"consec+secondary": 4.0, "short-range": 1.0, "nonconsec": 4.0, "medium": 5.0}) #{"consec": 3.0, "secondary": 3.0, "short-range": 2.0, "nonconsec": 2.0, "medium": 0.0}, {"consec+secondary": 3.0, "short-range": 2.0, "nonconsec": 2.0, "medium": 0.0}
	#mins = []
	#BPTI
	#sec_structs = ["helix,5,2,7", None, None, "sheet,0,1,7", None, "sheet,0,1,7", None, None, "helix,1,1,5", "helix,1,1,5"]
	#Insulin?
	'''for i in [5678, 5678910]: #range(1,11) + [12, 34, 56, 78, 910, 1234, 5678, 5678910]:
		mins.append(min_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments-test/seg" + str(i) + ".pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments/native_seg" + str(i) + ".pdb", dists=dists)) #sec_structs=sec_structs[i - 1]
	mins.append(min_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments-test/seg12345678910_back.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments/native.pdb", dists=dists))
	for m in mins:
		print m'''
	#batch_compare_charmm_sparc("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/casp-decoys", dists, "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/casp_charmm.txt")
	#score_structure_file("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments-test/seg12-3.pdb", distributions)
	#sparc_score_sequences("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Nonredundant/all_pdb_ids.txt", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/sequence_scores_charmm.txt", dists, charmm_scores=False)
	'''for path in os.listdir("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments"):
		if "native" in path or "pdb" not in path: continue
		print path, "======="
		score_structure_file(os.path.join("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments", path), distributions)'''
	'''for x in xrange(1, 4):
		for y in xrange(1, 4):
			for z in xrange(1, 4):
				print "Testing weights", x, y, z
				func(weights={ "consec": x, "secondary": x, "short-range": y, "nonconsec": 0.0, "medium": z }, base="segment_weight_test_sec/" + str(x) + str(y) + str(z))'''
	
	#Methods needed to calculate SPARC
	sparc_dir = "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC 3"
	#calculate_pdb_stats_network("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Nonredundant/all_pdb_ids.txt", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/expanded_data")
	#aggregate_networkdata("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/expanded_data", sparc_dir)
	#for path in os.listdir(sparc_dir):
	#	if os.path.isdir(os.path.join(sparc_dir, path)) and path != "medium":
	#		write_median_frequencies(os.path.join(sparc_dir, path))
	#format_hydrophobicity_dist(os.path.join(sparc_dir, "medium"), os.path.join(sparc_dir, "medium-2"))
	#os.system("rm -rf \"{}\"".format(os.path.join(sparc_dir, "medium")))
	#os.rename(os.path.join(sparc_dir, "medium-2"), os.path.join(sparc_dir, "medium"))
	#generate_permissible_sequences("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Nonredundant/all_pdb_ids.txt", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/permissible_sequences")
	aggregate_permissible_sequences("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/permissible_sequences", cutoff=0.005)
	
	#correlate_charmm_sparc("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/magainin_test_y.txt")
	#correlate_charmm_sparc("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/magainin_test_yz.txt")
	#analysis()
	#generate_structurepairs("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/small_pdb_ids original.txt", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/structure_pairs")