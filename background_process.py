from main import *
import gc

if __name__ == '__main__':
	'''root = "/Volumes/External Hard Drive/Science Fair 2014-15/tasser-decoys-2/"
	output = "/Volumes/External Hard Drive/Science Fair 2014-15/decoy-tm/"
	for struct in os.listdir(root):
		if not os.path.isdir(join(root, struct)): continue
		print struct
		calculate_tm_scores(join(root, struct), join(output, struct + ".txt"))'''
	#analyze_secondary_structure("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/all_pdb_ids.txt", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/secondary_structures", sample=1000)
	dists = load_dists(concurrent=True) #{"consec+secondary": 3.0, "short-range": 2.0, "nonconsec": 2.0, "medium": 0.0}
	#min1 = min_rmsd("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments-test/seg12-2.pdb", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments/native_seg12.pdb", dists=distributions)
	#score_structure_file("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments-test/seg12-3.pdb", distributions)
	sparc_score_sequences("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Nonredundant/all_pdb_ids_laptop.txt", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/sequence_scores_charmm.txt", dists, charmm_scores=True)
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
	#aggregate_permissible_sequences("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/permissible_sequences", cutoff=0.01)
	
	#correlate_charmm_sparc("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/magainin_test_y.txt")
	#correlate_charmm_sparc("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/magainin_test_yz.txt")
	#analysis()
	#generate_structurepairs("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/small_pdb_ids original.txt", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/structure_pairs")