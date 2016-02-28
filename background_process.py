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
	'''dists = load_dists({"consec+secondary": 4.0, "short-range": 1.0, "nonconsec": 4.0, "medium": 5.0}, concurrent=True, secondary=False)
	distributions = ["", "", "", ""]
	for dist in dists:
		if dist.identifier == "consec+secondary": distributions[0] = dist
		elif dist.identifier == "short-range": distributions[1] = dist
		elif dist.identifier == "nonconsec": distributions[2] = dist
		elif "medium" in dist.identifier: distributions[3] = dist'''
	#sparc_score_sequences("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Nonredundant/all_pdb_ids.txt", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/sequence_scores.txt", distributions)
	'''for path in os.listdir("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments"):
		if "native" in path or "pdb" not in path: continue
		print path, "======="
		score_structure_file(os.path.join("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Simulations/segments", path), distributions)'''
	'''for x in xrange(1, 4):
		for y in xrange(1, 4):
			for z in xrange(1, 4):
				print "Testing weights", x, y, z
				func(weights={ "consec": x, "secondary": x, "short-range": y, "nonconsec": 0.0, "medium": z }, base="segment_weight_test/" + str(x) + str(y) + str(z))'''
	#calculate_pdb_stats_network("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Nonredundant/all_pdb_ids.txt", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/expanded_data")
	#generate_permissible_sequences("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Nonredundant/all_pdb_ids.txt", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/permissible_sequences")
	aggregate_permissible_sequences("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/permissible_sequences", cutoff=0.1)
	#aggregate_networkdata("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/expanded_data", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Nonredundant/Combined C&S")
	#print total, omit, partial
	#correlate_charmm_sparc("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/magainin_test_xz.txt")
	#correlate_charmm_sparc("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/magainin_test_yz.txt")
	#analysis()
	#generate_structurepairs("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/small_pdb_ids original.txt", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/structure_pairs")
	'''basepath = "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC 2"
	for path in os.listdir(os.path.join(basepath, "secondary")):
		if ".txt" not in path: continue
		consecutive_allowed_alpha_zones(os.path.join(basepath, "secondary", path), os.path.join(basepath, "permissions-secondary", path))'''