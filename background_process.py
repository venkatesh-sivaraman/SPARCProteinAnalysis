from main import *
import gc

if __name__ == '__main__':
	'''root = "/Volumes/External Hard Drive/Science Fair 2014-15/tasser-decoys-2/"
	output = "/Volumes/External Hard Drive/Science Fair 2014-15/decoy-tm/"
	for struct in os.listdir(root):
		if not os.path.isdir(join(root, struct)): continue
		print struct
		calculate_tm_scores(join(root, struct), join(output, struct + ".txt"))'''
	analyze_secondary_structure("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/all_pdb_ids.txt", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/secondary_structures", sample=1000)
	#calculate_pdb_stats_network("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/all_pdb_ids_laptop.txt", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/block_data")
	#print total, omit, partial
	#correlate_charmm_sparc("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/magainin_test_xz.txt")
	#correlate_charmm_sparc("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/magainin_test_yz.txt")
	#analysis()
	#generate_structurepairs("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/small_pdb_ids original.txt", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/structure_pairs")
	'''basepath = "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/SPARC"
	for path in os.listdir(os.path.join(basepath, "consec")):
		if ".txt" not in path: continue
		consecutive_allowed_alpha_zones(os.path.join(basepath, "consec", path), os.path.join(basepath, "permissions", path))'''