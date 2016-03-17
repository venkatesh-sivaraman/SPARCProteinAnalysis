from main import *
import os, sys

if __name__ == '__main__':
	args = sys.argv[1:]
	input = None
	sparc_dir = None
	i = 0
	while i < len(args):
		if args[i].lower() == "-i":
			assert len(args) > i + 1, "Not enough arguments"
			input = args[i + 1]
			i += 2
		elif args[i].lower() == "-o":
			assert len(args) > i + 1, "Not enough arguments"
			sparc_dir = args[i + 1]
			i += 2
		else:
			assert False, "Unexpected command-line argument {}".format(args[i])

	data_dir = os.path.join(os.path.dirname(input), "sparc_data")
	print "Computing data from PDB structures..."
	calculate_pdb_stats_network(input, data_dir)
	print "Aggregating orientation and solvent data..."
	aggregate_networkdata(data_dir, sparc_dir)
	print "Formatting data..."
	for path in os.listdir(sparc_dir):
		if os.path.isdir(os.path.join(sparc_dir, path)) and path != "medium" and path != "permissible_sequences":
			write_median_frequencies(os.path.join(sparc_dir, path))
	format_hydrophobicity_dist(os.path.join(sparc_dir, "medium"), os.path.join(sparc_dir, "medium-2"))
	os.system("rm -rf \"{}\"".format(os.path.join(sparc_dir, "medium")))
	os.rename(os.path.join(sparc_dir, "medium-2"), os.path.join(sparc_dir, "medium"))
	#generate_permissible_sequences("/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/Nonredundant/all_pdb_ids.txt", "/Users/venkatesh-sivaraman/Documents/School/Science Fair/2016-proteins/permissible_sequences")
	print "Aggregating permissible sequence data..."
	aggregate_permissible_sequences(data_dir, os.path.join(sparc_dir, "permissible_sequences"), cutoff=0.01)
