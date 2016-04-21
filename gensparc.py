from main import *
import os, sys

if __name__ == '__main__':
	args = sys.argv[1:]
	input = None
	sparc_dir = None
	both = False
	calculate_interactions = False
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
		elif args[i].lower() == "-b":
			both = True
			calculate_interactions = False
			i += 1
		elif args[i].lower() == "-pi":
			calculate_interactions = True
			both = False
			i += 1
		else:
			assert False, "Unexpected command-line argument {}".format(args[i])

	data_dir = os.path.join(os.path.dirname(input), "sparc_data")

	#original_potential = os.path.join(os.path.dirname(os.path.realpath(__file__)), "potential")
	#permissions = AAPermissionsManager(os.path.join(original_potential, "permissions"), os.path.join(original_potential, "permissible_sequences", "all.txt"))
	#sec_struct_permissions = AASecondaryStructurePermissionsManager(os.path.join(original_potential, "permissible_sequences"))
	#cProfile.runctx('pdb_reference_state(input, data_dir, permissions, sec_struct_permissions)', {'pdb_reference_state': pdb_reference_state, 'input': input, 'data_dir': data_dir, 'permissions': permissions, 'sec_struct_permissions': sec_struct_permissions}, {})
	#pdb_reference_state(input, data_dir, permissions, sec_struct_permissions)
	#aggregate_consolidated_data(data_dir, sparc_dir)

	if calculate_interactions:
		print "Calculating the number of possible interactions between every amino acid type."
		num_possible_interacting(input, data_dir)
		aggregate_possible_interactions(data_dir, sparc_dir)
	else:
		'''print "Computing data from PDB structures..."
		if both:
			calculate_pdb_stats_network(input, data_dir, mode=both_orientations_mode)
		else:
			calculate_pdb_stats_network(input, data_dir, mode=default_network_mode)
		print "Aggregating orientation and solvent data..."
		#aggregate_consolidated_data(data_dir, sparc_dir)
		aggregate_networkdata(data_dir, sparc_dir, both=both)
		print "Formatting data..."'''
		for path in os.listdir(sparc_dir):
			if os.path.isdir(os.path.join(sparc_dir, path)) and path != "medium" and "permissible_sequences" not in path and path != "possible_interactions":
				print path
				write_median_frequencies(os.path.join(sparc_dir, path))
		'''if not both:
			format_hydrophobicity_dist(os.path.join(sparc_dir, "medium"), os.path.join(sparc_dir, "medium-2"))
			os.system("rm -rf \"{}\"".format(os.path.join(sparc_dir, "medium")))
			os.rename(os.path.join(sparc_dir, "medium-2"), os.path.join(sparc_dir, "medium"))
			#seq_dir = os.path.join(os.path.dirname(input), "sparc_permissible_sequences")
			#generate_permissible_sequences(input, seq_dir)
			print "Aggregating permissible sequence data..."
			aggregate_permissible_sequences(data_dir, os.path.join(sparc_dir, "permissible_sequences"), cutoff=0.01)
		print "Finished generating SPARC."'''