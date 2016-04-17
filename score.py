from main import *

if __name__ == '__main__':
	args = sys.argv[1:]
	terms = False
	sec_struct = False
	sparc_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "potential")
	inputs = []
	start = 0
	end = 0
	i = 0
	native = None
	while i < len(args):
		if args[i].lower() == "-terms":
			terms = True
			i += 1
		elif args[i].lower() == "-ss":
			sec_struct = True
			i += 1
		elif args[i].lower() == "-r":
			assert len(args) > i + 2, "Missing one or more range values"
			start = int(args[i + 1])
			end = int(args[i + 2])
			i += 3
		elif args[i].lower() == "-s":
			assert len(args) > i + 1, "Not enough arguments"
			sparc_dir = args[i + 1]
			i += 2
		elif args[i].lower() == "-n":
			assert len(args) > i + 1, "Not enough arguments"
			native = args[i + 1]
			i += 2
		else:
			inputs.append(args[i])
			i += 1
	
	distributions = load_dists(sparc_dir, secondary=True) #load_central_dist(sparc_dir, secondary=True)
	if terms:
		print "\nModel\tConsec\tSec\tShort\tLong\tSolvent\tTotal\n========================================================"
	else:
		print "\nModel\tSPARC Score\n========================="
	backup_sec_structs = []
	if native:
		inputs.insert(native, 0)
	for input in inputs:
		peptide = Polypeptide()
		inrange = None
		if input == native: inrange = range
		for modelno in peptide.iter_models(input, secondary_structure=True, range=inrange):
			print peptide.aminoacids
			if len(peptide.secondary_structures):
				backup_sec_structs = peptide.secondary_structures
			elif sec_struct:
				peptide.secondary_structures = backup_sec_structs
			
			if terms:
				spscore = ""
				sumscore = 0.0
				for d in distributions:
					subscore = d.score(peptide, peptide.aminoacids)
					sumscore += subscore
					spscore += "{:.2f}\t".format(subscore)
				print "{}\t{}{}".format(modelno, spscore, sumscore)
			else:
				print "{}\t{:.2f}".format(modelno, sum(d.score(peptide, peptide.aminoacids) for d in distributions))