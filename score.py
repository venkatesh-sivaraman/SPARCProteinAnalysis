from main import *

if __name__ == '__main__':
	args = sys.argv[1:]
	terms = False
	sec_struct = False
	inputs = []
	start = 0
	end = 0
	i = 0
	while i < len(args):
		if args[i].lower() == "-terms":
			terms = True
			i += 1
		elif args[i].lower() == "-ss":
			sec_struct = True
			i += 1
		elif args[i].lower() = "-r":
			assert len(args) > i + 2, "Missing one or more range values"
			start = int(args[i + 1])
			end = int(args[i + 2])
			i += 3
		else:
			inputs.append(args[i])
			i += 1
	
	sparc_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "potential")
	distributions = load_dists(sparc_dir, secondary=True)
	if terms:
		print "Model\tConsec\tSec\tShort\tLong\tSolvent\tTotal\n========================================================"
	else:
		print "Model\tSPARC Score\n========================="
	backup_sec_structs = []
	for input in inputs:
		peptide = Polypeptide()
		for modelno in peptide.iter_models(input, secondary_structure=True):
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
				print "{}\t{}".format(modelno, sum(d.score(peptide, peptide.aminoacids) for d in distributions))