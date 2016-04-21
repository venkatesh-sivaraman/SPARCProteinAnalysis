"""This module contains capabilities for evaluating models from SPARC simulations. Pass the word 'rmsd' or the word 'tmscore' to calculate that quantity.
	The -i parameter specifies either a file containing models or a directory in which all .pdb files will be evaluated.
	The -n parameter specifies a native structure file.
	The -d parameter specifies the path to the directive used by SPARC to generate the models. This parameter is necessary to examine segments of the native structure file; the files will be identified by the output names in the directive file."""

from main import *
from tmscore import *

rmsd_mode = "rmsd"
tmscore_mode = "tmscore"
sparc_mode = "sparc"

def sparc_validation(input, native, range, output=None, distributions=None):
	backup_sec_structs = []
	if output:
		outfile = open(output, "w")
	for i, path in enumerate([native, input]):
		peptide = Polypeptide()
		readrange = None
		if i == 0: readrange = range
		for modelno in peptide.iter_models(path, secondary_structure=True, range=readrange):
			if i == 0:
				print peptide.aminoacids, peptide.secondary_structures
				modelno = "model_orig.pdb"
			if len(peptide.secondary_structures):
				backup_sec_structs = peptide.secondary_structures
			else:
				peptide.secondary_structures = backup_sec_structs
		
			spscore = ""
			for d in distributions:
				subscore = d.score(peptide, peptide.aminoacids)
				spscore += "{:.4f},".format(subscore)
			spscore = spscore[:-1]
			if output:
				outfile.write("{};{}\n".format(modelno, spscore))
			else:
				print "{}\t{}".format(modelno, spscore)
	if output:
		outfile.close()

if __name__ == '__main__':
	args = sys.argv[1:]
	input = None
	native = None
	directive = None
	output = None
	mode = args[0]
	i = 1
	while i < len(args):
		if args[i].lower() == "-n":
			assert len(args) > i + 1, "Not enough arguments"
			native = args[i + 1]
			i += 2
		elif args[i].lower() == "-i":
			assert len(args) > i + 1, "Not enough arguments"
			input = args[i + 1]
			i += 2
		elif args[i].lower() == "-d":
			assert len(args) > i + 1, "Not enough arguments"
			directive = args[i + 1]
			i += 2
		elif args[i].lower() == "-o":
			assert len(args) > i + 1, "Not enough arguments"
			output = args[i + 1]
			i += 2
		else:
			assert False, "Unexpected command-line argument {}".format(args[i])

	input = os.path.realpath(input)
	native = os.path.realpath(native)

	distributions = None
	if mode == sparc_mode:
		sparc_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "potential")
		weights = {} #{ "consec": 4.0, "secondary": 4.0, "short_range": 1.0, "long_range": 4.0, "medium": 5.0 }
		distributions = load_dists(sparc_dir, secondary=True, weights=weights)

	if directive:
		seq = None
		runs = []
		with open(directive, "r") as file:
			lines = file.readlines()
			seq = lines[0].strip()
			del lines[0]
			processing_runs = False
			for line in lines:
				if len(line.strip()) == 0:
					if not processing_runs:
						processing_runs = True
					continue
				if processing_runs:
					if line[0] == "#": line = line[1:]
					comps = line.strip().split(";")
					runs.append([[y.strip() for y in x.split(",")] for x in comps])
		completed_outputs = []
		for run in runs:
			modelname = None
			if len(run[1]) > 1: # run[1] must be the input paths, and run[2] must be the output path name
				modelname = run[2][0]
				range = [[int(x) for x in run[0][0].split("-")][0], [int(x) for x in run[0][1].split("-")][1]]
			else:				# run[1] must be the output path name
				modelname = run[1][0]
				range = [int(x) for x in run[0][0].split("-")]
			if output:
				if os.path.exists(os.path.join(output, modelname[:-4] + ".txt")): continue
			if modelname in completed_outputs: continue
			if not os.path.exists(os.path.join(input, modelname)):
				print modelname, "does not exist, skipping"
				continue
			if mode == rmsd_mode:
				if output:
					print modelname, "-->", modelname[:-4] + ".txt"
					min_rmsd(os.path.join(input, modelname), native, range=range, output=os.path.join(output, modelname[:-4] + ".txt"))
				else:
					print modelname
					min_rmsd(os.path.join(input, modelname), native, range=range)
			elif mode == tmscore_mode:
				print modelname
				best_tm_score(os.path.join(input, modelname), native, range=range)
			elif mode == sparc_mode:
				if output:
					print modelname, "-->", modelname[:-4] + ".txt"
					sparc_validation(os.path.join(input, modelname), native, range, os.path.join(output, modelname[:-4] + ".txt"), distributions=distributions)
				else:
					print modelname
					sparc_validation(os.path.join(input, modelname), native, range, distributions=distributions)
			completed_outputs.append(modelname)
	else:
		if os.path.isdir(input):
			for path in os.listdir(input):
				if ".pdb" not in path or path[0] == "." or "native" in path: continue
				print path
				if mode == rmsd_mode:
					min_rmsd(os.path.join(input, path), native)
				elif mode == tmscore_mode:
					best_tm_score(os.path.join(input, path), native)
		else:
			if mode == rmsd_mode:
				min_rmsd(input, native)
			elif mode == tmscore_mode:
				best_tm_score(input, native)
