"""This module provides for testing SPARC on decoy sets. Running this module as __main__ supports decoy sets of the following type:
	INPUT
	--decoyset
	----decoy
	----decoy
	----decoy
	...
	--decoyset
	----decoy
	----decoy
	...
	To perform this decoy test, pass -i with a path to the INPUT directory and -o to a directory where the scores will be deposited. The native structure may be one of the decoys, or you may provide a parameter -n for the directory of native structures.
	To compute scores with the old version of SPARC, pass the old SPARC path after the argument -old.
	-s specifies the directory of SPARC, if not the default potential residing in the script directory."""

import sys, os
from main import *

def process_decoys_file((input, output, sparc_dir, nativepath, old)):
	if not os.path.isdir(input) or os.path.exists(output) or os.path.basename(input) == "doc": return
	protein_name = os.path.basename(input)
	print protein_name
	
	if old:
		dists_old = load_dists(old, concurrent=False, secondary=False)
		for d in dists_old: d.refstate = False
		dists_noref = load_dists(sparc_dir, concurrent=False, secondary=False)
		for d in dists_noref: d.refstate = False
		dists_yesref = load_dists(sparc_dir, concurrent=False, secondary=False)
		for d in dists_yesref: d.refstate = True
		distributions = dists_old + dists_noref + dists_yesref
	else:
		distributions = load_dists(sparc_dir, concurrent=False, secondary=False)

	paths = os.listdir(input)
	allpaths = [os.path.join(input, path) for path in paths]

	if "-" in protein_name:
		path = protein_name[:protein_name.find("-")] + ".pdb"
	else:
		path = protein_name + ".pdb"
	if nativepath: allpaths.append(os.path.join(nativepath, path))

	nativescores = None
	bounds = None
	peptide = Polypeptide()
	if nativepath:
		if os.path.exists(join(nativepath, path)):
			bounds, gaps, scores = sparc_scores_file(join(nativepath, path), distributions, retbounds=True, peptide=peptide) #, ignored_aas=gaps
			print path
			nativescores = scores
			if output and scores is not None:
				scorestr = ""
				for s in scores: scorestr += str(s) + ","
				scorestr = scorestr[:-1]
				with open(output, "w") as file:
					file.write("{}; {}\n".format(protein_name + "_orig.pdb", scorestr))
		else:
			print join(nativepath, path), "does not exist."
	else:
		scores = None
	for path in paths:
		if path == "list" or path == "rmsds": continue
		if True: #try:
			scores = sparc_scores_file(join(input, path), distributions, bounds=bounds, peptide=peptide) #, ignored_aas=gaps
		'''except Exception as e:
			print path, "exception ({})".format(e)
			continue'''
		if output and scores is not None:
			scorestr = ""
			for s in scores: scorestr += str(s) + ","
			scorestr = scorestr[:-1]
			with open(output, "a") as file:
				file.write("{}; {}\n".format(path, scorestr))
		del scores
		gc.collect()
	ret = ""
	print "Done"
	del paths
	del distributions[:]
	del distributions
	gc.collect()

def test_sparc(input, output, sparc_dir, natives=None, old=None):
	files = os.listdir(input)
	if not os.path.exists(output):
		os.mkdir(output)
	print len(files), "files"
	pool = multiprocessing.Pool(processes=2, maxtasksperchild=1)
	zipped = [(join(input, file), join(output, file + ".txt"), sparc_dir, natives, old) for file in files]
	map(process_decoys_file, zipped)
	pool.close()
	pool.join()
	print "done"

if __name__ == '__main__':
	args = sys.argv[1:]
	input = None
	output = None
	natives = None
	old = None
	sparc_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "potential")
	i = 0
	while i < len(args):
		if args[i].lower() == "-old":
			assert len(args) > i + 1, "Not enough arguments"
			old = args[i + 1]
			i += 2
		elif args[i].lower() == "-i":
			assert len(args) > i + 1, "Not enough arguments"
			input = args[i + 1]
			i += 2
		elif args[i].lower() == "-o":
			assert len(args) > i + 1, "Not enough arguments"
			output = args[i + 1]
			i += 2
		elif args[i].lower() == "-n":
			assert len(args) > i + 1, "Not enough arguments"
			natives = args[i + 1]
			i += 2
		elif args[i].lower() == "-s":
			assert len(args) > i + 1, "Not enough arguments"
			sparc_dir = args[i + 1]
			i += 2
		else:
			assert False, "Unexpected command-line argument {}".format(args[i])
	test_sparc(input, output, sparc_dir, natives, old)