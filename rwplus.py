import os, subprocess
from os.path import join
import multiprocessing

def process_decoys_file((input, output, nativepath)):
	if not os.path.isdir(input): return
	if os.path.exists(output): return
	if os.path.basename(input) == "doc": return
	protein_name = os.path.basename(input)
	print protein_name
	paths = os.listdir(input)
	allpaths = [os.path.join(input, path) for path in paths]

	if "-" in protein_name:
		path = protein_name[:protein_name.find("-")] + ".pdb"
	else:
		path = protein_name + ".pdb"
	if nativepath: allpaths.append(os.path.join(nativepath, path))

	if nativepath:
		if os.path.exists(join(nativepath, path)):
			out = subprocess.check_output("./calRWplus/calRWplus \"" + join(nativepath, path) + "\"", shell=True, stderr=fnull)
			score = out[out.find("=") + 1:].strip().split()[0]
			print "native scores:", score
			if output and score != 0.0:
				with open(output, "w") as file:
					file.write("{}; {}\n".format(protein_name + "_orig.pdb", score))
		else:
			print join(nativepath, path), "does not exist."
	for path in paths:
		if path == "list" or path == "rmsds": continue
		try:
			out = subprocess.check_output("./calRWplus/calRWplus \"" + join(input, path) + "\"", shell=True, stderr=fnull)
			score = out[out.find("=") + 1:].strip().split()[0]
		except:
			print path, "exception"
			continue
		if output and score != 0.0:
			with open(output, "a") as file:
				file.write("{}; {}\n".format(path, score))
		gc.collect()
	ret = ""
	print "Done"
	del paths
	del distributions
	gc.collect()


def test_sparc(input, output):
	files = os.listdir(input)
	print len(files), "files"
	if "casp11_seqs.txt" in files:
		files.remove("casp11_seqs.txt")
	pool = multiprocessing.Pool(processes=2, maxtasksperchild=1)
	zipped = [(join(input, file), join(output, file + ".txt"), "~/Documents/casp11.targets_unsplitted.release11242014") for file in reversed(files)]
	#print zipped
	pool.map(process_decoys_file, zipped)
	pool.close()
	pool.join()
	print "done"
