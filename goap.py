import os, subprocess
from os.path import join
import multiprocessing

def write_inp_file(inppath, path):
	with open(inppath, "w") as file:
		file.write("/home/ubuntu/Documents/goap-alone\n")
		file.write(path)

def process_decoys_file((input, output, nativepath)):
	if not os.path.isdir(input): return
	if os.path.exists(output): return
	if os.path.basename(input) == "doc": return
	protein_name = os.path.basename(input)
	print protein_name
	paths = os.listdir(input)
	allpaths = [os.path.join(input, path) for path in paths]
	inppath = "../goap.inp"

	if "-" in protein_name:
		path = protein_name[:protein_name.find("-")] + ".pdb"
	else:
		path = protein_name + ".pdb"
	if nativepath: allpaths.append(os.path.join(nativepath, path))

	if nativepath:
		if os.path.exists(join(nativepath, path)):
			write_inp_file(inppath, join(nativepath, path))
			out = subprocess.check_output("./goap<\"" + inppath + "\"", shell=True)
			score = out.split("\t")[1]
			print "native scores:", score
			if output and score != 0.0:
				with open(output, "w") as file:
					file.write("{}; {}\n".format(protein_name + "_orig.pdb", score))
		else:
			print join(nativepath, path), "does not exist."
	for path in paths:
		if path == "list" or path == "rmsds": continue
		try:
			write_inp_file(inppath, join(input, path))
			out = subprocess.check_output("./goap<\"" + inppath + "\"", shell=True)
		except subprocess.CalledProcessError as e:
			out = e.output
		try:
			score = out.split("\t")[1]
			print score
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
	zipped = [(join(input, file), join(output, file + ".txt"), "/home/ubuntu/Documents/casp11.targets_unsplitted.release11242014") for file in reversed(files)]
	#print zipped
	pool.map(process_decoys_file, zipped)
	pool.close()
	pool.join()
	print "done"

if __name__ == '__main__':
	test_sparc("PATH TO DECOYS", "/home/ubuntu/Documents/casp-goap")