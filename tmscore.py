"""This module calculates the TM-score using TMscore.jar, which is located in the same directory as this file."""

import subprocess
import os
from os.path import join
import shutil

def calculate_tm_score(native, candidate):
	string = "java -jar '/Users/venkatesh-sivaraman/Documents/Xcode Projects/PythonProteins/TMscore.jar' '{}' '{}'".format(candidate, native)
	output = subprocess.check_output(string, shell=True)
	lines = output.split("\n")
	scoreline = lines[14]
	num = float(scoreline.split("=")[1].split("(")[0])
	return num

def calculate_tm_scores(directory, output, natives=None):
	"""This function calculates the TM-scores for every structure in the directory you provide. The native structure used is the one that contains neither 'em' nor 'closc.'"""
	if os.path.isdir(directory):
		files = os.listdir(directory)
		nativepath = None
		if natives == None:
			for path in files:
				if ("em" not in path and "closc" not in path and "pdb" in path) or path == os.path.basename(directory) + ".pdb":
					nativepath = path
					break
		elif os.path.isdir(nativepath):
			nativepath = os.path.join(natives, os.path.basename(directory) + ".pdb")
		else:
			nativepath = natives
		assert nativepath is not None, "No native path found in {}".format(directory)
		print os.path.basename(directory), nativepath
		with open(output, "w") as file:
			for path in files:
				if ("T" not in path and "pdb" not in path): continue #or path.count("_") > 1
				score = calculate_tm_score(nativepath, join(directory, path))
				print path, score
				file.write("{},{}\n".format(path, score))
	else:
		nativepath = None
		if natives:
			if os.path.isdir(natives):
				nativepath = os.path.join(natives, os.path.basename(directory) + ".pdb")
			else:
				nativepath = natives
		assert nativepath is not None, "No native path found in {}".format(directory)
		tmpdir = join(os.path.dirname(directory), "sp_tmp")
		if not os.path.exists(tmpdir):
			os.mkdir(tmpdir)
		with open(output, "w") as outfile:
			tmpfile = open(os.path.join(tmpdir, "fun.pdb"), "w")
			with open(directory, "r") as file:
				modelno = 0
				for line in file:
					tmpfile.write(line)
					if "ENDMDL" in line:
						tmpfile.close()
						print "Writing model no", modelno
						score = calculate_tm_score(nativepath, join(tmpdir, "fun.pdb"))
						print modelno, score
						outfile.write("{},{}\n".format(modelno, score))
						tmpfile = open(os.path.join(tmpdir, "fun.pdb"), "w")
					elif "MODEL" in line:
						modelno = int(line.strip().split()[-1])
		shutil.rmtree(tmpdir)