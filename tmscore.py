"""This module calculates the TM-score using TMscore.jar, which is located in the same directory as this file."""

import subprocess
import os
from os.path import join

def calculate_tm_score(native, candidate):
	string = "java -jar '/Users/venkatesh-sivaraman/Documents/Xcode Projects/PythonProteins/TMscore.jar' '{}' '{}'".format(candidate, native)
	output = subprocess.check_output(string, shell=True)
	lines = output.split("\n")
	scoreline = lines[14]
	num = float(scoreline.split("=")[1].split("(")[0])
	return num

def calculate_tm_scores(directory, output):
	"""This function calculates the TM-scores for every structure in the directory you provide. The native structure used is the one that contains neither 'em' nor 'closc.'"""
	files = os.listdir(directory)
	nativepath = None
	for path in files:
		if "em" not in path and "closc" not in path and "pdb" in path:
			nativepath = path
			break
	assert nativepath is not None, "No native path found in {}".format(directory)
	with open(output, "w") as file:
		for path in files:
			if "pdb" not in path: continue
			score = calculate_tm_score(join(directory, nativepath), join(directory, path))
			print path, score
			file.write("{} {}\n".format(path, score))