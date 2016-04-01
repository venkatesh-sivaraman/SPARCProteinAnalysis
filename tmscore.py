"""This module calculates the TM-score using TMscore.jar, which is located in the same directory as this file."""

import subprocess
import os
from os.path import join
import shutil
from polypeptide import *

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

def best_tm_score(input, original_structure, range=None, dists=None, sec_structs=None, separate_scores=False):
	'''If you pass dists, the SPARC score will be logged next to the TM-score.
		If you pass range, the section of the original structure marked by range will be compared.'''
	tmpdir = os.path.join(os.path.dirname(input), "sp_tmp")
	if os.path.exists(tmpdir): shutil.rmtree(tmpdir)
	os.mkdir(tmpdir)
	max_tm = 0.0
	max_model = 0
	tm_sparc = []
	min_sparc = 0.0
	total_tm = 0.0
	num_tm = 0
	
	peptide = Polypeptide()
	if range:
		#Write the section of the original structure that corresponds to range
		peptide.read(original_structure)
		original_structure = os.path.join(tmpdir, "orig.pdb")
		with open(original_structure, "w") as file:
			file.write(peptide.pdb(range=range))

	with open(input, "r") as file:
		next_pdb = ""
		modelno = 0
		for line in file:
			next_pdb += line
			if "MODEL" in line:
				modelno = int(line[5:].strip())
			elif "ENDMDL" in line:
				with open(os.path.join(tmpdir, "model.pdb"), "w") as file:
					file.write(next_pdb)
				tm = calculate_tm_score(original_structure, join(tmpdir, "model.pdb"))
				if tm > 10000: continue
				if dists:
					try:
						peptide.read(os.path.join(tmpdir, "model.pdb"))
					except:
						next_pdb = ""
						del peptide.aminoacids[:]
						peptide.hashtable.clear()
						continue
					if sec_structs:
						if ',' in sec_structs:
							peptide.add_secondary_structures(sec_structs, format='csv')
						else:
							peptide.add_secondary_structures(sec_structs, format='pdb')
					if separate_scores:
						sp_score = ""
						for d in dists:
							sp_score += str(d.score(peptide, peptide.aminoacids)) + "\t"
						sp_score = sp_score[:-1]
						print str(modelno) + "\t" + str(tm) + "\t" + sp_score
					else:
						sp_score = sum(d.score(peptide, peptide.aminoacids) for d in dists)
						print str(modelno) + "\t" + str(tm) + "\t" + str(sp_score)
					
					tm_sparc.append((sp_score, tm))
					if tm > max_tm:
						min_sparc = sp_score
					del peptide.aminoacids[:]
					peptide.hashtable.clear()
				next_pdb = ""
				if tm > max_tm:
					max_tm = tm
					max_model = modelno
				total_tm += tm
				num_tm += 1
	shutil.rmtree(tmpdir)
	if len(tm_sparc):
		rank = len([x for x in tm_sparc if x[0] <= min_sparc])
		best_items = sorted(tm_sparc, key=lambda x: x[0])[:25]
		avg = sum(x[1] for x in best_items) / float(len(best_items))
		maximum = max(x[1] for x in best_items)
		minimum = min(x[1] for x in best_items)
		print "Best:", max_model, max_tm, ". Rank among structures:", rank, minimum, maximum, avg
		return (max_tm, min_sparc, sum(x[0] for x in tm_sparc) / float(len(tm_sparc)), rank, minimum, maximum, avg)
	else:
		if num_tm == 0: return 0.0
		print "Best:", max_model, max_tm, "avg:", total_tm / float(num_tm)
		return max_tm
