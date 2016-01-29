"""Comparison of CHARMM and SPARC"""

import subprocess
from subprocess import PIPE
from proteinmath import *
from polypeptide import *
import os
import gc
from os.path import join
from sparc_distribution import *
import urllib2
import time
import random
import shutil
import numpy

def run_charmm(pdbpath, mdpfile="../default.mdp", minimize=False, pairwise=False, bonded=True, solvate=False):
	"""Returns the CHARMM score of a protein located at pdbpath. The mdpfile is a parameter file needed for GROMACS to perform its preprocessing. The default is a file called default.mdp located in the directory outside the temp directory in which the pdb file is located.
		If minimize is True, a short MD simulation will be performed to minimize the energy. The resulting structure is deposited in the same directory as pdbpath, with filename final.pdb.
		Pairwise must be specified for the correct energy terms to be selected out of the gmx energy command.
		If bonded is True, the first item in the score tuple contains all bonded energy terms. If false, it contains all nonbonded energy terms."""
	origWD = os.getcwd() # remember our original working directory
	os.chdir(os.path.dirname(pdbpath))
	
	try:
		with open(os.devnull, "w") as fnull:
			#First, create a gmx file
			if solvate:
				water_model = "spc"
			else:
				water_model = "none"
			cmd1 = "rm \\#*\\#; gmx pdb2gmx -f " + os.path.basename(pdbpath) + " -ff charmm27 -water " + water_model + " -ignh"
			out = subprocess.check_output(cmd1, shell=True, stderr=fnull)
			
			#Then, edit the configuration to create a box
			cmd2 = "gmx editconf -f conf.gro -c -d 1.0"
			out = subprocess.check_output(cmd2, shell=True, stderr=fnull)
			
			mid_filenm = "out.gro"
			if solvate:
				#Solvate the protein
				cmd21 = "gmx solvate -cp out.gro -cs spc216.gro -o out_solv.gro -p topol.top"
				mid_filenm = "out_solv.gro"
				out = subprocess.check_output(cmd21, shell=True, stderr=fnull)

			#Now, run grompp to preprocess all the files
			cmd3 = "gmx grompp -c " + mid_filenm + " -p topol.top -f " + mdpfile
			out = subprocess.check_output(cmd3, shell=True, stderr=fnull)

			if minimize:
				print "Minimizing"
				#Run a short mdrun simulation to minimize the energy, then measure it
				cmd4 = "gmx mdrun -v -deffnm em -s topol.tpr -c final.gro"
				out = subprocess.check_output(cmd4, shell=True, stderr=fnull)
				#Convert the final structure back to PDB
				cmd41 = "gmx editconf -f final.gro -o final.pdb"
				out = subprocess.check_output(cmd41, shell=True, stderr=fnull)
				#subprocess.call("echo \"0\n0\n\" | gmx rms -s prot.pdb -f em.trr", shell=True, stderr=fnull)
			else:
				#Run an mdrun one-frame simulation to measure the energies
				cmd4 = "gmx mdrun -rerun out.gro -deffnm em -s topol.tpr"
				out = subprocess.check_output(cmd4, shell=True, stderr=fnull)

			#Finally, determine the energies using the energy function
			if bonded:
				elements = "1 2 3 4"
			else:
				elements = "6 7"
			if pairwise:
				cmd5 = "echo \"" + elements + " 9\n\n\" | gmx energy -f em.edr -s topol.tpr"
			else:
				cmd5 = "echo \"" + elements + " 10\n\n\" | gmx energy -f em.edr -s topol.tpr"
			out = subprocess.check_output(cmd5, shell=True, stderr=fnull)
	except:
		print "Exception"
		return None
		
	os.chdir(origWD) # get back to our original working directory
		
	gc.collect()
	if bonded and "Bond" in out:
		#Interpret the data and send it back
		lines = out[out.find("Bond"):].split("\n")
		comps = lines[0].split()
		bond_energy = float(comps[1])
		comps = lines[1].split()
		bond_energy += float(comps[1])
		comps = lines[2].split()
		bond_energy += float(comps[2])
		comps = lines[3].split()
		bond_energy += float(comps[2])
		comps = lines[4].split()
		nonbond_energy = float(comps[1])
		gc.collect()
		return (bond_energy, nonbond_energy)
	elif not bonded and "LJ-14" in out:
		lines = out[out.find("LJ-14"):].split("\n")
		comps = lines[0].split()
		bond_energy = float(comps[1])
		comps = lines[1].split()
		bond_energy += float(comps[1])
		comps = lines[2].split()
		nonbond_energy = float(comps[1])
		gc.collect()
		return (bond_energy, nonbond_energy)
	else:
		return None

def charmm_score(peptide, protpath, mdpfile="../default.mdp", aatags=[], minimize=False, pairwise=False, bonded=True, solvate=False):
	"""Obtains the bond and total CHARMM scores for the protein specified by peptide. Pass in an array of tag numbers to specify certain amino acids."""
	if len(aatags):
		newpeptide = Polypeptide(preaas=[peptide.aminoacids[i] for i in aatags])
		pdb = newpeptide.pdb()
		del newpeptide
	else:
		pdb = peptide.pdb()
	with open(protpath, "w") as file:
		file.write(pdb)
	cscores = run_charmm(protpath, mdpfile, minimize, pairwise, bonded, solvate)
	return cscores


def iter_residue_pair_pdbs(peptide, contactdist=5.0, bonded=False, ret_tags=False):
	"""Yields (idx, peptide) containing each pair of amino acids within contactdist in the peptide. The PDB string contains all the atoms if the peptide you specify was read using otheratoms."""
	idx = 0
	for i, aa1 in enumerate(peptide.aminoacids):
		for aa2 in peptide.nearby_aa(aa1, contactdist, consec=bonded):
			newpeptide = Polypeptide(preaas=[aa1, aa2])
			if ret_tags:
				print "TAGS:", aa1.tag, aa2.tag
				yield (idx, [aa1.tag, aa2.tag])
			else:
				yield (idx, newpeptide)
			idx += 1
			del newpeptide

def compare_charmm_sparc(pdbidfile, sparc_location, output, bonded=False, pairwise=True, idealized=True, minimize=False):
	"""Goes through each line in pdbidfile and evaluates CHARMM and SPARC on the structure. If idealized is True, then the program finds each combination of amino acid types and runs it through random permutations of their relative orientations. """
	
	#Load SPARC distribution
	sparc = SPARCBasicDistributionManager(sparc_location, bonded)
	
	tmpdir = os.path.join(os.path.dirname(output), "tmp_prot")
	if not os.path.exists(tmpdir):
		os.mkdir(tmpdir)
	
	with open(pdbidfile, 'r') as file:
		contents = file.readlines()
	peptide = Polypeptide()
	if idealized:
		done_pairs = []

	for pdbid in contents:
		pdbid = pdbid.strip()
		protpath = os.path.join(tmpdir, "prot.pdb")
		print "Processing " + pdbid + "..."
		try:
			response = urllib2.urlopen('http://www.rcsb.org/pdb/files/' + pdbid + '.pdb')
			peptide.read_file(response, otheratoms=True)
			response.close()
		except:
			print "==========================Omit %r" % pdbid
			del peptide.aminoacids[:]
			peptide.hashtable = None
			gc.collect()
			continue
		del response

		outputstring = ""
		if pairwise:
			if idealized and len(done_pairs) == AMINO_ACID_COUNT ** 2 / 2: break
			for i, pairpeptide in iter_residue_pair_pdbs(peptide, bonded=bonded):
				if idealized:
					cnt = 0
					tag1 = aacode(pairpeptide.aminoacids[0].type)
					tag2 = aacode(pairpeptide.aminoacids[1].type)
					if [tag1, tag2] in done_pairs or [tag2, tag1] in done_pairs: continue
					while cnt < 10:
						pairpeptide.aminoacids[0].acarbon = Point3D.zero()
						offset = Point3D(random.uniform(-3.0, 3.0), random.uniform(2.0, 6.0), random.uniform(-3.0, 3.0))
						pairpeptide.aminoacids[1].acarbon = pairpeptide.aminoacids[0].toglobal(offset)
						print tag1, tag2, offset

						#CHARMM
						score1 = charmm_score(pairpeptide, protpath, aatags=[0], minimize=minimize)
						score2 = charmm_score(pairpeptide, protpath, aatags=[1], minimize=minimize)
						score3 = charmm_score(pairpeptide, protpath, aatags=[], minimize=minimize)
						if not score1 or not score2 or not score3:
							outputstring += "ERROR" + pdbid + str(pairpeptide) + "\n"
							continue
						#SPARC
						spscore = sparc.score(pairpeptide, pairpeptide.aminoacids, onlyone=True)
						print spscore
						if spscore[0] != 0.0 or spscore[1] != 0.0:
							rel_loc = pairpeptide.aminoacids[0].tolocal(pairpeptide.aminoacids[1].acarbon)
							outputstring += "{},{},{},{},{},{},{},{},{},{},{}\n".format(pairpeptide.aminoacids[0].acarbon.distanceto(pairpeptide.aminoacids[1].acarbon), rel_loc.x, rel_loc.y, rel_loc.z, score1[1], score2[1], score3[1], spscore[0], spscore[1], spscore[2], spscore[3])
						cnt += 1
					time.sleep(15)
					done_pairs.append([tag1, tag2])
				else:
					#CHARMM
					score1 = charmm_score(pairpeptide, protpath, aatags=[0], minimize=minimize)
					score2 = charmm_score(pairpeptide, protpath, aatags=[1], minimize=minimize)
					score3 = charmm_score(pairpeptide, protpath, aatags=[], minimize=minimize)
					if not score1 or not score2 or not score3:
						outputstring += "ERROR" + pdbid + str(pairpeptide) + "\n"
						continue
					if minimize:
						newpeptide = Polypeptide()
						newpeptide.read(os.path.join(tmpdir, "final.pdb"))
					else:
						newpeptide = pairpeptide
					#SPARC
					spscore = sparc.score(newpeptide, newpeptide.aminoacids, onlyone=True)
					print spscore
					rel_loc = newpeptide.aminoacids[0].tolocal(newpeptide.aminoacids[1].acarbon)
					outputstring += "{},{},{},{},{},{},{},{},{},{},{}\n".format(newpeptide.aminoacids[0].acarbon.distanceto(newpeptide.aminoacids[1].acarbon), rel_loc.x, rel_loc.y, rel_loc.z, score1[1], score2[1], score3[1], spscore[0], spscore[1], spscore[2], spscore[3])
					if minimize:
						del newpeptide
				del pairpeptide
				gc.collect()
		else:
			#CHARMM
			score = charmm_score(peptide, protpath, minimize=minimize, pairwise=False)
			if not score:
				outputstring += "ERROR" + pdbid + "\n"
			else:
				#SPARC
				if minimize:
					newpeptide = Polypeptide()
					newpeptide.read(os.path.join(tmpdir, "final.pdb"))
				else:
					newpeptide = pairpeptide
				#SPARC
				spscore = sparc.score(newpeptide, newpeptide.aminoacids)
				print score[0], score[1], spscore
				rel_loc = newpeptide.aminoacids[0].tolocal(newpeptide.aminoacids[1].acarbon)
				outputstring += "{},{},{},{},{},{},{}\n".format(newpeptide.aminoacids[0].acarbon.distanceto(newpeptide.aminoacids[1].acarbon), rel_loc.x, rel_loc.y, rel_loc.z, score[0], score[1], spscore)
				if minimize:
					del newpeptide
			#time.sleep(30)

		with open(output, "a") as ofile:
			ofile.write(outputstring)
			outputstring = ""
		del peptide.aminoacids[:]
		peptide.hashtable = None
		gc.collect()
	del peptide

def iter_pz_variants(count, pz_alpha, *args):
	"""Yields `count` different points which fall in the position zone defined by pz_alpha (a point)."""
	for i in xrange(count):
		if len(args):
			pts = [Point3D(random.uniform(pz_alpha.x - 0.25, pz_alpha.x + 0.25), random.uniform(pz_alpha.y - 0.25, pz_alpha.y + 0.25), random.uniform(pz_alpha.z - 0.25, pz_alpha.z + 0.25))]
			for pz in args:
				pts.append(Point3D(random.uniform(pz.x - 0.25, pz.x + 0.25), random.uniform(pz.y - 0.25, pz.y + 0.25), random.uniform(pz.z - 0.25, pz.z + 0.25)))
			yield pts
		else:
			yield Point3D(random.uniform(pz_alpha.x - 0.25, pz_alpha.x + 0.25), random.uniform(pz_alpha.y - 0.25, pz_alpha.y + 0.25), random.uniform(pz_alpha.z - 0.25, pz_alpha.z + 0.25))

def compare_charmm_sparc_sp(structurepairs, sparc_location, output):
	"""Compares the CHARMM and SPARC scores based on the local structure pair files you provide. These are PDB models of two amino acids."""
	#Load SPARC distribution
	sparc = SPARCBasicDistributionManager(os.path.join(sparc_location, "consec"), True)
	sparc_nc = SPARCBasicDistributionManager(os.path.join(sparc_location, "nonconsec"), False)
	
	peptide = Polypeptide()
	
	tmpdir = os.path.join(os.path.dirname(output), "sp_tmp")
	if not os.path.exists(tmpdir): os.mkdir(tmpdir)
	protpath = os.path.join(tmpdir, "prot.pdb")
	for sp in os.listdir(structurepairs)[:2]:
		if ".pdb" not in sp: continue
		print "Processing " + sp + "..."
		try:
			peptide.read(os.path.join(structurepairs, sp), otheratoms=True)
		except:
			print "==========================Omit %r" % sp
			continue

		outputstring = ""
		cnt = 0
		fncomps = sp[:sp.find(".")].split("-")
		type1 = int(fncomps[0])
		type2 = int(fncomps[1])
		while cnt < 200:
			# zone1 expresses aa2's location in aa1's LCS, while zone2 expresses aa1's location in aa2's LCS.
			zone1 = Point3D.zero()
			keys = sparc.alpha_frequencies.keys()
			while sparc.alpha_frequencies[zone1][type1][type2] < 10 or zone1.y <= 0.0:
				zone1 = random.choice(keys)
			zone1 = Point3D(zone1.x + random.uniform(0.0, 0.9), zone1.y + random.uniform(0.0, 0.9), zone1.z + random.uniform(0.0, 0.9))
			zone2 = Point3D.zero()
			while sparc.alpha_frequencies[zone2][type2][type1] < 10 or zone2.y >= 0.0:
				zone2 = random.choice(keys)
			#Manipulate the locations of the two amino acids so that they respect zone1 and zone2.
			peptide.aminoacids[0].acarbon = Point3D.zero()
			peptide.aminoacids[1].acarbon = peptide.aminoacids[0].toglobal(zone1)
			axes = peptide.aminoacids[1].axes_for_zone(peptide.aminoacids[0].acarbon.subtract(peptide.aminoacids[1].acarbon), zone2, timeout=True)
			if not axes:
				print "Not valid orientation"
				continue
			print zone1, zone2, sparc.alpha_frequencies[zone1.floor()][type2][type1], sparc.alpha_frequencies[zone2][type1][type2]
			peptide.aminoacids[1].set_axes(*axes)
			print "Rel loc", peptide.aminoacids[1].tolocal(peptide.aminoacids[0].acarbon).floor()

			#CHARMM
			score1 = charmm_score(peptide, protpath, aatags=[0], mdpfile="../default_nonminimized.mdp", pairwise=True)
			score2 = charmm_score(peptide, protpath, aatags=[1], mdpfile="../default_nonminimized.mdp", pairwise=True)
			score3 = charmm_score(peptide, protpath, aatags=[], mdpfile="../default_nonminimized.mdp", pairwise=True)
			if not score1 or not score2 or not score3:
				outputstring += "ERROR" + sp + str(peptide) + "\n"
				continue
			'''newpeptide = Polypeptide()
			newpeptide.read(os.path.join(tmpdir, "final.pdb"))'''
			newpeptide = peptide
			#SPARC
			spscore1 = sparc.score(newpeptide, newpeptide.aminoacids, onlyone=True)
			spscore2 = sparc_nc.score(newpeptide, newpeptide.aminoacids)
			#rel_loc = newpeptide.aminoacids[0].tolocal(newpeptide.aminoacids[1].acarbon)
			distance = newpeptide.aminoacids[0].acarbon.distanceto(newpeptide.aminoacids[1].acarbon)
			outputstring += "{},{},{},{},{},{},{},{},{},{}\n".format(distance, score1[0], score1[1], score2[0], score2[1], score3[0], score3[1], spscore1[0], spscore1[1], spscore2) #'''(score3[1] + addlscores) / 26.0'''
			cnt += 1
			#del newpeptide
			gc.collect()

		with open(output, "a") as ofile:
			ofile.write(outputstring)
			outputstring = ""
		del peptide.aminoacids[:]
		peptide.hashtable = None
		gc.collect()
	del peptide
	shutil.rmtree(tmpdir)

#MARK: Magainin test

def protein_protein_energies(pdbfile, sparc_location, output):
	"""Takes the amino acids in the protein that you provide, duplicates them, and places the entire duplicated structure at various orientations with respect to the original structure, and evaluates both CHARMM and SPARC for each different orientation."""
	#Load SPARC distribution
	sparc = SPARCBasicDistributionManager(os.path.join(sparc_location, "consec"), True)
	sparc_nc = SPARCBasicDistributionManager(os.path.join(sparc_location, "nonconsec"), False)
	tmpdir = os.path.join(os.path.dirname(output), "sp_tmp")
	if not os.path.exists(tmpdir): os.mkdir(tmpdir)
	protpath = os.path.join(tmpdir, "prot.pdb")

	orig_dist = 0.0
	bchain_start = 0
	for dist in xrange(150):
		distx = 5.0 - dist / 5.0
		disty = -5.0 + dist / 5.0
		peptide = Polypeptide()
		peptide.read(pdbfile, otheratoms=True)
		inbchain = False
		for i in xrange(len(peptide.aminoacids)):
			if peptide.aminoacids[i].has_break:
				inbchain = True
				bchain_start = i + 1
				continue
			if inbchain:
				peptide.aminoacids[i].acarbon = peptide.aminoacids[i].acarbon.add(Point3D(0.0, distx, disty))
				if i == bchain_start:
					orig_dist = peptide.aminoacids[0].acarbon.distanceto(peptide.aminoacids[i].acarbon)
					print "Original distance:", orig_dist
		cscore = charmm_score(peptide, protpath, minimize=True, mdpfile="../default.mdp", bonded=False)
		print cscore
		if not cscore: continue
		newpeptide = Polypeptide()
		newpeptide.read(os.path.join(tmpdir, "final.pdb"))
		new_dist = newpeptide.aminoacids[0].acarbon.distanceto(newpeptide.aminoacids[bchain_start].acarbon)
		print "New:", new_dist
		spscore1 = sparc.score(newpeptide, newpeptide.aminoacids)
		spscore2 = sparc_nc.score(newpeptide, newpeptide.aminoacids)
		print spscore1, spscore2
		with open(output, "a") as ofile:
			ofile.write("{},{},{},{},{}\n".format(new_dist, cscore[0], cscore[1], spscore1, spscore2))
		del newpeptide
		del peptide
		gc.collect()

def batch_compare_charmm_sparc(input, sparc_location, output):
	"""This method goes through ALL structure files in input and calculates the CHARMM and SPARC scores for them, depositing the results in a text file located at output."""
	#Load SPARC distribution
	sparc = SPARCBasicDistributionManager(os.path.join(sparc_location, "consec"), True)
	sparc_nc = SPARCBasicDistributionManager(os.path.join(sparc_location, "nonconsec"), False)
	sparc_m = MediumDistributionManager(os.path.join(sparc_location, "medium"))
	tmpdir = os.path.join(os.path.dirname(output), "sp_tmp")
	if not os.path.exists(tmpdir): os.mkdir(tmpdir)
	protpath = os.path.join(tmpdir, "prot.pdb")
	
	for dir in os.listdir(input):
		if not os.path.isdir(os.path.join(input, dir)): continue
		paths = os.listdir(os.path.join(input, dir))
		for path in paths:
			if ".pdb" not in path: continue
			if random.uniform(0.0, 1.0) > 100.0 / len(paths): continue
			print path
			try:
				peptide = Polypeptide()
				peptide.read(os.path.join(input, dir, path), otheratoms=True)
				cscore = charmm_score(peptide, protpath, minimize=True, mdpfile="../default_solvated.mdp", bonded=False, solvate=True)
				print cscore
				if not cscore: continue
				newpeptide = Polypeptide()
				newpeptide.read(os.path.join(tmpdir, "final.pdb"))
				spscore1 = sparc.score(newpeptide, newpeptide.aminoacids)
				spscore2 = sparc_nc.score(newpeptide, newpeptide.aminoacids)
				spscore3 = sparc_m.score(newpeptide, newpeptide.aminoacids)
				with open(output, "a") as ofile:
					ofile.write("{},{},{},{},{},{}\n".format(path, cscore[0], cscore[1], spscore1, spscore2, spscore3))
				del newpeptide
				del peptide
				del spscore1, spscore2, spscore3, cscore
			except:
				print "Exception with the entire process for this path."
			gc.collect()
	print "done!"

#MARK: Testing individual energy terms

'''def bond_stretch(structurepairs, sparc_location, bonded=False, mindist=1.0, maxdist=9.5):
	"""Takes the structures found in structurepairs and shifts them to have various bond lengths between the carbon of the first amino acid and the nitrogen of the second one. The Ca, C, N, and Ca atoms will be constructed to be coplanar, with bond angles of 114 and 123 degrees, respectively."""
	#Load SPARC distribution
	sparc = SPARCBasicDistributionManager(sparc_location, bonded)
	
	peptide = Polypeptide()
	for sp in os.listdir(structurepairs):
		if ".pdb" not in sp: continue
		print "Processing " + sp + "..."
		try:
			peptide.read(os.path.join(structurepairs, sp), otheratoms=True)
		except:
			print "==========================Omit %r" % sp
			continue

		#First locate the two amino acids so that one is aligned with the GCS
		peptide.aminoacids[0].acarbon = Point3D.zero()
		peptide.aminoacids[0].set_axes(Point3D(1.0, 0.0, 0.0), Point3D(0.0, 1.0, 0.0), Point3D(0.0, 0.0, 1.0), normalized=True)
		#Now, set the second one so that its nitrogen atom is at 114 degrees from the carbon atom and 1.32 A away.

		del peptide.aminoacids[:]
		peptide.hashtable = None
		gc.collect()
	del peptide'''
