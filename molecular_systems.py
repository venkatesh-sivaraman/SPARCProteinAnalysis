"""This module provides the MolecularSystem class, which provides umbrella objects for systems of proteins (and later, possibly nucleic acids and other molecules)."""
from polypeptide import *

class MolecularSystem(object):
	"""Instances of the MolecularSystem class represent systems of one or more biological molecules. Currently, only proteins (of the Polypeptide class) are supported. This object currently contains simply a list of molecules."""
	def __init__(self, mols):
		self.molecules = mols

	def check_steric_clash(self, aa, mol, steric_cutoff=4.0, consec=2, mindiff=1, excluded=[]):
		"""Pass in an amino acid along with its parent molecule to determine if there are any steric clashes in the system. If you pass a value for mindiff, any steric clashes involving the same biomolecule must have at least mindiff separation to be counted. If you pass an array of aa tags for excluded, it's assumed you mean tags in the protein given at mol."""
		for molecule in self.molecules:
			if molecule == mol:
				if len(molecule.nearby_aa(aa, steric_cutoff, consec=consec, mindiff=mindiff, excluded=excluded)) > 0:
					return True
			else:
				if len(molecule.nearby_aa(aa.acarbon, steric_cutoff)) > 0:
					return True
		return False

	def nearby_aa(self, aa, mol, distance, consec=2, mindiff=1, excluded=[]):
		nearby = []
		for molecule in self.molecules:
			if molecule == mol:
				nearby += molecule.nearby_aa(aa, distance, consec=consec, mindiff=mindiff, excluded=excluded)
			else:
				nearby += molecule.nearby_aa(aa, distance)
		return nearby

	def center(self):
		"""Center the molecules at the origin."""
		center = Point3D.zero()
		ct = 0
		for mol in self.molecules:
			if not isinstance(mol, Polypeptide): continue #Need to change this when some other biomolecule comes up
			for aa in mol.aminoacids:
				center = center.add(aa.acarbon)
				ct += 1
		center = center.multiply(1.0 / ct)
		for mol in self.molecules:
			if not isinstance(mol, Polypeptide): continue
			for aa in mol.aminoacids:
				aa.acarbon = aa.acarbon.subtract(center)
				aa.nitrogen = aa.nitrogen.subtract(center)
				aa.carbon = aa.carbon.subtract(center)

	def pdb(self, modelno=1):
		"""Write out the PDB file with each polypeptide represented as a different chain."""
		total_pdb = "MODEL        %d\n" % modelno
		chain = 0
		for mol in self.molecules:
			if not isinstance(mol, Polypeptide): continue
			pdb = mol.pdb(modelno=modelno, chain=chain)
			total_pdb += pdb[pdb.find("\n") + 1 : pdb.find("ENDMDL")]
			chain += 1
		total_pdb += "ENDMDL\n"
		return total_pdb