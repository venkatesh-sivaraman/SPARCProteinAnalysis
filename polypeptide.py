from proteinmath import *
from aminoacids import *
from randomcoil import *
from secondary_structure import *
import os.path
import string, resource
import gc

class Polypeptide(object):
	'Represents a series of amino acids and provides for input from a file and further manipulation.'

	def __init__(self, preaas=None):
		self.hashtable = AAHashTable()
		self.aminoacids = []
		self.secondary_structures = []
		self.mass = 0.0
		if preaas:
			for i, aa in enumerate(preaas):
				aa = aa.withtag(i)
				self.hashtable.add(aa)
				self.aminoacids.append(aa)
				if aa:
					self.mass += aa.mass
	
	def randomcoil(self, seq, permissions=None, struct_permissions=None):
		del self.hashtable
		self.mass = 0.0
		self.hashtable = AAHashTable()
		if len(self.secondary_structures) and struct_permissions:
			self.aminoacids = generate_randomcoil(seq, permissions=permissions, secondary_structures=self.secondary_structures, struct_permissions=struct_permissions)
		else:
			self.aminoacids = generate_randomcoil(seq, permissions)
		for aa in self.aminoacids:
			if aa: self.mass += aa.mass
			self.hashtable.add(aa)
	
	def add_aa(self, aa):
		self.aminoacids.append(aa)
		self.hashtable.add(aa)
		self.mass = 0.0
		for aa in self.aminoacids:
			if aa: self.mass += aa.mass

	def add_aas(self, aas):
		self.mass = 0.0
		if not self.hashtable: self.hashtable = AAHashTable()
		for aa in aas:
			self.aminoacids.append(aa)
			self.hashtable.add(aa)
			if aa: self.mass += aa.mass

	def __str__(self):
		descriptions = [self.aminoacids[x].__str__() for x in range(len(self.aminoacids))]
		return "Polypeptide, %i amino acids:\n\t" % len(self.aminoacids) + str(descriptions)
	
	def nearby_aa(self, aa, distance, index=-1, consec=2, excluded=[], mindiff=1):
		"""Pass in an array of integers for excluded to specify any tags that should not be counted in the search. If you specify a mindiff and consec is False, this method will only return amino acids whose tag difference is greater than mindiff."""
		assert aa != None, "Must provide a valid amino acid object"
		#if index == -1:
		#	index = self.aminoacids.index(aa)
		ret = self.hashtable.nearby_aa(aa, distance)
		ret2 = []
		for ret_aa in ret:
			if ((consec == False and math.fabs(ret_aa.tag - aa.tag) > mindiff) or (consec == True and math.fabs(ret_aa.tag - aa.tag) == 1) or (consec == 2)) and ret_aa.tag not in excluded:
				if consec == True and aa.has_break and ret_aa.tag == aa.tag + 1: continue
				ret2.append(ret_aa)
		return ret2
	
	def xyz(self, highlight=None, escaped=True):
		if escaped is True:
			ret = "%d\\nNo comment\\n" % (len(self.aminoacids) * 4)
			for aa in self.aminoacids:
				name = "C"
				if highlight is not None and aa.tag in highlight: name = "Cl"
				ret += "%s\t%.4f\t%.4f\t%.4f\\n" % (name, aa.acarbon.x, aa.acarbon.y, aa.acarbon.z)
				h = aa.toglobal(Point3D(1.1, 0.0, math.acos(-1.0 / 3.0)).tocartesian())
				ret += "H\t%.4f\t%.4f\t%.4f\\n" % (h.x, h.y, h.z)
				ret += "C\t%.4f\t%.4f\t%.4f\\n" % (aa.carbon.x, aa.carbon.y, aa.carbon.z)
				ret += "N\t%.4f\t%.4f\t%.4f\\n" % (aa.nitrogen.x, aa.nitrogen.y, aa.nitrogen.z)
			return ret
		else:
			ret = "%d\nNo comment\n" % (len(self.aminoacids) * 4)
			for aa in self.aminoacids:
				name = "C"
				if highlight is not None and aa.tag in highlight: name = "Cl"
				ret += "%s\t%.4f\t%.4f\t%.4f\n" % (name, aa.acarbon.x, aa.acarbon.y, aa.acarbon.z)
				h = aa.toglobal(Point3D(1.1, 0.0, math.acos(-1.0 / 3.0)).tocartesian())
				ret += "H\t%.4f\t%.4f\t%.4f\n" % (h.x, h.y, h.z)
				ret += "C\t%.4f\t%.4f\t%.4f\n" % (aa.carbon.x, aa.carbon.y, aa.carbon.z)
				ret += "N\t%.4f\t%.4f\t%.4f\n" % (aa.nitrogen.x, aa.nitrogen.y, aa.nitrogen.z)
			return ret

	def pdb(self, modelno=1, chain=0):
		ret = "MODEL        %d\n" % modelno
		idx = 1
		chain_letters = string.ascii_uppercase
		for i, aa in enumerate(self.aminoacids):
			aa.nitrogen = aa.toglobal(Point3D(NITROGEN_BOND_LENGTH, math.pi + math.acos(-1.0 / 2.0) / 2.0, math.acos(-1.0 / 3.0)).tocartesian())
			aa.carbon = aa.toglobal(Point3D(CARBON_BOND_LENGTH, math.pi - math.acos(-1.0 / 2.0) / 2.0, math.acos(-1.0 / 3.0)).tocartesian())
			ret += "ATOM  {0:>5}  N   {1} {2}{3:>4}    {4:>8}{5:>8}{6:>8}\n".format(idx, aa.type, chain_letters[chain], i + 1, "%.3f" % aa.nitrogen.x, "%.3f" % aa.nitrogen.y, "%.3f" % aa.nitrogen.z)
			idx += 1
			ret += "ATOM  {0:>5}  CA  {1} {2}{3:>4}    {4:>8}{5:>8}{6:>8}\n".format(idx, aa.type, chain_letters[chain], i + 1, "%.3f" % aa.acarbon.x, "%.3f" % aa.acarbon.y, "%.3f" % aa.acarbon.z)
			idx += 1
			h = aa.toglobal(Point3D(1.1, 0.0, math.acos(-1.0 / 3.0)).tocartesian())
			ret += "ATOM  {0:>5}  H   {1} {2}{3:>4}    {4:>8}{5:>8}{6:>8}\n".format(idx, aa.type, chain_letters[chain], i + 1, "%.3f" % h.x, "%.3f" % h.y, "%.3f" % h.z)
			idx += 1
			ret += "ATOM  {0:>5}  C   {1} {2}{3:>4}    {4:>8}{5:>8}{6:>8}\n".format(idx, aa.type, chain_letters[chain], i + 1, "%.3f" % aa.carbon.x, "%.3f" % aa.carbon.y, "%.3f" % aa.carbon.z)
			idx += 1
			for atom, location in aa.otheratoms.iteritems():
				location = aa.toglobal(location)
				ret += "ATOM  {0:>5}  {1:<4}{2} {3}{4:>4}    {5:>8}{6:>8}{7:>8}\n".format(idx, atom.strip(), aa.type, chain_letters[chain], i + 1, "%.3f" % location.x, "%.3f" % location.y, "%.3f" % location.z)
				idx += 1
			if aa.has_break: chain += 1
		ret += "ENDMDL\n"
		return ret

	def rg(self):
		"""Calculates the radius of gyration of the polypeptide."""
		#First compute the center of mass
		cm = Point3D.zero()
		total_mass = 0.0
		for aa in self.aminoacids:
			if not aa: continue
			cm = cm.add(aa.acarbon.multiply(aa.mass))
			total_mass += aa.mass
		cm = cm.multiply(1.0 / total_mass)
		
		#Now compute Rg
		rg = 0.0
		for aa in self.aminoacids:
			if not aa: continue
			rg += aa.mass * (aa.acarbon.distanceto(cm) ** 2)
		rg = math.sqrt(rg / total_mass)
		print "Radius of gyration:", rg
		return rg

	def center(self):
		"""Center the polypeptide at the origin of the global coordinate system."""
		center = Point3D.zero()
		for aa in self.aminoacids:
			center = center.add(aa.acarbon)
		center = center.multiply(1.0 / len(self.aminoacids))
		for aa in self.aminoacids:
			aa.acarbon = aa.acarbon.subtract(center)
			aa.nitrogen = aa.nitrogen.subtract(center)
			aa.carbon = aa.carbon.subtract(center)

	def read_file(self, f, checkgaps=False, otheratoms=False, secondary_structure=False, fillgaps=False, fillends=False, cache_aas=None):
		"""Pass in a file or file-like object to read. Set checkgaps to True to return a list of missing amino acid indices (first in the tuple if necessary) if there is one or more gaps in the chain. Set otheratoms to True to add all atoms found in the PDB file to the amino acids (under the otheratoms array property of the amino acids)."""
		gaps = []
		current_aa = None
		current_chain = None
		current_seq = 0
		current_tag = -1
		minseq = 100000
		maxseq = -10000
		if secondary_structure:
			self.secondary_structures = []
			current_sheet = None
		
		for line in f:
			if line.find('ENDMDL') != -1:	#Only keep one model for each polypeptide
				if current_aa is not None and len(self.aminoacids) > 0 and current_aa.acarbon == Point3D.zero():
					self.aminoacids.pop()
				assert current_aa is None or current_aa.carbon != Point3D.zero() or current_aa.nitrogen != Point3D.zero(), "Cannot have an amino acid with no carbon or nitrogen location"
				if current_aa:
					current_aa.compute_coordinate_system_vectors()
					if current_aa.i == Point3D.zero() and len(self.aminoacids) > 0:
						self.aminoacids.pop()
				break
		
			if secondary_structure:
				if line.find('HELIX') == 0:
					start_seq = int(line[22:25])
					end_seq = int(line[34:37])
					type = line[39:41]
					if type[1] in "1234567890":
						type = int(type)
					else:
						type = int(type[0])
					chain_start = line[19]
					chain_end = line[31]
					if chain_start.upper() != "A" or chain_end.upper() != "A": continue
					self.secondary_structures.append(Helix(start_seq - 1, end_seq - 1, type))
				elif line.find('SHEET') == 0:
					sense = int(line[38:40])
					start_seq = int(line[23:26])
					end_seq = int(line[34:37])
					chain_start = line[21]
					chain_end = line[32]
					if chain_start.upper() != "A" or chain_end.upper() != "A": continue
					if not current_sheet or sense == 0:
						if current_sheet:
							self.secondary_structures.append(current_sheet)
						current_sheet = Sheet(start_seq - 1, end_seq - 1, sense)
					else:
						current_sheet.add_strand(start_seq - 1, end_seq - 1, sense)

			if line.find('ATOM') != 0: continue
			if line[17:20] == "SOL": break		#There are no more amino acids after the solvent molecules start

			seq = int(line[23:26])
			chain = line[21:22]
			if chain is not current_chain and current_aa:
				current_aa.has_break = True
			if chain is not current_chain or current_seq != seq:
				current_chain = chain
				if current_aa is not None and current_aa.acarbon == Point3D.zero() and len(self.aminoacids) > 0:
					self.aminoacids.pop()
				assert current_aa is None or current_aa.carbon != Point3D.zero() or current_aa.nitrogen != Point3D.zero(), "Cannot have an amino acid with no carbon or nitrogen location: %r" % current_aa
				if current_aa:
					#mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
					current_aa.compute_coordinate_system_vectors()
					#print "Mem check 2:", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss - mem
					#mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
					if current_aa.i == Point3D.zero() and len(self.aminoacids) > 0:
						self.aminoacids.pop()
						if fillgaps:
							self.aminoacids.append(None)
				if seq - current_seq > 1 and (current_seq != 0 or fillends):
					current_seq += 1
					if checkgaps:
						print "Gap between", current_seq, "and", seq
						for n in xrange(current_seq, seq): gaps.append(n)
					if fillgaps:
						while seq != current_seq:
							self.aminoacids.append(None)
							current_seq += 1
							current_tag += 1
				current_seq = seq
				if current_seq < minseq: minseq = current_seq
				if current_seq > maxseq: maxseq = current_seq
				current_tag += 1

				if cache_aas:
					current_aa = cache_aas[0]
					del cache_aas[0]
					current_aa.type = line[17:20]
					current_aa.tag = current_tag
				else:
					current_aa = AminoAcid(line[17:20], current_tag)
				self.aminoacids.append(current_aa)
	
			atom_code = line[13:17]
			replaced = atom_code.replace(" ", "")
			location = Point3D(float(line[30:38]), float(line[38:46]), float(line[46:54]))
			if replaced == "CA" or atom_code == "CA A":
				current_aa._acarbon = location
			elif replaced == "CB" or atom_code == "CB A":
				current_aa.sidechain = location
				if otheratoms:
					current_aa.add_other_atom(atom_code.split()[0], location)
			elif replaced == "C" or atom_code == "C  A":
				current_aa.carbon = location
			elif replaced == "N" or atom_code == "N  A":
				current_aa.nitrogen = location
			elif otheratoms and atom_code[-1] != "B" and atom_code[0] != "H" and line[12] != "H" and replaced != "OXT":
				current_aa.add_other_atom(atom_code.split()[0], location)
		
		if current_aa is not None and current_aa.acarbon == Point3D.zero() and len(self.aminoacids) > 0:
			self.aminoacids.pop()

		if self.hashtable:
			self.hashtable.clear()
		else:
			self.hashtable = AAHashTable()
		self.mass = 0.0
		for aa in self.aminoacids:
			if aa:
				self.mass += aa.mass
				self.hashtable.add(aa)

		if secondary_structure:
			if current_sheet:
				self.secondary_structures.append(current_sheet)

		if checkgaps is True:
			return (gaps, (minseq, maxseq))
		else:
			return (minseq, maxseq)

	def read(self, filepath, checkgaps=False, otheratoms=False, secondary_structure=False, fillgaps=False, fillends=False):
		"""Set checkgaps to True to return a True value (first in the tuple if necessary) if there is one or more gaps in the chain. Set otheratoms to True to add all atoms found in the PDB file to the amino acids (under the otheratoms array property of the amino acids). Set secondary_structure to True to populate the secondary_structures array of the protein with Helix and Sheet objects. Set fillgaps to True to fill the peptide amino acids array with None objects wherever there is a gap."""
		f = open(filepath, 'r')
		ret = self.read_file(f, checkgaps, otheratoms, secondary_structure, fillgaps)
		f.close()
		return ret
	
	def secondary_structure_aa(self, aa_idx):
		"""Returns the (structure, strand) of the protein that contains amino acid at the index provided, and None if there is no secondary structure defined there."""
		if isinstance(aa_idx, AminoAcid):
			return find_secondary_structure(self.secondary_structures, aa_idx.tag)
		return find_secondary_structure(self.secondary_structures, aa_idx)

	def add_secondary_structures(self, sec_struct_string, format='csv', range=None):
		"""Takes a string of secondary structures in PDB or CSV format, depending on the format option, and adds them to the secondary structure array. CSV format should be 'type', 'subtype', 'start', 'end'. For instance, helix,1,3,6 covers amino acids at indices 2-5 (values are 1-indexed). If you pass a tuple (start, end) for range, then the secondary structures will be shifted to only cover the amino acids within that range."""
		current_sheet = None
		self.secondary_structures = []
		if format == 'csv':
			for line in sec_struct_string.split('\n'):
				comps = line.split(',')
				start = int(comps[2]) - 1
				end = int(comps[3]) - 1
				if range:
					if end < range[0] or start > range[1]: continue
					start = max(0, start - range[0] + 1)
					end = min(end - range[0] + 1, range[1] - range[0])
				if comps[0].lower() == 'helix':
					self.secondary_structures.append(Helix(start, end, int(comps[1])))
				elif comps[0].lower() == 'sheet':
					if int(comps[1]) == 0:
						if current_sheet: self.secondary_structures.append(current_sheet)
						current_sheet = Sheet(start, end, int(comps[1]))
					else:
						current_sheet.add_strand(start, end, int(comps[1]))
				else: print "Unknown secondary structure"
		else:
			for line in sec_struct_string.split('\n'):
				if line.find('HELIX') == 0:
					start_seq = int(line[22:26]) - 1
					end_seq = int(line[34:38]) - 1
					if range:
						if end_seq < range[0] or start_seq > range[1]: continue
						start_seq = max(0, start_seq - range[0] + 1)
						end_seq = min(end_seq - range[0] + 1, range[1] - range[0])
					type = int(line[39:41])
					chain_start = line[19]
					chain_end = line[31]
					if chain_start.upper() != "A" or chain_end.upper() != "A": continue
					self.secondary_structures.append(Helix(start_seq, end_seq, type))
				elif line.find('SHEET') == 0:
					sense = int(line[38:40])
					start_seq = int(line[23:26]) - 1
					end_seq = int(line[34:37]) - 1
					if range:
						if end_seq < range[0] or start_seq > range[1]: continue
						start_seq = max(0, start_seq - range[0] + 1)
						end_seq = min(end_seq - range[0] + 1, range[1] - range[0])
					chain_start = line[21]
					chain_end = line[32]
					if chain_start.upper() != "A" or chain_end.upper() != "A": continue
					if not current_sheet or sense == 0:
						if current_sheet:
							self.secondary_structures.append(current_sheet)
						current_sheet = Sheet(start_seq, end_seq, sense)
					else:
						current_sheet.add_strand(start_seq, end_seq, sense)
		if current_sheet:
			self.secondary_structures.append(current_sheet)


