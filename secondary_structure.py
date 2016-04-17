secondary_struct_helix = "helix"
secondary_struct_sheet = "sheet"

class Strand(object):
	"""Defines a strand of a protein along with its metadata."""
	def __init__(self, start, end, ids=[]):
		self.start = start
		self.end = end
		self.identifiers = ids

	def __repr__(self):
		return "Strand {}: {}-{}".format(self.identifiers, self.start, self.end)

class SecondaryStructure(object):
	"""Defines a type of secondary structure and its geometry. Secondary structures are defined by one or more strands, which are represented by Strand objects referencing the start and end indices of the amino acids in the strand."""
	def __init__(self, t, strands=[]):
		self.type = t
		self.strands = strands

	def add_strand(self, start, end, identifier=None):
		"""Adds a strand in the protein with a new set of amino acids - e.g., a new component of a beta sheet."""
		self.strands.append(Strand(start, end, [identifier]))

	def __repr__(self):
		desc = "{{Secondary Structure ({}); ".format(self.type)
		for strand in self.strands:
			desc += "{}, {} thru {}; ".format(strand.identifiers, strand.start, strand.end)
		return desc[:-2] + "}}"

class Helix(SecondaryStructure):
	"""A subclass of SecondaryStructure that predefines a helix. Nothing about this class is unique except its initializer. Helices by default have only one strand."""
	def __init__(self, start=0, end=0, helix_type=0):
		if end - start: strands = [Strand(start, end, [helix_type])]
		else: strands = []
		super(Helix, self).__init__(secondary_struct_helix, strands)
		self.helix_type = helix_type

class Sheet(SecondaryStructure):
	"""A subclass of SecondaryStructure that represents a set of beta sheets. Nothing about this class is unique except its initializer. You initialize a beta sheet with one strand, then add more strands using the add_strand method."""
	def __init__(self, start=0, end=0, initial_type=0):
		"""Pass in the initial sheet type if desired, but it is most likely 0 because later strands will be parallel (1) or anti-parallel (-1) to the first strand."""
		if end - start: strands = [Strand(start, end, [initial_type])]
		else: strands = []
		super(Sheet, self).__init__(secondary_struct_sheet, strands)

def find_secondary_structure(secondary_structures, aa_idx):
	"""Returns the (structure, strand) of the protein that contains amino acid at the index provided, and None if there is no secondary structure defined there."""
	current_struct = None
	current_strand = None
	if len(secondary_structures) > 0:
		current_struct = None
		for secondary_struct in secondary_structures:
			for strand in secondary_struct.strands:
				if strand.start <= aa_idx and strand.end >= aa_idx:
					current_struct = secondary_struct
					current_strand = strand
					break
	if current_strand:
		return (current_struct, current_strand)
	return None

#MARK: - Built-in Secondary Structures

sec_structs = [	secondary_struct_helix + "1", secondary_struct_helix + "2",
			   secondary_struct_helix + "3", secondary_struct_helix + "4",
			   secondary_struct_helix + "5", secondary_struct_helix + "6",
			   secondary_struct_helix + "7", secondary_struct_helix + "8",
			   secondary_struct_helix + "9", secondary_struct_helix + "10",
			   secondary_struct_sheet + "0", secondary_struct_sheet + "1",
			   secondary_struct_sheet + "-1"]

def secondary_structures_dict(inner_value={}):
	global sec_structs
	sec_dict = {}
	for ss in sec_structs:
		sec_dict[ss] = inner_value
	return sec_dict

def is_valid_secondary_structure(struct_type):
	global sec_structs
	if not struct_type:
		return False
	if isinstance(struct_type, basestring):
		return struct_type in sec_structs
	if len(struct_type) > 1:
		return (struct_type[0].type + str(struct_type[1].identifiers[0])) in sec_structs

