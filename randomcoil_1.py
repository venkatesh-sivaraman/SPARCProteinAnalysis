def rms_speed(temperature):
	"""See Science Fair 2015 (12/14/14) for how this works."""
	return math.sqrt(temperature)

def simulate_random(peptide, temperature=298, time=0.001):
	"""Provide temperature in K and time in ms"""
	last_aa = None
	distance = rms_speed(temperature) * (time * 0.001) * (10 ** 5)
	for i, aa in enumerate(peptide.aminoacids):
		def randomizer():
			return aa.acarbon.add(Point3D(random.uniform(-distance, distance), random.uniform(-distance, distance), random.uniform(-distance, distance)))
		new_point = randomizer()
		while True:
			new_point = randomizer()
			if (last_aa == None or (new_point.distanceto(last_aa.acarbon) > 1.0 and new_point.distanceto(last_aa.acarbon) < 4.0)) and (i >= len(peptide.aminoacids) - 1 or (new_point.distanceto(peptide.aminoacids[i + 1].acarbon) > 1.0 and new_point.distanceto(peptide.aminoacids[i + 1].acarbon) < 4.0)): break
		aa.acarbon = new_point
		last_aa = aa