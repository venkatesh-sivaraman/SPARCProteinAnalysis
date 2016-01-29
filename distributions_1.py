def anglescore(self, protein, data):
	"""For frequency distributions, pass in one hypothetical amino acid for data."""
	return self.score(protein, [data])
	aa = data
	consec = 2
	if self.type == frequency_consec_disttype: consec = 1
	elif self.type == frequency_nonconsec_disttype: consec = 0
	nearby = protein.nearby_aa(aa, 10.0, consec=consec)
	score = 0.0
	print "Here?!"
	for aa2 in nearby:
		#Compute the distance from aa's x-axis and y-axis locations to aa2's alpha carbon.
		'''zone = aa2.tolocal(aa.acarbon).floor()
			xradius = zone.add(aa.i).distanceto(Point3D.zero())
			yradius = zone.add(aa.j).distanceto(Point3D.zero())
			if zone not in self.axis_regressions: continue
			coeffs = self.axis_regressions[zone]
			score -= math.log(max(coeffs[2] * (xradius - coeffs[0]) ** 2 + coeffs[1], 0.00001)) - math.log(0.0001)
			score -= math.log(max(coeffs[5] * (xradius - coeffs[3]) ** 2 + coeffs[4], 0.00001)) - math.log(0.0001)'''
		tag1 = aacode(aa.type)
		tag2 = aacode(aa2.type)
		zone = aa.tolocal(aa2.acarbon).floor()
		#try:
		#S = -ln(F/F0) if F > 0, +5 otherwise
		freq = self.alpha_frequency(tag1, tag2, zone)
		if freq > 0:
			subscore = -math.log(freq / self.median_frequencies[tag1][tag2])
		else:
			subscore = 5
		score += subscore
	
	return score
