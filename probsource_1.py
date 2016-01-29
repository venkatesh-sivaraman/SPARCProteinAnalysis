#Removed 2/6/15, indents not preserved

def iter_ensemble_alpha(self, aminoacids, proximity=1, step=1):
	"""This subclassed implementation adds optimization to the orientations of individual amino acids."""
		if len(aminoacids) == 1:
			for conformation in super(AAProbabilitySource, self).iter_ensemble_alpha(aminoacids, proximity, step):
				yield [self.optimize_orientation(aminoacids[0], conformation[0])]
	else:
		for conformation in super(AAProbabilitySource, self).iter_ensemble_alpha(aminoacids, proximity, step):
			yield [self.optimize_orientation(aminoacids[i], conformation[i]) for i in xrange(len(aminoacids))]

def optimize_orientation(self, aminoacid, conformation):
	"""Pass in an amino acid (position unchanged) and a conformation (position zone object). This function will determine the best orientation for the specified alpha zone, relative to the two consecutive amino acids, and return it in a new PositionZone object.
		This function is NOT used anymore. It is not effective because the correct alpha carbon location is not found often enough."""
			#See Notability - Science Fair Calculations.
			if aminoacid.tag > 0 and aminoacid.tag < len(self.protein.aminoacids) - 1:
			#The alpha zone of aminoacid, the C of the previous aa, and the N of the next aa form a triangle - therefore a plane. We need to find the bisector of the angle NAC.
			cprev = self.protein.aminoacids[aminoacid.tag - 1].toglobal(Point3D(1.54, math.acos(-1.0 / 3.0), math.acos(-1.0 / 3.0)).tocartesian()).subtract(aminoacid.acarbon)
			nnext = self.protein.aminoacids[aminoacid.tag + 1].toglobal(Point3D(1.5, 2 * math.pi - math.acos(-1.0 / 3.0), math.acos(-1.0 / 3.0)).tocartesian()).subtract(aminoacid.acarbon)
			if dotproduct(cprev, nnext) >= 0.99:
				#The vectors are coincident
				bisector = cprev
		elif dotproduct(cprev, nnext) <= -0.99:
			#The vectors point in opposite directions, so there are infinite bisectors - choose the bisector that lies in the plane of the existing C and N atoms
			c = aminoacid.carbon.subtract(aminoacid.acarbon)
				b1 = crossproduct(cprev, crossproduct(c, aminoacid.nitrogen.subtract(aminoacid.acarbon)))
				b2 = b1.multiply(-1.0)
				if b1.distanceto(c) < b2.distanceto(c):
					bisector = b1
			else:
				bisector = b2
			else:
				bisector = cprev.normalize().add(nnext.normalize())
#Now rotate the bisector around the plane defined by acarbon, cprev, and nnext in either direction.
#The axis of rotation is the normal vector to the plane.
axis = crossproduct(cprev, nnext).tolist()
	v1 = Point3D.list(np.dot(rotation_matrix(axis, math.acos(-1.0 / 3.0) / 2.0), bisector.tolist()).tolist())
		v2 = Point3D.list(np.dot(rotation_matrix(axis, -math.acos(-1.0 / 3.0) / 2.0), bisector.tolist()).tolist())
			if v1.distanceto(cprev) < v1.distanceto(nnext):
				carbon = v1
				nitrogen = v2
		else:
			carbon = v2
				nitrogen = v1
		elif aminoacid.tag > 0:
			#Position the nitrogen so that the alpha carbon, the carbon, and the right nitrogen are collinear.
			cprev = self.protein.aminoacids[aminoacid.tag - 1].toglobal(Point3D(1.54, math.acos(-1.0 / 3.0), math.acos(-1.0 / 3.0)).tocartesian()).subtract(aminoacid.acarbon)
			nitrogen = cprev.multiply(1.5 / cprev.magnitude())
			#There are infinite possible planes for the atoms to be located in if the points are collinear, so choose the carbon that's closest to the original one. To do that, rotate nitrogen by 109.5 degrees in the plane formed between nitrogen, alpha carbon, and the existing carbon.
			axis = crossproduct(nitrogen, aminoacid.carbon.subtract(aminoacid.acarbon)).tolist()
			v1 = Point3D.list(np.dot(rotation_matrix(axis, math.acos(-1.0 / 3.0)), nitrogen.tolist()).tolist())
			v2 = Point3D.list(np.dot(rotation_matrix(axis, -math.acos(-1.0 / 3.0)), nitrogen.tolist()).tolist())
			if v1.distanceto(aminoacid.carbon) < v2.distanceto(aminoacid.carbon):
				carbon = v1
			else:
				carbon = v2
elif aminoacid.tag < len(self.protein.aminoacids) - 1:
	#Position the carbon so that the alpha carbon, the carbon, and the right nitrogen are collinear.
	nnext = self.protein.aminoacids[aminoacid.tag + 1].toglobal(Point3D(1.5, 2 * math.pi - math.acos(-1.0 / 3.0), math.acos(-1.0 / 3.0)).tocartesian()).subtract(aminoacid.acarbon)
		carbon = nnext.multiply(1.54 / nnext.magnitude())
			#There are infinite possible planes for the atoms to be located in if the points are collinear, so choose the nitrogen that's closest to the original one. To do that, rotate nitrogen by 109.5 degrees in the plane formed between nitrogen, alpha carbon, and the existing carbon.
			axis = crossproduct(carbon, aminoacid.nitrogen.subtract(aminoacid.acarbon)).tolist()
			v1 = Point3D.list(np.dot(rotation_matrix(axis, math.acos(-1.0 / 3.0)), carbon.tolist()).tolist())
			v2 = Point3D.list(np.dot(rotation_matrix(axis, -math.acos(-1.0 / 3.0)), carbon.tolist()).tolist())
			if v1.distanceto(aminoacid.nitrogen) < v2.distanceto(aminoacid.nitrogen):
				nitrogen = v1
		else:
			nitrogen = v2
		
		axes = aminoacid.coordinate_axes(nitrogen, carbon)
		conformation.x_axis = axes[0]
		conformation.y_axis = axes[1]
		conformation.z_axis = axes[2]
					return conformation

