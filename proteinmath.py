import math
import numpy as np
from numpy import matlib
import random

def rotation_matrix(axis, theta):
	"""Return the rotation matrix associated with counterclockwise rotation about the given axis by theta radians."""
	axis = np.asarray(axis)
	theta = np.asarray(theta)
	axis = axis / math.sqrt(np.dot(axis, axis))
	a = math.cos(theta / 2)
	b, c, d = -axis * math.sin(theta / 2)
	aa, bb, cc, dd = a*a, b*b, c*c, d*d
	bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
	return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
					 [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
					 [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


class Point3D(object):
	'Represents a point or vector in 3D space'
	
	def _get_x(self):
		return self.__x
	def _set_x(self, value):
		self.__x = float(value)
		self.valid_mag = False
	x = property(_get_x, _set_x)
	
	def _get_y(self):
		return self.__y
	def _set_y(self, value):
		self.__y = float(value)
		self.valid_mag = False
	y = property(_get_y, _set_y)
	
	def _get_z(self):
		return self.__z
	def _set_z(self, value):
		self.__z = float(value)
		self.valid_mag = False
	z = property(_get_z, _set_z)
	

	def __init__(self, x, y, z):
		self.__x = 0.0
		self.__y = 0.0
		self.__z = 0.0
		self.x = float(x)
		self.y = float(y)
		self.z = float(z)
		self.valid_mag = False
				
	@classmethod
	def zero(cls):
		return Point3D(0.0, 0.0, 0.0)
	
	@classmethod
	def list(cls, array):
		return Point3D(array[0], array[1], array[2])
	
	def __eq__(self, other):
		if type(self) is type(other):
			return abs(self.x - other.x) <= 0.0001 and abs(self.y - other.y) <= 0.0001 and abs(self.z - other.z) <= 0.0001
		else:
			return False

	def __ne__(self, other):
		return not self.__eq__(other)
	
	def __hash__(self):
		return int(self.x) + int(self.y) * 20 + int(self.z) * 400
	
	def __str__(self):
		return "(%.3g, %.3g, %.3g)" % (float(self.x), float(self.y), float(self.z))
	
	def __repr__(self):
		return "(%.3g, %.3g, %.3g)" % (float(self.x), float(self.y), float(self.z))

	def nparray(self):
		return np.array(self.tolist())
	
	def tolist(self):
		return [self.x, self.y, self.z]
	
	def iteroffsets(self, maxoffset, step=1.0, startpt=None):
		if step == 0.0:
			yield self
			return
		for x in np.arange(-maxoffset, maxoffset + step, step):
			for y in np.arange(-maxoffset, maxoffset + step, step):
				for z in np.arange(-maxoffset, maxoffset + step, step):
					if startpt is not None and (x < startpt.x or (x == startpt.x and (y < startpt.y or (y == startpt.y and z < startpt.z)))): continue
					yield Point3D(self.x + x, self.y + y, self.z + z)

	def iteroffsetx(self, maxoffset, size, step=1.0, startpt=None, test=lambda x: True):
		pts = []
		for pt in self.iteroffsets(maxoffset, step, startpt):
			if len(pts) < size:
				if test(pt) is not False:
					pts.append(pt)
			else:
				yield pts
				del pts[:]

	def iter_randomoffsets(self, maxoffset, count=1):
		if maxoffset == 0.0:
			yield self
			return
		idx = 0
		while idx < count:
			r = random.uniform(0.0, maxoffset)
			theta = random.uniform(0.0, 2 * math.pi)
			phi = random.uniform(0.0, math.pi)
			yield self.add(Point3D(r, theta, phi).tocartesian())
			idx += 1

	def add(self, point):
		if (type(self) is not type(point)):
			print "Wrong type"
			return

		return Point3D(self.x + point.x, self.y + point.y, self.z + point.z)

	def add_spherical(self, point):
		if (type(self) is not type(point)):
			print "Wrong type"
			return
		
		return Point3D(self.x, self.y + point.y, self.z + point.z)

	def subtract(self, point):
		if (type(self) is not type(point)):
			print "Wrong type"
			return
			
		return Point3D(self.x - point.x, self.y - point.y, self.z - point.z)

	def subtract_spherical(self, point):
		if (type(self) is not type(point)):
			print "Wrong type"
			return
		
		return Point3D(self.x, self.y - point.y, self.z - point.z)

	def multiply(self, scalar):
		return Point3D(self.x * scalar, self.y * scalar, self.z * scalar)
	
	def floor(self, step=1.0):
		return Point3D(math.floor(self.x / step) * step, math.floor(self.y / step) * step, math.floor(self.z / step) * step)

	def floor_spherical(self, step=1.0):
		return Point3D(self.x, math.floor(self.y / step) * step, math.floor(self.z / step) * step)
	
	def distanceto(self, point):
		if (type(self) is not type(point)):
			print "Wrong type"
			return 0.0

		return math.sqrt((self.x - point.x) ** 2 + (self.y - point.y) ** 2 + (self.z - point.z) ** 2)

	def anglewith(self, vector):
		assert type(self) is type(vector), "Wrong type"
		if self.magnitude() * vector.magnitude() != 0:
			return math.acos(dotproduct(self, vector) / (self.magnitude() * vector.magnitude()))
		else:
			return float('nan')

	def _update_mag(self):
		if not hasattr(self, "_magfn"): self._magfn = math.sqrt
		self._mag = self._magfn(self.x ** 2 + self.y ** 2 + self.z ** 2)
		self.valid_mag = True

	def magnitude(self):
		if not self.valid_mag:
			self._update_mag()
		return self._mag

	def normalize(self):
		mag = self.magnitude()
		pt = Point3D(self.x / mag, self.y / mag, self.z / mag)
		return pt

	def tospherical(self):
		distance = self.magnitude()
		phi = math.acos(self.z / distance)
		if (phi >= 2 * math.pi):
			phi -= 2 * math.pi
		if (phi < 0.0):
			phi += 2 * math.pi
		if phi == 0.0:
			return Point3D(distance, 0.0, phi)
		theta = math.asin(max(min(self.y / (distance * math.sin(phi)), 1.0), -1.0))	#Within -pi/2 and pi/2
		if (self.x < 0.0):
			theta = math.pi - theta
		if (theta >= 2 * math.pi):
			theta -= 2 * math.pi
		if (theta < 0.0):
			theta += 2 * math.pi
		return Point3D(distance, theta, phi)

	def tocartesian(self):
		return Point3D(self.x * math.sin(self.z) * math.cos(self.y),
					   self.x * math.sin(self.z) * math.sin(self.y),
					   self.x * math.cos(self.z))

	def in_coordinate_system(self, u, v, w):
		'Computes the location of a point in a coordinate space defined by vectors u, v, and w.'
		mag = self.magnitude()
		# A = |P| [ (sqrt(2)/|P|) * P - w ]
		alpha = u.y * v.z - v.y * u.z
		beta = u.z * v.x - u.x * v.z
		gamma = u.x * v.y - v.x * u.y
		
		# A is the self vector FLATTENED onto the xy-plane.
		A = Point3D((beta * (self.x * w.y - self.y * w.x) + gamma * (self.x * w.z - self.z * w.x)) / (alpha * w.x + beta * w.y + gamma * w.z),
					(alpha * (self.y * w.x - self.x * w.y) + gamma * (self.y * w.z - self.z * w.y)) / (beta * w.y + alpha * w.x + gamma * w.z),
					(beta * (self.z * w.y - self.y * w.z) + alpha * (self.z * w.x - self.x * w.z)) / (gamma * w.z + beta * w.y + alpha * w.x))
		dot = dotproduct(u, A)
		# | wx	wy	 wz |
		# | ux	uy	 uz |
		# | Ax	Ay	 Az |
		det = w.x * u.y * A.z + w.y * u.z * A.x + w.z * u.x * A.y - w.z * u.y * A.x - w.y * u.x * A.z - w.x * u.z * A.y;
		theta = math.atan2(det, dot)
		if theta < 0.0:
			theta += 2 * math.pi
		if theta >= 2 * math.pi:
			theta -= 2 * math.pi
		
		if A.magnitude() <= 0.01:
			if dotproduct(w, self) > 0:
				phi = 0.0
			else:
				phi = math.pi
		elif mag * A.magnitude() != 0:
			aAngle = math.acos(max(min(dotproduct(A, self) / (mag * A.magnitude()), 1.0), -1.0))
			if dotproduct(w, self) > 0.0:				#Shows whether P points toward the new z-axis
				phi = math.pi / 2.0 - aAngle
			else:
				phi = math.pi / 2.0 + aAngle
		else:
			phi = math.pi / 2.0							#This means P and A coincide
		
		if phi < 0.0:
			phi += 2 * math.pi
		if phi >= 2 * math.pi:
			phi -= 2 * math.pi
		ret = Point3D(A.magnitude() * math.cos(theta),
					  A.magnitude() * math.sin(theta),
					  mag * math.cos(phi))
		assert math.fabs(ret.magnitude() - self.magnitude()) <= 0.1, "Incorrect in_coordinate_system calculation: %.2f, %.2f, %.2f (%r, %r)" % (mag, theta, phi, str(self), str(A))

		return ret

	def coordinate_system_for_transform(self, endpt):
		"""This function returns (i, j, k), 3 vectors that define a new coordinate system in which self can be expressed as endpt."""
		assert abs(endpt.magnitude() - self.magnitude()) <= 0.01, "Invalid rotation: {} and {} have different magnitudes".format(self, endpt)
		v = crossproduct(endpt.normalize(), self.normalize())
		s = v.magnitude()
		c = dotproduct(endpt.normalize(), self.normalize())
		skewsymmetric_cp = np.matrix( ((0, -v.z, v.y),
									   (v.z, 0, -v.x),
									   (-v.y, v.x, 0)) )
		rot = matlib.eye(3) + skewsymmetric_cp + skewsymmetric_cp * skewsymmetric_cp * ((1 - c) / (s ** 2))
		ptvectors = [endpt.normalize(), Point3D(1.0, 0.0, 0.0), Point3D(0.0, 1.0, 0.0), Point3D(0.0, 0.0, 1.0)]
		newpoints = np.dot([x.nparray() for x in ptvectors], rot.T).tolist()
		rotatedpoint = Point3D(newpoints[0][0], newpoints[0][1], newpoints[0][2]).multiply(self.magnitude())
		i = Point3D(newpoints[1][0], newpoints[1][1], newpoints[1][2]).normalize()
		j = Point3D(newpoints[2][0], newpoints[2][1], newpoints[2][2]).normalize()
		k = crossproduct(i, j).normalize()
		return (i, j, k)

	def rotation_matrix(self, new_vec):
		#Rotation - find the rotation matrix required to rotate orig_vec onto new_vec (see http://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d)
		v = crossproduct(self, new_vec)
		s = v.magnitude()
		c = dotproduct(self, new_vec)
		skewsymmetric_cp = np.matrix( ((0, -v.z, v.y),
									   (v.z, 0, -v.x),
									   (-v.y, v.x, 0)) )
		rot = matlib.eye(3) + skewsymmetric_cp + skewsymmetric_cp * skewsymmetric_cp * ((1 - c) / (s ** 2))
		return rot

	def rotate(self, rotation_matrix):
		new = rotation_matrix.dot(self.nparray()).tolist()
		return Point3D.list(new[0])

def midpoint(point1, point2):
	return Point3D((point1.x + point2.x) / 2.0, (point1.y + point2.y) / 2.0, (point1.z + point2.z) / 2.0)

def dotproduct(point1, point2):
	if (type(point1) is not type(point2)):
		print "Wrong type"
		return point1
	return point1.x * point2.x + point1.y * point2.y + point1.z * point2.z

def crossproduct(point1, point2):
	if (type(point1) is not type(point2)):
		raise TypeError("Cross product requires two Point3D objects")
	# | i	j	k |	i	j
	# | x1	y1	z1|	x1	y1
	# | x2	y2	z2|	x2	y2
	return Point3D(point1.y * point2.z - point1.z * point2.y,
				   point1.z * point2.x - point1.x * point2.z,
				   point1.x * point2.y - point1.y * point2.x)

# position zones

def zone_point(point):
	return (math.floor(point.x), math.floor(point.y), math.floor(point.z))

def zdistance(zone1, zone2):
	if (type(zone1) is not type(zone2)):
		print "Wrong type"
		return 0.0
		
	deltaX = math.fabs(zone1[0] - zone2[0])
	deltaY = math.fabs(zone1[1] - zone2[1])
	deltaZ = math.fabs(zone1[2] - zone2[2])
	return math.sqrt(deltaX ** 2 + deltaY ** 2 + deltaZ ** 2)

class PositionZone(object):

	def __init__(self, alpha=Point3D.zero(), x=Point3D.zero(), y=Point3D.zero(), z=Point3D.zero()):
		self.alpha_zone = alpha
		self.x_axis = x
		self.y_axis = y
		self.z_axis = z
		if self.alpha_zone != Point3D.zero() and self.x_axis != Point3D.zero() and self.y_axis != Point3D.zero() and self.z_axis != Point3D.zero():
			self.hash = self.calchash()
	
	def __str__(self):
		return "PZ {%r, %r, %r, %r}" % (str(self.alpha_zone), str(self.x_axis), str(self.y_axis), str(self.z_axis))

	def __repr__(self):
		return "PZ {%r, %r, %r, %r}" % (str(self.alpha_zone), str(self.x_axis), str(self.y_axis), str(self.z_axis))
	
	def calchash(self):
		alpha = int(self.alpha_zone.x) + int(self.alpha_zone.y) * 20 + int(self.alpha_zone.z) * 400
		x = int((20 ** 3) * 10 * (float(self.x_axis.x) + float(self.x_axis.y) * 20 + float(self.x_axis.z) * 400))
		y = int((20 ** 6) * 10 * (float(self.y_axis.x) + float(self.y_axis.y) * 20 + float(self.y_axis.z) * 400))
		z = int((20 ** 9) * 10 * (float(self.z_axis.x) + float(self.z_axis.y) * 20 + float(self.z_axis.z) * 400))
		return int(alpha + x + y + z)

	def __hash__(self):
		if not hasattr(self, "hash"): self.hash = self.calchash()
		return self.hash
	
	def __eq__(self, pz):
		if type(self) is not type(pz): return False
		return self.alpha_zone == pz.alpha_zone and self.x_axis == pz.x_axis and self.y_axis == pz.y_axis and self.z_axis == pz.z_axis
	
	def __ne__(self, pz):
		return not self.__eq__(pz)

	def stringforfile(self):
		return "%d, %d, %d; %.1f, %.1f, %.1f; %.1f, %.1f, %.1f; %.1f, %.1f, %.1f\n" % (self.alpha_zone.x, self.alpha_zone.y, self.alpha_zone.z, self.x_axis.x, self.x_axis.y, self.x_axis.z, self.y_axis.x, self.y_axis.y, self.y_axis.z, self.z_axis.x, self.z_axis.y, self.z_axis.z)
