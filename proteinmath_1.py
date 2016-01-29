import math

class Point3D(object):
	'Represents a point or vector in 3D space'
	
	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z
				
	@classmethod
	def zero(cls):
		return Point3D(0.0, 0.0, 0.0)
	
	def __eq__(self, other):
		if type(self) is type(other):
			return self.x == other.x and self.y == other.y and self.z == other.z
		else:
			return 0

	def __ne__(self, other):
		return not self.__eq__(other)
	
	def __str__(self):
		return "Point (%g, %g, %g)" % (self.x, self.y, self.z)

	def add(self, point):
		if (type(self) is not type(point)):
			print "Wrong type"
			return

		return Point3D(self.x + point.x, self.y + point.y, self.z + point.z)

	def subtract(self, point):
		if (type(self) is not type(point)):
			print "Wrong type"
			return
			
		return Point3D(self.x - point.x, self.y - point.y, self.z - point.z)

	def multiply(self, scalar):
		return Point3D(self.x * scalar, self.y * scalar, self.z * scalar)

	def distanceto(self, point):
		if (type(self) is not type(point)):
			print "Wrong type"
			return 0.0

		return math.sqrt((self.x - point.x) ** 2 + (self.y - point.y) ** 2 + (self.z - point.z) ** 2)

	def magnitude(self):
		return math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

	def normalize(self):
		mag = self.magnitude()
		return Point3D(self.x / mag, self.y / mag, self.z / mag)

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
		mag = self.magnitude();
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
		theta = math.atan2(det, dot);
		if theta < 0.0:
			theta += 2 * math.pi
		if theta >= 2 * math.pi:
			theta -= 2 * math.pi

		aAngle = math.acos(max(min(dotproduct(A, self) / (mag * A.magnitude()), 1.0), -1.0))
		if (aAngle == aAngle):
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
		if math.fabs(ret.magnitude() - self.magnitude()) >= 0.01:
			print "Problem"

		return ret


def midpoint(point1, point2):
	return Point3D((point1.x + point2.x) / 2.0, (point1.y + point2.y) / 2.0, (point1.z + point2.z) / 2.0)

def dotproduct(point1, point2):
	if (type(point1) is not type(point2)):
		print "Wrong type"
		return point1
	return point1.x * point2.x + point1.y * point2.y + point1.z * point2.z

# position zones

POSITION_ZONE_INTERVAL = 1.0
POSITION_ZONE_EXTENSION = 10

def zone_point(point):
	return ((POSITION_ZONE_EXTENSION * 2.0) * (POSITION_ZONE_EXTENSION * 2.0)) * (math.floor(point.z) + POSITION_ZONE_EXTENSION) + (POSITION_ZONE_EXTENSION * 2.0) * (math.floor(point.y) + POSITION_ZONE_EXTENSION) + (math.floor(point.x) + POSITION_ZONE_EXTENSION)

class PositionZone(object):

	def __init__(self, alpha, beta):
		self.alpha_zone = alpha
		self.beta_zone = beta
	
	def __str__(self):
		return "Position zone (%i, %i)" % (self.alpha_zone, self.beta_zone)

	def distanceto(self, zone):
		if (type(self) is not type(zone)):
			print "Wrong type"
			return 0.0
		
		deltaX = math.fabs((zone1 % (POSITION_ZONE_EXTENSION * 2)) - (zone2 % (POSITION_ZONE_EXTENSION * 2)))
		deltaY = math.fabs(math.floor((zone1 % ((POSITION_ZONE_EXTENSION * 2) * (POSITION_ZONE_EXTENSION * 2))) / (POSITION_ZONE_EXTENSION * 2.0)) - math.floor((zone2 % ((POSITION_ZONE_EXTENSION * 2) * (POSITION_ZONE_EXTENSION * 2))) / (POSITION_ZONE_EXTENSION * 2.0)))
		deltaZ = math.fabs(math.floor(zone1 / ((POSITION_ZONE_EXTENSION * 2.0) * (POSITION_ZONE_EXTENSION * 2.0))) - math.floor(zone2 / ((POSITION_ZONE_EXTENSION * 2.0) * (POSITION_ZONE_EXTENSION * 2.0))))
		return math.sqrt(deltaX ** 2 + deltaY ** 2 + deltaZ ** 2)

