from aminoacids import *
import random
from permissions import *

def random_axes(previous_aa):
	if previous_aa is None:
		i = Point3D(random.uniform(-1.0, 1.0), random.uniform(-1.0, 1.0), random.uniform(-1.0, 1.0))
		# Choose a random j such that i dot j = 0 -> ix * jx + iy * jy + iz * jz = 0
		j = Point3D(random.uniform(-1.0, 1.0), random.uniform(-1.0, 1.0), 0.0)
		j.z = -(i.x * j.x + i.y * j.y) / i.z
		j.normalize()
		k = crossproduct(i, j)
		return (i.normalize(), j, k.normalize())
	else:
		i = previous_aa.toglobal(Point3D(1.0, random.uniform(-math.pi * 0.1, math.pi * 0.1), random.uniform(math.pi * 0.45, math.pi * 0.55)).tocartesian()).subtract(previous_aa.acarbon)
		# Choose a random j such that i dot j = 0 -> ix * jx + iy * jy + iz * jz = 0
		diff = 0.1
		j = Point3D(random.uniform(previous_aa.j.x - diff, previous_aa.j.x + diff), random.uniform(previous_aa.j.y - diff, previous_aa.j.y + diff), 0.0)
		j.z = -(i.x * j.x + i.y * j.y) / i.z
		k = crossproduct(i, j)
		return (i.normalize(), j.normalize(), k.normalize())

def generate_randomcoil(sequence, permissions=None, steric_cutoff=3.0):
	'''This function generates a self-avoiding random walk in 3D space.'''
	current_pt = Point3D(random.uniform(-1.0, 1.0), random.uniform(-1.0, 1.0), random.uniform(-1.0, 1.0))
	last_aa = None
	allaa = []
	hashtable = AAHashTable()
	bondangle = math.acos(-1.0 / 3.0)
	for i, acode in enumerate(sequence):
		aminoacid = AminoAcid(aatypec(acode), tag=i,acarbon=current_pt)
		axes = []
		if last_aa is not None:
			#Choose a point to the right of the amino acid
			idx = 0
			while len(hashtable.nearby_aa(aminoacid, steric_cutoff, consec=False)) > 0 or last_aa.acarbon.distanceto(current_pt) <= 2.0:
				if idx == 100: return generate_randomcoil(sequence, permissions)
				#print i, str(current_pt), hashtable.nearby_aa(aminoacid, 2.0, consec=False)
				if permissions == None:
					current_pt = last_aa.toglobal(Point3D(random.uniform(3.0, 3.5), random.uniform(bondangle - 0.2, bondangle + 0.2), random.uniform(bondangle - 0.2, bondangle + 0.2)).tocartesian())
				else:
					candidates = permissions.allowed_conformations(aminoacid, last_aa)
					if not len(candidates):
						print "No permissible candidates"
						continue
					zone = random.choice(candidates)
					current_pt = zone.alpha_zone
					axes = (zone.x_axis, zone.y_axis, zone.z_axis)
				aminoacid.acarbon = current_pt
				idx += 1
		if permissions == None or len(axes) == 0:
			axes = random_axes(last_aa)
		aminoacid.set_axes(*axes, normalized=False)
		allaa.append(aminoacid)
		hashtable.add(aminoacid)
		last_aa = aminoacid
	del hashtable
	return allaa

def directed_randomwalk(startpt, endpt, numzones, steplength=1.0):
	"""This helper function generates a self-avoiding random walk from startpt to endpt and returns an array of position zones corresponding to each point along the way."""
	assert startpt.distanceto(endpt) < steplength * (numzones + 1), "Cannot make a random walk over {} angstroms with {} steps of length {}.".format(startpt.distanceto(endpt), numzones + 1, steplength)

	current_pt = startpt
	zones = []
	endpt = endpt.subtract(startpt)
	for idx in xrange(numzones):
		#Sample a point that is at distance steplength from current_pt, but within (numzones - idx) * steplength of endpt.
		#This can be modeled by the intersection of two spheres; first let's find the dimensions of this circular intersection. (See Notability, Science Fair Calculations)
		r1 = steplength
		r2 = (numzones - idx) * steplength
		t = (r1 ** 2 - r2 ** 2 + endpt.x ** 2 + endpt.y ** 2 + endpt.z ** 2) / (2 * (endpt.x ** 2 + endpt.y ** 2 + endpt.z ** 2))
		center = Point3D(endpt.x * t, endpt.y * t, endpt.z * t)
		print r1, r2, t, center
		#Something is incorrect here. I think it might be t above.
		assert r1 ** 2 - (t ** 2) * (endpt.x ** 2 + endpt.y ** 2 + endpt.z ** 2) > 0, "Incorrect math: endpoint {}".format(endpt)
		radius = math.sqrt(r1 ** 2 - (t ** 2) * (endpt.x ** 2 + endpt.y ** 2 + endpt.z ** 2))
		
		#Now choose a random point on the xy-plane.
		point2d = Point3D(random.uniform(0.0, radius), random.uniform(0.0, math.pi * 2.0), math.pi / 2.0).tocartesian()
		#Project this point onto the circular intersection region.
		K = (endpt.x * point2d.x + endpt.y * point2d.y) / (endpt.x ** 2 + endpt.y ** 2 + endpt.z ** 2) - t
		pt_projected = Point3D(point2d.x + K * endpt.x, point2d.y + K * endpt.y, K * endpt.z).add(Point3D(t * endpt.x, t * endpt.y, t * endpt.z))
		#Now turn this point into a vector with length steplength.
		pt_projected = pt_projected.multiply(steplength / pt_projected.magnitude())
		current_pt = current_pt.add(pt_projected)
		zones.append(PositionZone(current_pt))
		endpt = endpt.subtract(pt_projected)
	return zones