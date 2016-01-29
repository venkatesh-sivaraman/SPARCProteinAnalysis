from aminoacids import *
from proteinmath import *
import math
import random

def old_stuff():
	print str(aa)
	proximity_nonconsec = protein.nearby_aa(aa, distance, index=i, consec=False)
	proximity_consec = protein.nearby_aa(aa, distance, index=i, consec=True)
	probabilities = {} #[] for place_elements_sorted
	for dist in distributions:
		if dist.type == frequency_nonconsec_disttype:
			for aaprox in proximity_nonconsec:
				#place_elements_sorted(probabilities, dist.probabilities(aa, aaprox).items(), sorter=probsort)
				probabilities.update(dist.probabilities(aa, aaprox))
		elif dist.type == frequency_consec_disttype:
			for aaprox in proximity_consec:
				#place_elements_sorted(probabilities, dist.probabilities(aa, aaprox).items(), sorter=probsort)
				probabilities.update(dist.probabilities(aa, aaprox))
		elif dist.type == frequency_disttype:
			for aaprox in proximity_consec + proximity_nonconsec:
				#place_elements_sorted(probabilities, dist.probabilities(aa, aaprox).items(), sorter=probsort)
				probabilities.update(dist.probabilities(aa, aaprox))
	print len(probabilities)
	'''
		#This code is supposed to block together nearby scenarios and add their probabilities. But is it really worth it?
		hashtable = PointHashTable(probabilities)
		finalprobs = {}
		i = 0
		while len(probabilities) > 0:
		point = random.choice(probabilities.keys())
		nearby = hashtable.nearby_points(point, 1.0)
		finalprobs[point.floor()] = sum([probabilities[x] for x in nearby])
		for donepoint in nearby:
		probabilities.pop(donepoint)
		hashtable.remove(donepoint)
		print len(finalprobs)
		'''

simulation_frequency = 100 #iterations per microsecond

#arbitrarily allow 5 angstroms at 300 K
def farthest_distance_to_travel(temperature):
	return temperature / 60.0

def random_point_in_vicinity(point1, closepoint, proximity, maxdistance):
	#Give a random solution to the inequalities
	#sqrt((ret.x - closepoint.x)^2 + (ret.y - closepoint.y)^2 + (ret.z - closepoint.z)^2) < proximity
	#sqrt((ret.x - point1.x)^2 + (ret.y - point1.y)^2 + (ret.z - point1.z)^2) < maxdistance
	#by performing rejection sampling.
	random_point = Point3D(10000000, 10000000, 10000000)
	xmin = min(point1.x, closepoint.x - proximity)
	xmax = max(point1.x + maxdistance, closepoint.x)
	ymin = min(point1.y, closepoint.y - proximity)
	ymax = max(point1.y + maxdistance, closepoint.y)
	zmin = min(point1.z, closepoint.z - proximity)
	zmax = max(point1.z + maxdistance, closepoint.z)
	'''random_point.x ** 2 + random_point.y ** 2 + random_point.z ** 2 - 2 * (closepoint.x * random_point.x + closepoint.y * random_point.y + closepoint.z * random_point.z) < proximity ** 2 - (closepoint.x ** 2 + closepoint.y ** 2 + closepoint.z ** 2) && random_point.x ** 2 + random_point.y ** 2 + random_point.z ** 2 - 2 * (point1.x * random_point.x + point1.y * random_point.y + point1.z * random_point.z) < maxdistance ** 2 - (point1.x ** 2 + point1.y ** 2 + point1.z ** 2)'''
	while random_point.x ** 2 + random_point.y ** 2 + random_point.z ** 2 - 2 * (closepoint.x * random_point.x + closepoint.y * random_point.y + closepoint.z * random_point.z) >= proximity ** 2 - (closepoint.x ** 2 + closepoint.y ** 2 + closepoint.z ** 2) or random_point.x ** 2 + random_point.y ** 2 + random_point.z ** 2 - 2 * (point1.x * random_point.x + point1.y * random_point.y + point1.z * random_point.z) >= maxdistance ** 2 - (point1.x ** 2 + point1.y ** 2 + point1.z ** 2):
		#print "Recurring", random_point, closepoint, random_point.x ** 2 + random_point.y ** 2 + random_point.z ** 2 - 2 * (closepoint.x * random_point.x + closepoint.y * random_point.y + closepoint.z * random_point.z) + (closepoint.x ** 2 + closepoint.y ** 2 + closepoint.z ** 2), proximity ** 2
		random_point = Point3D(random.uniform(xmin, xmax), random.uniform(ymin, ymax), random.uniform(zmin, zmax))
	return random_point

def folding_iteration(aalist, temperature, r, scoresys):
	lastCarbon = Point3D.zero()
	for aa in aalist:
		distance = random.random() * farthest_distance_to_travel(temperature)
		theta = random.uniform(0, 2 * math.pi)
		phi = random.uniform(0, math.pi)
		new_point = random_point_in_vicinity(aa.acarbon, lastCarbon, 1.4, farthest_distance_to_travel(temperature)) #aa.acarbon.add(Point3D(distance, theta, phi).tocartesian())
		lastCarbon = new_point
		if (aalist.index(aa) == 0):
			print "Last carbon:", lastCarbon
		if (aalist.index(aa) == 1):
			print new_point

# aalist = array of amino acids, random = boolean indicating whether a score system will be used, scoresys = a score system object, duration = number of microseconds to simulate
def simulate_folding(aalist, temperature, random, scoresys, duration):
	for i in range(simulation_frequency * temperature):
		folding_iteration(aalist, temperature, random, scoresys)
