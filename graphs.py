"""This module runs as a separate target using Anaconda."""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.collections as col
import matplotlib.path as path
import numpy as np
import math
from proteinmath import *
from mpl_toolkits.mplot3d import Axes3D

def beginning():
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	plt.show()

def read_line_segment_data(input):
	data = []
	with open(input) as file:
		for line in file:
			line = line.replace(" ", "")
			comps = line.split(';')
			if len(comps) >= 2:
				x1, y1, z1 = (float(x) for x in comps[0].split(','))
				x2, y2, z2 = (float(x) for x in comps[1].split(','))
				#if abs(Point3D(x1, y1, z1).distanceto(Point3D.zero()) - 6.0) >= 1.0: continue
				if float(comps[2]) <= 100: continue
				x2 += x1; y2 += y1; z2 += z1
				data.append([[x1, x2], [y1, y2], [z1, z2]])
	return data

def read_quiver_data(input):
	data = [ [], [], [], [], [], [] ]
	with open(input) as file:
		for line in file:
			line = line.replace(" ", "")
			comps = line.split(';')
			if len(comps) >= 2:
				x1, y1, z1 = (float(x) for x in comps[0].split(','))
				x2, y2, z2 = (float(x) for x in comps[1].split(','))
				#if abs(Point3D(x1, y1, z1).distanceto(Point3D.zero()) - 6.0) >= 1.0: continue
				if float(comps[2]) <= 100: continue
				data[0].append(x1)
				data[1].append(y1)
				data[2].append(z1)
				data[3].append(x2)
				data[4].append(y2)
				data[5].append(z2)
	return data

def threed_lines(data):
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	'''for x in np.arange(0.0, 10.0, 1.0):
		for y in np.arange(0.0, 10.0, 1.0):
			ax.plot([x, x], [y, y], 'b-', zs=[math.sin(x) + math.sin(y), math.sin(x) + math.sin(y) + 1])'''
	for xs, ys, zs in data:
		ax.plot(xs, ys, 'b-', zs=zs)

	ax.set_xlabel('X')
	ax.set_xlim3d(-7, 7)
	ax.set_ylabel('Y')
	ax.set_ylim3d(-7, 7)
	ax.set_zlabel('Z')
	ax.set_zlim3d(-7, 7)

	plt.show()

def threed_quivers(data):
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	
	'''for x in np.arange(0.0, 10.0, 1.0):
		for y in np.arange(0.0, 10.0, 1.0):
		ax.plot([x, x], [y, y], 'b-', zs=[math.sin(x) + math.sin(y), math.sin(x) + math.sin(y) + 1])'''
	ax.quiver(*data)
	
	'''ax.set_xlabel('X')
	ax.set_xlim3d(-7, 7)
	ax.set_ylabel('Y')
	ax.set_ylim3d(-7, 7)
	ax.set_zlabel('Z')
	ax.set_zlim3d(-7, 7)'''

	plt.show()

if __name__ == '__main__':
	threed_quivers(read_quiver_data("/Volumes/External Hard Drive/Science Fair 2014-15/rankings/consec-axes/around-6-a-yaxes.txt"))