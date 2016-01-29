"""This module takes a method for plane regression from http://stackoverflow.com/questions/20699821/find-and-draw-regression-plane-to-a-set-of-points."""

import numpy as np
import scipy.optimize
import functools

def plane(x, y, params):
    return params[0] * x + params[1] * y + params[2]

def error(params, points):
	result = 0
	for (x,y,z) in points:
		plane_z = plane(x, y, params)
		diff = abs(plane_z - z)
		result += diff**2
	return result

def plane_regression(points):
	print points
	fun = functools.partial(error, points=points)
	params0 = [0, 0, 0]
	res = scipy.optimize.minimize(fun, params0)
	print fun(res.x)
	return res.x

import statsmodels.api as sm

def reg_m(y, x):
	
	ones = np.ones(len(x[0]))
	X = sm.add_constant(np.column_stack((x[0], ones)))
	for ele in x[1:]:
		X = sm.add_constant(np.column_stack((ele, X)))
	results = sm.OLS(y, X).fit()
	return results