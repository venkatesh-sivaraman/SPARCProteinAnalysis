import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.pyplot as plt

def fit_maxwell(input):
	print os.path.basename(input)
	data = []
	with open(input) as file:
		for line in file:
			data.append(float(line.strip()))
	#total = sum(data)
	#data = [d / total for d in data]
	#print data

	maxwell = stats.maxwell
	
	params = maxwell.fit(data, floc=0)
	print params
	d, p = stats.kstest(data, "maxwell", mode="asymp")
	print d, p
	norm = stats.norm

	params = norm.fit(data) #, floc=0)
	print params
	d, p = stats.kstest(data, "norm", mode="asymp")
	print d, p

	plt.hist(data, bins=50, normed=True, alpha=0.6, color='g')

	# Plot the PDF.
	xmin, xmax = plt.xlim()
	x = np.linspace(xmin, xmax, 100)
	p = norm.pdf(x, params[0], params[1])
	plt.plot(x, p, 'k', linewidth=2)
	title = "Fit results: mu = %.2f,  std = %.2f" % params
	plt.title(title)

	plt.show()