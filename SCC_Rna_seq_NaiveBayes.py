from sklearn import preprocessing
import numpy as np

def run():
	print "running Naive Bayes classification on Single Cell RNA Seq cata"

	x = np.array([
		[1., 2., -.5, 0., -2., 4., 1.2, .6],
		[2., 2., 1.5, 0., 0., .3, -.9, .3],
		[5., 2.1, -1.5, 0.2, -2.1, .1, .2, 1.6]
		])
	# x = ([[0.1, 1.5, -.5],
	# 	[15., 10., 2.]])
	
	print "Unprocessed data: "
	print x
	print ""

	x_scaled = preprocessing.scale(x)

	print "Processed data: "
	print x_scaled
	print ""