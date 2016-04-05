import numpy as np
from sklearn.naive_bayes import GaussianNB

gnb = GaussianNB()

# Source: http://scikit-learn.org/stable/modules/generated/sklearn.naive_bayes.GaussianNB.html#sklearn.naive_bayes.GaussianNB.fit
def fitTrainingData(training_data, nSamples):
	print "\nfitting training data to gaussian naive bayes classifier"

	# fit(X,Y)
		# x = 2D array of cells
		# y = 1D list of classifiers i.e. [group1, group1, group2, group2]

	# makes lists into np arrays
	training_dataNP = np.array(training_data)
	nSamplesNP = np.array(nSamples)

	# fit gaussian naive bayes according to samples and n features
	gnb.fit(training_dataNP, nSamplesNP)

def predictTestData(testing_data):
	print "\npredicting test data on gaussian naive bayes classifier"
	
	# make list into np array
	testing_dataNP = np.array(testing_data)

	# predict the values
	predicted = gnb.predict(testing_dataNP)

	return predicted



# a = np.array(a)

#  Naive Bayes Notes
# 	- P(outcome | known evidence)
#	- To predict an outcome given multiple evidence i.e. multiple cell types,
#	    we uncouple multiple pieces of evidence and teach each as independent (hence 'naive')
#		P(outcome | multiple evidence) = P(evidence1 | outcome) * P(evidence2 | outcome) * ... * P(evidenceN | outcome)
#		    *** scaled by P(multiple evidence)
#
#
#		P(outcome | evidence) = P(likelihood of evidence) * Prior prob. of outcome/ P(evidence)
#
#
#	- Gaussian: likelihood of features is assumed to gaussian
#
#	Sources: http://stackoverflow.com/questions/10059594/a-simple-explanation-of-naive-bayes-classification