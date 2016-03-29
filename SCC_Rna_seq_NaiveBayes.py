import numpy as np
from sklearn.naive_bayes import GaussianNB

# Source: http://scikit-learn.org/stable/modules/generated/sklearn.naive_bayes.GaussianNB.html#sklearn.naive_bayes.GaussianNB.fit
def fitTrainingData(training_data, nSamples):
	gnb = GaussianNB() # can I make this class/file wide?

	# fit gaussian naive bayes according to n samples and n features
	gnb.fit(training_data, nSamples)

def predictTestData(testing_data):
	# gnb.predit()



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