import sys
import numpy as np
from RNASeqData import RNASeqData
import preprocess
# import SCC_Rna_seq_NaiveBayes


if __name__ == '__main__':
	print "start"

	# initialize the data set class
	data = RNASeqData()
	
	# 
	data.setData(preprocess.run(data.rawDataFileName))

	data.getData()

# Resource: http://machinelearningmastery.com/get-your-hands-dirty-with-scikit-learn-now/
# Python for Java Programmers: http://python4java.necaiseweb.org/Fundamentals/TheBasics