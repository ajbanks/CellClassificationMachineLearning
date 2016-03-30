import sys
import numpy as np
from RNASeqData import RNASeqData
import preprocess
# import SCC_Rna_seq_NaiveBayes


if __name__ == '__main__':
	print "start"

	# initialize the data set class
	data = RNASeqData()
	
	# read raw RNA seq data into memory
	data.setRawData(preprocess.loadRawData(data.rawDataFileName))

	# read cell identifier annotations into memory
	data.setCellIdentifierAnnotations(preprocess.loadCellIdentifierAnnotations(data.annotationsFileName, 
		data.getNumCells()))

	# read molecule count annotations into memory
	data.setMoleculeCountAnnotations(preprocess.loadMoleculeCountAnnotations(data.annotationsFileName,
		data.getNumCells()))

	preprocess.downSampleByClusterSize(data.getRawData(), data.getCellIdentifierAnnotations())

# Resource: http://machinelearningmastery.com/get-your-hands-dirty-with-scikit-learn-now/
# Python for Java Programmers: http://python4java.necaiseweb.org/Fundamentals/TheBasics


