import sys
import numpy as np
from RNASeqData import RNASeqData
import preprocess
import guassianNB_RNASeq
import analysis


# Resource: http://machinelearningmastery.com/get-your-hands-dirty-with-scikit-learn-now/
# Python for Java Programmers: http://python4java.necaiseweb.org/Fundamentals/TheBasics

if __name__ == '__main__':
	print "start"

	# initialize the data set class
	data = RNASeqData()
	
	# read raw RNA seq data into memory
	data.setRawData(preprocess.loadRawData(data.rawDataFileName))

	# read cell identifier annotations into memory
	data.setCellIdentifierAnnotations(preprocess.loadCellIdentifierAnnotations(data.annotationsFileName, 
		data.getNumCellsRaw()))

	# read molecule count annotations into memory
	data.setMoleculeCountAnnotations(preprocess.loadMoleculeCountAnnotations(data.annotationsFileName,
		data.getNumCellsRaw()))

	# down sample the data by cluster size --> MAKE DOWN SAMPLING A CLA OPTION
	#	 i.e. scale all cluster size to the smallest cluster (by number of cells)
	# save down sampled data and random indices for accessing corresponding annotations
	downSampleClusterData, randIndices = preprocess.downSampleByClusterSize(data.getRawData(), 
		data.getCellIdentifierAnnotations())

	# add the data and random indices reference to the data class
	data.setDSClusterData(downSampleClusterData)
	data.setRandIndicesFromDS(randIndices)

	# down sample the data by the cell with the least number of molecules
	data.setDSCluster_MoleculeData(preprocess.downSampleByMoleculeCount(data.getDSClusterData(),
		data.getMoleculeCountAnnotations(), data.getRandIndices()))

	# DON'T NEED THIS ANYMORE
	# feature extraction --> find the avg type of each cluster
	# data.setFeatures(preprocess.featureExtraction(data.getDSCluster_MoleculeData(),
		# data.getCellIdentifierAnnotations(), data.getRandIndices()))

	# partition the data set into 70% training and 30% testing
	data.makeTrainingAndTestingData()

	# fit training data to gaussian nb
	guassianNB_RNASeq.fitTrainingData(data.getTrainingData(), data.getTargetValues())

	# predict values using gaussian nb
	guassianNB_predictionResults = guassianNB_RNASeq.predictTestData(data.getTestingData())
	
	# analyze results of guassian nb
	analysis.analyzeResults("Guassian Naive Bayes", guassianNB_predictionResults, data.getTestingDataTargetValues())


	print "\nexiting"
