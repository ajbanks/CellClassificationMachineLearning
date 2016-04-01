import sys
import numpy as np
from RNASeqData import RNASeqData
import preprocess
# import SCC_Rna_seq_NaiveBayes


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