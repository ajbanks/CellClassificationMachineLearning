import sys
import random

class RNASeqData(object):
	rawDataFileName = "GSE60361C13005Expression.txt"
	annotationsFileName = "expressionmRNAAnnotations.txt"

	# annotationsFileName = "";

	def __init__(self):
		print "\ninitializing RNASeqData"

	# initial raw data read from database file
	def setRawData(self, rawData):
		self.rawData = rawData

	# data that has been down sampled by cluster size
	def setDSClusterData(self, dsClusterData):
		self.dsClusterData = dsClusterData

	# data that has been down sampled by molecule count
	def setDSCluster_MoleculeData(self, dsCluster_MoleculeData):
		self.dsCluster_MoleculeData = dsCluster_MoleculeData

	# annotations that identify the group (1-9) of each cell
	def setCellIdentifierAnnotations(self, cellIdentifierAnnotations):
		self.cellIdentifierAnnotations = cellIdentifierAnnotations

	# annotations that list the molecule count for each cell
	def setMoleculeCountAnnotations(self, moleculeCountAnnotations):
		self.moleculeCountAnnotations = moleculeCountAnnotations

	def setRandIndicesFromDS(self, randIndices):
		# put the indices in ascending order
		self.randIndices = sorted(randIndices)

	def setFeatures(self, features):
		self.features = features

	def makeTrainingAndTestingData(self):
		print "\npartitioning data set - 70% training, 30% testing"
		# randomly selecte 70% of each cluster for training, 30% for training

		trainingData = []
		testingData = []

		# initialize lists to holds cells of each type
		type1 = []
		type2 = []
		type3 = []
		type4 = []
		type5 = []
		type6 = []
		type7 = []
		type8 = []
		type9 = []

		# initialize lists to hold parallel indices of cells of each type
		type1Idxs = []
		type2Idxs = []
		type3Idxs = []
		type4Idxs = []
		type5Idxs = []
		type6Idxs = []
		type7Idxs = []
		type8Idxs = []
		type9Idxs = []


		# iterate over random indices and add cells to corresponding lists
		iterator = 0
		for idx in self.randIndices:
			# extract type from cell identifier annotations
			_type = int(self.cellIdentifierAnnotations[idx])

			# random indices and down sample data are parallel, so add the corresponding cell to the list
			if _type == 1:
				type1.append(self.dsCluster_MoleculeData[iterator])
				type1Idxs.append(idx)
			elif _type == 2:
				type2.append(self.dsCluster_MoleculeData[iterator])
				type2Idxs.append(idx)
			elif _type == 3:
				type3.append(self.dsCluster_MoleculeData[iterator])
				type3Idxs.append(idx)
			elif _type == 4:
				type4.append(self.dsCluster_MoleculeData[iterator])
				type4Idxs.append(idx)
			elif _type == 5:
				type5.append(self.dsCluster_MoleculeData[iterator])
				type5Idxs.append(idx)
			elif _type == 6:
				type6.append(self.dsCluster_MoleculeData[iterator])
				type6Idxs.append(idx)
			elif _type == 7:
				type7.append(self.dsCluster_MoleculeData[iterator])
				type7Idxs.append(idx)
			elif _type == 8:
				type8.append(self.dsCluster_MoleculeData[iterator])
				type8Idxs.append(idx)
			elif _type == 9:
				type9.append(self.dsCluster_MoleculeData[iterator])
				type9Idxs.append(idx)

			iterator += 1

		# make sure all types have the same number of cells
		numCells = len(type1)
		if len(type2) != numCells or len(type3) != numCells or len(type4) != numCells or len(type5) != numCells or len(type6) != numCells or len(type7) != numCells or len(type8) != numCells or len(type9) != numCells:
			print "error: not all clusters have {numCells} cells".format(numCells=type1)
			return

		# randomly choose 70% of the indices of each cluster for training data, remaining for testing data

		# generate a list of indices 0-numCells to use to randomly select cells from each type
		indices = range(numCells)

		# find number of cells which constitue 70% of a cluster
		numTrainingCellsPerCluster = int(numCells * .7)
		
		# randomly select 70% of the indices of each type
		type1TrainingCellIdxs = random.sample(indices, numTrainingCellsPerCluster)
		type2TrainingCellIdxs = random.sample(indices, numTrainingCellsPerCluster)
		type3TrainingCellIdxs = random.sample(indices, numTrainingCellsPerCluster)
		type4TrainingCellIdxs = random.sample(indices, numTrainingCellsPerCluster)
		type5TrainingCellIdxs = random.sample(indices, numTrainingCellsPerCluster)
		type6TrainingCellIdxs = random.sample(indices, numTrainingCellsPerCluster)
		type7TrainingCellIdxs = random.sample(indices, numTrainingCellsPerCluster)
		type8TrainingCellIdxs = random.sample(indices, numTrainingCellsPerCluster)
		type9TrainingCellIdxs = random.sample(indices, numTrainingCellsPerCluster)

		# initialize type training lists for each cell (need these to keep parallel with indices)
		type1TrainingCells = []
		type2TrainingCells = []
		type3TrainingCells = []
		type4TrainingCells = []
		type5TrainingCells = []
		type6TrainingCells = []
		type7TrainingCells = []
		type8TrainingCells = []
		type9TrainingCells = []

		# initialize type training lists for each cell idx that correspond to annotation (need these to keep parallel with cells)
		type1TrainingCellIdxsAnn = []
		type2TrainingCellIdxsAnn = []
		type3TrainingCellIdxsAnn = []
		type4TrainingCellIdxsAnn = []
		type5TrainingCellIdxsAnn = []
		type6TrainingCellIdxsAnn = []
		type7TrainingCellIdxsAnn = []
		type8TrainingCellIdxsAnn = []
		type9TrainingCellIdxsAnn = []

		# initialize type testings list for each cell idx that correspond to annotations to use as answer key for analysis
		type1TestingCellIdxsAnn = []
		type2TestingCellIdxsAnn = []
		type3TestingCellIdxsAnn = []
		type4TestingCellIdxsAnn = []
		type5TestingCellIdxsAnn = []
		type6TestingCellIdxsAnn = []
		type7TestingCellIdxsAnn = []
		type8TestingCellIdxsAnn = []
		type9TestingCellIdxsAnn = []

		# iterate over cells in each type, add it to training or testing according to whether its index was randomly selected
		iterator = 0
		for cell in type1:
			if iterator in type1TrainingCellIdxs:
				type1TrainingCells.append(cell)
				type1TrainingCellIdxsAnn.append(type1Idxs[iterator])
			else:
				testingData.append(cell)
				type1TestingCellIdxsAnn.append(type1Idxs[iterator])

			iterator += 1

		iterator = 0
		for cell in type2:
			if iterator in type2TrainingCellIdxs:
				type2TrainingCells.append(cell)
				type2TrainingCellIdxsAnn.append(type2Idxs[iterator])
			else:
				testingData.append(cell)
				type2TestingCellIdxsAnn.append(type2Idxs[iterator])

			iterator += 1

		iterator = 0
		for cell in type3:
			if iterator in type3TrainingCellIdxs:
				type3TrainingCells.append(cell)
				type3TrainingCellIdxsAnn.append(type3Idxs[iterator])
			else:
				testingData.append(cell)
				type3TestingCellIdxsAnn.append(type3Idxs[iterator])

			iterator += 1

		iterator = 0
		for cell in type4:
			if iterator in type4TrainingCellIdxs:
				type4TrainingCells.append(cell)
				type4TrainingCellIdxsAnn.append(type4Idxs[iterator])
			else:
				testingData.append(cell)
				type4TestingCellIdxsAnn.append(type4Idxs[iterator])

			iterator += 1

		iterator = 0
		for cell in type5:
			if iterator in type5TrainingCellIdxs:
				type5TrainingCells.append(cell)
				type5TrainingCellIdxsAnn.append(type5Idxs[iterator])
			else:
				testingData.append(cell)
				type5TestingCellIdxsAnn.append(type5Idxs[iterator])

			iterator += 1

		iterator = 0
		for cell in type6:
			if iterator in type6TrainingCellIdxs:
				type6TrainingCells.append(cell)
				type6TrainingCellIdxsAnn.append(type6Idxs[iterator])
			else:
				testingData.append(cell)
				type6TestingCellIdxsAnn.append(type6Idxs[iterator])

			iterator += 1

		iterator = 0
		for cell in type7:
			if iterator in type7TrainingCellIdxs:
				type7TrainingCells.append(cell)
				type7TrainingCellIdxsAnn.append(type7Idxs[iterator])
			else:
				testingData.append(cell)
				type7TestingCellIdxsAnn.append(type7Idxs[iterator])

			iterator += 1

		iterator = 0
		for cell in type8:
			if iterator in type8TrainingCellIdxs:
				type8TrainingCells.append(cell)
				type8TrainingCellIdxsAnn.append(type8Idxs[iterator])
			else:
				testingData.append(cell)
				type8TestingCellIdxsAnn.append(type8Idxs[iterator])

			iterator += 1

		iterator = 0
		for cell in type9:
			if iterator in type9TrainingCellIdxs:
				type9TrainingCells.append(cell)
				type9TrainingCellIdxsAnn.append(type9Idxs[iterator])
			else:
				testingData.append(cell)
				type9TestingCellIdxsAnn.append(type9Idxs[iterator])

			iterator += 1

		# add the training cells and target values to class wide variables in parallel 1-9
		self.trainingData = []
		targetValuesIdxs = []

		if len(type1TrainingCells) != len(type1TrainingCellIdxsAnn):
			print "error: discrepancy between type 1 training cells and type 1 training cell indices"
		else:
			iterator = 0
			for cell in type1TrainingCells:
				self.trainingData.append(cell)
				targetValuesIdxs.append(type1TrainingCellIdxsAnn[iterator])
				iterator += 1

		if len(type2TrainingCells) != len(type2TrainingCellIdxsAnn):
			print "error: discrepancy between type 2 training cells and type 2 training cell indices"
		else:
			iterator = 0
			for cell in type2TrainingCells:
				self.trainingData.append(cell)
				targetValuesIdxs.append(type2TrainingCellIdxsAnn[iterator])
				iterator += 1

		if len(type3TrainingCells) != len(type3TrainingCellIdxsAnn):
			print "error: discrepancy between type 3 training cells and type 3 training cell indices"
		else:
			iterator = 0
			for cell in type3TrainingCells:
				self.trainingData.append(cell)
				targetValuesIdxs.append(type3TrainingCellIdxsAnn[iterator])
				iterator += 1

		if len(type4TrainingCells) != len(type4TrainingCellIdxsAnn):
			print "error: discrepancy between type 4 training cells and type 4 training cell indices"
		else:
			iterator = 0
			for cell in type4TrainingCells:
				self.trainingData.append(cell)
				targetValuesIdxs.append(type4TrainingCellIdxsAnn[iterator])
				iterator += 1

		if len(type5TrainingCells) != len(type5TrainingCellIdxsAnn):
			print "error: discrepancy between type 9 training cells and type 9 training cell indices"
		else:
			iterator = 0
			for cell in type5TrainingCells:
				self.trainingData.append(cell)
				targetValuesIdxs.append(type5TrainingCellIdxsAnn[iterator])
				iterator += 1

		if len(type6TrainingCells) != len(type6TrainingCellIdxsAnn):
			print "error: discrepancy between type 6 training cells and type 6 training cell indices"
		else:
			iterator = 0
			for cell in type6TrainingCells:
				self.trainingData.append(cell)
				targetValuesIdxs.append(type6TrainingCellIdxsAnn[iterator])
				iterator += 1

		if len(type7TrainingCells) != len(type7TrainingCellIdxsAnn):
			print "error: discrepancy between type 7 training cells and type 7 training cell indices"
		else:
			iterator = 0
			for cell in type7TrainingCells:
				self.trainingData.append(cell)
				targetValuesIdxs.append(type7TrainingCellIdxsAnn[iterator])
				iterator += 1

		if len(type8TrainingCells) != len(type8TrainingCellIdxsAnn):
			print "error: discrepancy between type 8 training cells and type 8 training cell indices"
		else:
			iterator = 0
			for cell in type8TrainingCells:
				self.trainingData.append(cell)
				targetValuesIdxs.append(type8TrainingCellIdxsAnn[iterator])
				iterator += 1

		if len(type9TrainingCells) != len(type9TrainingCellIdxsAnn):
			print "error: discrepancy between type 9 training cells and type 9 training cell indices"
		else:
			iterator = 0
			for cell in type9TrainingCells:
				self.trainingData.append(cell)
				targetValuesIdxs.append(type9TrainingCellIdxsAnn[iterator])
				iterator += 1


		# now that we have the indexes of the target values that correspond to the annotations, make a list of these identifiers
		self.targetValues = []
		for idx in targetValuesIdxs:
			self.targetValues.append(int(self.cellIdentifierAnnotations[idx]))

		# set the class wide variable testing data
		self.testingData = testingData

		# add the indices corresponding to identifier annotations for the testing data in parallel 1-9
		testingDataIdxsAnn = []

		for idx in type1TestingCellIdxsAnn:
			testingDataIdxsAnn.append(idx)

		for idx in type2TestingCellIdxsAnn:
			testingDataIdxsAnn.append(idx)

		for idx in type3TestingCellIdxsAnn:
			testingDataIdxsAnn.append(idx)

		for idx in type4TestingCellIdxsAnn:
			testingDataIdxsAnn.append(idx)

		for idx in type5TestingCellIdxsAnn:
			testingDataIdxsAnn.append(idx)

		for idx in type6TestingCellIdxsAnn:
			testingDataIdxsAnn.append(idx)

		for idx in type7TestingCellIdxsAnn:
			testingDataIdxsAnn.append(idx)

		for idx in type8TestingCellIdxsAnn:
			testingDataIdxsAnn.append(idx)

		for idx in type9TestingCellIdxsAnn:
			testingDataIdxsAnn.append(idx)

		# now that we have the indexes of the testing data target values that correspond to the annotations, make a list of these identifiers
		self.testingDataTargetValues = []
		for idx in testingDataIdxsAnn:
			self.testingDataTargetValues.append(int(self.cellIdentifierAnnotations[idx]))

		numTrainingDataCells = len(self.trainingData)
		numTestingDataCells = len(self.testingData)

		print "number training cells = {numTrainingDataCells}".format(numTrainingDataCells=numTrainingDataCells)
		print "number testing cells = {numTestingDataCells}".format(numTestingDataCells=numTestingDataCells)
		print "reference:"
		print "- total down sampled cells = {dsCells}".format(dsCells=len(self.dsCluster_MoleculeData))
		print "- down sampled cells * .7 = {dsCellsTraining} --> approx.".format(dsCellsTraining=int(len(self.dsCluster_MoleculeData)*.7))
		print "- down sampled cells * .3 = {dsCellsTesting} --> approx.".format(dsCellsTesting=int(len(self.dsCluster_MoleculeData)*.3))

	def getRandIndices(self):
		return self.randIndices

	def getRawData(self):
		return self.rawData

	def getDSClusterData(self):
		return self.dsClusterData

	def getDSCluster_MoleculeData(self):
		return self.dsCluster_MoleculeData

	def getTrainingData(self):
		return self.trainingData

	def getTestingData(self):
		return self.testingData

	def getFeatures(self):
		return self.features

	def getTargetValues(self):
		return self.targetValues

	def getTestingDataTargetValues(self):
		return self.testingDataTargetValues

	def getCellIdentifierAnnotations(self):
		return self.cellIdentifierAnnotations

	def getMoleculeCountAnnotations(self):
		return self.moleculeCountAnnotations

	def getNumCellsRaw(self):
		return len(self.rawData)

	def getNumCellsDSCluster(self):
		return len(self.dsClusterData)

	def getNumGenesRaw(self):
		return len(self.rawData[0])

	def getNumGenesDSCluster(self):
		return len(self.dsClusterData[0])

	def partitionData(self):
		print "partitioning data"
		# self.trainingData
		# self.testData