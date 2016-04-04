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


		# iterate over random indices and add cells to corresponding lists
		iterator = 0
		for idx in self.randIndices:
			# extract type from cell identifier annotations
			_type = int(self.cellIdentifierAnnotations[idx])

			# random indices and down sample data are parallel, so add the corresponding cell to the list
			if _type == 1:
				type1.append(self.dsCluster_MoleculeData[iterator])
			elif _type == 2:
				type2.append(self.dsCluster_MoleculeData[iterator])
			elif _type == 3:
				type3.append(self.dsCluster_MoleculeData[iterator])
			elif _type == 4:
				type4.append(self.dsCluster_MoleculeData[iterator])
			elif _type == 5:
				type5.append(self.dsCluster_MoleculeData[iterator])
			elif _type == 6:
				type6.append(self.dsCluster_MoleculeData[iterator])
			elif _type == 7:
				type7.append(self.dsCluster_MoleculeData[iterator])
			elif _type == 8:
				type8.append(self.dsCluster_MoleculeData[iterator])
			elif _type == 9:
				type9.append(self.dsCluster_MoleculeData[iterator])

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

		# iterate over cells in each type, add it to training or testing according to whether its index was randomly selected
		iterator = 0
		for cell in type1:
			if iterator in type1TrainingCellIdxs:
				trainingData.append(cell)
			else:
				testingData.append(cell)

			iterator += 1

		iterator = 0
		for cell in type2:
			if iterator in type2TrainingCellIdxs:
				trainingData.append(cell)
			else:
				testingData.append(cell)

			iterator += 1

		iterator = 0
		for cell in type3:
			if iterator in type3TrainingCellIdxs:
				trainingData.append(cell)
			else:
				testingData.append(cell)

			iterator += 1

		iterator = 0
		for cell in type4:
			if iterator in type4TrainingCellIdxs:
				trainingData.append(cell)
			else:
				testingData.append(cell)

			iterator += 1

		iterator = 0
		for cell in type5:
			if iterator in type5TrainingCellIdxs:
				trainingData.append(cell)
			else:
				testingData.append(cell)

			iterator += 1

		iterator = 0
		for cell in type6:
			if iterator in type6TrainingCellIdxs:
				trainingData.append(cell)
			else:
				testingData.append(cell)

			iterator += 1

		iterator = 0
		for cell in type7:
			if iterator in type7TrainingCellIdxs:
				trainingData.append(cell)
			else:
				testingData.append(cell)

			iterator += 1

		iterator = 0
		for cell in type8:
			if iterator in type8TrainingCellIdxs:
				trainingData.append(cell)
			else:
				testingData.append(cell)

			iterator += 1

		iterator = 0
		for cell in type9:
			if iterator in type9TrainingCellIdxs:
				trainingData.append(cell)
			else:
				testingData.append(cell)

			iterator += 1

		numTrainingDataCells = len(trainingData)
		numTestingDataCells = len(testingData)

		print "number training cells = {numTrainingDataCells}".format(numTrainingDataCells=numTrainingDataCells)
		print "number testing cells = {numTestingDataCells}".format(numTestingDataCells=numTestingDataCells)
		print "reference:"
		print "- total down sampled cells = {dsCells}".format(dsCells=len(self.dsCluster_MoleculeData))
		print "- down sampled cells * .7 = {dsCellsTraining} --> approx.".format(dsCellsTraining=int(len(self.dsCluster_MoleculeData)*.7))
		print "- down sampled cells * .3 = {dsCellsTesting} --> approx.".format(dsCellsTesting=int(len(self.dsCluster_MoleculeData)*.3))

		self.trainingData = trainingData
		self.testingData = testingData

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