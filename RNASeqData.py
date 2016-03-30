import sys

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

	# annotations that identify the group (1-9) of each cell
	def setCellIdentifierAnnotations(self, cellIdentifierAnnotations):
		self.cellIdentifierAnnotations = cellIdentifierAnnotations

	# annotations that list the molecule count for each cell
	def setMoleculeCountAnnotations(self, moleculeCountAnnotations):
		self.moleculeCountAnnotations = moleculeCountAnnotations

	def setRandIndicesFromDS(self, randIndices):
		self.randIndices = randIndices

	def getRandIndices(self):
		return self.randIndices

	def getRawData(self):
		return self.rawData

	def getDSClusterData(self):
		return self.dsClusterData

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