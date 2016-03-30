import sys

class RNASeqData(object):
	rawDataFileName = "GSE60361C13005Expression.txt"
	annotationsFileName = "expressionmRNAAnnotations.txt"

	# annotationsFileName = "";

	def __init__(self):
		print "\ninitializing RNASeqData"

	def setRawData(self, rawData):
		self.rawData = rawData

	def setCellIdentifierAnnotations(self, cellIdentifierAnnotations):
		self.cellIdentifierAnnotations = cellIdentifierAnnotations

	def setMoleculeCountAnnotations(self, moleculeCountAnnotations):
		self.moleculeCountAnnotations = moleculeCountAnnotations

	def getNumCells(self):
		return len(self.rawData)

	def getNumGenes(self):
		return len(self.rawData[0])

	def getRawData(self):
		return self.rawData

	def getCellIdentifierAnnotations(self):
		return self.cellIdentifierAnnotations

	def getMoleculeCountAnnotations(self):
		return self.moleculeCountAnnotations

	def partitionData(self):
		print "partitioning data"
		# self.trainingData
		# self.testData