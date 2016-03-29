import sys

class RNASeqData(object):
	rawDataFileName = "GSE60361C13005Expression.txt"

	# annotationsFileName = "";

	def __init__(self):
		print "initializing RNASeqData"

	def setData(self, data):
		self.data = data

	def getData(self):
		for row in data:
			print row

	def partitionData(self):
		print "partitioning data"
		# self.trainingData
		# self.testData