import sys
import xlrd
from sklearn import preprocessing
import numpy as np
from sklearn.cross_validation import train_test_split

def run(fileName):
	return loadDataSetFromFile(fileName)

def loadDataSetFromFile(fileName):
	print "reading data into memory"

	# the 9 cell type classifying genes
	Thy1 = "Thy1"
	Gad1 = "Gad1"
	Tbr1 = "Tbr1"
	Spink8 = "Spink8"
	Mbp = "Mbp"
	Aldoc = "Aldoc"
	Aif1 = "Aif1"
	Cldn5 = "Cldn5"
	Acta2 = "Acta2"

	# initialize list
	data = []

	cellIdx = 0
	geneIdx = 1 # initialize to 1 because the first entry in the list will hold the cell name
	genesRemoved = 0 # used to keep track of the reduction of genes
	with open(fileName) as ins:
		for line in ins:
			if cellIdx == 0: 
				# split the line into individual cell names
				cells = line.split("	")

				# remove 'cell_id'
				cells.pop(0)

				# initialize empty list for each cell, add cell name
				for cell in cells:
					data.append([])
					# data[cellIdx].append(cell) --> we don't want the cell names in our data
					cellIdx += 1
			else:
				# split the line into expressions of this gene to all cells
				geneExpressions = line.split("	")

				# remove the gene name
				gene = geneExpressions.pop(0)

				if gene == Thy1:
					print "Reading classifier Thy1"
				
				elif gene == Gad1:
					print "Processing classifier Gad1"
				
				elif gene == Tbr1:
					print "Processing classifier Tbr1"
				
				elif gene == Spink8:
					print "Processing classifier Spink8"
				
				elif gene == Mbp:
					print "Processing classifier Mbp"
				
				elif gene == Aldoc:
					print "Processing classifier Aldoc"
				
				elif gene == Aif1:
					print "Processing classifier Aif1"
				
				elif gene == Cldn5:
					print "Processing classifier Cldn5"
				
				elif gene == Acta2:
					print "Processing classifier Acta2"

				# initialize empty list to hold gene reads
				geneReads = []

				for expression in geneExpressions:
					geneReads.append(float(expression)) # cast string to float and append

				# sum the reads for the gene
				readTotal = np.sum(geneReads)
				if readTotal <= 25:
					genesRemoved+=1
				else: # if total reads > 25, add the expressions to their respective cells
					data.append(geneReads)
					geneIdx += 1 # move to next gene

	print "Added {numGenes} gene reads to {numCells} cells".format(numGenes=geneIdx, numCells=cellIdx)
	print "Reduced data dimensionality by removing {numGenesRemoved} genes with <= 25 total reads".format(numGenesRemoved=genesRemoved)

	return data

# normalize data
def preprocessData(data):
	print "preprocessing data"

	# create a numpy array from the list
	dataNP = np.array(data)

	data_normalized = preprocessing.scale(dataNP)

	partitionData(data)

# Preprocessing notes:
# 	- Test set has same format as training set
#	- Ideal is 70/30 partition
#	- Must partition randomly ***
