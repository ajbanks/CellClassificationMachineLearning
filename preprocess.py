import sys
import xlrd
from sklearn import preprocessing
import numpy as np

def run():
	return loadDataSetFromFile()

def loadDataSetFromFile():
	print "reading data into memory"

	# relative path to file
	fileName = "../Data/GSE60361C13005Expression.txt"

	# initialize list
	data = []

	with open(fileName) as ins:
		cellIdx = 0
		geneIdx = 1 # initialize to 1 because the first entry in the list will hold the cell name
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
				geneExpressions.pop(0)

				# add the expressions to their respective cells
				k = 0 # used to iterate across cells
				for expression in geneExpressions:
					data[k].append(expression)
					k += 1

				geneIdx += 1 # move to next gene i.e. row



	return preprocessData(data)

def preprocessData(data):
	print "preprocessing data"

	# create a numpy array from the list
	dataNP = np.array(data)

	data_normalized = preprocessing.scale(dataNP)

	return data_normalized