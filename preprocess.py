from __future__ import division
import sys
import xlrd
from sklearn import preprocessing
import numpy as np
import random
import math


# File: preprocess.py
#	This file does all of the preprocessing work before running classification. This file loads
#	raw data and annotations into memory. Next, this file is used to down dample by both
#	cluster size and molecule count.

def loadRawData(fileName):
	print "\nreading data into memory"

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
	genesUsed = 0 # tracks genes with > 25 total reads (start at 1 because we'll use this as idx for appending)
	genesRemoved = 0 # used to keep track of the reduction of genes i.e. those with less than 26 total reads
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
					print "Processing classifier Thy1"
				
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
					genesRemoved += 1
				else: # if total reads > 25, add the expressions to their respective cells
					# add each read to its cell
					if len(geneReads) == cellIdx: # make sure the number of reads matches the number of cells
						idx = 0 # used to iterate through cells
						for x in geneReads:
							data[idx].append(x) # append the read to the cell's list of gene reads
							idx+=1
						genesUsed += 1 # move to next gene (next row)
					else:
						print "{geneReads} - {cells}".format(geneReads=len(geneReads), cells=cellIdx) # error check

	print "Added {numGenes} gene reads to {numCells} cells".format(numGenes=genesUsed, numCells=cellIdx)
	print "Reduced data dimensionality by removing {numGenesRemoved} genes with <= 25 total reads".format(numGenesRemoved=genesRemoved)

	print "rows (cells) = {rows}".format(rows=len(data))
	print "cols (genes) = {cols}".format(cols=len(data[0]))

	return data

def loadCellIdentifierAnnotations(fileName, numCells):
	print "\nloading cell identifier annotations"
	with open(fileName) as ins:
		# ignore first line of ____ ?
		line_1 = ins.readline()

		line_2 = ins.readline()
		
		# split the line into individual cell identifiers
		identifiers = line_2.split("	")

		# remove whitespace and title
		identifiers.pop(0)
		identifiers.pop(0)

		if len(identifiers) == numCells:
			# remove newline identifier from final entry
			identifiers[len(identifiers)-1] = identifiers[len(identifiers)-1][:-len("\n")]

			return identifiers

		else: #debug statement
			print "Error loading annotations: num identifiers = {numIdentifiers}, num cells = {numCells}".format(numIdentifiers=len(identifiers), numCells=numCells)

def loadMoleculeCountAnnotations(fileName, numCells):
	print "\nloading molecule count annotations"
	with open(fileName) as ins:
		# ignore first lines of ___ and group numbers
		line_1 = ins.readline()
		line_2 = ins.readline()

		line_3 = ins.readline()

		moleculeCounts = line_3.split("	")

		# remove whitespace and title
		moleculeCounts.pop(0)
		moleculeCounts.pop(0)
		
		if len(moleculeCounts) == numCells:
			# remove newline identifier from final entry
			moleculeCounts[len(moleculeCounts)-1] = moleculeCounts[len(moleculeCounts)-1][:-len("\n")]

			return moleculeCounts
		else: # debug statement
			print "Error loading annotations: num molecule counts = {numMoleculeCounts}, num cells = {numCells}".format(numMoleculeCounts=len(moleculeCounts), numCells=numCells)

def downSampleByClusterSize(rawData, cellIdentifierAnnotations):
	print "\ndown sampling by cluster size"

	# initialize counters for the 9 types
	type1 = 0
	type2 = 0
	type3 = 0
	type4 = 0
	type5 = 0
	type6 = 0
	type7 = 0
	type8 = 0
	type9 = 0

	# initialize lists to hold indices of cells of a group
	type1Indices = []
	type2Indices = []
	type3Indices = []
	type4Indices = []
	type5Indices = []
	type6Indices = []
	type7Indices = []
	type8Indices = []
	type9Indices = []

	idx = 0
	for t in cellIdentifierAnnotations:
		_type = int(t)

		if _type == 1:
			type1Indices.append(idx)
			type1 += 1
		elif _type == 2:
			type2Indices.append(idx)
			type2 += 1
		elif _type == 3:
			type3Indices.append(idx)
			type3 += 1
		elif _type == 4:
			type4Indices.append(idx)
			type4 += 1
		elif _type == 5:
			type5Indices.append(idx)
			type5 += 1
		elif _type == 6:
			type6Indices.append(idx)
			type6 += 1
		elif _type == 7:
			type7Indices.append(idx)
			type7 += 1
		elif _type == 8:
			type8Indices.append(idx)
			type8 += 1
		elif _type == 9:
			type9Indices.append(idx)
			type9 += 1

		# increment the idx
		idx += 1


	print "type1 count = {type1Count}".format(type1Count=type1)
	print "type2 count = {type2Count}".format(type2Count=type2)
	print "type3 count = {type3Count}".format(type3Count=type3)
	print "type4 count = {type4Count}".format(type4Count=type4)
	print "type5 count = {type5Count}".format(type5Count=type5)
	print "type6 count = {type6Count}".format(type6Count=type6)
	print "type7 count = {type7Count}".format(type7Count=type7)
	print "type8 count = {type8Count}".format(type8Count=type8)
	print "type9 count = {type9Count}".format(type9Count=type9)

	# add all type counts to list to find smallest value
	typeCounts = [type1, type2, type3, type4, type5, type6, type7, type8, type9]

	minVal = min(typeCounts)

	print "minumum cluster cell count = {minClusterCount}".format(minClusterCount=minVal)

	# take minVal number of random cells from each cluster
	randType1Indices = random.sample(type1Indices, minVal)
	randType2Indices = random.sample(type2Indices, minVal)
	randType3Indices = random.sample(type3Indices, minVal)
	randType4Indices = random.sample(type4Indices, minVal)
	randType5Indices = random.sample(type5Indices, minVal)
	randType6Indices = random.sample(type6Indices, minVal)
	randType7Indices = random.sample(type7Indices, minVal)
	randType8Indices = random.sample(type8Indices, minVal)
	randType9Indices = random.sample(type9Indices, minVal)

	# make list of all randomly selected indices
	randIndices = []

	for randIdx in randType1Indices:
		randIndices.append(randIdx)

	for randIdx in randType2Indices:
		randIndices.append(randIdx)

	for randIdx in randType3Indices:
		randIndices.append(randIdx)

	for randIdx in randType4Indices:
		randIndices.append(randIdx)

	for randIdx in randType5Indices:
		randIndices.append(randIdx)

	for randIdx in randType6Indices:
		randIndices.append(randIdx)

	for randIdx in randType7Indices:
		randIndices.append(randIdx)

	for randIdx in randType8Indices:
		randIndices.append(randIdx)

	for randIdx in randType9Indices:
		randIndices.append(randIdx)

	# use the random indices to remove cells whose index is not in any random list 1-9 from the raw data set
	# BUG - SHOULD ONLY BE 234 CELLS LEFT AFTER BUT CURRENTLY 1540
	downSampleClusterData = []
	idx = 0
	for cell in rawData:
		if idx in randIndices: # if the current cell was not randomly selected, remove
			downSampleClusterData.append(cell)
		idx += 1 # move to next cell

	print "{numCells} total cells randomly selected, genes = {numGenes}".format(numCells=len(downSampleClusterData), 
		numGenes=len(downSampleClusterData[0]))

	return downSampleClusterData, randIndices

def downSampleByMoleculeCount(downSampleClusterData, moleculeCountAnnotations, randIndices):
	print "\ndown sampling by molecule count"

	# make sure the number of random indices matches the number of cells in the down sample cluster data
	if len(downSampleClusterData) != len(randIndices):
		print "error: discrepancy between number of cells in db and indices list"
		return;

	# make a list of the molecule count annotations corresponding to the randomly selected cells
	moleculeCountAnnsDS = []
	idx = 0
	for moleculeCount in moleculeCountAnnotations:
		if idx in randIndices:
			moleculeCountAnnsDS.append(int(moleculeCount))
		idx += 1

	# make sure that we have parallel annotations for all cells
	if len(downSampleClusterData) != len(moleculeCountAnnsDS):
		print "error: discrepancy between number of cells in db and molecule annotations"
		return;

	# find the smallest value in the list i.e. the cell with the least number of molecules. other values for reference
	minValue = min(moleculeCountAnnsDS)
	maxValue = max(moleculeCountAnnsDS)
	avgValue = sum(moleculeCountAnnsDS)/len(moleculeCountAnnsDS)

	print "down sampling all cells to minimum value = {minValue} molecules (max value = {maxValue}, avg value = {avgValue})".format(minValue=minValue,
	 maxValue=maxValue, avgValue=avgValue)

	# create data set down sampled by both cluster size and mnolecule count
	dsCluster_MoleculeData = []

	# iterate over all cells in down sample cluster data
	idx = 0
	for cell in downSampleClusterData:
		# initialize new list for genes down sampled by molecule count
		genesDS = []
		
		# find the proportion of this cell's molecule count to the smallest cell's molecule count
		prop = minValue/moleculeCountAnnsDS[idx]
		
		# iterate over all genes in the cell
		for gene in cell:

			# scale this gene to the proportion i.e. geneReadCount * proportion
			_gene = gene * prop

			# randomize a number between 0 and 1 to round up ( > .5) or down ( < .5)
			rand = random.random()

			# append the gene to the list
			genesDS.append(_gene)

			# if rand > .4999999999999:
			# 	genesDS.append(math.ceil(_gene))
			# else:
			# 	genesDS.append(math.floor(_gene))

		# add the cell's scaled gene data to the db
		dsCluster_MoleculeData.append(genesDS)

		# increment index
		idx += 1

	# make sure the number of cells and number of genes are correct
	if len(downSampleClusterData) == len(dsCluster_MoleculeData) and len(downSampleClusterData[0]) == len(dsCluster_MoleculeData[0]):
	 return dsCluster_MoleculeData






