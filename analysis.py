from __future__ import division
import sys
import os
import time
import math
import numpy as np


# File: analysis.py
#	After classification, this file is used to evaluate the performance of the classifier and write the results to an output file.
#	First, this file creates a confusion matrix and then computes the accuracy, sensitivity, specificity, MCC, and F1 Score at both
#	the class level and the global level across the supplied number of folds (10 folds for cross validation and 1 fold for 
#	non-cross validation). This class also uses a basic metric of merely counting the number of correct classifications that was used
#	initially to check the performance of the classifiers. After evaluating, the results are written to a file.

# global string variables for writing to file
tab = "	"
newLine = "\n"
types = ["Interneuron", "S1 Pyramidal", "CA1 Pyramidal", "Oligodendrocyte", "Microglia",
"Endothelial", "Astrocyte", "Ependymal", "Mural"]
GVs = "Global Evaluations"
headerBar = " " + tab + "Accuracy" + tab + "Sensitivity" + tab + "Specificity" + tab + "MCC" + tab + "F1"

def analyzeAndWriteToFile(classifier, predictions, answerKey, foldsEvaluations, k, runOption):
	topLevelDirName = "RNASeq_SingleCellClassification_Results"
	
	# directory names for each of the 4 combinations of down sampling and cross validation
	DSAndCVDirName = "DS & CV"
	DSAndNotCVDirName = "DS & !CV"
	NotDSAndCVDirName = "!DS & CV"
	NotDSAndNotCVDirName = "!DS & !CV"

	# date and time
	date = time.strftime("%m-%d-%Y")
	t = time.strftime("%I-%M-%S")

	# make file name
	clf = ""
	if classifier == "Radial Basis Function Support Vector Machine":
		clf = "rbfsvc"
	elif classifier == "Multi-Layer Perceptron (Neural Network)":
		clf = "nn"
	elif classifier.startswith("KNearestNeighbor Classifier_"):
		idx = classifier.index("_")
		clf = "{k}knn".format(k=classifier[idx+1:])
	elif classifier == "Random Forest Classifier":
		clf = "rf"


	fileName = clf + "_" + date + "_" + t + ".txt"

	print "\nAnalyzing results and writing to file {fileName}".format(fileName=fileName)
	
	# check if there is a 'results' directory in the project root
	if not os.path.exists(topLevelDirName):
		os.makedirs(topLevelDirName) # create the results directory heirarchy if not
		os.makedirs(topLevelDirName + "/" + DSAndCVDirName)
		os.makedirs(topLevelDirName + "/" + DSAndNotCVDirName)
		os.makedirs(topLevelDirName + "/" + NotDSAndCVDirName)
		os.makedirs(topLevelDirName + "/" + NotDSAndNotCVDirName)


	# create the file in the newly created directory based on the run options selected
	if runOption == 0:
		if not os.path.exists(topLevelDirName + "/" + DSAndCVDirName):
			os.makedirs(topLevelDirName + "/" + DSAndCVDirName)

		file = open(topLevelDirName + "/" + DSAndCVDirName + "/" + fileName, "w+")

	elif runOption == 1:
		if not os.path.exists(topLevelDirName + "/" + DSAndNotCVDirName):
			os.makedirs(topLevelDirName + "/" + DSAndNotCVDirName)

		file = open(topLevelDirName + "/" + DSAndNotCVDirName + "/" + fileName, "w+")

	elif runOption == 2:
		if not os.path.exists(topLevelDirName + "/" + NotDSAndCVDirName):
			os.makedirs(topLevelDirName + "/" + NotDSAndCVDirName)

		file = open(topLevelDirName + "/" + NotDSAndCVDirName + "/" + fileName, "w+")

	elif runOption == 3:
		if not os.path.exists(topLevelDirName + "/" + NotDSAndNotCVDirName):
			os.makedirs(topLevelDirName + "/" + NotDSAndNotCVDirName)

		file = open(topLevelDirName + "/" + NotDSAndNotCVDirName + "/" + fileName, "w+")


	print "** note: basic analysis only on final fold of cross validation"
	accuracies = analyzeResultsBasic(classifier, predictions, answerKey)
	
	sensitivities, specificities, mccs, f1Scores, globalVals = analyzeResultsRobust(foldsEvaluations, k)

	# make sure the results lists are the correct size
	if len(accuracies) != 9 and len(sensitivies) != 9 and len(specificities) != 9 and len(mccs) != 9 and len(f1Scores) != 9 and len(globalVals) != 5:
		print "error: results list are incorrect size - cannot write to file"
		file.close()
		return

	# write the header line to the file
	file.write(headerBar)

	# write the typed evaluations to the file
	typeIterator = 0
	while typeIterator < 9:
		file.write(formatLine(types[typeIterator], accuracies[typeIterator], sensitivities[typeIterator],
			specificities[typeIterator], mccs[typeIterator], f1Scores[typeIterator]))
		typeIterator += 1

	# write the global evaluations to the file
	file.write(formatLine(GVs, globalVals[0], globalVals[1], globalVals[2], globalVals[3], globalVals[4]))

	file.close()
	return


def formatLine(type, accuracy, sensitivity, specificity, mcc, f1score):
	line = newLine + type + tab + str(accuracy) + tab + str(sensitivity) + tab + str(specificity) + tab + str(mcc) + tab + str(f1score)

	return line

# this will calculate and store a 2D list with evaluations for each type in the order [accuracy, sensitivity, specificity, MCC, F1]
def calculateEvaluations(predictions, answerKey):
	# make sure lengths are equal
	if len(predictions) != len(answerKey):
		print "error: discrepancy between number of prediction results and answey keys"


	# confusionMatrix = [truePositives, falsePositives, falseNegatives, trueNegatives]

	foldEvaluations = []

	iterator = 1
	while iterator < 10: # confusion matrix for each type
		confusionMatrix = calculateConfusionMatrix(predictions, answerKey, iterator)
		foldEvaluation = calculateFoldEvaluations(confusionMatrix)
		foldEvaluations.append(foldEvaluation)
		iterator += 1

	return foldEvaluations

def calculateFoldEvaluations(confusionMatrix):

	# Reference:
	# confusionMatrix = [truePositives, falsePositives, falseNegatives, trueNegatives]
	# tp = 0
	# fp = 1
	# fn = 2
	# tn = 3

	# calculateAccurary(tp, tn, fp, fn)
	accuracy = calculateAccurary(confusionMatrix[0], confusionMatrix[3], confusionMatrix[1], confusionMatrix[2])

	# calculateSensitivity(tp, fn)
	sensitivity = calculateSensitivity(confusionMatrix[0], confusionMatrix[2])

	# calculateSpecificity(tn, fp)
	specificity = calculateSpecificity(confusionMatrix[3], confusionMatrix[1])

	# calculateMCC(tp, tn, fp, fn)
	mcc = calculateMCC(confusionMatrix[0], confusionMatrix[3], confusionMatrix[1], confusionMatrix[2])

	# calculateF1Score(tp, fp, fn):
	f1Score = calculateF1Score(confusionMatrix[0], confusionMatrix[1], confusionMatrix[2])

	foldEvaluation = [accuracy, sensitivity, specificity, mcc, f1Score]

	return foldEvaluation

def analyzeResultsRobust(foldsEvaluations, k):
	print "\nanalyzing evaluations for {k}-fold cross validation".format(k=k)

	if len(foldsEvaluations) != k and len(foldsEvaluations[0]) != 9 and len(foldsEvaluations[0][0]) != 5:
		print "error: 3D evaluations list is incorrect size"
		return

	# initialize values for avg accuracy, sensitivity, specificity, MCC, F1Score for each class and for all classes
	type1AccuracyTotal = 0
	type2AccuracyTotal = 0
	type3AccuracyTotal = 0
	type4AccuracyTotal = 0
	type5AccuracyTotal = 0
	type6AccuracyTotal = 0
	type7AccuracyTotal = 0
	type8AccuracyTotal = 0
	type9AccuracyTotal = 0

	type1SensitivityTotal = 0
	type2SensitivityTotal = 0
	type3SensitivityTotal = 0
	type4SensitivityTotal = 0
	type5SensitivityTotal = 0
	type6SensitivityTotal = 0
	type7SensitivityTotal = 0
	type8SensitivityTotal = 0
	type9SensitivityTotal = 0

	type1SpecificityTotal = 0
	type2SpecificityTotal = 0
	type3SpecificityTotal = 0
	type4SpecificityTotal = 0
	type5SpecificityTotal = 0
	type6SpecificityTotal = 0
	type7SpecificityTotal = 0
	type8SpecificityTotal = 0
	type9SpecificityTotal = 0

	type1MCCTotal = 0
	type2MCCTotal = 0
	type3MCCTotal = 0
	type4MCCTotal = 0
	type5MCCTotal = 0
	type6MCCTotal = 0
	type7MCCTotal = 0
	type8MCCTotal = 0
	type9MCCTotal = 0

	type1F1ScoreTotal = 0
	type2F1ScoreTotal = 0
	type3F1ScoreTotal = 0
	type4F1ScoreTotal = 0
	type5F1ScoreTotal = 0
	type6F1ScoreTotal = 0
	type7F1ScoreTotal = 0
	type8F1ScoreTotal = 0
	type9F1ScoreTotal = 0

	iterator = 0
	while iterator < k:
		#type 1 values
		type1AccuracyTotal += foldsEvaluations[iterator][0][0]
		type1SensitivityTotal += foldsEvaluations[iterator][0][1]
		type1SpecificityTotal += foldsEvaluations[iterator][0][2]
		type1MCCTotal += foldsEvaluations[iterator][0][3]
		type1F1ScoreTotal += foldsEvaluations[iterator][0][4]

		#type 2 values
		type2AccuracyTotal += foldsEvaluations[iterator][1][0]
		type2SensitivityTotal += foldsEvaluations[iterator][1][1]
		type2SpecificityTotal += foldsEvaluations[iterator][1][2]
		type2MCCTotal += foldsEvaluations[iterator][1][3]
		type2F1ScoreTotal += foldsEvaluations[iterator][1][4]

		#type 3 values
		type3AccuracyTotal += foldsEvaluations[iterator][2][0]
		type3SensitivityTotal += foldsEvaluations[iterator][2][1]
		type3SpecificityTotal += foldsEvaluations[iterator][2][2]
		type3MCCTotal += foldsEvaluations[iterator][2][3]
		type3F1ScoreTotal += foldsEvaluations[iterator][2][4]

		#type 4 values
		type4AccuracyTotal += foldsEvaluations[iterator][3][0]
		type4SensitivityTotal += foldsEvaluations[iterator][3][1]
		type4SpecificityTotal += foldsEvaluations[iterator][3][2]
		type4MCCTotal += foldsEvaluations[iterator][3][3]
		type4F1ScoreTotal += foldsEvaluations[iterator][3][4]

		#type 5 values
		type5AccuracyTotal += foldsEvaluations[iterator][4][0]
		type5SensitivityTotal += foldsEvaluations[iterator][4][1]
		type5SpecificityTotal += foldsEvaluations[iterator][4][2]
		type5MCCTotal += foldsEvaluations[iterator][4][3]
		type5F1ScoreTotal += foldsEvaluations[iterator][4][4]

		#type 6 values
		type6AccuracyTotal += foldsEvaluations[iterator][5][0]
		type6SensitivityTotal += foldsEvaluations[iterator][5][1]
		type6SpecificityTotal += foldsEvaluations[iterator][5][2]
		type6MCCTotal += foldsEvaluations[iterator][5][3]
		type6F1ScoreTotal += foldsEvaluations[iterator][5][4]

		#type 7 values
		type7AccuracyTotal += foldsEvaluations[iterator][6][0]
		type7SensitivityTotal += foldsEvaluations[iterator][6][1]
		type7SpecificityTotal += foldsEvaluations[iterator][6][2]
		type7MCCTotal += foldsEvaluations[iterator][6][3]
		type7F1ScoreTotal += foldsEvaluations[iterator][6][4]

		#type 8 values
		type8AccuracyTotal += foldsEvaluations[iterator][7][0]
		type8SensitivityTotal += foldsEvaluations[iterator][7][1]
		type8SpecificityTotal += foldsEvaluations[iterator][7][2]
		type8MCCTotal += foldsEvaluations[iterator][7][3]
		type8F1ScoreTotal += foldsEvaluations[iterator][7][4]

		#type 9 values
		type9AccuracyTotal += foldsEvaluations[iterator][8][0]
		type9SensitivityTotal += foldsEvaluations[iterator][8][1]
		type9SpecificityTotal += foldsEvaluations[iterator][8][2]
		type9MCCTotal += foldsEvaluations[iterator][8][3]
		type9F1ScoreTotal += foldsEvaluations[iterator][8][4]

		iterator += 1

	
	# find averages for each type
	avgAccuracyType1 = type1AccuracyTotal/k
	avgSensitivityType1 = type1SensitivityTotal/k
	avgSpecificityType1 = type1SpecificityTotal/k
	avgMccType1 = type1MCCTotal/k
	avgF1ScoreType1 = type1F1ScoreTotal/k

	avgAccuracyType2 = type2AccuracyTotal/k
	avgSensitivityType2 = type2SensitivityTotal/k
	avgSpecificityType2 = type2SpecificityTotal/k
	avgMccType2 = type2MCCTotal/k
	avgF1ScoreType2 = type2F1ScoreTotal/k

	avgAccuracyType3 = type3AccuracyTotal/k
	avgSensitivityType3 = type3SensitivityTotal/k
	avgSpecificityType3 = type3SpecificityTotal/k
	avgMccType3 = type3MCCTotal/k
	avgF1ScoreType3 = type3F1ScoreTotal/k

	avgAccuracyType4 = type4AccuracyTotal/k
	avgSensitivityType4 = type4SensitivityTotal/k
	avgSpecificityType4 = type4SpecificityTotal/k
	avgMccType4 = type4MCCTotal/k
	avgF1ScoreType4 = type4F1ScoreTotal/k

	avgAccuracyType5 = type5AccuracyTotal/k
	avgSensitivityType5 = type5SensitivityTotal/k
	avgSpecificityType5 = type5SpecificityTotal/k
	avgMccType5 = type5MCCTotal/k
	avgF1ScoreType5 = type5F1ScoreTotal/k

	avgAccuracyType6 = type6AccuracyTotal/k
	avgSensitivityType6 = type6SensitivityTotal/k
	avgSpecificityType6 = type6SpecificityTotal/k
	avgMccType6 = type6MCCTotal/k
	avgF1ScoreType6 = type6F1ScoreTotal/k

	avgAccuracyType7 = type7AccuracyTotal/k
	avgSensitivityType7 = type7SensitivityTotal/k
	avgSpecificityType7 = type7SpecificityTotal/k
	avgMccType7 = type7MCCTotal/k
	avgF1ScoreType7 = type7F1ScoreTotal/k

	avgAccuracyType8 = type8AccuracyTotal/k
	avgSensitivityType8 = type8SensitivityTotal/k
	avgSpecificityType8 = type8SpecificityTotal/k
	avgMccType8 = type8MCCTotal/k
	avgF1ScoreType8 = type8F1ScoreTotal/k

	avgAccuracyType9 = type9AccuracyTotal/k
	avgSensitivityType9 = type9SensitivityTotal/k
	avgSpecificityType9 = type9SpecificityTotal/k
	avgMccType9 = type9MCCTotal/k
	avgF1ScoreType9 = type9F1ScoreTotal/k

	# find global values
	accuracyTotal = avgAccuracyType1 + avgAccuracyType2 + avgAccuracyType3 + avgAccuracyType4 + avgAccuracyType5 + avgAccuracyType6 + avgAccuracyType7 + avgAccuracyType8 + avgAccuracyType9
	sensitivityTotal = avgSensitivityType1 + avgSensitivityType2 + avgSensitivityType3 + avgSensitivityType4 + avgSensitivityType5 + avgSensitivityType6 + avgSensitivityType7 + avgSensitivityType8 + avgSensitivityType9
	specificityTotal = avgSpecificityType1 + avgSpecificityType2 + avgSpecificityType3 + avgSpecificityType4 + avgSpecificityType5 + avgSpecificityType6 + avgSpecificityType7 + avgSpecificityType8 + avgSpecificityType9
	mccTotal = avgMccType1 + avgMccType2 + avgMccType3 + avgMccType4 + avgMccType5 + avgMccType6 + avgMccType7 + avgMccType8 + avgMccType9
	f1ScoreTotal = avgF1ScoreType1 + avgF1ScoreType2 + avgF1ScoreType3 + avgF1ScoreType4 + avgF1ScoreType5 + avgF1ScoreType6 + avgF1ScoreType7 + avgF1ScoreType8 + avgF1ScoreType9

	# find global values
	avgAccuracy = accuracyTotal/9
	avgSensitivity = sensitivityTotal/9
	avgSpecificity = specificityTotal/9
	avgMcc = mccTotal/9
	avgF1Score = f1ScoreTotal/9

	# print "Average accuracy per type: (BUG)"
	# print " - Type 1: {avgAccuracyType1}".format(avgAccuracyType1=avgAccuracyType1)
	# print " - Type 2: {avgAccuracyType2}".format(avgAccuracyType2=avgAccuracyType2)
	# print " - Type 3: {avgAccuracyType3}".format(avgAccuracyType3=avgAccuracyType3)
	# print " - Type 4: {avgAccuracyType4}".format(avgAccuracyType4=avgAccuracyType4)
	# print " - Type 5: {avgAccuracyType5}".format(avgAccuracyType5=avgAccuracyType5)
	# print " - Type 6: {avgAccuracyType6}".format(avgAccuracyType6=avgAccuracyType6)
	# print " - Type 7: {avgAccuracyType7}".format(avgAccuracyType7=avgAccuracyType7)
	# print " - Type 8: {avgAccuracyType8}".format(avgAccuracyType8=avgAccuracyType8)
	# print " - Type 9: {avgAccuracyType9}".format(avgAccuracyType9=avgAccuracyType9)

	print "\nAverage sensitivity per type: "
	print " - Type 1: {avgSensitivityType1}".format(avgSensitivityType1=avgSensitivityType1)
	print " - Type 2: {avgSensitivityType2}".format(avgSensitivityType2=avgSensitivityType2)
	print " - Type 3: {avgSensitivityType3}".format(avgSensitivityType3=avgSensitivityType3)
	print " - Type 4: {avgSensitivityType4}".format(avgSensitivityType4=avgSensitivityType4)
	print " - Type 5: {avgSensitivityType5}".format(avgSensitivityType5=avgSensitivityType5)
	print " - Type 6: {avgSensitivityType6}".format(avgSensitivityType6=avgSensitivityType6)
	print " - Type 7: {avgSensitivityType7}".format(avgSensitivityType7=avgSensitivityType7)
	print " - Type 8: {avgSensitivityType8}".format(avgSensitivityType8=avgSensitivityType8)
	print " - Type 9: {avgSensitivityType9}".format(avgSensitivityType9=avgSensitivityType9)

	print "\nAverage specificity per type: "
	print " - Type 1: {avgSpecificityType1}".format(avgSpecificityType1=avgSpecificityType1)
	print " - Type 2: {avgSpecificityType2}".format(avgSpecificityType2=avgSpecificityType2)
	print " - Type 3: {avgSpecificityType3}".format(avgSpecificityType3=avgSpecificityType3)
	print " - Type 4: {avgSpecificityType4}".format(avgSpecificityType4=avgSpecificityType4)
	print " - Type 5: {avgSpecificityType5}".format(avgSpecificityType5=avgSpecificityType5)
	print " - Type 6: {avgSpecificityType6}".format(avgSpecificityType6=avgSpecificityType6)
	print " - Type 7: {avgSpecificityType7}".format(avgSpecificityType7=avgSpecificityType7)
	print " - Type 8: {avgSpecificityType8}".format(avgSpecificityType8=avgSpecificityType8)
	print " - Type 9: {avgSpecificityType9}".format(avgSpecificityType9=avgSpecificityType9)

	print "\nAverage MCC per type: "
	print " - Type 1: {avgMccType1}".format(avgMccType1=avgMccType1)
	print " - Type 2: {avgMccType2}".format(avgMccType2=avgMccType2)
	print " - Type 3: {avgMccType3}".format(avgMccType3=avgMccType3)
	print " - Type 4: {avgMccType4}".format(avgMccType4=avgMccType4)
	print " - Type 5: {avgMccType5}".format(avgMccType5=avgMccType5)
	print " - Type 6: {avgMccType6}".format(avgMccType6=avgMccType6)
	print " - Type 7: {avgMccType7}".format(avgMccType7=avgMccType7)
	print " - Type 8: {avgMccType8}".format(avgMccType8=avgMccType8)
	print " - Type 9: {avgMccType9}".format(avgMccType9=avgMccType9)

	print "\nAverage F1 Score per type: "
	print " - Type 1: {avgF1ScoreType1}".format(avgF1ScoreType1=avgF1ScoreType1)
	print " - Type 2: {avgF1ScoreType2}".format(avgF1ScoreType2=avgF1ScoreType2)
	print " - Type 3: {avgF1ScoreType3}".format(avgF1ScoreType3=avgF1ScoreType3)
	print " - Type 4: {avgF1ScoreType4}".format(avgF1ScoreType4=avgF1ScoreType4)
	print " - Type 5: {avgF1ScoreType5}".format(avgF1ScoreType5=avgF1ScoreType5)
	print " - Type 6: {avgF1ScoreType6}".format(avgF1ScoreType6=avgF1ScoreType6)
	print " - Type 7: {avgF1ScoreType7}".format(avgF1ScoreType7=avgF1ScoreType7)
	print " - Type 8: {avgF1ScoreType8}".format(avgF1ScoreType8=avgF1ScoreType8)
	print " - Type 9: {avgF1ScoreType9}".format(avgF1ScoreType9=avgF1ScoreType9)

	print "\n Average global values: "
	print " - Accuracy: {avgAccuracy}".format(avgAccuracy=avgAccuracy)
	print " - Sensitivity: {avgSensitivity}".format(avgSensitivity=avgSensitivity)
	print " - Specificity: {avgSpecificity}".format(avgSpecificity=avgSpecificity)
	print " - MCC: {avgMcc}".format(avgMcc=avgMcc)
	print " - F1 Score: {avgF1Score}".format(avgF1Score=avgF1Score)

	# make lists of each category to return to calling method (for write to file)
	sensitivities = [avgSensitivityType1, avgSensitivityType2, avgSensitivityType3, avgSensitivityType4, avgSensitivityType5, 
	avgSensitivityType6, avgSensitivityType7, avgSensitivityType8, avgSensitivityType9]

	specificities = [avgSpecificityType1, avgSpecificityType2, avgSpecificityType3, avgSpecificityType4, avgSpecificityType5, 
	avgSpecificityType6, avgSpecificityType7, avgSpecificityType8, avgSpecificityType9]

	mccs = [avgMccType1, avgMccType2, avgMccType3, avgMccType4, avgMccType5, avgMccType6, avgMccType7, avgMccType8, avgMccType9]

	f1Scores = [avgF1ScoreType1, avgF1ScoreType2, avgF1ScoreType3, avgF1ScoreType4, avgF1ScoreType5, avgF1ScoreType6, 
	avgF1ScoreType7, avgF1ScoreType8, avgF1ScoreType9]

	globalVals = [avgAccuracy, avgSensitivity, avgSpecificity, avgMcc, avgF1Score]

	return sensitivities, specificities, mccs, f1Scores, globalVals


def analyzeResultsBasic(classifier, predictions, answerKey):
	print "\nanalyzing results for {classifier} classifier using basic metric".format(classifier=classifier)

	# make sure lengths are equal
	if len(predictions) != len(answerKey):
		print "error: discrepancy between number of prediction results and answey keys"
		return

	# initialize lists for each identifier for prediction results and answer key
	type1ResultsCount = 0
	type1KeyCounts = 0

	type2ResultsCount = 0
	type2KeyCounts = 0

	type3ResultsCount = 0
	type3KeyCounts = 0

	type4ResultsCount = 0
	type4KeyCounts = 0

	type5ResultsCount = 0
	type5KeyCounts = 0

	type6ResultsCount = 0
	type6KeyCounts = 0

	type7ResultsCount = 0
	type7KeyCounts = 0

	type8ResultsCount = 0
	type8KeyCounts = 0

	type9ResultsCount = 0
	type9KeyCounts = 0

	# find the counts for each type in the predicion results
	for _type in predictions:
		if _type == 1:
			type1ResultsCount += 1
		elif _type == 2:
			type2ResultsCount += 1
		elif _type == 3:
			type3ResultsCount += 1
		elif _type == 4:
			type4ResultsCount += 1
		elif _type == 5:
			type5ResultsCount += 1
		elif _type == 6:
			type6ResultsCount += 1
		elif _type == 7:
			type7ResultsCount += 1
		elif _type == 8:
			type8ResultsCount += 1
		elif _type == 9:
			type9ResultsCount += 1

	# find the counts for each type in the answer key
	for _type in answerKey:
		if _type == 1:
			type1KeyCounts += 1
		elif _type == 2:
			type2KeyCounts += 1
		elif _type == 3:
			type3KeyCounts += 1
		elif _type == 4:
			type4KeyCounts += 1
		elif _type == 5:
			type5KeyCounts += 1
		elif _type == 6:
			type6KeyCounts += 1
		elif _type == 7:
			type7KeyCounts += 1
		elif _type == 8:
			type8KeyCounts += 1
		elif _type == 9:
			type9KeyCounts += 1

	# we need to establish a metric to penalize over classification i.e. classifying 9 cells as type 1 when only 8 cells are type 1
	# we'll find the proportion of 1 cell/all cells of type and subtract this with each over classified cell
	type1Proportion = calculateProportion(1, type1KeyCounts)
	type2Proportion = calculateProportion(1, type2KeyCounts)
	type3Proportion = calculateProportion(1, type3KeyCounts)
	type4Proportion = calculateProportion(1, type4KeyCounts)
	type5Proportion = calculateProportion(1, type5KeyCounts)
	type6Proportion = calculateProportion(1, type6KeyCounts)
	type7Proportion = calculateProportion(1, type7KeyCounts)
	type8Proportion = calculateProportion(1, type8KeyCounts)
	type9Proportion = calculateProportion(1, type9KeyCounts)
	
	# type9Proportion = 1/type9KeyCounts

	# check each result for over classification, and substract accordingly. else, set by proportion

	type1Results = -1.
	if type1ResultsCount > type1KeyCounts:
		type1Results = 1.0 - (type1Proportion*(type1ResultsCount - type1KeyCounts))
	else:
		type1Results = calculateProportion(type1ResultsCount, type1KeyCounts)

	type2Results = -1.
	if type2ResultsCount > type2KeyCounts:
		type2Results = 1.0 - (type2Proportion*(type2ResultsCount - type2KeyCounts))
	else:
		type2Results = calculateProportion(type2ResultsCount, type2KeyCounts)

	type3Results = -1.
	if type3ResultsCount > type3KeyCounts:
		type3Results = 1.0 - (type3Proportion*(type3ResultsCount - type3KeyCounts))
	else:
		type3Results = calculateProportion(type3ResultsCount, type3KeyCounts)

	type4Results = -1.
	if type4ResultsCount > type4KeyCounts:
		type4Results = 1.0 - (type4Proportion*(type4ResultsCount - type4KeyCounts))
	else:
		type4Results = calculateProportion(type4ResultsCount, type4KeyCounts)

	type5Results = -1.
	if type5ResultsCount > type5KeyCounts:
		type5Results = 1.0 - (type5Proportion*(type5ResultsCount - type5KeyCounts))
	else:
		type5Results = calculateProportion(type5ResultsCount, type5KeyCounts)

	type6Results = -1.
	if type6ResultsCount > type6KeyCounts:
		type6Results = 1.0 - (type6Proportion*(type6ResultsCount - type6KeyCounts))
	else:
		type6Results = calculateProportion(type6ResultsCount, type6KeyCounts)

	type7Results = -1.
	if type7ResultsCount > type7KeyCounts:
		type7Results = 1.0 - (type7Proportion*(type7ResultsCount - type7KeyCounts))
	else:
		type7Results = calculateProportion(type7ResultsCount, type7KeyCounts)

	type8Results = -1.
	if type8ResultsCount > type8KeyCounts:
		type8Results = 1.0 - (type8Proportion*(type8ResultsCount - type8KeyCounts))
	else:
		type8Results = calculateProportion(type8ResultsCount, type8KeyCounts)

	type9Results = -1.
	if type9ResultsCount > type9KeyCounts:
		type9Results = 1.0 - (type9Proportion*(type9ResultsCount - type9KeyCounts))
	else:
		type9Results = calculateProportion(type9ResultsCount, type9KeyCounts)


	# compute total prediction results
	totalPredictionResults = calculateTotalAccuracy(type1Results, type2Results, type3Results, type4Results, type5Results,
		type6Results, type7Results, type8Results, type9Results)

	print "Accuracy determined by:"
	print "if number of predictions > total known cells of type:"
	print "    accuracy = 1.0 - ((1/number_cells_type) * (number of predictions - known number of cells)"
	print "else:"
	print "    accuracy = Number_of_Type_Predicted / Actual_Number_of_Type"
	print "Type 1 accuracy (predictions: {type1ResultsCount}, known: {type1KeyCounts}) = {type1Results}".format(type1ResultsCount=type1ResultsCount, 
		type1KeyCounts=type1KeyCounts, type1Results=type1Results)
	print "Type 2 accuracy (predictions: {type2ResultsCount}, known: {type2KeyCounts}) = {type2Results}".format(type2ResultsCount=type2ResultsCount, 
		type2KeyCounts=type2KeyCounts, type2Results=type2Results)
	print "Type 3 accuracy (predictions: {type3ResultsCount}, known: {type3KeyCounts}) = {type3Results}".format(type3ResultsCount=type3ResultsCount, 
		type3KeyCounts=type3KeyCounts, type3Results=type3Results)
	print "Type 4 accuracy (predictions: {type4ResultsCount}, known: {type4KeyCounts}) = {type4Results}".format(type4ResultsCount=type4ResultsCount, 
		type4KeyCounts=type4KeyCounts, type4Results=type4Results)
	print "Type 5 accuracy (predictions: {type5ResultsCount}, known: {type5KeyCounts}) = {type5Results}".format(type5ResultsCount=type5ResultsCount, 
		type5KeyCounts=type5KeyCounts, type5Results=type5Results)
	print "Type 6 accuracy (predictions: {type6ResultsCount}, known: {type6KeyCounts}) = {type6Results}".format(type6ResultsCount=type6ResultsCount, 
		type6KeyCounts=type6KeyCounts, type6Results=type6Results)
	print "Type 7 accuracy (predictions: {type7ResultsCount}, known: {type7KeyCounts}) = {type7Results}".format(type7ResultsCount=type7ResultsCount, 
		type7KeyCounts=type7KeyCounts, type7Results=type7Results)
	print "Type 8 accuracy (predictions: {type8ResultsCount}, known: {type8KeyCounts}) = {type8Results}".format(type8ResultsCount=type8ResultsCount, 
		type8KeyCounts=type8KeyCounts, type8Results=type8Results)
	print "Type 9 accuracy (predictions: {type9ResultsCount}, known: {type9KeyCounts}) = {type9Results}".format(type9ResultsCount=type9ResultsCount, 
		type9KeyCounts=type9KeyCounts, type9Results=type9Results)
	print "Total accuracy (non_zero_prediction_results/num_non_zero_values) = {totalPredictionResults}".format(totalPredictionResults=totalPredictionResults)

	# make list of basic accuracies
	accuracies = [type1Results, type2Results, type3Results, type4Results, type5Results, type6Results, type7Results, type8Results, type9Results]

	return accuracies

def calculateProportion(n, total):
	if total < 1:
		return -1

	return n/total

def calculateTotalAccuracy(type1Results, type2Results, type3Results, type4Results, type5Results, type6Results, type7Results
	, type8Results, type9Results):

	total = 0

	posValues = 0

	if type1Results != -1:
		total += type1Results
		posValues += 1
	if type2Results != -1:
		total += type2Results
		posValues += 1
	if type3Results != -1:
		total += type3Results
		posValues += 1
	if type4Results != -1:
		total += type4Results
		posValues += 1
	if type5Results != -1:
		total += type5Results
		posValues += 1
	if type6Results != -1:
		total += type6Results
		posValues += 1
	if type7Results != -1:
		total += type7Results
		posValues += 1
	if type8Results != -1:
		total += type8Results
		posValues += 1
	if type9Results != -1:
		total += type9Results
		posValues += 1

	return total/posValues


def calculateConfusionMatrix(predictions, answerKey, _type):
	# compute the confusion matrix (true positive, false negative, false positive, true negative) for the given type
	if len(predictions) != len(answerKey):
		print "error: discrepancy between number of predictions and answer key"

	size = len(predictions)

	truePositives = 0
	falsePositives = 0
	falseNegatives = 0
	trueNegatives = 0

	iterator = 0
	while iterator < size:
		prediction = predictions[iterator]
		key = answerKey[iterator]

		if key == _type:
			# if prediction matches, true positive / else, false negative
			if prediction == key:
				truePositives += 1
			else:
				falseNegatives += 1

		else: # if key is "other" type
			# if prediction matches, true negative ("other" classified as same "other") / else, false positive ("other" classified as type)
			if prediction == key:
				trueNegatives += 1
			else:
				falsePositives += 1

		iterator += 1

	# print "tp, fp, fn, tn = {tp}, {fp}, {fn}, {tn}".format(tp=truePositives, fp=falsePositives, fn=falseNegatives, tn=trueNegatives)

	confusionMatrix = [truePositives, falsePositives, falseNegatives, trueNegatives]

	return confusionMatrix

def calculateAccurary(tp, tn, fp, fn):
	numerator = (tp + tn)
	denominator = (tp + fp + fn + tn)
	
	if denominator > 0:
		return numerator / denominator

	return 0.

def calculateSensitivity(tp, fn):
	numerator = tp
	denominator = tp + fn

	if denominator > 0:
		return numerator / denominator

	return 0.

def calculateSpecificity(tn, fp):
	numerator = tn
	denominator = tn + fp

	if denominator > 0:
		return numerator / denominator

	return 0.

def calculateMCC(tp, tn, fp, fn):
	numerator = ((tp * tn) - (fp * fn))
	denominator = math.sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))

	if denominator > 0:
		return numerator / denominator

	return 0.

def calculateF1Score(tp, fp, fn):
	numerator = (2 * tp)
	denominator = ((2 * tp) + fp + fn)
	
	if denominator > 0:
		return numerator / denominator

	return 0.
