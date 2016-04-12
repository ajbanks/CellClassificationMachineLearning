from __future__ import division
import sys
import math
import numpy as np

def analyzeFoldResults(predictions, answerKey):
	# make sure lengths are equal
	if len(predictions) != len(answerKey):
		print "error: discrepancy between number of prediction results and answey keys"

	type1ConfusionMatrix = calculateConfusionMatrix(predictions, answerKey, 1)

	print type1ConfusionMatrix

def analyzeCrossValidationAccuracyResults(accuracyResults, k=10):
	print "calculating average accuracy for {k}-fold cross validation".format(k=k)
	avg = np.mean(accuracyResults)
	print "Prediction results of 10-fold cross validation: {avg}".format(avg=avg)

def analyzeResults(classifier, predictions, answerKey):
	print "\nanalyzing results for {classifier} classifier".format(classifier=classifier)

	# make sure lengths are equal
	if len(predictions) != len(answerKey):
		print "error: discrepancy between number of prediction results and answey keys"

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

	confusionMatrix = [truePositives, falsePositives, falseNegatives, trueNegatives]

	return confusionMatrix

def calculateAccurary(tp, tn, fp, fn):
	print "calculating accuracy = (#_true_positives + #_true_negatives) / (#_true_positives + #_false_positives + #_false_negatives + #_true_negatives)"

	accuracy = (tp + tn) / (tp + fp + fn + tn)

	return accuracy

def calculateSensitivity(tp, fn):
	print "calculating sensitivity = #_true_positives / (#_true_positives + #_false_negatives)"

	sensitivty = tp / (tp + fn)

	return sensitivity

def calculateSpecificity(tn, fp):
	print "calculating specificity = #_true_negatives / (#_true_negatives + #_false_positives)"

	specificity = tn / (tn + fp)

	return specificity

def calculateMCC(tp, tn, fp, fn):
	print "calculating MCC (Matthew's Correlation Coefficient) = ((#_true_positives * #_true_negatives) - (#_false_positives * #_false_negatives)) / sqroot((#_true_positives + #_false_positives)*(#_true_positives + #_false_negatives)*(#_true_negatives + #_false_positives)*(#_true_negatives + #_false_negatives))"

	MCC = ((tp * tn) - (fp * fn)) / math.sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))

	return MCC

def calculateF1Score(tp, fp, fn):
	print "calculating F1 Score = (2 * #_true_positives) / ((2 * #_true_positives) + #_false_positives + #_false_negatives)"

	F1Score = (2 * tp) / ((2 * tp) + fp + fn)

	return F1Score













# ---------------------------------- BASIC METRIC ACCURACY FOR CROSS VALIDATION -------------------------------
# def analyzeFoldResults(predictions, answerKey):
# 	# make sure lengths are equal
# 	if len(predictions) != len(answerKey):
# 		print "error: discrepancy between number of prediction results and answey keys"

# 	# initialize lists for each identifier for prediction results and answer key
# 	type1ResultsCount = 0
# 	type1KeyCounts = 0

# 	type2ResultsCount = 0
# 	type2KeyCounts = 0

# 	type3ResultsCount = 0
# 	type3KeyCounts = 0

# 	type4ResultsCount = 0
# 	type4KeyCounts = 0

# 	type5ResultsCount = 0
# 	type5KeyCounts = 0

# 	type6ResultsCount = 0
# 	type6KeyCounts = 0

# 	type7ResultsCount = 0
# 	type7KeyCounts = 0

# 	type8ResultsCount = 0
# 	type8KeyCounts = 0

# 	type9ResultsCount = 0
# 	type9KeyCounts = 0

# 	# find the counts for each type in the predicion results
# 	for _type in predictions:
# 		if _type == 1:
# 			type1ResultsCount += 1
# 		elif _type == 2:
# 			type2ResultsCount += 1
# 		elif _type == 3:
# 			type3ResultsCount += 1
# 		elif _type == 4:
# 			type4ResultsCount += 1
# 		elif _type == 5:
# 			type5ResultsCount += 1
# 		elif _type == 6:
# 			type6ResultsCount += 1
# 		elif _type == 7:
# 			type7ResultsCount += 1
# 		elif _type == 8:
# 			type8ResultsCount += 1
# 		elif _type == 9:
# 			type9ResultsCount += 1

# 	# find the counts for each type in the answer key
# 	for _type in answerKey:
# 		if _type == 1:
# 			type1KeyCounts += 1
# 		elif _type == 2:
# 			type2KeyCounts += 1
# 		elif _type == 3:
# 			type3KeyCounts += 1
# 		elif _type == 4:
# 			type4KeyCounts += 1
# 		elif _type == 5:
# 			type5KeyCounts += 1
# 		elif _type == 6:
# 			type6KeyCounts += 1
# 		elif _type == 7:
# 			type7KeyCounts += 1
# 		elif _type == 8:
# 			type8KeyCounts += 1
# 		elif _type == 9:
# 			type9KeyCounts += 1

# 	# we need to establish a metric to penalize over classification i.e. classifying 9 cells as type 1 when only 8 cells are type 1
# 	# we'll find the proportion of 1 cell/all cells of type and subtract this with each over classified cell
# 	type1Proportion = calculateProportion(1, type1KeyCounts)
# 	type2Proportion = calculateProportion(1, type2KeyCounts)
# 	type3Proportion = calculateProportion(1, type3KeyCounts)
# 	type4Proportion = calculateProportion(1, type4KeyCounts)
# 	type5Proportion = calculateProportion(1, type5KeyCounts)
# 	type6Proportion = calculateProportion(1, type6KeyCounts)
# 	type7Proportion = calculateProportion(1, type7KeyCounts)
# 	type8Proportion = calculateProportion(1, type8KeyCounts)
# 	type9Proportion = calculateProportion(1, type9KeyCounts)

# 	# check each result for over classification, and substract accordingly. else, set by proportion
# 	type1Results = -1.
# 	if type1ResultsCount > type1KeyCounts:
# 		type1Results = 1.0 - (type1Proportion*(type1ResultsCount - type1KeyCounts))
# 	else:
# 		type1Results = calculateProportion(type1ResultsCount, type1KeyCounts)

# 	type2Results = -1.
# 	if type2ResultsCount > type2KeyCounts:
# 		type2Results = 1.0 - (type2Proportion*(type2ResultsCount - type2KeyCounts))
# 	else:
# 		type2Results = calculateProportion(type2ResultsCount, type2KeyCounts)

# 	type3Results = -1.
# 	if type3ResultsCount > type3KeyCounts:
# 		type3Results = 1.0 - (type3Proportion*(type3ResultsCount - type3KeyCounts))
# 	else:
# 		type3Results = calculateProportion(type3ResultsCount, type3KeyCounts)

# 	type4Results = -1.
# 	if type4ResultsCount > type4KeyCounts:
# 		type4Results = 1.0 - (type4Proportion*(type4ResultsCount - type4KeyCounts))
# 	else:
# 		type4Results = calculateProportion(type4ResultsCount, type4KeyCounts)

# 	type5Results = -1.
# 	if type5ResultsCount > type5KeyCounts:
# 		type5Results = 1.0 - (type5Proportion*(type5ResultsCount - type5KeyCounts))
# 	else:
# 		type5Results = calculateProportion(type5ResultsCount, type5KeyCounts)

# 	type6Results = -1.
# 	if type6ResultsCount > type6KeyCounts:
# 		type6Results = 1.0 - (type6Proportion*(type6ResultsCount - type6KeyCounts))
# 	else:
# 		type6Results = calculateProportion(type6ResultsCount, type6KeyCounts)

# 	type7Results = -1.
# 	if type7ResultsCount > type7KeyCounts:
# 		type7Results = 1.0 - (type7Proportion*(type7ResultsCount - type7KeyCounts))
# 	else:
# 		type7Results = calculateProportion(type7ResultsCount, type7KeyCounts)

# 	type8Results = -1.
# 	if type8ResultsCount > type8KeyCounts:
# 		type8Results = 1.0 - (type8Proportion*(type8ResultsCount - type8KeyCounts))
# 	else:
# 		type8Results = calculateProportion(type8ResultsCount, type8KeyCounts)

# 	type9Results = -1.
# 	if type9ResultsCount > type9KeyCounts:
# 		type9Results = 1.0 - (type9Proportion*(type9ResultsCount - type9KeyCounts))
# 	else:
# 		type9Results = calculateProportion(type9ResultsCount, type9KeyCounts)


# 	# compute total prediction results
# 	totalPredictionResults = calculateTotalAccuracy(type1Results, type2Results, type3Results, type4Results, type5Results,
# 		type6Results, type7Results, type8Results, type9Results)

# 	return totalPredictionResults