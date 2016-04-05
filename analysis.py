from __future__ import division
import sys

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

	# compute prediction results for each type
	type1Results = type1ResultsCount/type1KeyCounts
	type2Results = type2ResultsCount/type2KeyCounts
	type3Results = type3ResultsCount/type3KeyCounts
	type4Results = type4ResultsCount/type4KeyCounts
	type5Results = type5ResultsCount/type5KeyCounts
	type6Results = type6ResultsCount/type6KeyCounts
	type7Results = type7ResultsCount/type7KeyCounts
	type8Results = type8ResultsCount/type8KeyCounts
	type9Results = type9ResultsCount/type9KeyCounts

	# compute total prediction results
	allPredictionResults = type1Results + type2Results + type3Results + type4Results + type5Results + type6Results + type7Results
	+ type8Results + type9Results

	totalPredictionResults = allPredictionResults/9

	print "Accuracy determing by = Number_of_Type_Predicted / Actual_Number_of_Type"
	print "Type 1 accuracy - {type1ResultsCount}/{type1KeyCounts} = {type1Results}".format(type1ResultsCount=type1ResultsCount, 
		type1KeyCounts=type1KeyCounts, type1Results=type1Results)
	print "Type 2 accuracy - {type2ResultsCount}/{type2KeyCounts} = {type2Results}".format(type2ResultsCount=type2ResultsCount, 
		type2KeyCounts=type2KeyCounts, type2Results=type2Results)
	print "Type 3 accuracy - {type3ResultsCount}/{type3KeyCounts} = {type3Results}".format(type3ResultsCount=type3ResultsCount, 
		type3KeyCounts=type3KeyCounts, type3Results=type3Results)
	print "Type 4 accuracy - {type4ResultsCount}/{type4KeyCounts} = {type4Results}".format(type4ResultsCount=type4ResultsCount, 
		type4KeyCounts=type4KeyCounts, type4Results=type4Results)
	print "Type 5 accuracy - {type5ResultsCount}/{type5KeyCounts} = {type5Results}".format(type5ResultsCount=type5ResultsCount, 
		type5KeyCounts=type5KeyCounts, type5Results=type5Results)
	print "Type 6 accuracy - {type6ResultsCount}/{type6KeyCounts} = {type6Results}".format(type6ResultsCount=type6ResultsCount, 
		type6KeyCounts=type6KeyCounts, type6Results=type6Results)
	print "Type 7 accuracy - {type7ResultsCount}/{type7KeyCounts} = {type7Results}".format(type7ResultsCount=type7ResultsCount, 
		type7KeyCounts=type7KeyCounts, type7Results=type7Results)
	print "Type 8 accuracy - {type8ResultsCount}/{type8KeyCounts} = {type8Results}".format(type8ResultsCount=type8ResultsCount, 
		type8KeyCounts=type8KeyCounts, type8Results=type8Results)
	print "Type 9 accuracy - {type9ResultsCount}/{type9KeyCounts} = {type9Results}".format(type9ResultsCount=type9ResultsCount, 
		type9KeyCounts=type9KeyCounts, type9Results=type9Results)
	print "Total accuracy (all prediction results/9) = {totalPredictionResults}".format(totalPredictionResults=totalPredictionResults)




