import sys
import time
import numpy as np
from RNASeqData import RNASeqData
import preprocess
import guassianNB_RNASeq
import rbfSVC_RNASeq
import analysis


# Resource: http://machinelearningmastery.com/get-your-hands-dirty-with-scikit-learn-now/
# Python for Java Programmers: http://python4java.necaiseweb.org/Fundamentals/TheBasics

# run with down sampling and cross validation: python main.py GSE60361C13005Expression.txt expressionmRNAAnnotations.txt 1 1
# run with down sampling and without cross validation: python main.py GSE60361C13005Expression.txt expressionmRNAAnnotations.txt 1 0
# run without downsampling and with cross validation: python main.py GSE60361C13005Expression.txt expressionmRNAAnnotations.txt 0 1
# run without downsampling and without cross validation: python main.py GSE60361C13005Expression.txt expressionmRNAAnnotations.txt 0 0

if __name__ == '__main__':
	t0 = time.clock()
	print "start"
	
	# check for correct number of args
	if len(sys.argv) != 5:
		print "Usage: python main.py <raw_data_file> <annotations_file> <down sample? --> 0,1> <cross validate? --> 0,1>"
		sys.exit(0)

	raw_data_file = sys.argv[1]
	annotations_file = sys.argv[2]
	downSampleFlag = False
	crossValidateFlag = False
	if sys.argv[3] == "1":
		downSampleFlag = True
	if sys.argv[4] == "1":
		crossValidateFlag = True


	print "Using:"
	print " - raw data: {raw}".format(raw=raw_data_file)
	print " - annotations: {ann}".format(ann=annotations_file)
	if downSampleFlag:
		print "** Down sampling enabled **"
	else:
		print "** Down sampling disabled **"
	if crossValidateFlag:
		print "** Cross validation enabled **"
	else:
		print "** Cross validation disabled **"


	# initialize the data set class
	data = RNASeqData(raw_data_file, annotations_file)
	
	# read raw RNA seq data into memory
	data.setRawData(preprocess.loadRawData(data.getRawDataFileName()))

	# read cell identifier annotations into memory
	data.setCellIdentifierAnnotations(preprocess.loadCellIdentifierAnnotations(data.getAnnotationsFileName(), 
		data.getNumCellsRaw()))

	# read molecule count annotations into memory
	data.setMoleculeCountAnnotations(preprocess.loadMoleculeCountAnnotations(data.getAnnotationsFileName(),
		data.getNumCellsRaw()))

	if downSampleFlag:
		# down sample the data by cluster size --> MAKE DOWN SAMPLING A CLA OPTION
		#	 i.e. scale all cluster size to the smallest cluster (by number of cells)
		# save down sampled data and random indices for accessing corresponding annotations
		downSampleClusterData, randIndices = preprocess.downSampleByClusterSize(data.getRawData(), 
			data.getCellIdentifierAnnotations())

		# add the data and random indices reference to the data class
		data.setDSClusterData(downSampleClusterData)
		data.setRandIndicesFromDS(randIndices)

		# down sample the data by the cell with the least number of molecules
		data.setDSCluster_MoleculeData(preprocess.downSampleByMoleculeCount(data.getDSClusterData(),
			data.getMoleculeCountAnnotations(), data.getRandIndices()))

		if crossValidateFlag:
			# make 10-fold cross validation data
			data.makeCrossValidationTrainingAndTestingData(downSampleFlag)

			folds = data.getFolds()

			foldsKey = data.getFoldsKey()

			# make sure the data is parallel
			if len(folds) != len(foldsKey) or len(folds[0]) != len(foldsKey[0]):
				print "error: folds and folds key are not parallel data sets"
				sys.exit(0)

			print "\nfitting training data and predicting test data on GNBC and RBFSVC"
			iterator = 0 # we'll use this to iterate through folds and use each as the training data
			foldsEvaluationsGNB = []
			foldsEvaluationsRBFSVC = []
			while iterator < 10:
				testingData = folds[iterator]
				testingDataKey = foldsKey[iterator]

				# make 2D arrays of training cells and keys
				trainingFolds = []
				trainingKeys = []
				i = 0
				while i < 10:
					if i != iterator:
						for cell in folds[i]:
							trainingFolds.append(cell)
						for key in foldsKey[i]:
							trainingKeys.append(key)
					i += 1

				# ***************** GUASSIAN NB *****************
				# fit the 9 training folds using nb classifier
				guassianNB_RNASeq.fitTrainingData(trainingFolds, trainingKeys)

				# make predictions on the 1 testing fold using nb classifier
				guassianNB_predictionResults = guassianNB_RNASeq.predictTestData(testingData)

				# add the accuracies for this fold to accuracies list
				foldsEvaluationsGNB.append(analysis.calculateEvaluations(guassianNB_predictionResults, testingDataKey))
				# ***************** END GUASSIAN NB *****************


				# ***************** RBF SVC *****************
				# fit the 9 training folds using rbf svc
				rbfSVC_RNASeq.fitTrainingData(trainingFolds, trainingKeys)

				# make predictions on the 1 testing fold using rbf svc
				rbfSVC_predictionResults = rbfSVC_RNASeq.predictTestData(testingData)

				# add the accuracies for this fold to accuracies list
				foldsEvaluationsRBFSVC.append(analysis.calculateEvaluations(rbfSVC_predictionResults, testingDataKey))
				# ***************** END RBF SVC *****************


				# increment iterator to process the next fold as testing data
				iterator += 1

			# ***************** GUASSIAN NB *****************
			# analysis with basic metric
			print "** note: basic analysis only on final fold of cross validation"
			analysis.analyzeResultsBasic("Guassian Naive Bayes", guassianNB_predictionResults, testingDataKey)

			# analysis with robust evaluations
			analysis.analyzeResultsRobust(foldsEvaluationsGNB, 10)
			# ***************** END GUASSIAN NB *****************


			# ***************** RBF SVC *****************
			# analysis with basic metric
			analysis.analyzeAndWriteToFile("Radial Basis Function Support Vector Machine", rbfSVC_predictionResults, testingDataKey, foldsEvaluationsRBFSVC, 10)
			# analysis.analyzeResultsBasic("Radial Basis Function Support Vector Machine", rbfSVC_predictionResults, testingDataKey)

			# # analysis with robust evaluations
			# analysis.analyzeResultsRobust(foldsEvaluationsRBFSVC, 10)
			# ***************** END RBF SVC *****************

		else:
			# partition the down sampled data set into 70% training and 30% testing
			data.makeDSTrainingAndTestingData()

			# ***************** GUASSIAN NB *****************
			# fit down sampled training data to gaussian nb
			print "\nfitting training data to gaussian naive bayes classifier"
			guassianNB_RNASeq.fitTrainingData(data.getDSTrainingData(), data.getDSTargetValues())

			# predict values using gaussian nb with down sampled data
			print "predicting test data on gaussian naive bayes classifier"
			guassianNB_predictionResults = guassianNB_RNASeq.predictTestData(data.getDSTestingData())
		
			# analyze results of guassian nb with basic metric
			analysis.analyzeResultsBasic("Guassian Naive Bayes", guassianNB_predictionResults, data.getDSTestingDataTargetValues())

			# analyze results using robust evaluations
			foldsEvaluationsGNB = [] # single fold list but we still need to use a 3D list
			foldsEvaluationsGNB.append(analysis.calculateEvaluations(guassianNB_predictionResults, data.getDSTestingDataTargetValues()))
			analysis.analyzeResultsRobust(foldsEvaluationsGNB, 1)
			# ***************** END GUASSIAN NB *****************



			# ***************** RBF SVC *****************
			# fit down sampled training data to radial basis function support vector machine
			print "\n\nfitting training data to radial basis function support vector machine"
			rbfSVC_RNASeq.fitTrainingData(data.getDSTrainingData(), data.getDSTargetValues())

			# predict values using rbf support vector machine with down sampled data
			print "predicting test data on radial basis function support vector machine"
			rbfSVC_predictionResults = rbfSVC_RNASeq.predictTestData(data.getDSTestingData())

			# # analyze results of rbf support vector machine with basic metric
			# analysis.analyzeResultsBasic("Radial Basis Function Support Vector Machine", rbfSVC_predictionResults, data.getDSTestingDataTargetValues())

			# analyze results using robust evaluations
			foldsEvaluationsRBFSVC = [] # single fold list but we still need to use a 3D list
			foldsEvaluationsRBFSVC.append(analysis.calculateEvaluations(rbfSVC_predictionResults, data.getDSTestingDataTargetValues()))
			# analysis.analyzeResultsRobust(foldsEvaluationsRBFSVC, 1)

			analysis.analyzeAndWriteToFile("Radial Basis Function Support Vector Machine", rbfSVC_predictionResults, data.getDSTestingDataTargetValues(), foldsEvaluationsRBFSVC, 1)
			# ***************** END RBF SVC *****************


	else:
		if crossValidateFlag:
			# make 10-fold cross validation data
			data.makeCrossValidationTrainingAndTestingData(downSampleFlag)

			folds = data.getFolds()

			foldsKey = data.getFoldsKey()

			# make sure the data is parallel
			if len(folds) != len(foldsKey) or len(folds[0]) != len(foldsKey[0]):
				print "error: folds and folds key are not parallel data sets"
				sys.exit(0)

			print "\nfitting training data and predicting test data on GNBC and RBFSVC"
			iterator = 0 # we'll use this to iterate through folds and use each as the training data
			foldsEvaluationsGNB = []
			foldsEvaluationsRBFSVC = []
			while iterator < 10:
				testingData = folds[iterator]
				testingDataKey = foldsKey[iterator]

				# make 2D arrays of training cells and keys
				trainingFolds = []
				trainingKeys = []
				i = 0
				while i < 10:
					if i != iterator:
						for cell in folds[i]:
							trainingFolds.append(cell)
						for key in foldsKey[i]:
							trainingKeys.append(key)
					i += 1

				# ***************** GUASSIAN NB *****************
				# fit the 9 training folds using nb classifier
				guassianNB_RNASeq.fitTrainingData(trainingFolds, trainingKeys)

				# make predictions on the 1 testing fold using nb classifier
				guassianNB_predictionResults = guassianNB_RNASeq.predictTestData(testingData)

				# add the accuracies for this fold to accuracies list
				foldsEvaluationsGNB.append(analysis.calculateEvaluations(guassianNB_predictionResults, testingDataKey))
				# ***************** END GUASSIAN NB *****************


				# ***************** RBF SVC *****************
				# fit the 9 training folds using nb classifier
				rbfSVC_RNASeq.fitTrainingData(trainingFolds, trainingKeys)

				# make predictions on the 1 testing fold using nb classifier
				rbfSVC_predictionResults = rbfSVC_RNASeq.predictTestData(testingData)

				# add the accuracies for this fold to accuracies list
				foldsEvaluationsRBFSVC.append(analysis.calculateEvaluations(rbfSVC_predictionResults, testingDataKey))
				# ***************** END RBF SVC *****************

				# increment iterator to process the next fold as testing data
				iterator += 1

			# ***************** GUASSIAN NB *****************
			# analysis with basic metric
			print "** note: basic analysis only on final fold of cross validation"
			analysis.analyzeResultsBasic("Guassian Naive Bayes", guassianNB_predictionResults, testingDataKey)

			# analysis with robust evaluations
			analysis.analyzeResultsRobust(foldsEvaluations, 10)
			# ***************** END GUASSIAN NB *****************



			# ***************** RBF SVC *****************
			
			analysis.analyzeAndWriteToFile("Radial Basis Function Support Vector Machine", rbfSVC_predictionResults, testingDataKey, foldsEvaluationsRBFSVC, 10)

			# # analysis with basic metric
			# print "** note: basic analysis only on final fold of cross validation"
			# analysis.analyzeResultsBasic("Radial Basis Function Support Vector Machine", rbfSVC_predictionResults, testingDataKey)

			# # analysis with robust evaluations
			# analysis.analyzeResultsRobust(foldsEvaluationsRBFSVC, 10)
			# ***************** END RBF SVC *****************


		else:
			# partition the data set into 70% training and 30% testing
			data.makeTrainingAndTestingData()

			# ***************** GUASSIAN NB *****************
			# fit training data to gaussian nb
			print "\nfitting training data to guassian naive bayes classifier"
			guassianNB_RNASeq.fitTrainingData(data.getTrainingData(), data.getTrainingDataTargetValues())

			# predict values using gaussian nb with down sampled data
			print "predicting test data on guassian naive bayes classifier"
			guassianNB_predictionResults = guassianNB_RNASeq.predictTestData(data.getTestingData())
		
			# analyze results of guassian nb on down sampled data using basic metric
			analysis.analyzeResultsBasic("Guassian Naive Bayes", guassianNB_predictionResults, data.getTestingDataTargetValues())

			# analyze results using robust evaluations
			foldsEvaluationsGNB = [] # single fold list but we still need to use a 3D list
			foldsEvaluationsGNB.append(analysis.calculateEvaluations(guassianNB_predictionResults, data.getTestingDataTargetValues()))
			analysis.analyzeResultsRobust(foldsEvaluationsGNB, 1)
			# ***************** END GUASSIAN NB *****************



			# ***************** RBF SVC *****************
			# fit training data to gaussian nb
			print "\nfitting training data to radial basis function support vector machine"
			rbfSVC_RNASeq.fitTrainingData(data.getTrainingData(), data.getTrainingDataTargetValues())

			# predict values using gaussian nb with down sampled data
			print "predicting test data on radial basis function support vector machine"
			rbfSVC_predictionResults = rbfSVC_RNASeq.predictTestData(data.getTestingData())

			# analyze results of guassian nb on down sampled data using basic metric
			# analysis.analyzeResultsBasic("Radial Basis Function Support Vector Machine", rbfSVC_predictionResults, data.getTestingDataTargetValues())

			# analyze results using robust evaluations
			foldsEvaluationsRBFSVC = [] # single fold list but we still need to use a 3D list
			foldsEvaluationsRBFSVC.append(analysis.calculateEvaluations(rbfSVC_predictionResults, data.getTestingDataTargetValues()))
			# analysis.analyzeResultsRobust(foldsEvaluationsRBFSVC, 1)

			analysis.analyzeAndWriteToFile("Radial Basis Function Support Vector Machine", rbfSVC_predictionResults, data.getTestingDataTargetValues(), foldsEvaluationsRBFSVC, 1)

			# ***************** END RBF SVC *****************



	

	print "\nprogram execution: {t} seconds".format(t=time.clock()-t0)
	print "exiting"
