import sys
import preprocess
import SCC_Rna_seq_NaiveBayes


if __name__ == '__main__':
	print "start"

	data = preprocess.run()

	# print len(data[0]) --> should be 19972 i.e. 19972 genes

	# for d in data[0]:
	# 	print d

	# SCC_Rna_seq_NaiveBayes.run()

	# FIGURE OUT THIS WARNING: /Users/bradenkatzman/anaconda2/lib/python2.7/site-packages/sklearn/utils/validation.py:420: DataConversionWarning: Data with input dtype |S5 was converted to float64 by the scale function.
 # warnings.warn(msg, DataConversionWarning)