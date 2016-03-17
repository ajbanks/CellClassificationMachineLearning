import sys
import numpy as np
import preprocess
import SCC_Rna_seq_NaiveBayes


if __name__ == '__main__':
	print "start"

	data = preprocess.run()

	for x in np.nditer(data):
		print x

	# print len(data[0]) --> should be 19972 i.e. 19972 genes

	# for d in data[0]:
	# 	print d

	# SCC_Rna_seq_NaiveBayes.run()