from sklearn.neural_network import MLPClassifier
from tensorflow import keras
import tensorflow as tf
from hyperlib.nn.layers import lin_hyp
from hyperlib.nn.optimizers import rsgd
from hyperlib.manifold import poincare
import numpy as np
import datetime

# File: neuralNetwork_RNASeq.py
#	This file defines the Multi-Layer Perceptron (Neural Network). It fits the training data and the
#	samples to the classifier. Then, it takes training data and makes predictions, returning the 
#	results of the predictions.
hyp_model = tf.keras.models.Sequential([
	lin_hyp.LinearHyperbolic(10, poincare.Poincare(), 1, activation='relu'),
  lin_hyp.LinearHyperbolic(10, poincare.Poincare(), 1, activation='softmax')
])

log_dir = "logs/fit/" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=log_dir, histogram_freq=1)
sgd = keras.optimizers.SGD(learning_rate=0.1)
rsgd = rsgd.RSGD(learning_rate=1)

hyp_model.compile(
    optimizer=rsgd,
    loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True),
    metrics=[tf.keras.metrics.SparseCategoricalAccuracy()],
)


euclidean_model = tf.keras.models.Sequential([
	tf.keras.layers.Dense(10, activation='relu'),
  tf.keras.layers.Dense(10, activation='softmax')
])

euclidean_model.compile(
    optimizer=sgd,
    loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True),
    metrics=[tf.keras.metrics.SparseCategoricalAccuracy()],
)

mlpClf = MLPClassifier(hidden_layer_sizes=[100], max_iter=200)
# activation='relu', alpha=1e-05,
#        batch_size='auto', beta_1=0.9, beta_2=0.999, early_stopping=False,
#        epsilon=1e-08, hidden_layer_sizes=(100,), learning_rate='constant',
#        learning_rate_init=0.001, max_iter=200, momentum=0.9,
#        nesterovs_momentum=True, power_t=0.5, random_state=1, shuffle=True,
#        tol=0.0001, validation_fraction=0.1, verbose=False,
#        warm_start=False)

# Source: http://scikit-learn.org/dev/modules/neural_networks_supervised.html
def fitTrainingData(training_data, nSamples):
	# fit(X,Y)
		# x = 2D array of cells
		# y = 1D list of classifiers i.e. [group1, group1, group2, group2]

	# makes lists into np arrays
	print(np.array(training_data)[0])
	training_dataNP = np.array(training_data)
	nSamplesNP = np.array(nSamples)

	# mlpClf.fit(training_dataNP, nSamplesNP)

	print(nSamplesNP)
	print("fit hyperbolic model")
	hyp_model.fit(
		training_dataNP, nSamplesNP,
		epochs=40
	)

def predictTestData(testing_data):
	# make list into np array
	testing_dataNP = np.array(testing_data)

	# predict the values
	# predicted = mlpClf.predict(testing_dataNP)


	raw_predicted = hyp_model.predict(testing_dataNP)
	predicted = []
	for pred in raw_predicted:
		result = np.where(pred == np.amax(pred))
		predicted.append(result[0][0])

	print("predicted hyperbolic model")
	print(predicted)
	return predicted