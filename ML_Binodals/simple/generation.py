import tensorflow               as tf
import tensorflow.keras.backend as K
import numpy                    as np

from tensorflow              import keras
from keras                   import layers
from sklearn.model_selection import train_test_split

import time
import os
import argparse


parser = argparse.ArgumentParser(description="Run a neural network on my test system.")
parser.add_argument("--mat", action='store_true', dest='mat', help="Create outputs in the shape of a matrix.", default=False)
parser.add_argument("--vec", action='store_true', dest='vec', help="Create outputs in the shape of a vector.", default=False)
args   = parser.parse_args()

def create_mat(x1, x2):
	M = np.zeros((2,2))
	noise = np.random.normal(0, 0.005, 4)
	M[0][0] = np.sin(2*x1 )+noise[0]
	M[0][1] = np.cos(x1+x2)+noise[1]
	M[1][0] = np.sin(x1+x2)+noise[2]
	M[1][1] = np.cos(2*x2 )+noise[3]

	return M

def create_vec(x1, x2):
	M = np.zeros(4)
	noise = np.random.normal(0, 0.005, 4)
	M[0] = np.sin(2*x1 )+noise[0]
	M[1] = np.cos(x1+x2)+noise[1]
	M[2] = np.sin(x1+x2)+noise[2]
	M[3] = np.cos(2*x2 )+noise[3]

	return M



if __name__=="__main__":
	start = time.time()

	print(f"tensorflow version = {tf.__version__}")

	# generated all the inputs and outputs
	inputs  = []
	outputs = []
	for x1 in np.linspace(0.1, 1, 20):
		for x2 in np.linspace(0.1, 1, 20):
			if args.mat:
				M = create_mat(x1, x2)
			elif args.vec:
				M = create_vec(x1, x2)
			else:
				print ("There is an issue. Bad argument provided.")
				exit ()
			inputs.append(np.array([x1,x2]))
			outputs.append(M)

	inputs  = np.array(inputs )
	outputs = np.array(outputs)

	

	train_input, test_input, train_output, test_output = train_test_split(inputs, outputs, test_size=0.2, random_state=42)
	tf_inputs  = tf.convert_to_tensor(train_input,  dtype=tf.float32)

	print (f"tf_inputs = {tf_inputs.shape}")

	input_layer  = keras.Input(shape=(2,))
	x            = layers.Dense(64, activation="relu" )(input_layer)
	x            = layers.Dense(64, activation="relu" )(x)
	output_layer = layers.Dense(4, activation="linear")(x)
	if args.mat:
		output_layer = layers.Reshape((2,2), activation="linear")(output_layer)

	def custom_loss(y_true, y_pred):
	    
		print (f"y_true.shape = {y_true.shape}")
		print (f"type(y_true) = {type(y_true)}")
		mse_loss  = tf.reduce_mean(tf.square(y_true - y_pred))

		# sin_2x1      = tf.sin(tf_inputs[:,0])
		# sum_sin_x1x2 = tf.cos(tf_inputs[:,0]+tf_inputs[:,1])
		# sum_cos_x1x2 = tf.sin(tf_inputs[:,0]+tf_inputs[:,1])
		# cos_2x2      = tf.cos(tf_inputs[:,1])

		# trig_loss = y_pred[:,0]/sin_2x1 + y_pred[:,1]/sum_cos_x1x2 - y_pred[:,2]/sum_sin_x1x2 - y_pred[:,3]/cos_2x2
		# trig_loss = tf.reduce_mean(tf.abs(trig_loss))

		return mse_loss # + trig_loss


	# setting up the model in totality
	model = keras.Model(inputs=input_layer, outputs=output_layer)

	# compiling the model
	model.compile (optimizer="adam", loss=custom_loss, metrics=[keras.metrics.MeanSquaredError()])

	# allow it to fit
	history = model.fit (train_input, train_output, epochs=512, validation_split=0.2)
	print(f"MSE loss: {history.history['loss'][-1]}")
	print(f"val loss: {history.history['val_loss'][-1]}")

	results = model.evaluate(test_input, test_output)
	print (results)
	stop = time.time()
	print (f"Time for computation is {stop-start}")
