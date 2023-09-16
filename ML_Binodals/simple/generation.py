import tensorflow as tf
import numpy      as np
import os

from tensorflow              import keras
from keras                   import layers
from sklearn.model_selection import train_test_split

print(f"tensorflow version = {tf.__version__}")

def create_mat(x1, x2):
	M = np.zeros((2,2))
	noise = np.random.normal(0, 0.005)
	M[0][0] = np.sin(2*x1 )+noise
	M[0][1] = np.cos(x1+x2)+noise
	M[1][0] = np.sin(x1+x2)+noise
	M[1][1] = np.cos(2*x2 )+noise

	return M

def create_vec(x1, x2):
	M = np.zeros(4)
	noise = np.random.normal(0, 0.005)
	M[0] = np.sin(2*x1 )+noise
	M[1] = np.cos(x1+x2)+noise
	M[2] = np.sin(x1+x2)+noise
	M[3] = np.cos(2*x2 )+noise

	return M



if __name__=="__main__":

	# generated all the inputs and outputs
	inputs  = []
	outputs = []
	for x1 in np.linspace(0.1, 1, 20):
		for x2 in np.linspace(0.1, 1, 20):
			M = create_vec(x1, x2)
			inputs.append(np.array([x1,x2]))
			outputs.append(M)


	inputs  = np.array(inputs)
	outputs = np.array(outputs)


	train_input, test_input, train_output, test_output = train_test_split(inputs, outputs, test_size=0.2, random_state=42)

	input_layer  = keras.Input(shape=(2,))
	x            = layers.Dense(4, activation="relu")(input_layer)
	output_layer = layers.Dense(4, activation="linear")(x)


	# setting up the model in totality
	model = keras.Model(inputs=input_layer, outputs=output_layer)

	# compiling the model
	model.compile (optimizer="adam", loss="mse", metrics=[keras.metrics.MeanSquaredError()])

	# allow it to fit
	history = model.fit (train_input, train_output, epochs=512, validation_split=0.2)
	print(history.history)

	results = model.evaluate(test_input, test_output)
	print (results)