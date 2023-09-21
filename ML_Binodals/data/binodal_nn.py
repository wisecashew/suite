import tensorflow        as tf
import numpy             as np
import pandas            as pd
import matplotlib.pyplot as plt
import mpltern
import os
import re

from tensorflow              import keras
from keras                   import layers
from tensorflow.keras.losses import Loss
from sklearn.model_selection import train_test_split

print(f"tensorflow version = {tf.__version__}")

if __name__=="__main__":

	# get all the trimmed binodals 
	current_directory = os.getcwd()
	all_files         = os.listdir(current_directory)
	trimmed_files     = [file for file in all_files if (file.startswith("trim") and file.endswith(".binodal"))]

	# go through all the trimmed files and get them as dataframes
	the_output_binodals = []
	for tfile in trimmed_files:
		df = pd.read_csv(tfile, sep='\|', names=["phi_s_top","phi_p_top","phi_c_top","phi_s_bot","phi_p_bot","phi_c_bot"], skiprows=1, engine='python')
		data = df.values.reshape(40000, 6)

		# convery to numpy object
		data = np.array(data, dtype=np.float64)
		the_output_binodals.append(data)

	the_outputs = np.array(the_output_binodals)

	# go through all the trimmed files and get the inputs as vectors
	# this will be a vector that looks like [vs, vp, vc, chi_sc, chi_pc, chi_pc]

	# define a regex that captures numerical (positive and negative) decimal values
	pattern = r'(?:vs|vc|vp|chisc|chips|chipc)_([-+]?\d+\.\d+)'
	the_inputs = []
	for tstring in trimmed_files:
	
		# find all matches
		matches = re.findall(pattern, tstring)
		the_inputs.append(np.array(matches, dtype=np.float64))
	
	the_inputs = np.array(the_inputs)

	# I now have my inputs and outputs well-defined
	# make sure they are the same size (i.e. for every input, there should be a corresponding output)

	if (the_inputs.shape[0] == the_outputs.shape[0]):
		print(f"Number of inputs = {the_inputs.shape[0]}, number of outputs = {the_outputs.shape[0]}")
		print("Inputs and outputs have one-to-one correspondence! This is good. Moving on to training the model.\n")

	else:
		print(f"Number of inputs = {the_inputs.shape[0]}, number of outputs = {the_outputs.shape[0]}")
		print("There is a problem.Inputs and outputs do not have a one-to-one correspondence.")
		exit ()

	# get the total size of the dataset
	total_samples = the_inputs.shape[0]

	# define a split
	train_inputs, test_inputs, train_outputs, test_outputs = train_test_split(the_inputs, the_outputs, test_size=0.2, random_state=42)

	# start creating the structure of the neural network
	inputs  = keras.Input(shape=(6,))
	x       = layers.Dense(64, activation="relu")(inputs)
	x       = layers.Dense(128, activation="relu")(x)
	x       = layers.Dense(256, activation="relu")(x)
	x       = layers.Dense(40000*6, activation="linear")(x)
	outputs = layers.Reshape((40000, 6))(x)

	# now that we have the skeleton of our neural network we need to define a custom loss function
	# we will try to faithfully recreate the ones in the training set, but also 
	# tack on the loss that comes with getting the chemical potentials wrong
	mu_a = lambda phi_a, phi_b, phi_c, vs, vc, vp, chi_sc, chi_ps, chi_pc: np.log(phi_a) + 1 - phi_a \
				- vs/vp * phi_b - vs/vc * (phi_c) + vs * (phi_b**2 * chi_ps + (phi_c)**2 * \
				chi_sc + phi_b * (phi_c) * (chi_ps + chi_sc - chi_pc) ) 

	mu_b = lambda phi_a, phi_b, phi_c, vs, vc, vp, chi_sc, chi_ps, chi_pc: np.log(phi_b) + 1 - phi_b \
				- vp/vs * phi_a - vp/vc * (phi_c) + vp * (phi_a**2 * chi_ps + (phi_c)**2 * \
				chi_pc + phi_a * (phi_c) * (chi_ps + chi_pc - chi_sc) )

	mu_c = lambda phi_a, phi_b, phi_c, vs, vc, vp, chi_sc, chi_ps, chi_pc: np.log(phi_c) + 1 - phi_c \
				- vc/vs * phi_a - vc/vp * phi_b + vc * (phi_a**2 * chi_sc + phi_b**2 * \
				chi_pc + phi_a * phi_b * (chi_sc + chi_pc - chi_ps) )

	delta_mu_a = lambda pa1, pb1, pc1, pa2, pb2, pc2, vs, vc, vp, chi_sc, chi_ps, chi_pc: \
				np.abs(mu_a(pa1, pb1, pc1, vs, vc, vp, chi_sc, chi_ps, chi_pc) - \
				mu_a(pa2, pb2, pc2, vs, vc, vp, chi_sc, chi_ps, chi_pc))

	delta_mu_b = lambda pa1, pb1, pc1, pa2, pb2, pc2, vs, vc, vp, chi_sc, chi_ps, chi_pc: \
				np.abs(mu_b(pa1, pb1, pc1, vs, vc, vp, chi_sc, chi_ps, chi_pc) - \
				mu_b(pa2, pb2, pc2, vs, vc, vp, chi_sc, chi_ps, chi_pc))

	delta_mu_c = lambda pa1, pb1, pc1, pa2, pb2, pc2, vs, vc, vp, chi_sc, chi_ps, chi_pc: \
				np.abs(mu_c(pa1, pb1, pc1, vs, vc, vp, chi_sc, chi_ps, chi_pc) - \
				mu_c(pa2, pb2, pc2, vs, vc, vp, chi_sc, chi_ps, chi_pc))

	# set up the inputs and outputs
	model = keras.Model(inputs=inputs, outputs=outputs)

	# compile the thing before we can fit things to it
	model.compile (optimizer='adam', loss=keras.losses.MeanSquaredError(), metrics=[keras.metrics.MeanSquaredError()])

	# allow it to fit
	history = model.fit(train_inputs, train_outputs, epochs=500, validation_split=0.1)

	# check how well the model does on 
	model.save ("binodal_maker.h5")


