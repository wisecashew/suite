import tensorflow               as tf
import tensorflow.keras.backend as K
import numpy                    as np

from tensorflow              import keras
from keras                   import layers
from sklearn.model_selection import train_test_split
from keras.layers import Lambda

import time
import os
import argparse


def create_vec(x1, x2):
    M = np.zeros(4)
    noise = np.random.normal(0, 0.005, 4)
    M[0] = np.sin(2*x1 ) + noise[0]
    M[1] = np.cos(x1+x2) + noise[1]
    M[2] = np.sin(x1+x2) + noise[2]
    M[3] = np.cos(2*x2 ) + noise[3]

    return M

if __name__=="__main__":

    start = time.time()
    print(f"tensorflow version = {tf.__version__}")
    # generated all the inputs and outputs
    inputs  = []
    outputs = []
    for x1 in np.linspace(0.1, 1, 20):
        for x2 in np.linspace(0.1, 1, 20):
            M = create_vec(x1, x2)
            inputs.append(np.array([x1,x2]))
            outputs.append(M)

    inputs  = np.array(inputs )
    outputs = np.array(outputs)

    # inputs[:,0] = (inputs[:,0] - np.mean(inputs[:,0]))/np.std(inputs[:,0])
    # inputs[:,1] = (inputs[:,1] - np.mean(inputs[:,1]))/np.std(inputs[:,1])

    train_input, test_input, train_output, test_output = train_test_split(inputs, outputs, test_size=0.2, random_state=42)

    input_layer  = keras.Input(shape=(2,))
    x            = layers.Dense(64  , activation="sigmoid" )(input_layer)
    output_layer = layers.Dense(4, activation="linear", name="output_layer")(x)

    # another output, minimize its difference to 0.
    sin_2x1      = tf.sin(input_layer[:, 0])
    sum_sin_x1x2 = tf.cos(input_layer[:, 0]+input_layer[:, 1])
    sum_cos_x1x2 = tf.sin(input_layer[:, 0]+input_layer[:, 1])
    cos_2x2      = tf.cos(input_layer[:, 1])
    trig         = tf.divide(output_layer[:, 0], sin_2x1) + tf.divide(output_layer[:, 1], sum_cos_x1x2) - tf.divide(output_layer[:, 2], sum_sin_x1x2) - tf.divide(output_layer[:, 3], cos_2x2)
    trig_output  = Lambda(lambda x: tf.reduce_mean(tf.abs(x)), name="trig_output")(trig)

    # setting up the model in totality
    model = keras.Model(inputs=input_layer, outputs=[output_layer, trig_output])

    # compiling the model
    model.compile(optimizer="adam", 
                loss={"output_layer": "mse", "trig_output": "mse"},
                metrics={"output_layer": ["mae"],
                           "trig_output": ["mae"]},
                loss_weights={'output_layer': 1.0, 
                "trig_output": 0.2})

    # allow it to fit
    history = model.fit (train_input, [train_output, np.zeros((len(train_input),))], epochs=1000, validation_split=0.2)
    print(f"MSE loss: {history.history['loss'][-1]}")
    print(f"val loss: {history.history['val_loss'][-1]}")

    results = model.evaluate(test_input, test_output)
    print(model.predict([[0.5,0.5], [0.1, 0.9]]))
    print (results)
    stop = time.time()
    print (f"Time for computation is {stop-start}")
