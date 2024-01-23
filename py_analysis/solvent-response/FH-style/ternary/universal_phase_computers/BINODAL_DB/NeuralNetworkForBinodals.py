import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import torch
from torch import nn
from torch.utils.data import Dataset, DataLoader
from torch.utils.data import random_split
from torchvision import datasets
from torchvision.transforms import ToTensor

import argparse
parser = argparse.ArgumentParser(description='Create a skeleton solution for the binodal. This is a memory-intensive computation.'         )
parser.add_argument('-w_fn',    metavar='w_fn',     dest='w_fn',    type=float, action='store', help='enter weight of the default loss function.', default=0.5)
parser.add_argument('-w_split', metavar='w_split',  dest='w_split', type=float, action='store', help='enter weight of my loss function.',          default=0.5)
args = parser.parse_args() 

# get cpu, gpu or mps device for training
device = (
    "cuda" if torch.cuda.is_available()
    else "mps"
    if torch.backends.mps.is_available()
    else "cpu")

print(f"Using {device} device.", flush=True)

class CustomDataset(Dataset):
	def __init__(self, dataframe, input_cols, output_cols):
		self.dataframe   = dataframe
		self.input_cols  = input_cols
		self.output_cols = output_cols

	def __len__(self):
		return len(self.dataframe)

	def __getitem__(self, idx):
		inputs    = torch.tensor(self.dataframe.to_numpy()[idx, self.input_cols] , dtype=torch.float32)
		pred_y    = torch.tensor(self.dataframe.to_numpy()[idx, self.output_cols], dtype=torch.float32)
		pred_mask = torch.sign  (torch.tensor(self.dataframe.to_numpy()[idx, self.output_cols]))
		return inputs, pred_y, pred_mask

class NeuralNetwork(nn.Module):
	def __init__(self):
		super().__init__()
		self.flatten = nn.Flatten()
		self.stack   = nn.Sequential(
			nn.Linear(8, 16),
			nn.ReLU(),
			nn.Linear(16,16),
			nn.ReLU(),
			nn.Linear(16, 6),
			nn.Sigmoid()
			)

	def forward(self, x):
		vols  = self.stack(x)
		masks = torch.sigmoid(vols)
		return vols, masks

model = NeuralNetwork().to(device)

# to train a model, you need a loss function and an optimizer
learning_rate = 1e-3
batch_size    = 256
w_fn          = args.w_fn   /(args.w_fn+args.w_split)
w_split       = args.w_split/(args.w_fn+args.w_split)

# define an optimizer
optimizer     = torch.optim.Adam(model.parameters(), lr=learning_rate)

# define losses
# the mask is a bunch of ones and zeros
# so this is a binary cross entropy loss 
mask_loss = nn.BCELoss()

def phase_split_loss(pred_fracs, pred_mask, y):
	pred = pred_fracs * pred_mask
	loss = torch.mean((pred - y)**2)
	return loss

# define the training loop 
def train(dataloader, model, mask_loss, phase_split_loss, optimizer):
	train_loss = []
	size = len(dataloader)
	model.train()
	for batch, (X, fracs, mask) in enumerate(dataloader):

		# compute prediction error
		pred_fracs, pred_mask  = model(X)
		split_loss   = phase_split_loss(pred_fracs, pred_mask, fracs)
		pred_mask    = pred_mask.float()
		mask         = mask.float()
		m_loss       = mask_loss(pred_mask, mask)
		loss         = w_split*split_loss + w_fn*m_loss

		# backprop
		loss.backward()
		optimizer.step()
		optimizer.zero_grad()
		train_loss.append(loss.item())

		if batch % 100 == 0:
			loss    = loss.item()
			current = (batch+1)*len(X)
			print(f"batch loss: {loss:>7f} [{current:>5d}/{size*len(X):>5d}]", flush=True)

	return np.mean(train_loss)

# define the test loop
def test(dataloader, model, mask_loss, phase_split_loss):
	# set the model to evaluation mode -- important for batch normalization and dropout layers
	# unnecessary in this situation but added for best practices 
	model.eval()
	size        = len(dataloader.dataset)
	num_batches = len(dataloader)
	test_loss   = 0

	# Evaluating the model with torch.no_grad() ensures that no gradients are computed during test mode
	# also serves to reduce unnecessary gradient computations and memory usage for tensors with requires_grad=True
	with torch.no_grad():
		for X, fracs, mask in dataloader:
			pred_fracs, pred_mask = model(X)
			split_loss = phase_split_loss(pred_fracs, pred_mask, fracs)
			m_loss     = mask_loss(pred_mask, mask)
			test_loss += w_split*split_loss + w_fn*m_loss 

	test_loss /= num_batches
	print(f"Test error:\n Accuracy: avg loss: {test_loss:>8f}\n", flush=True)
	return test_loss 

# being actual computation
df      = pd.read_csv("chips_2.4-chipc_2.4-chisc_2.4-vs_1-vc_1-vp_1.db", sep="\|", engine="python")
dataset = CustomDataset(df, [0, 1, 2, 3, 4, 5, 6, 7], [8, 9, 10, 11, 12, 13])
train_size = int(0.8 * len(dataset))
test_size  = len(dataset) - train_size

train_dataset, test_dataset = random_split(dataset, [train_size, test_size])

dataloader_train = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
dataloader_test  = DataLoader(test_dataset , batch_size=batch_size, shuffle=True)

print(f"Training dataset size: {len(dataloader_train)}", flush=True)
print(f"Testing dataset size:  {len(dataloader_test)} ", flush=True)

epochs = 10
epoch_train_loss = []
epoch_test_loss  = [] 

for t in range(epochs):
	# define the losses
	train_loss = train(dataloader_train, model, mask_loss, phase_split_loss, optimizer)
	test_loss  = test (dataloader_test , model, mask_loss, phase_split_loss)

	# accumulate all the epochal losses
	epoch_train_loss.append(train_loss)
	epoch_test_loss.append(test_loss)


print(f"Done with training and testing the neural network!", flush=True)
print(f"Making image...", end=' ', flush=True)
fig = plt.figure(figsize=(2,2))
ax  = plt.axes()
ax.plot(epoch_train_loss, range(epochs), c='darkred',   marker='o', mec='k', lw=1, label="training loss")
ax.plot(epoch_test_loss,  range(epochs), c='steelblue', marker='^', mec='k', lw=1, label="testing loss" )
ax.legend()
ax.grid()
fig.savefig("epochal_losses.png", dpi=1200, bbox_inches="tight")
print(f"done!", flush=True)

torch.save(model.state_dict(), "model.pth")
print(f"Saved the binodal solver in \"model.pth\".", flush=True)





