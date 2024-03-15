#!/home/satyend/.conda/envs/ML/bin/python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import mpltern
import torch
from torch import nn
from torch.utils.data import Dataset, DataLoader
from torch.utils.data import random_split
import time

import argparse
parser = argparse.ArgumentParser(description='Creates a neural network to predict how a system phase separates.')
parser.add_argument('--label',             dest='label',            type=str,   action='store', help="Enter style of neural network: \"with_splits\", \"with_weights_sep\", \"with_weights_comb\"")
parser.add_argument('--database',          dest='db',               type=str,   action='store', help='enter address of the database.')
parser.add_argument('--epochs',            dest='ep',               type=int,   action='store', help='enter number of epochs.')
parser.add_argument('--w_comp',            dest='w_comp',           type=float, action='store', help="weight for the comp MSE loss.",    default=0.5)
parser.add_argument('--w_wts',             dest='w_wts',            type=float, action='store', help="weight for the weights MSE loss.", default=0.5)
parser.add_argument('--w_split',           dest='w_split',          type=float, action='store', help="weight for the splits MSE loss.",  default=0.5)
parser.add_argument('--w_cross',           dest='w_cross',          type=float, action='store', help="weight for the entropy loss.",     default=0.5)
parser.add_argument('--w_unity',           dest='w_unity',          type=float, action='store', help="weight for the unity loss.",       default=0)
parser.add_argument('--w_mu',              dest='w_mu',             type=float, action='store', help="weight for the mu loss.",          default=0)
parser.add_argument('--prefix',            dest='prefix',           type=str,   action='store', help='enter prefix for the image.',      default="")
parser.add_argument('--dump-to',           dest='dumpfile',         type=str,   action='store', help='enter the address where to dump serialized neural network.')
parser.add_argument('--get-from',          dest='getfile',          type=str,   action='store', help='enter the address where the serialized neural network is.')
args = parser.parse_args() 

start = time.time()

def make_figures(epochs, epoch_train_losses, epoch_valid_losses):
	# create the grand loss figure
	fig_loss = plt.figure(num=1, figsize=(4,4))
	ax_loss  = plt.axes()

	# plot the loss plots 
	ax_loss.plot(range(epochs), epoch_train_losses[:,0],  c='darkred', marker='o', mec='k', lw=1, label="training loss")
	ax_loss.plot(range(epochs), epoch_valid_losses[:,0],  c='tomato',  marker='^', mec='k', lw=1, label="validation loss")
	ax_loss.legend()
	ax_loss.grid()
	fig_loss.savefig(args.prefix+"_total_loss", dpi=1200, bbox_inches="tight")

	# create the mse loss figure
	fig_mse = plt.figure(num=2, figsize=(4,4))
	ax_mse  = plt.axes() 

	# plot the mse plots
	ax_mse.plot(range(epochs), epoch_train_losses[:,1],  c='darkred', marker='o', mec='k', lw=1, label="mse training loss")
	ax_mse.plot(range(epochs), epoch_valid_losses[:,1],  c='tomato',  marker='^', mec='k', lw=1, label="mse validation loss")
	ax_mse.legend()
	ax_mse.grid()
	fig_mse.savefig(args.prefix+"_mse_loss", dpi=1200, bbox_inches="tight")

	# create the entropy loss figure
	fig_entropy = plt.figure(num=3, figsize=(4,4))
	ax_entropy  = plt.axes()

	# plot the entropy plots
	ax_entropy.plot(range(epochs), epoch_train_losses[:,2],  c='darkred', marker='o', mec='k', lw=1, label="label entropy training loss")
	ax_entropy.plot(range(epochs), epoch_valid_losses[:,2],  c='tomato',  marker='^', mec='k', lw=1, label="label entropy validation loss")
	ax_entropy.legend()
	ax_entropy.grid()
	fig_entropy.savefig(args.prefix+"_entropy_loss", dpi=1200, bbox_inches="tight")

	# create the accuracy curve
	fig_unity = plt.figure(num=4, figsize=(4,4))
	ax_unity  = plt.axes() 

	# create the accuracy plots
	ax_unity.plot(range(epochs), epoch_train_losses[:,3],  c='darkred',   marker='o', mec='k', lw=1, label="training  unity loss")
	ax_unity.plot(range(epochs), epoch_valid_losses[:,3],  c='tomato',   marker='^', mec='k', lw=1,  label="validation unity loss")
	ax_unity.legend()
	ax_unity.grid()
	fig_unity.savefig(args.prefix+"_unityloss", dpi=1200, bbox_inches="tight")

	# create the unity loss curve 
	fig_mu = plt.figure(num=5, figsize=(4,4))
	ax_mu  = plt.axes()

	ax_mu.plot(range(epochs), epoch_train_losses[:,4], c='darkred', marker='o', mec='k', lw=1,  label="training mu loss")
	ax_mu.plot(range(epochs), epoch_valid_losses [:,4], c='tomato',  marker='^', mec='k', lw=1, label="validation mu loss")
	ax_mu.legend()
	ax_mu.grid()
	fig_mu.savefig(args.prefix+"_mu_loss", dpi=1200, bbox_inches="tight")

	# create the chemical potential loss curve
	fig_acc = plt.figure(num=6, figsize=(4,4))
	ax_acc  = plt.axes()

	ax_acc.plot(range(epochs), epoch_train_losses[:,-1], c='darkred', marker='o', mec='k', lw=1, label="training accuracy loss")
	ax_acc.plot(range(epochs), epoch_valid_losses[:,-1], c='tomato',  marker='^', mec='k', lw=1, label="validation accuracy loss" )
	ax_acc.legend()
	ax_acc.grid()
	fig_acc.savefig(args.prefix+"_accuracy", dpi=1200, bbox_inches="tight")


	if args.label == "with_weights_sep":
		fig_weights = plt.figure(num=7, figsize=(4,4))
		ax_weights  = plt.axes()

		ax_weights.plot(range(epochs), epoch_train_losses[:,5], c='darkred', marker='o', mec='k', lw=1, label="training mse weights loss")
		ax_weights.plot(range(epochs), epoch_valid_losses[:,5], c='tomato',  marker='o', mec='k', lw=1, label="validation mse weights loss")

		fig_weights.savefig(args.prefix + "_weights_loss", dpi=1200, bbox_inches="tight")

		fig_split = plt.figure(num=8, figsize=(4,4))
		ax_split  = plt.axes()

		ax_split.plot(range(epochs), epoch_train_losses[:,6], c='darkred', marker='o', mec='k', lw=1, label="training mse splits loss")
		ax_split.plot(range(epochs), epoch_valid_losses[:,6], c='tomato',  marker='o', mec='k', lw=1, label="validation mse splits loss")

		fig_split.savefig(args.prefix + "_splits_loss", dpi=1200, bbox_inches="tight")

	elif args.label == "with_weights_comb":
		fig_split = plt.figure(num=7, figsize=(4,4))
		ax_split  = plt.axes()

		ax_split.plot(range(epochs), epoch_train_losses[:,5], c='darkred', marker='o', mec='k', lw=1, label="training mse splits loss")
		ax_split.plot(range(epochs), epoch_valid_losses[:,5], c='tomato',  marker='o', mec='k', lw=1, label="validation mse splits loss")

		fig_split.savefig(args.prefix + "_splits_loss", dpi=1200, bbox_inches="tight")
		

	return

# get cpu, gpu or mps device for training
device = (
    "cuda" if torch.cuda.is_available()
    else "mps"
    if torch.backends.mps.is_available()
    else "cpu")

print(f"Using {device} device.", flush=True)

# define custom dataset
if args.label == "with_splits":

	# define the chemical potential loss
	def mu_compute(X, ttensor):
		phi_s = ttensor[:,0] 
		phi_p = ttensor[:,1] 
		phi_c = 1-phi_s-phi_p 

		mu_s = torch.log(phi_s) + 1 - phi_s - X[:,0]/X[:,2] * phi_p - \
			X[:,0]/X[:,1] * (phi_c) + X[:,0] * (phi_p**2 * X[:,4] + \
			(phi_c)**2 * X[:,3] + phi_p * phi_c \
			* (X[:,4] + X[:,3] - X[:,5])) 

		mu_p = torch.log(phi_p) + 1 - phi_p - X[:,2]/X[:,0] * phi_s - X[:,2]/X[:,1] * phi_c \
			+ X[:,2] * (phi_s**2 * X[:,4] + phi_c**2 * X[:,5] + phi_s * phi_c * \
			(X[:,4] + X[:,5] - X[:,3]))

		mu_c = torch.log(phi_c) + 1 - phi_c - X[:,1]/X[:,0] * phi_s - \
			X[:,1]/X[:,2] * phi_p + X[:,1] * (phi_s**2 * X[:,3] + \
			phi_p**2 * X[:,5] + phi_s * phi_p * \
			(X[:,3] + X[:,5] - X[:,4]))
		
		# if you are going nan, there is nothing you can do 
		mu_s = torch.nan_to_num(mu_s, nan=0.0, posinf=0.0, neginf=0.0)
		mu_p = torch.nan_to_num(mu_p, nan=0.0, posinf=0.0, neginf=0.0)
		mu_c = torch.nan_to_num(mu_c, nan=0.0, posinf=0.0, neginf=0.0)

		return mu_s, mu_p, mu_c

	# define the chemical potential loss function
	def mu_loss_fn(X, pred_labels, pred_composition):

		# get the labels
		_, argmax_pred_labels  = torch.max(pred_labels, dim=1)

		# mask_true_pred = (argmax_true_labels == argmax_pred_labels).unsqueeze(1)
		mask_for_one_phase    = (argmax_pred_labels == 0).unsqueeze(1).float()
		mask_for_two_phases   = (argmax_pred_labels == 1).unsqueeze(1).float()

		# if you predict just the one phase, you have to make sure there is no chemical 
		# potential loss

		# if you predict two phases, you have to check if the two phases have the same
		# chemical potential 

		pred_composition_phase_one   = pred_composition[:, 0:2].clone()
		pred_composition_phase_two   = pred_composition[:, 3:5].clone()
		pred_composition_phase_three = pred_composition[:, 6:8].clone()

		# compute chemical potentials 
		mu_s1, mu_p1, mu_c1 = mu_compute(X, pred_composition_phase_one)
		mu_s2, mu_p2, mu_c2 = mu_compute(X, pred_composition_phase_two)
		mu_s3, mu_p3, mu_c3 = mu_compute(X, pred_composition_phase_three)

		# calculate the differences in chemical potential between first phase and second phase 
		delta_mu_first_second = ((mu_s1 - mu_s2)**2 + (mu_p1 - mu_p2)**2 + (mu_c1 - mu_c2)**2)
		delta_mu_first_third  = ((mu_s1 - mu_s3)**2 + (mu_p1 - mu_p3)**2 + (mu_c1 - mu_c3)**2)
		delta_mu_second_third = ((mu_s2 - mu_s3)**2 + (mu_p2 - mu_p3)**2 + (mu_c2 - mu_c3)**2)

		# if you have only one phase, ignore all the computation
		delta_mu_first_second *= (1-mask_for_one_phase)
		delta_mu_first_third  *= (1-mask_for_one_phase)
		delta_mu_second_third *= (1-mask_for_one_phase)

		#if you have two phases, only count the delta_mu between the first and second phase 
		delta_mu_first_third  *= (1-mask_for_two_phases)
		delta_mu_second_third *= (1-mask_for_two_phases)

		# get every delta_mu together and compute the mean 
		delta_mu_total = torch.mean(delta_mu_first_second + delta_mu_first_third + delta_mu_second_third)
		
		return delta_mu_total

	# define the unity loss
	def unity_loss_fn(pred_composition):
		# I need to calculate penalties for all three compositions
		# Calculate the penalty if the sum of elements is greater than 1
		phi_c1 = torch.sum(pred_composition[:,:2], dim=1) - 1
		phi_c1 = torch.clamp(phi_c1, min=0)     # Ensure penalty is non-negative
		phi_c1 = torch.mean(phi_c1)             # Apply penalty factor

		phi_c2 = torch.sum(pred_composition[:,2:4], dim=1) - 1
		phi_c2 = torch.clamp(phi_c2, min=0)     # Ensure penalty is non-negative
		phi_c2 = torch.mean(phi_c2)             # Apply penalty factor

		phi_c3 = torch.sum(pred_composition[:,4:6], dim=1) - 1
		phi_c3 = torch.clamp(phi_c3, min=0)     # Ensure penalty is non-negative
		phi_c3 = torch.mean(phi_c3)             # Apply penalty factor

		return phi_c1 + phi_c2 + phi_c3

	# define the loss for the compositions
	comp_loss_fn  = nn.MSELoss()

	# define the loss for the labels
	label_loss_fn = nn.CrossEntropyLoss()

	# get the weights
	W       = args.w_comp + args.w_cross + args.w_unity + args.w_mu 
	w_comp  = args.w_comp  / W
	w_cross = args.w_cross / W
	w_unity = args.w_unity / W
	w_mu    = args.w_mu    / W

	# define the col indices
	input_col_idx       = [0, 1, 2, 3, 4, 5, 6, 7]
	label_col_idx       = [8, 9, 10]
	composition_col_idx = [11, 12, 13, 14, 15, 16]

	# define the dataset
	class CustomDataset(Dataset):
		def __init__(self, dataframe, input_cols, label_cols, comp_cols):
			self.dataframe   = dataframe
			self.input_cols  = input_cols
			self.label_cols  = label_cols
			self.comp_cols   = comp_cols

		def __len__(self):
			return len(self.dataframe)

		def __getitem__(self, idx):
			inputs      = torch.tensor(self.dataframe.to_numpy()[idx, self.input_cols], dtype=torch.float32)
			label       = torch.tensor(self.dataframe.to_numpy()[idx, self.label_cols], dtype=torch.float32)
			composition = torch.tensor(self.dataframe.to_numpy()[idx, self.comp_cols],  dtype=torch.float32)
			return inputs, label, composition

	# define neural network
	class NeuralNetwork(nn.Module):
		def __init__(self):
			super().__init__()
			self.flatten    = nn.Flatten()
			self.classifier = nn.Sequential(
				nn.Linear(8, 16),
				nn.ReLU(),
				nn.Linear(16, 16),
				nn.ReLU(),
				nn.Linear(16, 16),
				nn.ReLU(),
				nn.Linear(16, 3),
				)
			self.splitter = nn.Sequential(
				nn.Linear(8, 16),
				nn.ReLU(),
				nn.Linear(16, 16),
				nn.ReLU(),
				nn.Linear(16, 16),
				nn.ReLU(),
				nn.Linear(16,6),
				nn.Sigmoid()
			)

		def forward(self, x):
			# send the inputs (x) to device
			x          = x.to(device)

			# predict labels and composition
			pred_label   = self.classifier(x)
			pred_comp    = self.splitter(x)

			# get the arguments
			_, argmax_label = torch.max(pred_label, 1)

			# start getting masks 
			mask_for_one_phase  = (argmax_label == 0).unsqueeze(1).float()
			mask_for_two_phases = (argmax_label == 1).unsqueeze(1).float()

			# i do not need to make a mask_for_three_phases because 
			# when three phases have been predicted, there is no reason to "zero-out" anything

			# create a clone
			pred_comp_clone = pred_comp.clone()

			# for all those compositions where the only one phase is predicted, everything from column index 3 needs to be zero'd out
			pred_comp_clone[:, 3:] *= (1 - mask_for_one_phase)

			# if you only have one phase, the number in column index 2 (weight of the phase) has to be 1
			pred_comp_clone[:, 2]   = torch.where(mask_for_one_phase == 1, torch.tensor(1.0).to(device), pred_comp_clone[:, 2])

			# if you have two phases, then everything from column index 6 needs to be zero'd out 
			pred_comp_clone[:, 6:] *= (1 - mask_for_two_phases)

			return pred_label, pred_comp
			
	# define the training loop 
	def train(dataloader, model, optimizer):
		# rev up the model to being training
		model.train()

		# define the empty loss lists
		training_loss = []
		comp_mse_loss = []
		entropy_loss  = []
		unity_loss    = []
		mu_loss       = []
		accuracy      = []
		correct       = []
		total         = []
		current       = 0

		# get the size of the dataloader
		size = len(dataloader)

		for batch, (X, labels, composition) in enumerate(dataloader):

			# get X, fracs and mask on device (explicitly)
			X           = X.to(device)
			labels      = labels.to(device)
			composition = composition.to(device)

			# get the predictions
			pred_labels, pred_composition = model(X)

			# get the loss from the labels 
			label_loss    = label_loss_fn(pred_labels, labels)
			comp_loss     = comp_loss_fn(pred_composition, composition)
			u_loss        = unity_loss_fn(pred_composition)
			chem_pot_loss = mu_loss_fn(X, labels.clone(), pred_labels, pred_composition)

			# get the individual losses
			comp_mse_loss.append(comp_loss.item())
			entropy_loss.append(label_loss.item()) 
			unity_loss.append(u_loss.item())
			mu_loss.append(chem_pot_loss.item())

			# evaluate total loss
			loss = w_comp * comp_loss + w_cross * label_loss + w_unity * u_loss + w_mu * chem_pot_loss

			# backprop
			loss.backward()
			optimizer.step()
			optimizer.zero_grad()

			# update loss lists 
			training_loss.append(loss.item())

			# count how many right answers have been given
			_, argmax_predlabel = torch.max(pred_labels, 1)
			_, argmax_label     = torch.max(labels, 1)
			correct.append((argmax_predlabel==argmax_label).sum().item())
			total.append(len(X))
			accuracy.append(correct[-1]/total[-1])

			# update current
			current += len(X)
			
			if batch % 20 == 0:
				loss     = loss.item()
				# print(f"pred composition = {pred_composition[-5:]}")
				print(f"Current comp MSE = {comp_loss.item()}; " +\
		  			  f"current cross entropy = {label_loss.item()}; " + \
					  f"current unity loss = {unity_loss[-1]}; " + \
					  f"current chem pot loss = {mu_loss[-1]}; " + \
					  f"current net accuracy = {sum(correct)/sum(total)}.", flush=True)
				print(f"batch #{batch+1} loss: {loss:>7f} [processed: {current:>5d}/{size*len(X):>5d}]", flush=True)
			
		return [np.mean(training_loss), np.mean(comp_mse_loss), np.mean(entropy_loss), np.mean(unity_loss), np.mean(mu_loss), np.mean(accuracy)]

	# define the test loop
	def test(dataloader, model):
		# set the model to evaluation mode -- important for batch normalization and dropout layers
		# unnecessary in this situation but added for best practices 
		model.eval()

		# set up all the other numbers to keep track of
		testing_loss     = []
		comp_mse_loss    = []
		entropy_loss     = []
		unity_loss       = []
		mu_loss          = []
		correct          = []
		total            = []
		accuracy         = []

		# Evaluating the model with torch.no_grad() ensures that no gradients are computed during test mode
		# also serves to reduce unnecessary gradient computations and memory usage for tensors with requires_grad=True
		with torch.no_grad():
			for batch, (X, labels, composition) in enumerate(dataloader):
				# get all the objects on to device
				X           = X.to(device)
				labels      = labels.to(device)
				composition = composition.to(device)

				# start running the model
				pred_labels, pred_composition = model(X)

				# get the loss from the labels
				label_loss    = label_loss_fn(pred_labels, labels)
				comp_loss     = comp_loss_fn (pred_composition, composition)
				u_loss        = unity_loss_fn(pred_composition)
				chem_pot_loss = mu_loss_fn(X, labels.clone(), pred_labels, pred_composition)			

				# get the individual losses
				comp_mse_loss.append(comp_loss.item())
				entropy_loss.append(label_loss.item()) 
				unity_loss.append(u_loss.item())
				mu_loss.append(chem_pot_loss.item())

				# evaluated total loss
				loss = w_cross * label_loss + w_comp * comp_loss + w_unity * u_loss + w_mu * chem_pot_loss

				# count how many right answers have been given
				_, argmax_predlabel = torch.max(pred_labels,   1)
				_, argmax_label     = torch.max(labels, 1)
				correct.append((argmax_predlabel==argmax_label).sum().item())
				total.append(len(X))

				# update loss lists
				testing_loss.append(loss.item())
				accuracy.append(correct[-1]/total[-1])
		
		print(f"total correct predictions = {sum(correct)} out of total = {sum(total)}.")
		return [np.mean(testing_loss), np.mean(comp_mse_loss), np.mean(entropy_loss), np.mean(unity_loss), np.mean(mu_loss), np.mean(accuracy)]

elif args.label == "with_weights_sep":
	
	# define the chemical potential loss
	def mu_test(X, ttensor):
		phi_s = ttensor[:,0] # torch.clamp(ttensor[:,0], min=0.001, max=0.999)
		phi_p = ttensor[:,1] # torch.clamp(ttensor[:,1], min=0.001, max=0.999)
		phi_c = 1-phi_s-phi_p # torch.clamp(1-ttensor[:,0]-ttensor[:,1], min=0.001, max=0.999)

		mu_s = torch.log(phi_s) + 1 - phi_s - X[:,0]/X[:,2] * phi_p - \
			X[:,0]/X[:,1] * (phi_c) + X[:,0] * (phi_p**2 * X[:,4] + \
			(phi_c)**2 * X[:,3] + phi_p * phi_c \
			* (X[:,4] + X[:,3] - X[:,5])) 

		mu_p = torch.log(phi_p) + 1 - phi_p - X[:,2]/X[:,0] * phi_s - X[:,2]/X[:,1] * phi_c \
			+ X[:,2] * (phi_s**2 * X[:,4] + phi_c**2 * X[:,5] + phi_s * phi_c * \
			(X[:,4] + X[:,5] - X[:,3]))

		mu_c = torch.log(phi_c) + 1 - phi_c - X[:,1]/X[:,0] * phi_s - \
			X[:,1]/X[:,2] * phi_p + X[:,1] * (phi_s**2 * X[:,3] + \
			phi_p**2 * X[:,5] + phi_s * phi_p * \
			(X[:,3] + X[:,5] - X[:,4]))
		
		mu_s = torch.nan_to_num(mu_s, nan=0.0, posinf=0.0, neginf=0.0)
		mu_p = torch.nan_to_num(mu_p, nan=0.0, posinf=0.0, neginf=0.0)
		mu_c = torch.nan_to_num(mu_c, nan=0.0, posinf=0.0, neginf=0.0)

		return mu_s, mu_p, mu_c

	# define the loss function
	def mu_loss_fn(X, true_labels, pred_labels, pred_composition):
		mu_loss = torch.tensor(0.0, device=device)

		_, argmax_true_labels  = torch.max(true_labels, dim=1)
		_, argmax_pred_labels  = torch.max(pred_labels, dim=1)

		mask_true_pred = (argmax_true_labels == argmax_pred_labels).unsqueeze(1)
		mask_for_1 = (argmax_true_labels == 1).unsqueeze(1)

		# print(f"true_label = {true_labels}, pred_label = {pred_labels}")
		# print(f"pred_composition = {pred_composition}")
		if torch.logical_and(mask_for_1, mask_true_pred).any():
			mask_for_1     = mask_for_1.float()
			mask_true_pred = mask_true_pred.float()
			mask_combined  = mask_for_1 * mask_true_pred

			# Clone pred_composition before modifying it
			pred_composition_1 = pred_composition[:, 0:2].clone()
			pred_composition_2 = pred_composition[:, 2:4].clone()

			pred_composition_1 *= mask_combined
			pred_composition_2 *= mask_combined

			pred_composition_1 += 0.3 * (1 - mask_combined)
			pred_composition_2 += 0.3 * (1 - mask_combined)

			mu_s1, mu_p1, mu_c1 = mu_test(X, pred_composition_1)
			mu_s2, mu_p2, mu_c2 = mu_test(X, pred_composition_2)

			mu_diff = torch.mean(((mu_s1 - mu_s2) ** 2 + (mu_p1 - mu_p2) ** 2 + (mu_c1 - mu_c2) ** 2) * mask_for_1)
			mu_loss += mu_diff

		mask_for_2 = (argmax_true_labels == 2).unsqueeze(1)
		if mask_for_2.any():
			mask_for_2 = mask_for_2.float()
			mu_s1, mu_p1, mu_c1 = mu_test(X, pred_composition[:, 0:2])
			mu_s2, mu_p2, mu_c2 = mu_test(X, pred_composition[:, 2:4])
			mu_s3, mu_p3, mu_c3 = mu_test(X, pred_composition[:, 4:6])

			mu_diff = torch.mean(((mu_s1 - mu_s2) ** 2   + (mu_p1 - mu_p2) ** 2 + (mu_c1 - mu_c2) ** 2 +
									(mu_s2 - mu_s3) ** 2 + (mu_s1 - mu_s3) ** 2 + (mu_p1 - mu_p3) ** 2 +
									(mu_p2 - mu_p3) ** 2 + (mu_c1 - mu_c3) ** 2 + (mu_c2 - mu_c3) ** 2) * mask_for_2)
			mu_loss += mu_diff

		return mu_loss

	# define the unity loss
	def unity_loss_fn(pred_composition):
		# I need to calculate penalties for all three compositions
		# Calculate the penalty if the sum of elements is greater than 1
		phi_c1 = torch.sum(pred_composition[:,:2], dim=1) - 1
		phi_c1 = torch.clamp(phi_c1, min=0)     # Ensure penalty is non-negative
		phi_c1 = torch.mean(phi_c1)             # Apply penalty factor

		phi_c2 = torch.sum(pred_composition[:,2:4], dim=1) - 1
		phi_c2 = torch.clamp(phi_c2, min=0)     # Ensure penalty is non-negative
		phi_c2 = torch.mean(phi_c2)             # Apply penalty factor

		phi_c3 = torch.sum(pred_composition[:,4:6], dim=1) - 1
		phi_c3 = torch.clamp(phi_c3, min=0)     # Ensure penalty is non-negative
		phi_c3 = torch.mean(phi_c3)             # Apply penalty factor

		return phi_c1 + phi_c2 + phi_c3

	# define the split_loss 
	def split_loss_fn(X, pred_composition, pred_weights):
		final_split = pred_composition[:,0:2] * pred_weights[:,0:1] + pred_composition[:,2:4] * pred_weights[:,1:2] + pred_composition[:,4:6] * pred_weights[:,2:3]
		split_loss  = torch.mean(torch.sum((X[:,-2:]-final_split)**2, axis=1))
		return split_loss

	# define the loss for the compositions
	comp_loss_fn   = nn.MSELoss()

	# define the loss for the labels
	label_loss_fn  = nn.CrossEntropyLoss()

	# define the loss for the weights 
	weight_loss_fn = nn.MSELoss()

	# define the dataset
	class CustomDataset(Dataset):
		def __init__(self, dataframe, input_cols, label_cols, comp_cols, weight_cols):
			self.dataframe   = dataframe
			self.input_cols  = input_cols
			self.label_cols  = label_cols
			self.comp_cols   = comp_cols
			self.weight_cols = weight_cols

		def __len__(self):
			return len(self.dataframe)

		def __getitem__(self, idx):
			inputs      = torch.tensor(self.dataframe.to_numpy()[idx, self.input_cols],  dtype=torch.float32)
			label       = torch.tensor(self.dataframe.to_numpy()[idx, self.label_cols],  dtype=torch.float32)
			composition = torch.tensor(self.dataframe.to_numpy()[idx, self.comp_cols],   dtype=torch.float32)
			weights     = torch.tensor(self.dataframe.to_numpy()[idx, self.weight_cols], dtype=torch.float32)
			return inputs, label, composition, weights

	# define neural network
	class NeuralNetwork(nn.Module):
		def __init__(self):
			super().__init__()
			self.flatten    = nn.Flatten()
			self.classifier = nn.Sequential(
				nn.Linear(8, 16),
				nn.ReLU(),
				nn.Linear(16, 16),
				nn.ReLU(),
				nn.Linear(16, 16),
				nn.ReLU(),
				nn.Linear(16, 3),
				)
			self.splitter = nn.Sequential(
				nn.Linear(8, 16),
				nn.ReLU(),
				nn.Linear(16, 16),
				nn.ReLU(),
				nn.Linear(16, 16),
				nn.ReLU(),
				nn.Linear(16,6),
				nn.Sigmoid()
			)
			self.weights = nn.Sequential(
				nn.Linear(8, 16),
				nn.ReLU(),
				nn.Linear(16, 16),
				nn.ReLU(),
				nn.Linear(16, 16),
				nn.ReLU(),
				nn.Linear(16,3),
				nn.Sigmoid()
			)

		def forward(self, x):
			# send the inputs (x) to device
			x          = x.to(device)

			# predict labels and composition
			pred_label   = self.classifier(x)
			pred_comp    = self.splitter(x)
			pred_weights = self.weights(x)

			# get the arguments
			_, argmax_label = torch.max(pred_label, 1)

			# start getting masks
			mask_for_0 = (argmax_label==0).unsqueeze(1) 
			if mask_for_0.any():
				mask_for_0                = 1 - mask_for_0.float()
				pred_comp_clone           = pred_comp.clone() # create a clone of pred_comps
				pred_weights_clone        = pred_weights.clone()
				pred_comp_clone[:,2:]    *= mask_for_0
				pred_weights_clone[:,1:] *= mask_for_0
				pred_weights_clone[:,0]   = 1
				pred_comp                 = pred_comp_clone
				pred_weights              = pred_weights_clone
			
			mask_for_1 = (argmax_label==1).unsqueeze(1)
			if mask_for_1.any():
				mask_for_1                = 1 - mask_for_1.float()
				pred_comp_clone           = pred_comp.clone()
				pred_weights_clone        = pred_weights.clone()
				pred_comp_clone[:,4:]    *= mask_for_1
				pred_weights_clone[:,2:] *= mask_for_1
				pred_comp                 = pred_comp_clone
				pred_weights              = pred_weights_clone
						
			return pred_label, pred_comp, pred_weights

	# get the weights
	W       = args.w_wts + args.w_comp + args.w_cross + args.w_unity + args.w_mu + args.w_split
	w_wts   = args.w_wts   / W
	w_comp  = args.w_comp  / W 
	w_cross = args.w_cross / W
	w_unity = args.w_unity / W 
	w_mu    = args.w_mu    / W 
	w_split = args.w_split / W

	# define the training loop 
	def train(dataloader, model, optimizer):
		# rev up the model to being training
		model.train()

		# define the empty loss lists
		training_loss  = []
		comp_mse_loss  = []
		wts_mse_loss   = []
		split_mse_loss = []
		entropy_loss   = []
		unity_loss     = []
		mu_loss        = []
		accuracy       = []
		correct        = []
		total          = []
		current        = 0

		# get the size of the dataloader
		size = len(dataloader)

		for batch, (X, labels, composition, weights) in enumerate(dataloader):

			# get X, fracs and mask on device (explicitly)
			X           = X.to(device)
			labels      = labels.to(device)
			composition = composition.to(device)
			weights     = weights.to(device)

			# get the predictions
			pred_labels, pred_composition, pred_weights = model(X)

			# get the loss from the labels 
			weights_loss  = weight_loss_fn(pred_weights, weights)
			split_loss    = split_loss_fn (X, pred_composition, pred_weights)
			label_loss    = label_loss_fn(pred_labels, labels)
			comp_loss     = comp_loss_fn(pred_composition, composition)
			u_loss        = unity_loss_fn(pred_composition)
			chem_pot_loss = mu_loss_fn(X, labels.clone(), pred_labels, pred_composition)

			# get the individual losses
			comp_mse_loss.append (comp_loss.item())
			wts_mse_loss.append  (weights_loss.item())
			split_mse_loss.append(split_loss.item())
			entropy_loss.append  (label_loss.item())
			unity_loss.append    (u_loss.item())
			mu_loss.append       (chem_pot_loss.item())

			# evaluate total loss
			loss = w_comp * comp_loss + w_wts * weights_loss + w_split * split_loss + w_cross * label_loss + w_unity * u_loss + w_mu * chem_pot_loss

			# backprop
			loss.backward()
			optimizer.step()
			optimizer.zero_grad()

			# update loss lists 
			training_loss.append(loss.item())

			# count how many right answers have been given
			_, argmax_predlabel = torch.max(pred_labels, 1)
			_, argmax_label     = torch.max(labels, 1)
			correct.append((argmax_predlabel==argmax_label).sum().item())
			total.append(len(X))
			accuracy.append(correct[-1]/total[-1])

			# update current
			current += len(X)
			
			if batch % 20 == 0:
				loss     = loss.item()
				print(f"Current comp loss = {comp_loss.item()}; "       + \
					  f"current cross entropy = {label_loss.item()}; "      + \
					  f"current unity loss = {unity_loss[-1]}; "         + \
					  f"current chem pot loss = {mu_loss[-1]};"             + \
					  f"current wt loss = {weights_loss.item()}; "    + \
					  f"current split loss = {split_loss.item()}; "      + \
					  f"current net accuracy = {sum(correct)/sum(total)}; ", flush=True)
				print(f"batch #{batch+1} loss: {loss:>7f} [processed: {current:>5d}/{size*len(X):>5d}]", flush=True)
			
		return [np.mean(training_loss), np.mean(comp_mse_loss), np.mean(entropy_loss), np.mean(unity_loss), np.mean(mu_loss), np.mean(wts_mse_loss), np.mean(split_mse_loss), np.mean(accuracy)]

	# define the test loop
	def test(dataloader, model):
		# set the model to evaluation mode -- important for batch normalization and dropout layers
		# unnecessary in this situation but added for best practices 
		model.eval()

		# set up all the other numbers to keep track of
		testing_loss     = []
		comp_mse_loss    = []
		wts_mse_loss     = []
		split_mse_loss   = []
		entropy_loss     = []
		unity_loss       = []
		mu_loss          = []
		correct          = []
		total            = []
		accuracy         = []

		# Evaluating the model with torch.no_grad() ensures that no gradients are computed during test mode
		# also serves to reduce unnecessary gradient computations and memory usage for tensors with requires_grad=True
		with torch.no_grad():
			for batch, (X, labels, composition, weights) in enumerate(dataloader):
				# get all the objects on to device
				X           = X.to(device)
				labels      = labels.to(device)
				composition = composition.to(device)
				weights     = weights.to(device)

				# start running the model
				pred_labels, pred_composition, pred_weights = model(X)

				# get the loss from the labels
				label_loss    = label_loss_fn(pred_labels, labels)
				comp_loss     = comp_loss_fn (pred_composition, composition)
				weights_loss  = weight_loss_fn(pred_weights, weights)
				u_loss        = unity_loss_fn(pred_composition)
				chem_pot_loss = mu_loss_fn(X, labels.clone(), pred_labels, pred_composition)			
				split_loss    = split_loss_fn(X, pred_composition, pred_weights)

				# get the individual losses
				wts_mse_loss.append(weights_loss.item())
				comp_mse_loss.append(comp_loss.item())
				split_mse_loss.append(split_loss.item())
				entropy_loss.append(label_loss.item()) 
				unity_loss.append(u_loss.item())
				mu_loss.append(chem_pot_loss.item())
				

				# evaluated total loss
				loss = w_cross * label_loss + w_comp * comp_loss + w_split * split_loss + w_wts * weights_loss + w_unity * u_loss + w_mu * chem_pot_loss

				# count how many right answers have been given
				_, argmax_predlabel = torch.max(pred_labels,   1)
				_, argmax_label     = torch.max(labels, 1)
				correct.append((argmax_predlabel==argmax_label).sum().item())
				total.append(len(X))

				# update loss lists
				testing_loss.append(loss.item())
				accuracy.append(correct[-1]/total[-1])
		
		print(f"total correct predictions = {sum(correct)} out of total = {sum(total)}.")
		return [np.mean(testing_loss), np.mean(comp_mse_loss), np.mean(entropy_loss), np.mean(unity_loss), np.mean(mu_loss), np.mean(wts_mse_loss), np.mean(split_mse_loss), np.mean(accuracy)]

	# define the col indices
	input_col_idx       = [0, 1, 2, 3, 4, 5, 6, 7]
	label_col_idx       = [8, 9, 10]
	composition_col_idx = [11, 12, 13, 14, 15, 16]
	weight_col_idx      = [17, 18, 19]

elif args.label == "with_weights_comb":
	
	# define the chemical potential loss
	def mu_compute(X, ttensor):
		phi_s = ttensor[:,0] 
		phi_p = ttensor[:,1] 
		phi_c = 1-phi_s-phi_p 

		mu_s = torch.log(phi_s) + 1 - phi_s - X[:,0]/X[:,2] * phi_p - \
			X[:,0]/X[:,1] * (phi_c) + X[:,0] * (phi_p**2 * X[:,4] + \
			(phi_c)**2 * X[:,3] + phi_p * phi_c \
			* (X[:,4] + X[:,3] - X[:,5])) 

		mu_p = torch.log(phi_p) + 1 - phi_p - X[:,2]/X[:,0] * phi_s - X[:,2]/X[:,1] * phi_c \
			+ X[:,2] * (phi_s**2 * X[:,4] + phi_c**2 * X[:,5] + phi_s * phi_c * \
			(X[:,4] + X[:,5] - X[:,3]))

		mu_c = torch.log(phi_c) + 1 - phi_c - X[:,1]/X[:,0] * phi_s - \
			X[:,1]/X[:,2] * phi_p + X[:,1] * (phi_s**2 * X[:,3] + \
			phi_p**2 * X[:,5] + phi_s * phi_p * \
			(X[:,3] + X[:,5] - X[:,4]))
		
		mu_s = torch.nan_to_num(mu_s, nan=0.0, posinf=0.0, neginf=0.0)
		mu_p = torch.nan_to_num(mu_p, nan=0.0, posinf=0.0, neginf=0.0)
		mu_c = torch.nan_to_num(mu_c, nan=0.0, posinf=0.0, neginf=0.0)

		return mu_s, mu_p, mu_c

	# define the chemical potential loss function
	def mu_loss_fn(X, pred_labels, pred_composition):

		# get the labels
		_, argmax_pred_labels  = torch.max(pred_labels, dim=1)

		# mask_true_pred = (argmax_true_labels == argmax_pred_labels).unsqueeze(1)
		mask_for_one_phase    = (argmax_pred_labels == 0).unsqueeze(1).float()
		mask_for_two_phases   = (argmax_pred_labels == 1).unsqueeze(1).float()

		# if you predict just the one phase, you have to make sure there is no chemical 
		# potential loss

		# if you predict two phases, you have to check if the two phases have the same
		# chemical potential 

		pred_composition_phase_one   = pred_composition[:, 0:2].clone()
		pred_composition_phase_two   = pred_composition[:, 3:5].clone()
		pred_composition_phase_three = pred_composition[:, 6:8].clone()

		# compute chemical potentials 
		mu_s1, mu_p1, mu_c1 = mu_compute(X, pred_composition_phase_one)
		mu_s2, mu_p2, mu_c2 = mu_compute(X, pred_composition_phase_two)
		mu_s3, mu_p3, mu_c3 = mu_compute(X, pred_composition_phase_three)

		# calculate the differences in chemical potential between first phase and second phase 
		delta_mu_first_second = ((mu_s1 - mu_s2)**2 + (mu_p1 - mu_p2)**2 + (mu_c1 - mu_c2)**2)
		delta_mu_first_third  = ((mu_s1 - mu_s3)**2 + (mu_p1 - mu_p3)**2 + (mu_c1 - mu_c3)**2)
		delta_mu_second_third = ((mu_s2 - mu_s3)**2 + (mu_p2 - mu_p3)**2 + (mu_c2 - mu_c3)**2)

		# if you have only one phase, ignore all the computation
		delta_mu_first_second *= (1-mask_for_one_phase)
		delta_mu_first_third  *= (1-mask_for_one_phase)
		delta_mu_second_third *= (1-mask_for_one_phase)

		#if you have two phases, only count the delta_mu between the first and second phase 
		delta_mu_first_third  *= (1-mask_for_two_phases)
		delta_mu_second_third *= (1-mask_for_two_phases)

		# get every delta_mu together and compute the mean 
		delta_mu_total = torch.mean(delta_mu_first_second + delta_mu_first_third + delta_mu_second_third)
		
		return delta_mu_total

	# define the unity loss
	def unity_loss_fn(pred_composition):
		# I need to calculate penalties for all three compositions
		# Calculate the penalty if the sum of elements is greater than 1
		phi_c1 = torch.sum  (pred_composition[:,0:2], dim=1) - 1
		phi_c1 = torch.clamp(phi_c1, min=0)     # Ensure penalty is non-negative
		phi_c1 = torch.mean (phi_c1)             # Apply penalty factor

		phi_c2 = torch.sum  (pred_composition[:,3:5], dim=1) - 1
		phi_c2 = torch.clamp(phi_c2, min=0)     # Ensure penalty is non-negative
		phi_c2 = torch.mean (phi_c2)             # Apply penalty factor

		phi_c3 = torch.sum  (pred_composition[:,6:8], dim=1) - 1
		phi_c3 = torch.clamp(phi_c3, min=0)     # Ensure penalty is non-negative
		phi_c3 = torch.mean (phi_c3)             # Apply penalty factor

		return phi_c1 + phi_c2 + phi_c3

	# define the split_loss 
	def split_loss_fn(X, pred_composition):
		final_split = pred_composition[:,0:2] * pred_composition[:,2:3] + pred_composition[:,3:5] * pred_composition[:,5:6] + pred_composition[:,6:8] * pred_composition[:,8:9]
		split_loss  = torch.mean(torch.sum((X[:,-2:]-final_split)**2, axis=1))
		return split_loss

	# define the loss for the compositions
	comp_loss_fn  = nn.MSELoss()

	# define the loss for the labels
	label_loss_fn = nn.CrossEntropyLoss()

	# define the dataset
	class CustomDataset(Dataset):
		def __init__(self, dataframe, input_cols, label_cols, comp_cols):
			self.dataframe   = dataframe
			self.input_cols  = input_cols
			self.label_cols  = label_cols
			self.comp_cols   = comp_cols

		def __len__(self):
			return len(self.dataframe)

		def __getitem__(self, idx):
			inputs      = torch.tensor(self.dataframe.to_numpy()[idx, self.input_cols],  dtype=torch.float32)
			label       = torch.tensor(self.dataframe.to_numpy()[idx, self.label_cols],  dtype=torch.float32)
			composition = torch.tensor(self.dataframe.to_numpy()[idx, self.comp_cols],   dtype=torch.float32)
			return inputs, label, composition

	# define neural network
	class NeuralNetwork(nn.Module):
		def __init__(self):
			super().__init__()
			self.flatten    = nn.Flatten()
			self.classifier = nn.Sequential(
				nn.Linear(8, 16),
				nn.ReLU(),
				nn.Linear(16, 16),
				nn.ReLU(),
				nn.Linear(16, 16),
				nn.ReLU(),
				nn.Linear(16, 3),
				)
			self.splitter = nn.Sequential(
				nn.Linear(8, 16),
				nn.ReLU(),
				nn.Linear(16, 16),
				nn.ReLU(),
				nn.Linear(16, 16),
				nn.ReLU(),
				nn.Linear(16,9),
				nn.Sigmoid()
			)

		def forward(self, x):
			# send the inputs (x) to device
			x          = x.to(device)

			# predict labels and composition
			pred_label   = self.classifier(x)
			pred_comp    = self.splitter(x)

			# get the arguments
			_, argmax_label = torch.max(pred_label, 1)

			# start getting masks 
			mask_for_one_phase  = (argmax_label == 0).unsqueeze(1).float()
			mask_for_two_phases = (argmax_label == 1).unsqueeze(1).float()

			# i do not need to make a mask_for_three_phases because 
			# when three phases have been predicted, there is no reason to "zero-out" anything

			# create a clone
			pred_comp_clone = pred_comp.clone()

			# for all those compositions where the only one phase is predicted, everything from column index 3 needs to be zero'd out
			pred_comp_clone[:, 3:] *= (1 - mask_for_one_phase)

			# if you only have one phase, the number in column index 2 (weight of the phase) has to be 1
			pred_comp_clone[:, 2]   = torch.where(mask_for_one_phase == 1, torch.tensor(1.0).to(device), pred_comp_clone[:, 2])

			# if you have two phases, then everything from column index 6 needs to be zero'd out 
			pred_comp_clone[:, 6:] *= (1 - mask_for_two_phases)

			return pred_label, pred_comp

	# get the weights
	W       = args.w_comp + args.w_cross + args.w_unity + args.w_mu + args.w_split
	w_comp  = args.w_comp  / W 
	w_cross = args.w_cross / W
	w_unity = args.w_unity / W 
	w_mu    = args.w_mu    / W 
	w_split = args.w_split / W

	# define the training loop 
	def train(dataloader, model, optimizer):
		# rev up the model to being training
		model.train()

		# define the empty loss lists
		training_loss  = []
		comp_mse_loss  = []
		split_mse_loss = []
		entropy_loss   = []
		unity_loss     = []
		mu_loss        = []
		accuracy       = []
		correct        = []
		total          = []
		current        = 0

		# get the size of the dataloader
		size = len(dataloader)

		for batch, (X, labels, composition) in enumerate(dataloader):

			# get X, fracs and mask on device (explicitly)
			X           = X.to(device)
			labels      = labels.to(device)
			composition = composition.to(device)

			# get the predictions
			pred_labels, pred_composition = model(X)

			# get the loss from the labels 
			label_loss    = label_loss_fn(pred_labels, labels)
			comp_loss     = comp_loss_fn(pred_composition, composition)
			u_loss        = unity_loss_fn(pred_composition)
			chem_pot_loss = mu_loss_fn(X, labels.clone(), pred_labels, pred_composition)
			split_loss    = split_loss_fn(X, pred_composition) 

			# get the individual losses
			comp_mse_loss.append(comp_loss.item())
			entropy_loss.append(label_loss.item()) 
			unity_loss.append(u_loss.item())
			mu_loss.append(chem_pot_loss.item())
			split_mse_loss.append(split_loss.item())

			# evaluate total loss
			loss = w_comp * comp_loss + w_cross * label_loss + w_split * split_loss + w_unity * u_loss + w_mu * chem_pot_loss

			# backprop
			loss.backward()
			optimizer.step()
			optimizer.zero_grad()

			# update loss lists 
			training_loss.append(loss.item())

			# count how many right answers have been given
			_, argmax_predlabel = torch.max(pred_labels, 1)
			_, argmax_label     = torch.max(labels, 1)
			correct.append((argmax_predlabel==argmax_label).sum().item())
			total.append(len(X))
			accuracy.append(correct[-1]/total[-1])

			# update current
			current += len(X)

			if batch % 20 == 0:
				loss     = loss.item()
				print(f"Current comp loss = {comp_loss.item()};"      + \
		  			  f"current cross entropy = {label_loss.item()};" + \
					  f"current unity loss = {unity_loss[-1]}; "      + \
					  f"current split loss = {split_mse_loss[-1]}; "      + \
					  f"current chem pot loss = {mu_loss[-1]};", end=' ', flush=True)
				print(f"current net accuracy = {sum(correct)/sum(total)};", end=' ', flush=True)
				print(f"batch #{batch+1} loss: {loss:>7f} [processed: {current:>5d}/{size*len(X):>5d}]", flush=True)
			
		return [np.mean(training_loss), np.mean(comp_mse_loss), np.mean(entropy_loss), np.mean(unity_loss), np.mean(mu_loss), np.mean(split_mse_loss), np.mean(accuracy)]

	# define the test loop
	def test(dataloader, model):
		# set the model to evaluation mode -- important for batch normalization and dropout layers
		# unnecessary in this situation but added for best practices 
		model.eval()

		# set up all the other numbers to keep track of
		testing_loss     = []
		comp_mse_loss    = []
		split_mse_loss   = []
		entropy_loss     = []
		unity_loss       = []
		mu_loss          = []
		correct          = []
		total            = []
		accuracy         = []

		# Evaluating the model with torch.no_grad() ensures that no gradients are computed during test mode
		# also serves to reduce unnecessary gradient computations and memory usage for tensors with requires_grad=True
		with torch.no_grad():
			for batch, (X, labels, composition) in enumerate(dataloader):
				# get all the objects on to device
				X           = X.to(device)
				labels      = labels.to(device)
				composition = composition.to(device)

				# start running the model
				pred_labels, pred_composition = model(X)

				# get the loss from the labels
				label_loss    = label_loss_fn(pred_labels, labels)
				comp_loss     = comp_loss_fn (pred_composition, composition)
				u_loss        = unity_loss_fn(pred_composition)
				chem_pot_loss = mu_loss_fn(X, labels.clone(), pred_labels, pred_composition)			
				split_loss    = split_loss_fn(X, pred_composition)

				# get the individual losses
				comp_mse_loss.append(comp_loss.item())
				entropy_loss.append(label_loss.item()) 
				unity_loss.append(u_loss.item())
				mu_loss.append(chem_pot_loss.item())
				split_mse_loss.append(split_loss.item())

				# evaluated total loss
				loss = w_cross * label_loss + w_split * split_loss + w_comp * comp_loss + w_unity * u_loss + w_mu * chem_pot_loss

				# count how many right answers have been given
				_, argmax_predlabel = torch.max(pred_labels,   1)
				_, argmax_label     = torch.max(labels, 1)
				correct.append((argmax_predlabel==argmax_label).sum().item())
				total.append(len(X))

				# update loss lists
				testing_loss.append(loss.item())
				accuracy.append(correct[-1]/total[-1])
		
		print(f"total correct predictions = {sum(correct)} out of total = {sum(total)}.")
		return [np.mean(testing_loss), np.mean(comp_mse_loss), np.mean(entropy_loss), np.mean(unity_loss), np.mean(mu_loss), np.mean(accuracy)]

	# define the col indices
	input_col_idx       = [0, 1, 2, 3, 4, 5, 6, 7]
	label_col_idx       = [8, 9, 10]
	composition_col_idx = [11, 12, 17, 13, 14, 18, 15, 16, 19]

if __name__=="__main__":

	start = time.time()

	if args.label not in ["with_splits", "with_weights_sep", "with_weights_comb"]:
		print(f"Incorrect label provided. Exiting.")
		exit()

	# define the model
	model     = NeuralNetwork().to(device)
	try:
		model.load_state_dict(torch.load(args.getfile, map_location=torch.device(device)))
	except:
		print(f"Starting training afresh.", flush=True)

	# define an optimizer
	learning_rate = 1e-3
	optimizer     = torch.optim.Adam(model.parameters(), lr=learning_rate)

	# begin extracting the database
	df      = pd.read_csv(args.db, sep="\|", engine="python")

	# set everything for the neural network with splits
	if args.label == "with_splits":
		# set up the dataset
		dataset = CustomDataset(df, input_col_idx, label_col_idx, composition_col_idx)

		# define the training size and testing size
		train_size = int(0.8 * len(dataset))
		valid_size  = len(dataset) - train_size

		# define other hyperparameters
		batch_size    = train_size//100 if train_size // 100 > 0 else 1

		# split the data set into training and testing
		train_dataset, valid_dataset = random_split(dataset, [train_size, valid_size], generator=torch.Generator().manual_seed(42))

		# cast them as DataLoader objects
		dataloader_train = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
		dataloader_valid = DataLoader(valid_dataset, batch_size=batch_size, shuffle=True)

		# print some things...
		print(f"Training dataset size: {len(dataloader_train)}", flush=True)
		print(f"Testing dataset size:  {len(dataloader_valid)}", flush=True)

		# define some epochal containers
		epochs             = args.ep
		epoch_train_losses = []
		epoch_valid_losses = []

		print(f"Begin training loop!", flush=True)
		for t in range(epochs):
			print(f"epoch = {t}...", flush=True)

			# begin training the model
			training_losses   = train(dataloader_train, model, optimizer)
			validation_losses = test (dataloader_valid, model)

			# accumulate all the epochal losses
			epoch_train_losses.append(training_losses)
			epoch_valid_losses.append(validation_losses)
			print(f"==================================================", flush=True)

		# convert everything to a numpy array
		epoch_train_losses = np.array(epoch_train_losses)
		epoch_valid_losses  = np.array(epoch_valid_losses)

		# create all the figures 
		make_figures(epochs, epoch_train_losses, epoch_valid_losses)

	# set everything for the neural network with weights, and involves 3 neural networks
	elif args.label == "with_weights_sep":
		# set up the dataset
		dataset = CustomDataset(df, input_col_idx, label_col_idx, composition_col_idx, weight_col_idx)

		# define the training size and testing size
		train_size = int(0.8 * len(dataset))
		valid_size  = len(dataset) - train_size

		# define other hyperparameters
		batch_size    = train_size//100 if train_size // 100 > 0 else 1

		# split the data set into training and testing
		train_dataset, valid_dataset = random_split(dataset, [train_size, valid_size], generator=torch.Generator().manual_seed(42))

		# cast them as DataLoader objects
		dataloader_train = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
		dataloader_valid = DataLoader(valid_dataset, batch_size=batch_size, shuffle=True)

		# print some things...
		print(f"Training dataset size: {len(dataloader_train)}", flush=True)
		print(f"Testing dataset size:  {len(dataloader_valid)}", flush=True)

		# define some epochal containers
		epochs             = args.ep
		epoch_train_losses = []
		epoch_valid_losses = []

		print(f"Begin training loop!", flush=True)
		for t in range(epochs):
			print(f"epoch = {t}...", flush=True)

			# begin training the model
			training_losses   = train(dataloader_train, model, optimizer)
			validation_losses = test (dataloader_valid, model)

			# accumulate all the epochal losses
			epoch_train_losses.append(training_losses  )
			epoch_valid_losses.append(validation_losses)
			print(f"==================================================", flush=True)

		# convert everything to a numpy array
		epoch_train_losses  = np.array(epoch_train_losses)
		epoch_valid_losses  = np.array(epoch_valid_losses)

		# create all the figures 
		make_figures(epochs, epoch_train_losses, epoch_valid_losses)

	# set everything for the neural network with weights, and involves the standard two enural networks, with the weights as the part of the output
	elif args.label == "with_weights_comb":
		# set up the dataset
		dataset = CustomDataset(df, input_col_idx, label_col_idx, composition_col_idx)

		# define the training size and testing size
		train_size = int(0.8 * len(dataset))
		valid_size  = len(dataset) - train_size

		# define other hyperparameters
		batch_size    = train_size//100 if train_size // 100 > 0 else 1

		# split the data set into training and testing
		train_dataset, valid_dataset = random_split(dataset, [train_size, valid_size], generator=torch.Generator().manual_seed(42))

		# cast them as DataLoader objects
		dataloader_train = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
		dataloader_valid = DataLoader(valid_dataset, batch_size=batch_size, shuffle=True)

		# print some things...
		print(f"Training dataset size: {len(dataloader_train)}", flush=True)
		print(f"Testing dataset size:  {len(dataloader_valid)}", flush=True)

		# define some epochal containers
		epochs             = args.ep
		epoch_train_losses = []
		epoch_valid_losses = []

		print(f"Begin training loop!", flush=True)
		for t in range(epochs):
			print(f"epoch = {t}...", flush=True)

			# begin training the model
			training_losses   = train(dataloader_train, model, optimizer)
			validation_losses = test (dataloader_valid, model)

			# accumulate all the epochal losses
			epoch_train_losses.append(training_losses)
			epoch_valid_losses.append(validation_losses)
			print(f"==================================================", flush=True)

		# convert everything to a numpy array
		epoch_train_losses = np.array(epoch_train_losses)
		epoch_valid_losses = np.array(epoch_valid_losses)

		# create all the figures 
		make_figures(epochs, epoch_train_losses, epoch_valid_losses)

	# done with training the model
	print(f"Saving the model...", flush=True, end=' ')
	torch.save(model.state_dict(), args.dumpfile)
	print(f"Saved!", flush=True)
	
	stop = time.time()
	print(f"Time for model training is {stop-start} seconds.", flush=True)
