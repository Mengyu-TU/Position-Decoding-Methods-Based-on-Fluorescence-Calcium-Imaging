#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 14:06:22 2019

@author: mengyu
"""
## This script uses Bayesian hidden Markov model (HMM) developed to decode animals positions
## For details of the HMM model please refer to:
## S.W. et al (2016) "A Bayesian nonparamatric approach fro uncovering rat hippocampal population codes during spatial navigation"
## https://github.com/slinderman/pyhsmm_spiketrains

## The output Z_infer_mat is the inferred state at each time step of the testing data
## The outputs are saved and should be passed to demo_part2_HMM_Decoding.m in MATLBA for calculating
## the decoding errors and visualization purposes

import math

import numpy as np
import matplotlib.pyplot as plt

import scipy.io

import os
os.chdir('Calcium_Demo/HMM_package/pyhsmm_spiketrains ') # Must direct to the folder "pyhsmm_spiketrains" to load the HMM models

from pybasicbayes.util.text import progprint_xrange
from pyhsmm.basic.distributions import PoissonDuration

import pyhsmm_spiketrains.models
reload(pyhsmm_spiketrains.models)

# Set the seed
seed = 0
print "setting seed to ", seed
np.random.seed(seed)

os.chdir('Calcium_Demo/data')   # Direct to the folder with data

# Load data set for training and testing
file_name = 'FilteredMPP4Decoding4HMM.mat'
data = scipy.io.loadmat(file_name)

FilteredMPP = data['FilteredMPP4Decoding'] 

T_total = len(FilteredMPP[:,1])         # Dimension: Number Time Steps x Number Neurons

# Training Data
T_train = math.floor(T_total * 0.9);    # Use 90 % of data for training
S_train = FilteredMPP[:int(T_train),:]  
S_train = np.int64(S_train)

# Testing Data
S_test =  FilteredMPP[int(T_train):,:]  # Data format: column: Number of neurons                 
S_test = np.int64(S_test)

T_test = len(S_test[:,0])
T_train = len(S_train[:,0])

NumNeuron = len(S_train[0,:])           # No. of neuron ensembles
N_iter = 500                            # Number of iterations of Gibbs sampling

max_K = 50   # Number of inferred states in the HMM
mc_iter = 1  # No. of  random selection of neurons to be done
NumNeuronsUsed = 200; 
NeuronSelectedIndex = np.zeros([NumNeuronsUsed,mc_iter])

PosteriorProbTest = np.zeros([T_test,max_K,mc_iter])
N_used_vec = np.zeros([mc_iter])
Z_inf_mat = np.zeros([mc_iter,T_train])


for j in range(mc_iter):                    
    # From all neurons, randomly select 'NumNeuronsUsed' neurons                              
    idxAll = np.random.permutation(NumNeuron)  
    idx = idxAll [:NumNeuronsUsed]
    NeuronSelectedIndex[:,j] = idx
    
    S_test_PARTIAL = S_test[:,idx]
    S_test_PARTIAL = S_test_PARTIAL.copy(order='C')
    S_train_PARTIAL = S_train[:,idx]                     # Data format: column: Number of neurons
    S_train_PARTIAL = S_train_PARTIAL.copy(order='C')
    
    test_hmm = pyhsmm_spiketrains.models.PoissonHMM(N = NumNeuronsUsed, K = max_K)

    test_hmm.add_data(S_train_PARTIAL)
    
    # Fit the test model with Gibbs sampling
    for itr in progprint_xrange(N_iter): 
        test_hmm.resample_model()        
        test_hmm.resample_obs_hypers_hmc()
        
    # Get the inferred state sequence
    test_hmm.relabel_by_usage()
    N_used_inf = len(list(test_hmm.used_states))
    print(N_used_inf)
    N_used_vec[j]=N_used_inf
        
    # Generate state sequence
    Z_inf = test_hmm.stateseqs[0]
    Z_inf_mat[j,:] = Z_inf.transpose()
    
    # Compute the marginal distribution over states at each of the test time steps
    PosteriorProbTest[:,:,j] = test_hmm.heldout_state_marginals(S_test_PARTIAL)
                        
# Save data
scipy.io.savemat('CaIm_exp_inf_1216201901.mat',{'Z_inf_mat':Z_inf_mat, \
                         'N_used_vec':N_used_vec,\
                         'PosteriorProbTest':PosteriorProbTest})


