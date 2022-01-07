# Position-Decoding-Methods-Based-on-Fluorescence-Calcium-Imaging
Demonstration codes for: Tu M, Zhao R, Adler A, Gan W-B, Chen ZS. 2020. Efficient Position Decoding Methods Based on Fluorescence Calcium Imaging in the Mouse Hippocampus. Neural Comput. 32(6):1144–1167. doi:10.1162/neco_a_01281. https://www.mitpressjournals.org/doi/abs/10.1162/neco_a_01281.

I am editing the README file. Adding some more details about the project description.

# Required packages: 
1. CNMF_E-master (Available at https://github.com/zhoupc/CNMF_E): MATLAB package required for demo_Simulation.m
2. pyhsmm-spiketrains (Available at https://github.com/slinderman/pyhsmm_spiketrains): Python package required for demo_part1_HMM_Decoding.py
3. Circular Statistics Toolbox (Available at https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics): MATLAB toolbox required for demo_OLE_Decoding.m

# Demonstration include:

demo_Simulation.m: Simulate the calcium fluorescence traces with second order autoregressive model. At the end, 2 figures will be plotted: the fluorescence calcium traces and the true spikes, inferred spikes from using spike deconvolution, same as Figure 6B and 6C respectively.

demo_Maximum_Likelihood_Decoding.m: Use maximum likelihood estimator to decode positions with filtered MPP data and plot the inferred and true trajectories.

demo_OLE_Decoding.m: Use optimal linear decoder to decode positions with filtered MPP data and plot the inferred and true trajectories.

demo_part1_HMM_Decoding.py: Use hidden Markov model to uncover latent structures in filtered MPP data and save the inferred hidden states for further decoding in demo_part2_HMM_Decoding.m 

demo_part2_HMM_Decoding.m: Map inferred states to positions and plots the inferred and true positions.

demo_Marked_Point_Process.m: Illustration of how the marked MPP is generated from raw fluorescence trace 
