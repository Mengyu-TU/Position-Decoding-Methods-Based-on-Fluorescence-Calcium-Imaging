%% This script demonstrates how to extract filtered Marked Point Process (MPP) from raw df/f 
clear
clc

addpath('data')

%% Load data
load('data/Raw_DF_F_Two_Photon.mat','RawDF_F')
% RawDF_F: Dimension: Number of Neurons x Number of TimeSteps
NumNeuron = size(RawDF_F,1);

%% Identify and extract peaks from df/f
FluorscePeak = zeros(size(RawDF_F));

for iNeuron = 1:NumNeuron
    y = RawDF_F(iNeuron,:);
    MinPeak = 0; % Minimum peak height
    MinProm = 0;  % Minimum peak prominence
    [pks,locs] = findpeaks(y,'MinPeakHeight',MinPeak,'MinPeakProminence',MinProm); % MATLAB built-in function 
    
    if ~isempty(locs)
        FluorscePeak(iNeuron,locs) = pks;
    end
    
end

%% Convolution
FilterMPP = zeros(size(FluorscePeak));
h1 = [0.1, 0.2, 0.2, 0, 0];  % Optmized filter for this particular RawDF_F
h = h1/sum(h1);              % Normalize so that energy conserved

for iNeuron = 1: NumNeuron
    x = FluorscePeak(iNeuron,:);
    y = conv(x,h,'same');    % MATLAB built-in function 
    FilterMPP(iNeuron,:) = y;
end

%% Illustration: plotting filtered MPP and Raw df/f with a randomly selected neurons
NeuronID =  randi(NumNeuron);
Interval = 1:100;

Raw = RawDF_F(NeuronID,Interval);
Peak = FluorscePeak(NeuronID,Interval);
MPP = FilterMPP(NeuronID,Interval);

MPP_location = find(MPP ~= 0);
MPP_magnitutde = MPP(MPP_location);

peak_location = find(Peak ~= 0);
peak_magnitude = Peak(peak_location);

% Plot
figure;hold on;
plot(Raw,'.-','LineWidth',0.8);
plot(peak_location,peak_magnitude,'*b','MarkerSize',9)
h = stem(MPP_location,MPP_magnitutde,'r','LineWidth',1);
set(h, 'Marker', 'none');
legend('Raw DF/F','Peak Identified','Filtered MPP')
xlabel('Time Steps')
title(['Neuron ',num2str(NeuronID)])





