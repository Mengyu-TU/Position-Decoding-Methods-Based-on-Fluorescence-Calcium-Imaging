%% Optimal Linear Estimation (OLE) Decoder using filtered MPP data generated with Dataset 1 Mouse 1 
clear
clc

addpath(genpath('OLE_package'))
addpath('data')

%% load data for decoding
load('FilteredMPP4Decoding4OLE.mat','FilteredMPP4Decoding4OLE')
Data4Decoding = FilteredMPP4Decoding4OLE;
NumNeuron = size(Data4Decoding,2);

%% load true position
load('Position4DecodingOLE.mat','position','trackLen')
% Varible: position: true positions of the animal
% trackLen: full track length 108 cm

%% Convert position to angles
position2angle = (position/trackLen *2*pi)'; % Convert to 1xNumTSteps

%% Set parameters for OLE decoder
nFold = 10;     % 10 fold cross-validation
NumBasis = 150; % Number of basis functions
kappa = 50;     % Coefficient kappa in the basis function

%% OLE decoding and calculate the mean decoding errors
num_iter = 1;

Median_error = zeros(num_iter,1);
Median_std_of_error = zeros(num_iter,1);
Mean_error = zeros(num_iter,1);
Mean_std_of_error = zeros(num_iter,1);

for iIter = 1 :num_iter
    %% Optimal linear estimation
    % Input Dimension: NumTSteps x NumNeuron
    [probability,pvec,max_position] = OLE_CrossValidate_MY(Data4Decoding,...
        NumBasis,kappa,position2angle,nFold);
    
    %% Circular error
    err =  circ_dist(     position2angle, (max_position/length(pvec)*2*pi))*trackLen/pi;
    err_cir = err;
    
    %% Bootstrap to get the median error
    bootstat = bootstrp(500,'mean',abs(err_cir));
    mean_error = mean(bootstat);
    mean_std_error = std(bootstat);
    
    Mean_error(iIter) = mean_error;
    Mean_std_of_error(iIter) = mean_std_error;
    
    bootstat = bootstrp(500,'median',abs(err_cir));
    median_error = mean(bootstat);
    median_std_error = std(bootstat);
    
    Median_error(iIter) = median_error;
    Median_std_of_error(iIter) = median_std_error;
end

%% Plot inferred and true trajectories
inferred_pos = max_position/length(pvec)*2*pi*trackLen /2/pi;
figure; plot(position,'.-','LineWidth',2,'MarkerSize',10);hold on;
plot(inferred_pos,'.-','LineWidth',1,'MarkerSize',10);
xlabel('Time Steps')
ylabel('X [cm]')
legend('True trajectory','Inferred trajectory with filtered MPP')

%% Get the median and mean error for iterated decoding
Median = mean(Median_error,1)
MedianSTD = mean(Median_std_of_error,1)


