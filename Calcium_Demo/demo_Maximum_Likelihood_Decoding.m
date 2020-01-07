%% This demo script uses maximum likelihood decoding to infer the animal's position
%% Refer to "Hippocampal replay of extended experience" for Maximum Likelihood Estimator. (Davidson TJ1, Kloosterman F, Wilson MA. 2009)
clear
clc

addpath('data')
addpath(genpath('Maximum_Likelihood_functions')

%% Load Data
% Load filtered MPP for decoding
load('FilteredMPP4Decoding4MLE.mat')
Data4Decoding = SpikeMag';   % Transpose to Dimension: Number TimeSteps x Number Neurons

% Position
load('Position4DecodingMLE.mat')

%% Settings for maximum likelihood decoding
temporal_bin_size = 0.5; % second

ZSCORE = 1; % Z SCORE: 1; Don't Z SCORE: 0 ;
nFold = 10; % the number of folds of cross-validation to be done

%% Z Score Lap by Lap for Every timestep across neurons and Spatial
if ZSCORE
    %% Temporal Z score
    tempData = zeros(size(Data4Decoding));
    for iLap = 1: NumLap
        % Dimension Data4Decoding: NumTSteps x NumNeuron
        startIndex = NumTStepperLap*(iLap - 1) + 1;
        endIndex = NumTStepperLap*iLap;
        
        % Z score
        PartialData = Data4Decoding(startIndex:endIndex ,:);
        tempData(startIndex:endIndex ,:) = zscore(PartialData,0,1); % if dim = 1, then zscore along the columns of X,
    end
    %% SPATIAL Z score
    tempData2 = zscore(tempData,0,2);                               % if dim = 2, then zscore along the rows of X.
else % Not Z score
    tempData2 = Data4Decoding;
end

%% Resample to Poisson Distribution
PoissonMean = 5;
Data4Decoding = resample2poission_CZ(tempData2,PoissonMean,0);

______________________________________________________________________________________________________
%% Maximum likelihood decoding: cross-validated
%% Choose the midpt as the state
dx = 2.5; % cm
xgrid = linspace(min(position),max(position),ceil((max(position)-min(position))/dx)); % discretize the arena
posBin = xgrid(1:end-1) + mean(diff(xgrid))/2; % xgrid(2:end); %

%% Setting for construct place field
NORMALIZE = 0; % 0: don't normalize place field; 1: normalized place field
PLOT = 0; % 0: don't plot place field; 1: plot place field

%% Maximum likelihood: cross-validated
Z_test_whole = zeros(size(position));

m.cv = getCVidx(size(Data4Decoding,1),nFold);
for i=1:m.cv.nfoldcv
    TrainData = Data4Decoding(m.cv.tr{i},:);
    TestData = Data4Decoding(m.cv.ts{i},:);
    TrainPos = position(m.cv.tr{i});
    TestPos = position(m.cv.ts{i});
    
    %% Construct Place Field
    [~,lamda,sorted_firing_smoothed,xgrid] = PlaceFields_func_MY(TrainPos,TrainData,...
        dx,temporal_bin_size,NumLap,NORMALIZE,PLOT);  % lamda: PosBin x NumNeuron
    
    %% Maximum likelihood
    Z_test = maximum_likelihood_func_MY(TestData,posBin,lamda,temporal_bin_size);
    Z_test_whole(m.cv.ts{i}) = Z_test;
    
end

%% Median Error
Z_index = discretize(position,xgrid);
Z_true = zeros(1,length(Z_index));
for i  = 1:length(Z_index)
    Z_true(i) = (xgrid(Z_index(i) + 1) + xgrid(Z_index(i)))/2;
end
error = sqrt( (Z_test_whole' - Z_true).^2);
bootstat = bootstrp(500,'median',abs(error));
median_std_error = std(bootstat);
median_error = mean(bootstat)

%% Plot true and inferred trajectories
figure;plot(Z_true,'.-');hold on;plot(Z_test_whole,'*-')
legend('True','Inferred')
ylabel('X [cm]')
xlabel('Time Step')


