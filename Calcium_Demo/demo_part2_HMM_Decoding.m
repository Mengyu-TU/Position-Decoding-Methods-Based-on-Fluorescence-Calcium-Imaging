%% This script continues from the Python "demo_part1_HMM_Decoding.py"
%% It maps the inferred states to positions and plots the inferred and true positions
%% and calculates the median error with bootstrapping
clear
clc

addpath('data')

%% Load decoding outcome
DecodingDate = '12162019';
DecodingIndex = '01';  
data = sprintf('CaIm_exp_inf_%s%s.mat',DecodingDate,DecodingIndex);

% Variables: 
% Z_inf_mat: Inferred states; Dimension: Iteration x Number of Training Time Steps
% PosteriorProbTest: Posterior Probablity for every testing time steps;
% Dimension:
load(data,'Z_inf_mat','PosteriorProbTest','N_used_vec');

NTStepsTrain = length(Z_inf_mat(1,:));
NTStepsTest =  length(PosteriorProbTest(:,1,1));

%% Load true position
load('Position4DecodingHMM.mat','position','trackLen')

if NTStepsTrain + NTStepsTest == length(position)
    position_train = position(1:NTStepsTrain);
    position_test = position(NTStepsTrain+1:end);
else
    error('Wrong position loaded!/HMM Run Wrong')
end

%% Preallocate variables for storing errors and inferred positions
NO_iteration = length(Z_inf_mat(:,1));
decoding_error = zeros(NO_iteration,NTStepsTrain);
infer_xposition_train = zeros(NO_iteration,NTStepsTrain);
infer_xposition_test = zeros(NO_iteration,NTStepsTest);
error_train = zeros(NO_iteration,NTStepsTrain);
error_test = zeros(NO_iteration,NTStepsTest);

%% Add 1 to inferred states because inferred state number in Python starts from 0
if min(min(Z_inf_mat)) == 0
    Z_inf_mat_matlab = Z_inf_mat(:,:,:) + 1;    
else
    Z_inf_mat_matlab = Z_inf_mat(:,:,:);
end

%% Map states to positions and calculate errors
for iteration = 1 : NO_iteration
    %% Map inferred latent states to positions 
    M = max(N_used_vec(iteration)); % maximal inferred state no.
    
    Z_inf = Z_inf_mat_matlab(iteration,:);
    MAP_inf_state_pos = zeros(1,M);
    
    for i = 1:M
        index = find(Z_inf == i);
        x_pos_corr = position_train(index);
        
        median_x_pos = median(x_pos_corr);  
        MAP_inf_state_pos(:,i) = median_x_pos;
        
    end % End of finding correspondence between inferred states and positions
    
    %% Find the inferred positions for the testing data
    for h = 1:NTStepsTest
        infer_xposition_test(iteration,h) = PosteriorProbTest(h,1:M,iteration)/sum(PosteriorProbTest(h,1:M,iteration))...
            *(MAP_inf_state_pos)';
    end
    
    %% Error at Every Time Step of Testing Data    
    error = infer_xposition_test(iteration,:) - position_test;
    error_temp = error;
    % Track is closed, largest error can't exceed half of the track length
    error(abs(error_temp) > 0.5 *trackLen) = trackLen - abs(error_temp(abs(error_temp) > 0.5 *trackLen));
    error_test(iteration,:) = (error.^2 ).^0.5;
    
end

%% Bootstrap to obtain the median errors
IterBootstraps = 20;
Median = zeros(IterBootstraps,1);
Median_STD = zeros(IterBootstraps,1);

for iIter = 1:IterBootstraps
    bootstat = bootstrp(500,'median',error_test(:));
    median_error = mean(bootstat);
    std_error = std(bootstat);
    Median(iIter) = median_error;
    Median_STD(iIter) = std_error;
end

FinalMedian = mean(Median)
FinalMedian_STD = mean(Median_STD)

%% Plot the inferred positions and true positions of testing data
figure;plot(infer_xposition_test(1,:),'*-','MarkerSize',15,'LineWidth',1);
hold on;plot(position_test,'.-','MarkerSize',20,'LineWidth',2);legend('Inferred Trajectory','True Trajectory')
xlabel('Time Steps','FontWeight','bold')
ylabel('X [cm]','FontWeight','bold')
set(gca,'FontWeight','bold','FontSize',15)
set(gcf, 'Position',  [300, 300, 800, 350])
