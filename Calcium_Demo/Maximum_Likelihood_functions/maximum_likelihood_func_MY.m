function Z_test = maximum_likelihood_func_MY(TestData,posBin,lamda,temporal_bin_size)
 
num_TStep = length(TestData(:,1));
num_posBin = length(posBin);

logll = zeros(num_TStep,num_posBin); % Size: TimeSteps x PostionBin
prob = zeros(num_TStep,num_posBin);  % Size: TimeSteps x PositionBin
Z_test = zeros(num_TStep,1);         % Size: TimeSteps x 1
% probAtT= zeros(num_TStep,1);         % Size: TimeSteps x 1

% lamda: Size PositionBin x Neuron; Add small number to avoid ln0
Lam = lamda + (10^-8) * ones(size(lamda));

for TStep = 1: num_TStep
    TestDataAtTStep = TestData(TStep,:); % 1 x NumNeuron
    
    %% Compute the log likelihood of all position bins for the current time step
    logll(TStep,:) = TestDataAtTStep * log(Lam' * temporal_bin_size) - (temporal_bin_size*sum(Lam,2))'- ...
        sum(log(factorial(TestDataAtTStep')));
    
    %% Posterior probability for all position bins for the current time step
    prob(TStep,:) = exp(logll(TStep,:)) / sum(exp(logll(TStep,:)));
%     prob(TStep,:) = exp(logll(TStep,1:end-1)) / sum(exp(logll(TStep,1:end-1)));

    %% Choose the position bin with the largest probability as the position for the current time step
    [~,Index] = max(prob(TStep,:));
    %     Index = find(prob(TStep,:) == max(prob(TStep,:)),1,'first');
    Z_test(TStep) = posBin(Index);
    %     probAtT(TStep) = prob(TStep,Index) / sum(prob(TStep,:));
end

end