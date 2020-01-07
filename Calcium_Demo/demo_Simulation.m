%% This script demonstrates the procedures for simulating fluorescence calcium traces with plots for illustration
%% Package required before running the demo: CNMF_E-master (https://github.com/zhoupc/CNMF_E)
clear
clc

addpath(genpath('Simulation_functions'))

%% Set up
seed = 99;
rng(seed);
framerate = 20;       % Hz

%% Construct Positions
trackLen = 100;       % cm
NumLap = 20;
NumTStepperLap = 100;
TotalTStep = NumLap*NumTStepperLap;

PosPerLap = linspace(0,trackLen,NumTStepperLap);
temp = (repmat(PosPerLap,NumLap))';  % Input dimension for function: NumTSteps x NumNeurons
position = temp(:,1);

%% Construct Tuning Curves
NumNeuron = 50;
max_tuning_pos = linspace(0,trackLen,NumNeuron); % Peak firing positions for each neuron
std_tuning = 8;        % cm

firingrate_max = 5;   % Hz
firingrate_timeseries = zeros(NumNeuron,TotalTStep);

for iNeuron = 1:NumNeuron
    firingrate_timeseries(iNeuron,:) = firingrate_max*...
        exp(-0.5*( (position - max_tuning_pos(iNeuron))/std_tuning).^2);
end

%% Spikes
trueSpikes = zeros(NumNeuron,TotalTStep);
% dt = 1
for iTStep = 1:TotalTStep
    for iNeuron = 1:NumNeuron
        trueSpikes(iNeuron,iTStep) = (rand < firingrate_timeseries(iNeuron,iTStep)/framerate);
    end
end % Ground True Spikes

%% AR 2: Spikes to Calcium Traces
trueConcentration = double(trueSpikes);
AR_order = 2;
b = 0;
noise_std = 0.3;              % Variance of Gaussian Noise to the fluorescence

% Generate AR coeffiencients:
% Credit to Wei X et.al.,2019. A zero-inflated gamma model for post-deconvolved calcium imaging traces
temporal_bin = 1/framerate * 1000;   % milisecond
decay_cst = 400;                     % decay constant: milisecond
tau_r = 1;
e_ac = [-1/decay_cst, -1/tau_r];      % if AR(1) then tau_r=0
lambda_1 = exp(e_ac(1)*temporal_bin); % how much decay in one time bin
lambda_2 = exp(e_ac(2)*temporal_bin);
AR_coeff1 = lambda_1 + lambda_2;
AR_coeff2 = -1*lambda_1*lambda_2;

AR_coeff = [AR_coeff2;AR_coeff1;1];

% Calculate the calcium concentration based on AR 2
for t = (AR_order + 1): TotalTStep
    trueConcentration(:, t) = trueConcentration(:, (t - AR_order):t) * AR_coeff;   % Concentration
end

Fluorescence = b + trueConcentration + noise_std * randn(NumNeuron, TotalTStep);

%% Plot Fluorescence Calcium Traces: Figure 6B
randiNeuron = [50 42 30 10 1]; 
NumTraces = 5;
time_interval = (1:1000)/framerate;

figure;hold on; set(gcf,'Position',[235  235  845  376])
for iTrace = 1:NumTraces
    plot(time_interval,Fluorescence(randiNeuron(iTrace),1:1000) + 4*(iTrace-1)*ones(size(time_interval)),...
        '-','LineWidth',1.5)
end
title('Simulated Fluorescence Calcium Traces')
xlabel('Time [s]')

%% Plot One Neuron's Fluorescence Calcium Trace, True Spikes, Inferred Spikes
%% Peak of Calcium Trace Identified, Filtered MPP: Figure 6C
NeuronID =  5; % randi(NumNeuron)
y = Fluorescence(NeuronID,:);

MinProm = 0.00;
peak_min = max(y)*0.30; 
[pks,locs] = findpeaks(y,'MinPeakHeight',peak_min,'MinPeakProminence',MinProm);

[cUnit, sUnit, optionsUnit] = deconvolveCa(y, 'ar1','constrained','optimize_b', 'optimize_pars');
FluorsceMag = zeros(size(y));
FluorsceMag(locs) = pks;
h1 = [0.1, 0.2, 0.4, 0, 0];
h = h1/sum(h1);
x = FluorsceMag;
FilterMPP = conv(x,h,'same');

figure; set(gcf,'Position',[37,75,1196, 551]) 
subplot(5,1,1);
plot(y,'.-','LineWidth',0.8);
xlim([0,450])
set(gca,'FontSize',14,'FontWeight','bold');
set(gca,'XTick',[])
set(gca,'YTick',[])
legend('Raw Fluorescence Trace')

subplot(5,1,2);
h2 = stem(trueSpikes(NeuronID,:),'k','LineWidth',1);
set(h2, 'Marker', 'none');
ylim([0,2])
xlim([0,450])
set(gca,'FontSize',14,'FontWeight','bold');
set(gca,'XTick',[])
set(gca,'YTick',[])
legend('True Spikes')

subplot(5,1,3);
h = stem(sUnit,'r','LineWidth',1);
set(h, 'Marker', 'none');
xlim([0,450])
set(gca,'FontSize',14,'FontWeight','bold');
set(gca,'XTick',[])
set(gca,'YTick',[])
legend('Inferred Spikes from Spike Deconvolution')

subplot(5,1,4);
h3 = stem(locs,pks,'b','LineWidth',1);
set(h3, 'Marker', 'none');
xlim([0,450])
set(gca,'FontSize',14,'FontWeight','bold');
set(gca,'XTick',[])
set(gca,'YTick',[])
legend('Peak Identified')

subplot(5,1,5)
h4 = stem(FilterMPP,'cyan','LineWidth',1);
set(h4, 'Marker', 'none');
xlim([0,450])
set(gca,'FontSize',12,'FontWeight','bold');
set(gca,'XTick',[0:200:400])
set(gca,'YTick',[])
legend('Filtered MPP')
xlabel('Time Steps')

