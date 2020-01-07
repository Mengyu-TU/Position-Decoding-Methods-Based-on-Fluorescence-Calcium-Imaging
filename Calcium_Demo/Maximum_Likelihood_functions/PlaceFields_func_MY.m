function [sorted_rank,firing_smoothed,firing_sorted_smoothed,xgrid] = PlaceFields_func_MY(position,spike,...
    dx,temporal_bin_size,NumLap,NORMALIZE,PLOT)
% spike: Dimension: NumTStep x NumNeuron
NumNeuron = size(spike,2);
NumTSteps =  size(spike,1);

%% Prepare to plot place fields
xgrid = linspace(min(position),max(position),ceil((max(position)-min(position))/dx)); % discretize the arena

% Firing rate matrix
firing_rate_time_series = spike/temporal_bin_size;

firingrate = zeros(numel(xgrid)-1,NumNeuron);
numvisits = NumLap * ones(numel(xgrid)-1,1);   % number of times agent visits bin

for iunit=1:NumNeuron
    
    for t=1:NumTSteps
        % Find the firing rate for ith unit at the position evaluated at
        if position(t)>=min(xgrid) && position(t)<=max(xgrid)
            ixposition = find(position(t) < xgrid ,1,'first') - 1;% if xposition < min(xgrid) ixposition = 0, if xposition > max(xgrid) ixposition is an empty matrix
            
%             numvisits(ixposition) = numvisits(ixposition) + 1;
            firingrate(ixposition,iunit) = firingrate(ixposition,iunit)...
                + firing_rate_time_series(t,iunit);
        end% time step 't'
        
    end
    
end  % for iunit = 1:numh

%% Smooth the 1D place fields
firing_smoothed = zeros(size(firingrate));
for iunit = 1:NumNeuron
    firing_smoothed(:,iunit) = smooth_ratemap_Colgin(firingrate(:,iunit)./numvisits);
end

%% Sort firing rate matrix according to peak position of each neuron and save the ranking
pos_rank = zeros(1,NumNeuron);
for q = 1:NumNeuron
    temp = find(firing_smoothed(:,q)== max(firing_smoothed(:,q)),1,'first'); %max by default: 'omitnan'
    pos_rank(q) = temp(1);
end

rank = [pos_rank;1:NumNeuron];
temp = (sortrows(rank'))';
sorted_rank = [1:NumNeuron;temp(2,:)]; % Row 2: Sorted peak firing position: from small X to big

%% Plot normalized smoothed 1D place fields
sorted_firing_mat_smoothed = zeros(size(firing_smoothed));
for iUnit = 1: NumNeuron
    sorted_firing_mat_smoothed(:,iUnit) = firing_smoothed(:,sorted_rank(2,iUnit));
end

firing_sorted_smoothed = zeros(size(sorted_firing_mat_smoothed));
for i = 1:NumNeuron
    if NORMALIZE
        firing_sorted_smoothed(:,i) = sorted_firing_mat_smoothed(:,i)/max(sorted_firing_mat_smoothed(:,i));
    else
        firing_sorted_smoothed(:,i) = sorted_firing_mat_smoothed(:,i);
    end
end

if PLOT
    imagesc(xgrid(2:end),1:NumNeuron,firing_sorted_smoothed');
    colormap(flipud(gray));colorbar;
end

end