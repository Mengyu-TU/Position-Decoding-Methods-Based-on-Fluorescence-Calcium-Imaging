

% Reference: Spatially distributed local fileds in the hippocampus encode
% rat position, Science (2014) 

% Supplementary Information: Methods

% regress LFP features to spatial position 

clear;close all;clc
addpath(genpath('circstat'))
load([sprintf('d001_b%04i_lite.mat',100)])





% absAct       38118 x 63 (absolute amplitude) 
% Act_demod    38118 x 63 (complex value)
% pos          38118 x 2
% rpos         38118 x 1
% spk          38118 x 85
% v            38118 x 2

% p       1 x 21  ???
% cT      1 x 85  ???(value: 1/2/3)

trackLen = 250; % cm

% velocity threshold to choose RUN episodes

% Use only data where speed>5% max
vel = [0; diff(pos(:,1))];
spd = smooth(sqrt(vel.^2));  % moving average with span 5


sidx = find(spd>(max(spd)*0.15));  % 19203
%sidx = find(spd>0.005);  % 10795 (5 cm/s) 

pos = pos(sidx,:);    % 19203 x 1
rpos = rpos(sidx,:);  % 19203 x 1, range [0, 2]

% Observed Signals to use for decoding
yAll=cell(0);
yAll{1} = zscore(spk);
yAll{2} = zscore(absAct);
yAll{3} = zscore([real(Act_demod) imag(Act_demod)]);    
yAll{4} = zscore([spk real(Act_demod) imag(Act_demod)]);
mnames={'Spikes','abs(A)','Real(A_d)+Imag(A_d)','Joint'};



% Decoding: OLE
for model_num=1:length(yAll)
    y = yAll{model_num}(sidx,:);
    
    % Generate basis...
    figure(3)
    [basis,Xrbf] = get1Dbasis('vonmises',75,rpos*pi,400);
    % basis.n: 75   
    % basis.m: 1 x 75
    % basis.s: 400
    % basis.rbf_basis: 256 x 75
    % Xrbf:  19203 x 75
    
    
    pvec = linspace(0,2*pi,256); % circular phase 
    [tmp,dbasis] = get1Dbasis('vonmises',basis.n,pvec,basis.s);
    
    xlabel('Position')
    title('Basis')
    set(gca,'fontsize',20);
    
    % Cross-validated predictions
    m.cv = getCVidx(size(y,1),10);
    f=[];
    for i=1:m.cv.nfoldcv
        % least squared estimation 
        w = y(m.cv.tr{i},:) \ Xrbf(m.cv.tr{i},:);  % training
        f(m.cv.ts{i},:) = y(m.cv.ts{i},:) * (w*dbasis');  % testing
    end
    
    % Plotting and error calcs...
    xl=[4000 8000]; % range to show predictions
    f = zscore(f')';
    decoderSummary  % figure illustration 
  
    % Save results
    errAll(:,model_num) = err;
    errpAll(:,model_num) = errp;
    errmAll(:,model_num) = errm;
    % non-directional error
    [err_nd_All(:,model_num),mini] = (min([abs(errmAll(:,model_num)) abs(errpAll(:,model_num))]'));
    pos_hat(:,model_num) = maxpost'/length(pvec)*2*pi;
end

%% Median error in cm with s.e.
for model_num=1:length(yAll)
    bootstat = bootstrp(500,'median',abs(err_nd_All(:,model_num)));
    [mean(bootstat) std(bootstat)]
end

%% Hinton diagram confusion matrices

figure 
edges = linspace(0,2*pi,64);
e{1}=edges;
e{2}=edges;
ex = edges+mean(diff(edges))/2;
[mx,my] = meshgrid(ex,ex);

pvec =[1 3 4]; 
for i=1:length(pvec)
    subplot(1,3,i)
    [n,c]=hist3([rpos*pi pos_hat(:,pvec(i))],'Edges',e);
    nnorm = bsxfun(@rdivide,n,sum(n));
    idx = nnorm>10e-3;
    scatter(mx(idx)/pi*trackLen-trackLen,my(idx)/pi*trackLen-trackLen,nnorm(idx)*100,'filled','sk')
    axis image
    xlim([-1 1]*trackLen)
    ylim([-1 1]*trackLen)
    set(gca,'TickDir','out')
    set(gca,'XTick',[-trackLen 0 trackLen])
    set(gca,'YTick',[-trackLen 0 trackLen])
    title(mnames{pvec(i)})
end