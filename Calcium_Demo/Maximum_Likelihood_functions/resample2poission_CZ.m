function resampled = resample2poission_CZ(data,lambda,seed)
    
    if nargin == 0
        load('NewData_run.mat');
        data = rmuaa; lambda= 10; seed = 0;
    end
    
    num = numel(data);
    
    %% remove zero value
    reshaped = reshape(data,[],1);
    lowthr = prctile(reshaped(reshaped>0),0.5);
    
    idx = find(reshaped>lowthr);
    reshaped(reshaped<=lowthr) = 0;
    
    %% random number
    rng(seed); 
    rndnum = poissrnd(lambda,length(idx),1);
    
    [~,s1] = sort(rndnum);
    [~,s2] = sort(reshaped(idx));
    
    resampled = zeros(num,1);
    resampled(idx(s2)) = rndnum(s1);
    
    resampled = reshape(resampled,size(data,1),size(data,2));
end