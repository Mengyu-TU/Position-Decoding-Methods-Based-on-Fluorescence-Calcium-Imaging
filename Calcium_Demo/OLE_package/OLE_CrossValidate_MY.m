function [probability,pvec,max_pos] = OLE_CrossValidate_MY(Data4Decoding,NumBasis,kappa,position2angle,nFold)
%% This function uses Optimal Linear Estimation to decode positions
%% Additional package required: OLE package

% Input: 
% Data4Decoding: Training and Testing Data for Decoding; Dimension: NumTSteps x NumNeurons
% NumBasis: number of equally-spaced von Mises functions
% kappa: vonmises base coefficient 
% angularPos: position converted to angle
% nFold: how many folds of cross validation

% Output:
% probability: 
% max_pos: decoded position

% Xrbf (basis function for position at each TStep (NumTStep x NumBasis)
% Call format: [b,Xrbf] = get1Dbasis(btype,NumBasis,X,kappa)
[basis,Xrbf] = get1Dbasis('vonmises',NumBasis,position2angle,kappa);

% Circular phase
pvec = linspace(0,2*pi,256);
[~,dbasis] = get1Dbasis('vonmises',basis.n,pvec,basis.s);
probability=[];

% Cross-validated predictions
if logical(nFold)
    m.cv = getCVidx(size(Data4Decoding,1),nFold);
    
    for i=1:m.cv.nfoldcv
        % least squared estimation
        w = Data4Decoding(m.cv.tr{i},:) \ Xrbf(m.cv.tr{i},:);  % training
        probability(m.cv.ts{i},:) = Data4Decoding(m.cv.ts{i},:) * (w*dbasis');  %#ok<*AGROW> % testing
    end
    
else % Train only; NO Testing
    w = Data4Decoding \ Xrbf;  % training
    probability = Data4Decoding * (w*dbasis');
end

[~,max_pos] = max(probability'); %#ok<*UDIM>

end