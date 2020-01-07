% Get indices for cross-validation...

function cv = getCVidx(numel,nfoldcv,isShuff)

cv.n = numel;
cv.nfoldcv = nfoldcv;
% Evenly spaced...
break_points = floor(linspace(0,numel,nfoldcv+1));

if nargin<3 || ~isShuff
    ridx=1:numel;
else
    ridx=randperm(numel);
end
for i=1:nfoldcv
    cv.tr{i} = ridx([1:break_points(i) (break_points(i+1)+1):numel]);
    cv.ts{i} = ridx([(break_points(i)+1):break_points(i+1)]);
end