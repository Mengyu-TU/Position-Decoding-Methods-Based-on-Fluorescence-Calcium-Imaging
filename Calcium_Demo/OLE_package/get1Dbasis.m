
function [b,Xrbf] = get1Dbasis(btype,n,X,s)

btype=lower(btype);
b.n=n;
b.type=btype;
if strcmp(btype,'gaussian')
    b.m = linspace(0,1,n);
    if nargin>3
        b.s=s;
    else
        b.s = mean(diff(b.m));
    end
    b.x0=linspace(0,1,256);
    b.rbf_basis = get1Dbasis_gaussian(b.x0,b.m,b.s);
    if nargin>2
        Xrbf = get1Dbasis_gaussian(X,b.m,b.s);
    end
elseif strcmp(btype,'vonmises')
    b.m = linspace(0,2*pi-2*pi/n,n);
    if nargin>3
        b.s=s;
    else
        b.s = 100;
    end
    b.x0=linspace(0,2*pi,256);
    b.rbf_basis = get1Dbasis_vonmises(b.x0,b.m,b.s);  
    if nargin>2
        Xrbf = get1Dbasis_vonmises(X,b.m,b.s);
    end
elseif strcmp(btype,'fourier')
    b.x0=linspace(0,2*pi,256);
    b.rbf_basis = get1Dbasis_fourier(b.x0,n);
    if nargin>2
        Xrbf = get1Dbasis_fourier(X,n);
    end
elseif strcmp(btype,'hist')
    b.x0=linspace(0,2*pi,256);
    b.rbf_basis = get1Dbasis_hist(b.x0,n);
    if nargin>2
        Xrbf = get1Dbasis_hist(X,n);
    end
end

% clf % clear current figure window
% plot(b.x0,b.rbf_basis);


function Xrbf = get1Dbasis_gaussian(X,rbfm,rbfs)

if nargin<3,
    rbfs = ones(1,rbfn)*mean(diff(rbfm));
end

Xrbf=[];
for i=1:length(rbfm)
    Xrbf(:,i) = exp(-(X-rbfm(i)).^2/rbfs.^2/2)/rbfs/sqrt(2*pi);
end

function Xrbf = get1Dbasis_vonmises(X,rbfm,rbfs)

Xrbf=[];
for i=1:length(rbfm)
    Xrbf(:,i) = exp(rbfs*cos(X-rbfm(i)));
    Xrbf(:,i) = Xrbf(:,i)/exp(rbfs);
end

function Xrbf = get1Dbasis_fourier(X,rbfn)

Xrbf=[];
for i=1:rbfn
    Xrbf(:,2*i-1) = sin(X*i);
    Xrbf(:,2*i)   = cos(X*i);
end

function Xrbf = get1Dbasis_hist(X,n)

Xrbf=[];
edges = linspace(0,2*pi,n+1);
for i=1:n
    Xrbf(:,i) = (X>edges(i)) & (X<edges(i+1));
end