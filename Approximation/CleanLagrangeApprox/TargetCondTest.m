%% Specify a target function
syms x y;
f = 1./(1 + 4*(x.^2 + y.^2));  
f = matlabFunction(f); clear x y;

%% Load up the interp node set
st = load('DiskPoissonNodesClustered.mat');
xi = st.fullintnodes{4};
xb = st.bdrynodes{4};
x = [xi;xb];

%% Load up the evaluation node set
st = load('DiskPoissonNodesLarge.mat');
xie = st.fullintnodes{7};
xbe = st.bdrynodes{7};
xe = [xie;xbe];

%% Specify a kernel
rbf = @(e,r) (r.^6.*(35.*r.^2 - 88.*r + 56*spones(r)))./3;

%% Set a target condition number and an initial shape parameter
kt = 1e8;

%% Specify the function to root find on: log ( cond(rbf(ep,r))/target) = 0
%% Also create a kd-tree
tree = KDTreeSearcher(x);
options.TolX = 1e-4;
ep_func = @(ep,x,xc,tc) log( condest(rbf(ep,DistanceMatrixCSRBFwt(x,xc,ep,tree))) / tc );

%% Find the shape parameter
guess = 2; %a guess for epsilon
ep= fzero(ep_func,guess,options,x,x,kt);

%% Solve for the interpolation coefficients
r = DistanceMatrixCSRBFwt(x,x,ep,tree);
A = rbf(ep,r);
f_x = f(x(:,1),x(:,2));
c = A\f_x;

%% Evaluate the interpolant somewhere
re = DistanceMatrixCSRBFwt(xe,x,ep,tree);
Ae = rbf(ep,re);
s_e = Ae*c;

%% Measure errors
f_xe = f(xe(:,1),xe(:,2));
norm(s_e - f_xe)/norm(f_xe)