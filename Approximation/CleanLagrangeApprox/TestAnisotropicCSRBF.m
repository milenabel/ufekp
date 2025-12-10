%% We'll test anisotropic interpolation with CSRBFs

%% Load up a point set
load('DiskPoissonNodesClustered.mat');
xi = fullintnodes{1};
xb = bdrynodes{1};
x = [xi;xb];
y = bump_cinf(x(:,1), x(:,2), 0.3, 0.1, 0.5,5);
xie = fullintnodes{6};
xbe = bdrynodes{6};
xe = [xie;xbe];
ye = bump_cinf(xe(:,1),xe(:,2), 0.3, 0.1, 0.5,5);

%% Do standard CSRBF + poly to fit this
tree = KDTreeSearcher(x);
ell = floor(0.8*nthroot(size(x,1),size(x,2)));
alph = 0;%Legendre
%rbf = @(e,r) -r.^4.*(4.*r - 5*spones(r));
rbf = @(e,r) (r.^6.*(35.*r.^2 - 88.*r + 56*spones(r)))./3; %C4
%rbf = @(e,r) r.^8.*(66*spones(r) - 154*r + 121*r.^2 - 32*r.^3); %C6

ep = 2;
[el2,elinf,a_time,e_time,c_poly,cond_num,c_rbf, sparsity] = CSRBFGen(x,y,ell,xe,alph,rbf,ep,tree,ye);

%% Now do a version where the RBF sees warped points but the polynomial sees the standard ones
ep2 = 2;
[el2_2,elinf_2,a_time_2,e_time_2,c_poly_2,cond_num_2,c_rbf_2,sparsity_2] = CSRBFGenWarp(x,y,ell,xe,alph,rbf,ep,ye);