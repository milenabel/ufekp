%% Test our method on manifolds. The polynomial matrix should be
%% rank-deficient on any algebraic Riemannian manifold.
%% We will tackle this with a rank-revealing QR factorization of
%% matrices that involve P.

%% Numbers of points on the torus for a convergence study
n = [5,10,15,20,25,30,35,40];

%% Number of eval points
ne = 50;
xe = computeHexagonTorusNodes(ne);

%% Target function
rb = 1/3;
f = @(x) exp( (x(:,1) + x(:,2) + x(:,3)));
%f = @(x) exp( (x(:,1) + x(:,2) + x(:,3)).^2);
%f = @(x) (acos(x(:,1))<rb).* (1 + cos(pi*acos(x(:,1))/rb))/2; %C^1

%% Evaluate this target at the eval points
f_xe = f(xe);
f_xe(isnan(f_xe)) = 0;

%% Define a Wendland RBF
rbf = @(e,r) -r.^4.*(4.*r - 5*spones(r)); %C2
%rbf = @(e,r) (r.^6.*(35.*r.^2 - 88.*r + 56*spones(r)))./3; %C4
%rbf = @(e,r) r.^8.*(66*spones(r) - 154*r + 121*r.^2 - 32*r.^3); %C6

%% Convergence study
for k=1:length(n)
    %% Generate the interpolation nodes
    x = computeHexagonTorusNodes(n(k));
    Ns(k,1) = length(x);

    %% Get a tree
    tree = KDTreeSearcher(x);

    %% Evaluate the target function there
    f_x = f(x);  
    f_x(isnan(f_x)) = 0;

    %% Get the polynomial degree. We'll use ambient polynomials in R^3
    ell = floor(0.8*nthroot(length(x),3));   

    %% Do a CSRBF + poly fit that uses a rank-revealing QR
   [el2(k,1),elinf(k,1),a_time(k,1),e_time(k,1),~,cond_num(k,1),~, sparsity(k,1)] = CSRBFGenManifold(x,f_x,ell,xe,0,rbf,5,tree,f_xe);   
end