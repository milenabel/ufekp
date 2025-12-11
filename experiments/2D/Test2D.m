%% Interpolation in 2D domains using CSRBFs combined with
%%Jacobi polynomials. 


%% Spatial dimension
dim = 2;

%% Load up the node set
%% Get evaluation nodes
st = load('DiskPoissonNodesLarge.mat');
xe =  [st.fullintnodes{7}; st.bdrynodes{7}];
st = load('DiskPoissonNodesClustered.mat');
clear p N pxe;

start_nodes = 1;
end_nodes = size(st.fullintnodes,1);

% Sparsity storage initialization 
sparsity_fs1 = zeros(end_nodes, 3);

% Polynomial degree storage initialization
ell_poly = zeros(end_nodes, 1);         % For pure polynomial interpolation
ell_diag = zeros(end_nodes, 3);         % For diagonal approximation (3 smoothness levels)
ell_fs1 = zeros(end_nodes, 1);          % For fixed support method


% Support and shape parameters initialization
all_eps_fs = zeros(1, 1);  % 3 smoothness levels × 3 K_targets
supports_fs = zeros(1, 1);  

%% This partially determines the number of polynomial terms
%slightly smaller fac helps with increasing accuracy for both rough and
%smooth functions. There seems to be a notion of the "best" polynomial
%degree for a given number of nodes
fac = 0.8;

%% Setup anonymous function for function interpolation true solution
%%We'll use Trefethen's function gallery from his Spectral methods book
%%with 2d and 3d analogues
syms x y;    
f = (x.^2 + y.^2).^(3/2);                function_name = 'xy_p_2d';
% f = exp( ((x + y).^(2))/0.2 );           function_name = 'exp_p_2d';
dfx = diff(f,x); dfy = diff(f,y);
f = matlabFunction(f);
if dim==1
    dfx = matlabFunction(dfx);
    clear x;
elseif dim==2
    dfx = matlabFunction(dfx);
    dfy = matlabFunction(dfy);
    clear x y;
elseif dim==3
    dfx = matlabFunction(dfx);
    dfy = matlabFunction(dfy);
    dfz = matlabFunction(dfz);
    clear x y z;
end
%% All possible tests:
%% 1. Different smoothness for CSRBF
%% 2. Different shape parameter strats
%% 3. Different interp techniques

%% Get the standard polynomial least squares results
for k=start_nodes:end_nodes
    xi = st.fullintnodes{k};
    xb = st.bdrynodes{k};
    x = [xi;xb];
    alph = 0; %legendre
    ell = floor(fac*nthroot(length(x),dim));     
    ell = max([ell,1]); 
    ell_poly(k,1) = ell;
    if dim==1
        y = f(x(:,1));
        ye_true = f(xe(:,1));
    elseif dim==2
        y = f(x(:,1),x(:,2));
        ye_true = f(xe(:,1),xe(:,2));
    elseif dim==3
        y = f(x(:,1),x(:,2),x(:,3));
        ye_true = f(xe(:,1),xe(:,2),xe(:,3));
    end    
    [el2_poly(k,1),elinf_poly(k,1),a_time_poly(k,1),e_time_poly(k,1),c_poly{k,1}] = PLS(x,y,ell,xe,alph,ye_true);   
    sN(k,1) = nthroot(length(x),dim);
end


%% Next, get the diagonal approximation results 
smoothness=1;
%% Wendland C2 in 3d, pd in all lower dimensions
rbf = @(e,r) -r.^4.*(4.*r - 5*spones(r));
drbfor = @(e,r) 20.*e.^2.*r.^3;
lrbf = @(e,r) -20.*e.^2.*r.^2.*(r - 3*spones(r));

for k=start_nodes:end_nodes
    xi = st.fullintnodes{k};
    xb = st.bdrynodes{k};
    x = [xi;xb];
    alph = 0; %legendre
    ell = floor(fac*nthroot(length(x),dim));     
    ell = max([ell,1]); 
    ell_diag(k,smoothness) = ell;

    y = f(x(:,1),x(:,2));
    ye_true = f(xe(:,1),xe(:,2));  
    tree = KDTreeSearcher(x);
    [el2_diag(k,smoothness),elinf_diag(k,smoothness),a_time_diag(k,smoothness),e_time_diag(k,smoothness),c_poly_diag{k,smoothness}] = CSRBFDiag(x,y,ell,xe,alph,rbf,tree,ye_true);    
end

%% Now find shape parameters that induce a target condition number on the finest node set
%% Use those on the coarser node sets, and it looks like the supports are increasing
smoothness=1 ;
%% Wendland C2 in 3d, pd in all lower dimensions
rbf = @(e,r) -r.^4.*(4.*r - 5*spones(r));
drbfor = @(e,r) 20.*e.^2.*r.^3;
lrbf = @(e,r) -20.*e.^2.*r.^2.*(r - 3*spones(r));

for k=start_nodes:end_nodes
    xi = st.fullintnodes{k};
    xb = st.bdrynodes{k};
    x = [xi;xb];
    alph = 0; %legendre
    ell = floor(fac*nthroot(length(x),dim));     
    ell = max([ell,1]); 
    ell_fs1(k,smoothness) = ell;  % For K=1e12
    y = f(x(:,1),x(:,2));
    ye_true = f(xe(:,1),xe(:,2));

    tree = KDTreeSearcher(x);

    % For K_target = 1e12 (j=1)
    ep1 = 10;
    [el2_fs1(k,smoothness), elinf_fs1(k,smoothness), a_time_fs1(k,smoothness), e_time_fs1(k,smoothness), c_poly_fs1{k,smoothness}, cond_fs1(k,smoothness), ~, sparsity_fs1(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep1,tree,ye_true);
end


%% Save results 
% Timestamp for uniqueness
timestamp = datestr(datetime('now'), 'yyyyMMdd_HHmmss');
results_dir = fullfile('results/', sprintf('%s', function_name),'/high');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end
results_filename = fullfile(results_dir, sprintf('results_%s.mat', function_name));


% Save everything
save(results_filename, ...
    'el2_poly', 'elinf_poly', 'a_time_poly', 'e_time_poly', 'c_poly', 'ell_poly', ...
    'el2_diag', 'elinf_diag', 'a_time_diag', 'e_time_diag', 'c_poly_diag', 'ell_diag', ...
    'el2_fs1', 'elinf_fs1', 'a_time_fs1', 'e_time_fs1', 'c_poly_fs1', 'cond_fs1', 'sparsity_fs1', 'ell_fs1', ...
    'sN', 'dim', 'function_name', 'timestamp');