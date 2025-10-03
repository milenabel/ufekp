%% Interpolation in 1D/2D/3D domains using CSRBFs combined with
%%Jacobi polynomials. 


%% Spatial dimension
dim = 3;


%% Load up the node set
if dim==1
    N = 2.^(2:8); N = N';
    for k=1:length(N)
        X = chebspace2(-1,1,N(k));
        xi = X(2:end-1,:);
        xb = [X(1,:); X(end,:)];
        st.fullintnodes{k,1} = xi;
        st.bdrynodes{k,1} = xb;
    end
    clear xb xi k;
    xe = linspace(-1,1,2^14).';
elseif dim==2
    %% Get evaluation nodes
    st = load('DiskPoissonNodesLarge.mat');
    xe =  [st.fullintnodes{7}; st.bdrynodes{7}];
    st = load('DiskPoissonNodesClustered.mat');
elseif dim==3
     st = load('SpherePoissonNodesLarge.mat');
     xe = [st.fullintnodes{1}; st.bdrynodes{1}];
     st = load('SpherePoissonNodesClustered.mat');
end
start_nodes = 1;
end_nodes = size(st.fullintnodes,1);

% Sparsity storage initialization 
sparsity_fs1 = zeros(end_nodes, 1);

% Polynomial degree storage initialization
ell_poly = zeros(end_nodes, 1);          % For pure polynomial interpolation
ell_diag = zeros(end_nodes, 3);         % For diagonal approximation (3 smoothness levels)
ell_fs1 = zeros(end_nodes, 3);          % For fixed support method, K=1e2

% Support and shape parameters initialization
all_eps_fs = zeros(1, 1);  % 1 smoothness level × 1 K_target
supports_fs = zeros(1, 1);  

%% This partially determines the number of polynomial terms
%slightly smaller fac helps with increasing accuracy for both rough and
%smooth functions. There seems to be a notion of the "best" polynomial
%degree for a given number of nodes
if dim==1
    fac = 1.0;
elseif dim==2
    % fac = 0.8;
    fac = 1.0;
elseif dim==3
    fac = 1.0;
end

%% Setup anonymous function for function interpolation true solution
%%We'll use Trefethen's function gallery from his Spectral methods book
%%with 2d and 3d analogues
if dim==1
    syms x;       
    %f = abs(x);                function_name = 'abs_1d';
    %f = exp(-x.^(-2));     function_name = 'exp_1d';
    f = 1./(1 + 25*x.^2);  function_name = 'rk_1d';
    %f = x.^(10);           function_name = 'poly_1d';
    dfx = diff(f,x);
elseif dim==2
    syms x y;    
    %f = abs(x-0.5).^3* abs(y+0.4).^3;              %function_name = 'abs_2d';       
    %f = exp(-x.^(-2)).*exp(-y.^(-2));        function_name = 'exp_2d';
    %f = 1./(1 + 9*(x.^2 + y.^2));           function_name = 'rk_2d';
    %f = exp(-10*((x-.3).^2+y.^2));           function_name = 'exp10inv_2d';
    %f = x.^45 .* y.^35;                        function_name = 'poly_2d';
    %f = exp(x + y);                          function_name = 'exp_xy_2d';

    f = (x.^2 + y.^2).^(3/2);                function_name = 'xy_p_2d';
    %f = exp( ((x + y).^(2))/0.1 );           function_name = 'exp_p_2d';
    dfx = diff(f,x); dfy = diff(f,y);
elseif dim==3
    syms x y z;      
    %f = abs(x - 0.1).^3.*abs(y + 0.2).^3.*abs(z - 0.3).^3;              function_name = 'abs_3d';
    %f = exp(-10*((x-.3).^(-2)+y.^(-2) + z.^(-2)));    function_name = 'exp10inv_3d'
    %f = exp(-10*((x-.3).^2+y.^2 + z.^2));             function_name = 'exp10_3d';
    %f = exp( x + y + z);
    % f = 1./(1 + 9*(x.^2 + y.^2 + z.^2));             function_name = 'rk_3d';
    %f = x.^(4).*y.^(2).*z.^(2);                       function_name = 'poly_3d';

    %f = (x.^2 + y.^2 + z^2).^(3/2);                 function_name = 'xy_p_3d';
    f = exp( ((x + y + z).^(2))/0.8 );              function_name = 'exp_p_3d';
    dfx = diff(f,x); dfy = diff(f,y); dfz = diff(f,z);
end
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

%% Get the standard polynomial least squares stuff out of the way
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


%% Next, get the diagonal approximation stuff out of the way for different smoothnesses
smoothness=1;
%% Wendland C2 in 3d, pd in all lower dimensions
rbf = @(e,r) -r.^4.*(4.*r - 5*spones(r));
drbfor = @(e,r) 20.*e.^2.*r.^3;
if dim==1
    lrbf = @(e,r) -20.*e.^2.*r.^2.*(2.*r - 3*spones(r));
elseif dim==2
    lrbf = @(e,r) -20.*e.^2.*r.^2.*(r - 3*spones(r));
elseif dim==3
    lrbf = @(e,r) 60.*e.^2.*r.^2;
end

for k=start_nodes:end_nodes
    xi = st.fullintnodes{k};
    xb = st.bdrynodes{k};
    x = [xi;xb];
    alph = 0; %legendre
    ell = floor(fac*nthroot(length(x),dim));    
    ell = max([ell,1]); 
    ell_diag(k,smoothness) = ell;

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
    tree = KDTreeSearcher(x);
    [el2_diag(k,smoothness),elinf_diag(k,smoothness),a_time_diag(k,smoothness),e_time_diag(k,smoothness),c_poly_diag{k,smoothness}] = CSRBFDiag(x,y,ell,xe,alph,rbf,tree,ye_true);    
end

smoothness=1;
rbf = @(e,r) -r.^4.*(4.*r - 5*spones(r));
for k=start_nodes:end_nodes
    xi = st.fullintnodes{k};
    xb = st.bdrynodes{k};
    x = [xi;xb];
    alph = 0; %legendre
    ell = floor(fac*nthroot(length(x),dim));
    ell = max([ell,1]); 
    ell_fs1(k,smoothness) = ell;  % For K=1e2
    y = f(x(:,1),x(:,2),x(:,3));
    ye_true = f(xe(:,1),xe(:,2),xe(:,3));

    tree = KDTreeSearcher(x);

    ep1 = 5;
    [el2_fs1(k,smoothness), elinf_fs1(k,smoothness), a_time_fs1(k,smoothness), e_time_fs1(k,smoothness), c_poly_fs1{k,smoothness}, cond_fs1(k,smoothness), ~, sparsity_fs1(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep1,tree,ye_true);
end


%% Save results and clear
timestamp = datestr(datetime('now'), 'yyyyMMdd_HHmmss');

results_dir = fullfile('results/', sprintf('%s', function_name),'/high');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end
results_filename = fullfile(results_dir, sprintf('results_%s.mat', function_name));
save( results_filename, ...
    'el2_poly', 'elinf_poly', 'a_time_poly', 'e_time_poly', 'c_poly', 'ell_poly', ...
    'el2_diag', 'elinf_diag', 'a_time_diag', 'e_time_diag', 'c_poly_diag', 'ell_diag', ...
    'el2_fs1', 'elinf_fs1', 'a_time_fs1', 'e_time_fs1', 'c_poly_fs1', 'cond_fs1', 'sparsity_fs1', 'ell_fs1', ...
    'all_eps_fs', 'supports_fs', ...
    'sN', 'dim', 'function_name', 'timestamp');
