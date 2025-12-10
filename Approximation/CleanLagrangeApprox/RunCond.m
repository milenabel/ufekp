%% This file includes different conditional strategies: FS strategy, effective, and evaluation condition numbers for 2D.

%% Interpolation in 1D/2D/3D domains using CSRBFs combined with
%%Jacobi polynomials. 


%% Spatial dimension
dim = 2;

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
    % clear N xb xi k;
    clear xb xi k;
    xe = linspace(-1,1,2^14).';
elseif dim==2
    %% Get evaluation nodes
    st = load('DiskPoissonNodesLarge.mat');
    xe =  [st.fullintnodes{7}; st.bdrynodes{7}];

    %st = load('DiskPoissonNodes.mat');
    st = load('DiskPoissonNodesClustered.mat');
    %st = load('StarPoissonNodesClustered.mat');
    %st = load('DomainPoissonNodes2.mat');
%     pxe = haltonset(dim,'Skip',1e3,'Leap',1e2); %Quasi-random node set
%     pxe = scramble(pxe,'RR2');
%     N = (4:5:60).^2; N = N';
%     for k=1:length(N)
%         p = net(pxe,N(k));
%         p = p*2 - 1;
%         st.fullintnodes{k,1} = p;    
%     end
%     xe = net(pxe,10000);
%     xe = xe*2 - 1;
    clear p N pxe;
elseif dim==3
     st = load('SpherePoissonNodesLarge.mat');
     xe = [st.fullintnodes{2}; st.bdrynodes{2}];
% 
     st = load('SpherePoissonNodesClustered.mat');

%    st = load('RBCPoissonNodesClustered.mat');
%    xe = [st.fullintnodes{7}; st.bdrynodes{7}];
%    st = load('BumpySpherePoissonNodes.mat');
%    st = load('BumpySpherePoissonNodesClustered.mat');
%     pxe = haltonset(dim,'Skip',1e3,'Leap',1e2); %Quasi-random node set
%     pxe = scramble(pxe,'RR2');
%     N = (4:5:20).^3; N = N';
%     for k=1:length(N)
%         p = net(pxe,N(k));
%         p = p*2 - 1;
%         st.fullintnodes{k,1} = p;    
%     end
%     xe = net(pxe,10000);
%     xe = xe*2 - 1;
%     clear p N pxe;
end
start_nodes = 1;
end_nodes = size(st.fullintnodes,1);

% Sparsity storage initialization 
sparsity_fs1 = zeros(end_nodes, 3);
sparsity_vs1 = zeros(end_nodes, 3);

% Stability diagnostics (FS, smoothness=1 only)
kappa_eff_fs1 = nan(end_nodes, 1);
lam1_max_fs1 = nan(end_nodes, 1);
lam1_med_fs1 = nan(end_nodes, 1);
lam2_max_fs1 = nan(end_nodes, 1);
lam2_med_fs1 = nan(end_nodes, 1);

% Polynomial degree storage initialization
ell_poly = zeros(end_nodes, 1);         % For pure polynomial interpolation
ell_diag = zeros(end_nodes, 3);         % For diagonal approximation (3 smoothness levels)
ell_fs1 = zeros(end_nodes, 1);          % For fixed support method
ell_vs1 = zeros(end_nodes, 3);          % For variable support method, K=1e12


% Support and shape parameters initialization
all_eps_fs = zeros(1, 1);  
supports_fs = zeros(1, 1);  
all_eps_vs = zeros(1, 1);   
supports_vs = zeros(1, 1);

% Diagnostic-driven FS variants (choose on finest grid, reuse)
% FS-κeff
ell_fs_eff = zeros(end_nodes,1);
el2_fs_eff = nan(end_nodes,1);
elinf_fs_eff = nan(end_nodes,1);
a_time_fs_eff = nan(end_nodes,1);
e_time_fs_eff = nan(end_nodes,1);
cond_fs_eff = nan(end_nodes,1); % achieved κ_eff
eps_fs_eff_chosen = NaN; % chosen on finest, reused

% FS-evalcond (minimize median ||L||_2)
ell_fs_eval = zeros(end_nodes,1);
el2_fs_eval = nan(end_nodes,1);
elinf_fs_eval = nan(end_nodes,1);
a_time_fs_eval = nan(end_nodes,1);
e_time_fs_eval = nan(end_nodes,1);
evalcond_med_fs_eval = nan(end_nodes,1); % achieved median ||L||_2
eps_fs_eval_chosen = NaN; % chosen on finest, reused


%% This partially determines the number of polynomial terms
%slightly smaller fac helps with increasing accuracy for both rough and
%smooth functions. There seems to be a notion of the "best" polynomial
%degree for a given number of nodes
if dim==1
    fac = 1.0;
elseif dim==2
    % fac = 0.8;
    fac = 0.8;
elseif dim==3
    fac = 1.0;
end

%% Setup anonymous function for function interpolation true solution
%%We'll use Trefethen's function gallery from his Spectral methods book
%%with 2d and 3d analogues
if dim==1
    syms x;       
    %f = abs(x);            function_name = 'abs_1d';
    %f = exp(-x.^(-2));     function_name = 'exp_1d';
    f = 1./(1 + 25*x.^2);  function_name = 'rk_1d';
    %f = x.^(10);           function_name = 'poly_1d';
    dfx = diff(f,x);
elseif dim==2
    syms x y;    
    %f = abs(x-0.5).^3* abs(y+0.4).^3;       function_name = 'abs_2d';       
    %f = exp(-x.^(-2)).*exp(-y.^(-2));       function_name = 'exp_2d';
    %f = 1./(1 + 9*(x.^2 + y.^2));           function_name = 'rk_2d';
    %f = exp(-10*((x-.3).^2+y.^2));          function_name = 'exp10inv_2d';
    %f = x.^45 .* y.^35;                     function_name = 'poly_2d';
    %f = exp(x + y);                         function_name = 'exp_xy_2d';

    f = (x.^2 + y.^2).^(3/2);                function_name = 'xy_p_2d';
    % f = exp( ((x + y).^(2))/0.2 );         function_name = 'exp_p_2d';
    dfx = diff(f,x); dfy = diff(f,y);
elseif dim==3
    syms x y z;      
    %f = abs(x - 0.1).^3.*abs(y + 0.2).^3.*abs(z - 0.3).^3; function_name = 'abs_3d';
    %f = exp(-10*((x-.3).^(-2)+y.^(-2) + z.^(-2)));    function_name = 'exp10inv_3d'
    %f = exp(-10*((x-.3).^2+y.^2 + z.^2));             function_name = 'exp10_3d';
    %f = exp( x + y + z);
    % f = 1./(1 + 9*(x.^2 + y.^2 + z.^2));             function_name = 'rk_3d';
    %f = x.^(4).*y.^(2).*z.^(2);                       function_name = 'poly_3d';

    f = (x.^2 + y.^2 + z^2).^(3/2);                   function_name = 'xy_p_3d';
    %f = exp( ((x + y + z).^(2))/0.1 );                function_name = 'exp_p_3d';
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


%% Now find shape parameters that induce a target condition number (FS) on the finest node set
%% Use those on the coarser node sets (fixed shape). Add two diagnostic-driven FS variants.
%% One choice of smoothnesses
smoothness=1 ;
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

% One-time chosen epsilons for the two diagnostic-driven variants (decided on finest grid)
eps_fs_eff_chosen  = NaN;   % FS-keff (target kappa_eff = 1e12)
eps_fs_eval_chosen = NaN;   % FS-eval (min median ||L||_2)

% DESCEND so we hit the finest level first (to choose epsilons there)
for k=end_nodes:-1:start_nodes
    xi = st.fullintnodes{k};
    xb = st.bdrynodes{k};
    x = [xi;xb];
    alph = 0; %legendre
    ell = floor(fac*nthroot(length(x),dim));     
    ell = max([ell,1]); 
    ell_fs1(k,smoothness)  = ell; % For K=1e12
    ell_fs_eff(k, 1) = ell; % FS-keff series
    ell_fs_eval(k, 1) = ell; % FS-eval series
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

    % Precompute polynomial basis on scaled nodes once per level (shared)
    xmax = max(x,[],1); 
    xmin = min(x,[],1);
    c1 = 0.5*(xmax + xmin); 
    c2 = 2./(xmax - xmin);
    recurrence = @(N) jacobi_recurrence(N, alph, alph);
    idxs = total_degree_indices(size(x,2), ell);
    P_scaled = mpoly_eval((x - c1).*c2, idxs, recurrence); % n×m

    % FS baseline series: for K_target = 1e12 (j=1)
    ep1 = 10;
    [el2_fs1(k,smoothness), elinf_fs1(k,smoothness), a_time_fs1(k,smoothness), e_time_fs1(k,smoothness), c_poly_fs1{k,smoothness}, cond_fs1(k,smoothness), ~, sparsity_fs1(k,smoothness)] = ...
        CSRBFGen(x,y,ell,xe,alph,rbf,ep1,tree,ye_true);

    % (A) FS-keff: choose eps on finest grid so kappa_eff ≈ 1e12, reuse elsewhere
    if isnan(eps_fs_eff_chosen)
        target_keff = 1e12; % match the "high" setting
        eps_bracket = [3, 15];
        eps_fs_eff_chosen = find_eps_for_kappa_eff(x, tree, rbf, P_scaled, target_keff, eps_bracket);
        fprintf('Chosen eps_fs_eff (kappa_eff≈1e12) on finest grid: %.4g\n', eps_fs_eff_chosen);
    end
    ep_eff = eps_fs_eff_chosen;
    [el2_fs_eff(k,1), elinf_fs_eff(k,1), a_time_fs_eff(k,1), e_time_fs_eff(k,1), ~] = CSRBFGen(x, y, ell, xe, alph, rbf, ep_eff, tree, ye_true);

    % Record achieved κ_eff for reporting (and clear temporaries)
    K_eff = rbf(ep_eff, DistanceMatrixCSRBFwt(x, x, ep_eff, tree));
    cond_fs_eff(k,1) = rbf_effective_cond(K_eff, P_scaled);
    clear K_eff

    % (B) FS-eval: choose eps on finest grid to minimize median ||L||_2, reuse elsewhere
    if isnan(eps_fs_eval_chosen)
        eps_bounds = [3, 15];
        eps_fs_eval_chosen = find_eps_min_evalcond(x, xe, tree, rbf, P_scaled, c1, c2, idxs, eps_bounds, recurrence);
        fprintf('Chosen eps_fs_eval (min median ||L||_2) on finest grid: %.4g\n', eps_fs_eval_chosen);
    end
    ep_eval = eps_fs_eval_chosen;
    [el2_fs_eval(k,1), elinf_fs_eval(k,1), a_time_fs_eval(k,1), e_time_fs_eval(k,1), ~] = CSRBFGen(x, y, ell, xe, alph, rbf, ep_eval, tree, ye_true);

    % store achieved median ||L||_2 for reporting (diagnostics only)
    K_eval = rbf(ep_eval, DistanceMatrixCSRBFwt(x, x, ep_eval, tree));
    Aaug_eval = [K_eval, P_scaled; P_scaled', sparse(size(P_scaled,2), size(P_scaled,2))];
    Dqr_eval = decomposition(Aaug_eval,'qr');
    solveTt_eval = @(rhs) Dqr_eval \ rhs;
    Kx_fun_eval = @(xstar) rbf(ep_eval, DistanceMatrixCSRBFwt(xstar, x, ep_eval, tree)).';
    Px_fun = @(xstar) mpoly_eval((xstar - c1).*c2, idxs, recurrence);

    qstep_eval = max(1, floor(size(xe,1)/2000)); % smaller subsample to reduce memory
    Xsample_eval = xe(1:qstep_eval:end, :);

    [~, ~, ~, med2_tmp] = eval_condition_numbers(solveTt_eval, Kx_fun_eval, Px_fun, Xsample_eval);
    evalcond_med_fs_eval(k,1) = med2_tmp;

    clear K_eval Aaug_eval Dqr_eval solveTt_eval Kx_fun_eval Xsample_eval

    % Diagnostics for the original FS baseline (ep1), placed last to lower peak memory
    K_ep1 = rbf(ep1, DistanceMatrixCSRBFwt(x, x, ep1, tree));
    kappa_eff_fs1(k,1) = rbf_effective_cond(K_ep1, P_scaled);

    Aaug_ep1 = [K_ep1, P_scaled; P_scaled', sparse(size(P_scaled,2), size(P_scaled,2))];
    Dqr_ep1 = decomposition(Aaug_ep1, 'qr');
    solveTt_ep1 = @(rhs) Dqr_ep1 \ rhs;
    Kx_fun_ep1 = @(xstar) rbf(ep1, DistanceMatrixCSRBFwt(xstar, x, ep1, tree)).';
    Px_fun = @(xstar) mpoly_eval((xstar - c1).*c2, idxs, recurrence);

    qstep = max(1, floor(size(xe,1)/2000)); % smaller subsample for further memory relief
    Xsample = xe(1:qstep:end, :);

    [lam1_max_fs1(k,1), lam1_med_fs1(k,1), lam2_max_fs1(k,1), lam2_med_fs1(k,1)] = eval_condition_numbers(solveTt_ep1, Kx_fun_ep1, Px_fun, Xsample);

    clear K_ep1 Aaug_ep1 Dqr_ep1 solveTt_ep1 Kx_fun_ep1 Xsample Px_fun P_scaled c1 c2 idxs recurrence
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
    'kappa_eff_fs1','lam1_max_fs1','lam1_med_fs1','lam2_max_fs1','lam2_med_fs1', ...
    'el2_fs_eff','elinf_fs_eff','a_time_fs_eff','e_time_fs_eff','cond_fs_eff','ell_fs_eff','eps_fs_eff_chosen', ...
    'el2_fs_eval','elinf_fs_eval','a_time_fs_eval','e_time_fs_eval','evalcond_med_fs_eval','ell_fs_eval','eps_fs_eval_chosen', ...
    'sN', 'dim', 'function_name', 'timestamp');

function eps_opt = find_eps_for_kappa_eff(x, tree, rbf, P_scaled, target_keff, eps_bracket)
% Root-find in log-space for kappa_eff ≈ target_keff.
    target_log = log10(target_keff);

    function v = mismatch(ep)
        K = rbf(ep, DistanceMatrixCSRBFwt(x, x, ep, tree));
        ke = rbf_effective_cond(K, P_scaled);
        v = log10(ke) - target_log;
    end

    a = eps_bracket(1); 
    b = eps_bracket(2);
    grid = linspace(a,b,8);
    vals = arrayfun(@(e) mismatch(e), grid);
    id = find(sign(vals(1:end-1)).*sign(vals(2:end))<=0, 1);
    if ~isempty(id)
        eps_opt = fzero(@mismatch, [grid(id), grid(id+1)]);
    else
        % fallback: minimize absolute mismatch
        eps_opt = fminbnd(@(e) abs(mismatch(e)), a, b);
    end
end

function eps_opt = find_eps_min_evalcond(x, xe, tree, rbf, P_scaled, c1, c2, idxs, eps_bounds, recurrence)
% Minimize the median ||L(x*)||_2 over a subsample of Xe.
    qstep = max(1, floor(size(xe,1)/2000));
    Xsample = xe(1:qstep:end, :);

    function m = median_L2(ep)
        K = rbf(ep, DistanceMatrixCSRBFwt(x, x, ep, tree));
        Aaug = [K, P_scaled; P_scaled', sparse(size(P_scaled,2), size(P_scaled,2))];
        Dqr  = decomposition(Aaug,'qr');
        solveTt = @(rhs) Dqr \ rhs;
        Kx_fun = @(xstar) rbf(ep, DistanceMatrixCSRBFwt(xstar, x, ep, tree)).';
        Px_fun = @(xstar) mpoly_eval((xstar - c1).*c2, idxs, recurrence);
        [~, ~, ~, med2] = eval_condition_numbers(solveTt, Kx_fun, Px_fun, Xsample);
        m = med2;
    end

    a = eps_bounds(1); 
    b = eps_bounds(2);
    try
        eps_opt = fminbnd(@median_L2, a, b);
    catch
        grid = linspace(a,b,12);
        vals = arrayfun(@median_L2, grid);
        [~,ii] = min(vals);
        eps_opt = grid(ii);
    end
end