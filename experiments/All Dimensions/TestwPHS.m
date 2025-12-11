% TestwPHS.m
% Compute global PHS + polynomial results on the same setups 
% and save/append them into results_<function_name>.mat

clear; clc;

%% Spatial dimension
dim = 1; % 1, 2, or 3
phsdeg = 5; % manually change the curretn PHS degree to record results correctly 
assert(ismember(phsdeg,[5 7 9]), 'phsdeg must be one of {5,7,9}.'); % (current picks: 5, 7, 9)
m_fixed = (phsdeg - 1)/2; % m = 2, 3, or 4
assert(m_fixed == round(m_fixed), 'phsdeg must be odd so (phsdeg-1)/2 is integer.');
m_fixed = round(m_fixed);

fac = (dim==2) * 1.0 + (dim~=2) * 1.0; % keep consistent with other runs
alph = 0; % Legendre (Jacobi alpha=0), matches PLS

%% Load the same node sets
if dim==1
    N = 2.^(2:8); N = N';
    for k=1:length(N)
        X = chebspace2(-1,1,N(k));
        xi = X(2:end-1,:); xb = [X(1,:); X(end,:)];
        st.fullintnodes{k,1} = xi;
        st.bdrynodes{k,1}    = xb;
    end
    clear N xb xi k;
    xe = linspace(-1,1,2^14).';
elseif dim==2
    % Evaluation nodes (finest)
    st_e = load('DiskPoissonNodesLarge.mat');
    xe = [st_e.fullintnodes{7}; st_e.bdrynodes{7}]; 
    % Training nodes
    st = load('DiskPoissonNodesClustered.mat');
    clear p N pxe;
elseif dim==3
    st_e = load('SpherePoissonNodesLarge.mat');
    xe = [st_e.fullintnodes{2}; st_e.bdrynodes{2}];
    st = load('SpherePoissonNodesClustered.mat');
end

start_nodes = 1;
end_nodes = size(st.fullintnodes,1);
if dim == 3
    end_nodes = min(end_nodes, 5);
end

% Choose the same target function you used in the corresponding run 
if dim==1
    syms x;       
    f = abs(x);            function_name = 'abs_1d';
    % f = 1./(1 + 25*x.^2);  function_name = 'rk_1d';
    dfx = diff(f,x);
elseif dim==2
    syms x y;    
    f = (x.^2 + y.^2).^(3/2);                function_name = 'xy_p_2d';
    % f = exp( ((x + y).^(2))/0.2 );         function_name = 'exp_p_2d';
    dfx = diff(f,x); dfy = diff(f,y);
elseif dim==3
    syms x y z;      
    f = (x.^2 + y.^2 + z^2).^(3/2);                   function_name = 'xy_p_3d';
    % f = exp( ((x + y + z).^(2))/0.1 );                function_name = 'exp_p_3d';
    dfx = diff(f,x); dfy = diff(f,y); dfz = diff(f,z);
end
f = matlabFunction(f); % true function for errors

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

%% Allocate outputs for PHS + poly 
el2_phs = nan(end_nodes,1);
elinf_phs = nan(end_nodes,1);
a_time_phs = nan(end_nodes,1);
e_time_phs = nan(end_nodes,1);
ell_phs = nan(end_nodes,1);
m_phs_all = nan(end_nodes,1);
c_phs = cell(end_nodes,1); % store coefficients struct

sN = nan(end_nodes,1);

%% Main sweep over N 
for k = start_nodes:end_nodes
    xi = st.fullintnodes{k};
    xb = st.bdrynodes{k};
    x = [xi; xb];

    % samples
    if dim==1
        y = f(x(:,1));
        ye_true = f(xe(:,1));
    elseif dim==2
        y = f(x(:,1), x(:,2));
        ye_true = f(xe(:,1), xe(:,2));
    else % dim==3
        y = f(x(:,1), x(:,2), x(:,3));
        ye_true = f(xe(:,1), xe(:,2), xe(:,3));
    end

    % polynomial total degree 
    ell = floor(fac * nthroot(size(x,1), dim));
    ell = max(ell, 1);
    ell_phs(k) = ell;

    % choose PHS order m. Safe default: m = max(1, floor(ell))
    % (ensures ℓ ≥ m for unisolvency with r^(2m+1))
    m_phs =  m_fixed;
    m_phs_all(k) = m_phs;

    % Run the PHS + poly solver
    [el2_phs(k), elinf_phs(k), a_time_phs(k), e_time_phs(k), c_phs{k}] = ...
        PHS_poly(x, y, ell, xe, alph, ye_true, m_phs);

    sN(k) = nthroot(size(x,1), dim);
    fprintf('k=%d | N=%d | ell=%d | m=%d | relL2=%.3e | Linf=%.3e\n', ...
            k, size(x,1), ell, m_phs, el2_phs(k), elinf_phs(k));
end

%% Save results 
timestamp = datestr(datetime('now'), 'yyyyMMdd_HHmmss');
results_dir = fullfile('results/', sprintf('%s', function_name), '/high');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
results_filename = fullfile(results_dir, sprintf('results_%s_phsp_%d.mat', function_name, phsdeg));

% If results_<function_name>.mat already exists from CSRBF/PLS runs,
% append the new PHS fields. Otherwise, create the file.
if exist(results_filename, 'file')
    save(results_filename, ...
        'el2_phs','elinf_phs','a_time_phs','e_time_phs','c_phs', ...
        'ell_phs','m_phs_all','phsdeg','m_fixed','sN','dim','function_name','timestamp', ...
        '-append');
else
    save(results_filename, ...
        'el2_phs','elinf_phs','a_time_phs','e_time_phs','c_phs', ...
        'ell_phs','m_phs_all','phsdeg','m_fixed','sN','dim','function_name','timestamp');
end

disp(['Saved PHS+poly results to: ', results_filename]);
