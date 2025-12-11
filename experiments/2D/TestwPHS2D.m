% TestwPHS.m for 2D
% Compute global PHS + polynomial results on the same setups 
% and save/append them into results_<function_name>.mat

clear; clc;

%% Spatial dimension
dim = 2; % 2
phsdeg = 9; % manually change the curretn PHS degree to record results correctly 
assert(ismember(phsdeg,[5 7 9]), 'phsdeg must be one of {5,7,9}.'); % (current picks: 5, 7, 9)
m_fixed = (phsdeg - 1)/2; % m = 2, 3, or 4
assert(m_fixed == round(m_fixed), 'phsdeg must be odd so (phsdeg-1)/2 is integer.');
m_fixed = round(m_fixed);

fac = (dim==2) * 1.0 + (dim~=2) * 1.0; % keep consistent with other runs
alph = 0; % Legendre (Jacobi alpha=0), matches PLS

%% Load the same node sets
% Evaluation nodes (finest)
st_e = load('DiskPoissonNodesLarge.mat');
xe = [st_e.fullintnodes{7}; st_e.bdrynodes{7}]; 

% Training nodes
st = load('DiskPoissonNodesClustered.mat');
clear p N pxe;

start_nodes = 1;
end_nodes = size(st.fullintnodes,1);
if dim == 3
    end_nodes = min(end_nodes, 5);
end

% Choose the same target function you used in the corresponding run 
syms x y;    
f = (x.^2 + y.^2).^(3/2);                function_name = 'xy_p_2d';
% f = exp( ((x + y).^(2))/0.2 );         function_name = 'exp_p_2d';
dfx = diff(f,x); dfy = diff(f,y);
f = matlabFunction(f); % true function for errors

dfx = matlabFunction(dfx);
dfy = matlabFunction(dfy);
clear x y;

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
    y = f(x(:,1), x(:,2));
    ye_true = f(xe(:,1), xe(:,2));

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
