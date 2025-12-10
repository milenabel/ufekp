function DriverRun2DAbl()
%% Interpolation in 2D using CSRBFs + Jacobi polynomials
%% Ablation: boundary ("edge") effects vs interior-only evaluations
%% FS-only (fixed ε). No keff/eval strategies, no diagnostics.

%% Setup
dim = 2;  % fixed for this ablation
assert(dim==2,'DriverRun2DAbl is for dim=2 only.');

% Load nodes (same as your working 2D setup)
stL = load('DiskPoissonNodesLarge.mat');
xe_full = [stL.fullintnodes{7}; stL.bdrynodes{7}]; % evaluation set (full)
clear stL
st = load('DiskPoissonNodesClustered.mat'); % training sets (levels)

start_nodes = 1;
end_nodes = size(st.fullintnodes,1);

% Interior-only evaluation mask on xe_full: r <= 1 - deltaE
deltaE = 0.10; % "how far" from boundary we call interior
rf = sqrt(sum(xe_full.^2,2));
xe_int = xe_full(rf <= (1 - deltaE), :);

% Polynomial degree scaling factor
fac = 0.8;

% RBF: Wendland C2 (your stable choice)
rbf = @(e,r) -r.^4.*(4.*r - 5*spones(r));

% Baseline FS (fixed ε)
eps_baseline = 10;

%% Target functions (in the requested order)
tauBL = 0.05;
funcs(1).f = @(x,y) exp(-(1 - sqrt(x.^2 + y.^2))./tauBL);
funcs(1).nm = 'edge_bl_exp';

funcs(2).f = @(x,y) max(1 - sqrt(x.^2 + y.^2), 0).^0.25;
funcs(2).nm = 'edge_pow';

funcs(3).f = @(x,y) exp(x + y);
funcs(3).nm = 'mild_exp_xy';

funcs(4).f = @(x,y) 1./(1 + 9*(x.^2 + y.^2));
funcs(4).nm = 'mild_runge_r';

funcs(5).f = @(x,y) (x.^2 + y.^2).^(3/2);
funcs(5).nm = 'interior_cusp';

%% Sweep over target functions (save one .mat per function)
for jf = 1:numel(funcs)
    fhandle = funcs(jf).f;
    function_name = funcs(jf).nm;

    nfLevels = end_nodes - start_nodes + 1;
    el2_full_fs_base = nan(nfLevels,1);
    elinf_full_fs_base = nan(nfLevels,1);
    el2_int_fs_base = nan(nfLevels,1);
    elinf_int_fs_base = nan(nfLevels,1);

    lev = 0;
    for k = start_nodes:end_nodes
        lev = lev + 1;

        xi = st.fullintnodes{k};
        xb = st.bdrynodes{k};
        x = [xi; xb];

        % polynomial degree and KD-tree
        ell = max(1, floor(fac * nthroot(size(x,1), dim)));
        alph = 0; % Legendre
        tree = KDTreeSearcher(x);

        % Training values
        y = fhandle(x(:,1), x(:,2));

        % Full evaluation set
        ye_full = fhandle(xe_full(:,1), xe_full(:,2));
        [el2_full_fs_base(lev,1), elinf_full_fs_base(lev,1)] = ...
            call_CSRBFGen_errors_only(x, y, ell, xe_full, alph, rbf, eps_baseline, tree, ye_full);

        % Interior-only evaluation set
        if ~isempty(xe_int)
            ye_int = fhandle(xe_int(:,1), xe_int(:,2));
            [el2_int_fs_base(lev,1), elinf_int_fs_base(lev,1)] = ...
                call_CSRBFGen_errors_only(x, y, ell, xe_int, alph, rbf, eps_baseline, tree, ye_int);
        end
    end

    % (Optional convenience) edge gain = E_full / E_int
    edge_gain_fs_base = nan(size(el2_full_fs_base));
    ok = isfinite(el2_full_fs_base) & isfinite(el2_int_fs_base) & (el2_int_fs_base > 0);
    edge_gain_fs_base(ok) = el2_full_fs_base(ok) ./ el2_int_fs_base(ok);

    %% Save per-function results
    timestamp = datestr(datetime('now'), 'yyyyMMdd_HHmmss');
    results_dir = fullfile('results_ablation', function_name, 'edge_vs_interior');
    if ~exist(results_dir, 'dir'), mkdir(results_dir); end

    results_filename = fullfile(results_dir, sprintf('abl2d_%s_%s.mat', function_name, timestamp));
    save(results_filename, ...
        'function_name','deltaE','eps_baseline', ...
        'el2_full_fs_base','elinf_full_fs_base','el2_int_fs_base','elinf_int_fs_base', ...
        'edge_gain_fs_base', ...
        'start_nodes','end_nodes','dim','fac');

    fprintf('Saved FS-only ablation for %-15s -> %s\n', function_name, results_filename);
end

fprintf('\nAblation finished. Files are under results_ablation/<func>/edge_vs_interior/.\n');
end % DriverRun2DAbl


%% helper (thin wrapper)
function [el2, elinf] = call_CSRBFGen_errors_only(x, y, ell, xe, alph, rbf, ep, tree, ye_true)
% Call CSRBFGen with the provided evaluation set and ye_true,
% but only return L2 / Linf errors to keep the driver focused.
[el2, elinf] = CSRBFGen(x, y, ell, xe, alph, rbf, ep, tree, ye_true);
end
