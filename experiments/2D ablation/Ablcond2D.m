function Ablcond2D()
%% Interpolation in 2D using CSRBFs + Jacobi polynomials
%% Ablation: boundary ("edge") effects vs interior-only evaluations
%% No changes to CSRBFGen.m required.

%% Setup
dim = 2;  % fixed for this ablation
assert(dim==2,'Ablcond2D is for dim=2 only.');

% Load nodes 
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

% RBF: Wendland C2 
rbf = @(e,r) -r.^4.*(4.*r - 5*spones(r));

% Storage (per function, per level)
nF = 5; % number of target functions below
nfLevels = end_nodes - start_nodes + 1;

el2_full_fs_base = nan(nfLevels, nF);
elinf_full_fs_base = nan(nfLevels, nF);
el2_int_fs_base = nan(nfLevels, nF);
elinf_int_fs_base = nan(nfLevels, nF);

el2_full_fs_keff = nan(nfLevels, nF);
elinf_full_fs_keff = nan(nfLevels, nF);
el2_int_fs_keff = nan(nfLevels, nF);
elinf_int_fs_keff = nan(nfLevels, nF);

el2_full_fs_eval = nan(nfLevels, nF);
elinf_full_fs_eval = nan(nfLevels, nF);
el2_int_fs_eval = nan(nfLevels, nF);
elinf_int_fs_eval = nan(nfLevels, nF);

kappa_eff_keff = nan(nfLevels, nF);    % achieved κ_eff (diagnostic; full system)
medianL2_eval = nan(nfLevels, nF);     % achieved median ||L||_2 (diagnostic)

eps_chosen_keff = nan(1, nF);          % chosen ε on finest (reused)
eps_chosen_eval = nan(1, nF);          % chosen ε on finest (reused)

% For reproducible baseline FS (fixed ε)
eps_baseline = 10;                        

%% Target functions 
% Edge-heavy (funny at boundary):
% 1) boundary layer: exp(-(1-r)/tau), tau small
% 2) singular-ish slope at boundary: (1 - r)^{1/4}_+  (C^0 at boundary, steep gradient)
% Edge-mild (not boundary-focused):
% 3) analytic everywhere: exp(x+y)
% 4) mild rational: 1/(1 + 9 r^2)
% 5) center-cusp (interior feature, not edge): (x^2+y^2)^{3/2}

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

%% Sweep over target functions 
for jf = 1:nF
    fhandle = funcs(jf).f;
    function_name = funcs(jf).nm;

    % Build eps choices on the finest grid ONCE per function (then reuse)
    k_fin = end_nodes;
    xfin = [st.fullintnodes{k_fin}; st.bdrynodes{k_fin}];
    ellf = max(1, floor(fac * nthroot(size(xfin,1), dim)));
    [Pfin, idxsfin, recfin, c1fin, c2fin] = buildP(xfin, ellf, 0);
    treef = KDTreeSearcher(xfin);

    % Pick ε for FS-κeff and FS-eval on the finest grid, once per function
    if isnan(eps_chosen_keff(1,jf))
        eps_chosen_keff(1,jf) = find_eps_for_kappa_eff( xfin, treef, rbf, Pfin, 1e12, [3, 15] );
    end
    if isnan(eps_chosen_eval(1,jf))
        eps_chosen_eval(1,jf) = find_eps_min_evalcond( xfin, xe_full, treef, rbf, Pfin, c1fin, c2fin, idxsfin, [3, 15], recfin );
    end

    % Reuse across all levels; evaluate on full vs interior xe
    lev = 0;
    for k = start_nodes:end_nodes
        lev = lev + 1;

        xi = st.fullintnodes{k};
        xb = st.bdrynodes{k};
        x  = [xi; xb];

        ell = max(1, floor(fac * nthroot(size(x,1), dim)));
        alph = 0; % Legendre
        tree = KDTreeSearcher(x);

        % Training values
        y = fhandle(x(:,1), x(:,2));

        % FULL evaluation set 
        ye_full = fhandle(xe_full(:,1), xe_full(:,2));

        % Baseline FS (fixed eps)
        [el2_full_fs_base(lev,jf), elinf_full_fs_base(lev,jf)] = ...
            call_CSRBFGen_errors_only(x, y, ell, xe_full, alph, rbf, eps_baseline, tree, ye_full);

        % FS-κeff (reuse eps from finest)
        epK = eps_chosen_keff(1,jf);
        [el2_full_fs_keff(lev,jf), elinf_full_fs_keff(lev,jf)] = ...
            call_CSRBFGen_errors_only(x, y, ell, xe_full, alph, rbf, epK, tree, ye_full);

        % κ_eff diagnostic (for reporting): uses the training system (not xe)
        P_scaled = buildP_only(x, ell, alph);
        K_eff    = rbf(epK, DistanceMatrixCSRBFwt(x, x, epK, tree));
        kappa_eff_keff(lev,jf) = rbf_effective_cond(K_eff, P_scaled);

        % FS-evalcond (reuse eps from finest)
        epE = eps_chosen_eval(1,jf);
        [el2_full_fs_eval(lev,jf), elinf_full_fs_eval(lev,jf)] = ...
            call_CSRBFGen_errors_only(x, y, ell, xe_full, alph, rbf, epE, tree, ye_full);

        % median ||L||_2 diagnostic (using a subsample of xe_full)
        medianL2_eval(lev,jf) = compute_medianL2(x, xe_full, alph, rbf, epE, tree, P_scaled);

        % INTERIOR-ONLY evaluation set
        if isempty(xe_int)
            warning('xe_int is empty for deltaE=%.3f, skipping interior-only metrics.', deltaE);
        else
            ye_int = fhandle(xe_int(:,1), xe_int(:,2));

            % Baseline FS (fixed eps)
            [el2_int_fs_base(lev,jf), elinf_int_fs_base(lev,jf)] = ...
                call_CSRBFGen_errors_only(x, y, ell, xe_int, alph, rbf, eps_baseline, tree, ye_int);

            % FS-κeff
            [el2_int_fs_keff(lev,jf), elinf_int_fs_keff(lev,jf)] = ...
                call_CSRBFGen_errors_only(x, y, ell, xe_int, alph, rbf, epK, tree, ye_int);

            % FS-evalcond
            [el2_int_fs_eval(lev,jf), elinf_int_fs_eval(lev,jf)] = ...
                call_CSRBFGen_errors_only(x, y, ell, xe_int, alph, rbf, epE, tree, ye_int);
        end
    end

    %% Save per-function results (keeps files small & tidy)
    timestamp = datestr(datetime('now'), 'yyyyMMdd_HHmmss');
    results_dir = fullfile('results_ablation/', function_name, '/edge_vs_interior');
    if ~exist(results_dir, 'dir'), mkdir(results_dir); end

    results_filename = fullfile(results_dir, sprintf('abl2d_%s_%s.mat', function_name, timestamp));

    save(results_filename, ...
        'function_name', 'deltaE', ...
        'eps_baseline', 'eps_chosen_keff', 'eps_chosen_eval', ...
        'el2_full_fs_base','elinf_full_fs_base','el2_int_fs_base','elinf_int_fs_base', ...
        'el2_full_fs_keff','elinf_full_fs_keff','el2_int_fs_keff','elinf_int_fs_keff', ...
        'el2_full_fs_eval','elinf_full_fs_eval','el2_int_fs_eval','elinf_int_fs_eval', ...
        'kappa_eff_keff','medianL2_eval', ...
        'start_nodes','end_nodes','dim','fac');
    fprintf('Saved ablation results for %s -> %s\n', function_name, results_filename);
end

fprintf('\nAblation finished. Files are under results_ablation/<func>/edge_vs_interior/.\n');

end % Ablcond2D


%% Small helpers (self-contained)
function [P_scaled, idxs, recurrence, c1, c2] = buildP(x, ell, alph)
    xmax = max(x,[],1); xmin = min(x,[],1);
    c1 = 0.5*(xmax + xmin); 
    c2 = 2./(xmax - xmin);
    recurrence = @(N) jacobi_recurrence(N, alph, alph);
    idxs = total_degree_indices(size(x,2), ell);
    P_scaled = mpoly_eval((x - c1).*c2, idxs, recurrence);
end

function P_scaled = buildP_only(x, ell, alph)
    % version that only returns P (for diagnostics), to keep memory low
    xmax = max(x,[],1); xmin = min(x,[],1);
    c1 = 0.5*(xmax + xmin); 
    c2 = 2./(xmax - xmin);
    recurrence = @(N) jacobi_recurrence(N, alph, alph);
    idxs = total_degree_indices(size(x,2), ell);
    P_scaled = mpoly_eval((x - c1).*c2, idxs, recurrence);
end

function [el2, elinf] = call_CSRBFGen_errors_only(x, y, ell, xe, alph, rbf, ep, tree, ye_true)
    % Wrapper: call CSRBFGen with the provided evaluation set and ye_true,
    % but only return the L2 / Linf errors (to keep the driver focused).
    [el2, elinf] = CSRBFGen(x, y, ell, xe, alph, rbf, ep, tree, ye_true);
end

function eps_opt = find_eps_for_kappa_eff(x, tree, rbf, P_scaled, target_keff, eps_bracket)
% Root-find epsilon so that κ_eff(A_aug) ≈ target_keff, robust to no sign-change.
    target_log = log10(target_keff);
    function v = mismatch(ep)
        K = rbf(ep, DistanceMatrixCSRBFwt(x, x, ep, tree));
        v  = log10(rbf_effective_cond(K, P_scaled)) - target_log;
    end
    a = eps_bracket(1); b = eps_bracket(2);
    grid = linspace(a,b,8);
    vals = arrayfun(@(e) mismatch(e), grid);
    id = find(sign(vals(1:end-1)).*sign(vals(2:end))<=0, 1);
    if ~isempty(id)
        eps_opt = fzero(@mismatch, [grid(id), grid(id+1)]);
    else
        eps_opt = fminbnd(@(e) abs(mismatch(e)), a, b);
    end
end

function eps_opt = find_eps_min_evalcond(x, xe, tree, rbf, P_scaled, c1, c2, idxs, eps_bounds, recurrence)
% Minimize the median ||L(x*)||_2 over a subsample of Xe to pick epsilon.
    qstep = max(1, floor(size(xe,1)/2000));
    Xsample = xe(1:qstep:end, :);
    function m = median_L2(ep)
        K = rbf(ep, DistanceMatrixCSRBFwt(x, x, ep, tree));
        Aaug = [K, P_scaled; P_scaled', sparse(size(P_scaled,2), size(P_scaled,2))];
        Dqr  = decomposition(Aaug,'qr');
        solveTt = @(rhs) Dqr \ rhs;
        Kx_fun  = @(xstar) rbf(ep, DistanceMatrixCSRBFwt(xstar, x, ep, tree)).';
        Px_fun  = @(xstar) mpoly_eval((xstar - c1).*c2, idxs, recurrence);
        [~, ~, ~, med2] = eval_condition_numbers(solveTt, Kx_fun, Px_fun, Xsample);
        m = med2;
    end
    a = eps_bounds(1); b = eps_bounds(2);
    try
        eps_opt = fminbnd(@median_L2, a, b);
    catch
        grid = linspace(a,b,12);
        vals = arrayfun(@median_L2, grid);
        [~,ii] = min(vals);
        eps_opt = grid(ii);
    end
end

function med2 = compute_medianL2(x, xe, alph, rbf, ep, tree, P_scaled)
% Compute median ||L(x*)||_2 (diagnostic) on a subsample of xe.
    K = rbf(ep, DistanceMatrixCSRBFwt(x, x, ep, tree));
    Aaug = [K, P_scaled; P_scaled', sparse(size(P_scaled,2), size(P_scaled,2))];
    Dqr = decomposition(Aaug,'qr');
    solveTt = @(rhs) Dqr \ rhs;
    Kx_fun = @(xstar) rbf(ep, DistanceMatrixCSRBFwt(xstar, x, ep, tree)).';
    % Rebuild scaling for Px(x*)
    xmax = max(x,[],1); xmin = min(x,[],1);
    c1 = 0.5*(xmax + xmin); 
    c2 = 2./(xmax - xmin);
    recurrence = @(N) jacobi_recurrence(N, alph, alph);
    ell = size(P_scaled,2); %#ok<NASGU>  % not used directly
    idxs = total_degree_indices(size(x,2), infer_degree_from_cols(size(P_scaled,2), size(x,2)));
    Px_fun = @(xstar) mpoly_eval((xstar - c1).*c2, idxs, recurrence);
    qstep = max(1, floor(size(xe,1)/2000));
    Xsample = xe(1:qstep:end, :);
    [~, ~, ~, med2] = eval_condition_numbers(solveTt, Kx_fun, Px_fun, Xsample);
end

function ell = infer_degree_from_cols(mcols, d)
% crude inversion of total-degree basis size C(d+ell, d)
    % Solve C(d+ell, d) = mcols for ell (small d)
    ell = 0;
    while nchoosek(d+ell, d) < mcols
        ell = ell + 1;
        if ell > 2000, break; end
    end
end