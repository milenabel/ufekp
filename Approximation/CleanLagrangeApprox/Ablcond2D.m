function Ablcond2D()
%% Interpolation in 2D using CSRBFs + Jacobi polynomials
%% Ablation: boundary ("edge") effects vs interior-only evaluations
%% No changes to CSRBFGen.m required.

%% Setup
dim = 2;  % fixed for this ablation
assert(dim==2,'Ablcond2D is for dim=2 only.');

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

% RBF: Wendland C2 (same as your stable choice)
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

%% Plotting (separate file is needed)
% function plot_edgeabl(fs_variant, func_names, which_series, base_dir)
% % One figure per FS variant, functions in the legend.
% % fs_variant: 'base' | 'keff' | 'eval'
% % which_series: 'full' | 'interior' | 'edge_gain'
% % base_dir: root results folder (default 'results_ablation')
% 
%     if nargin < 4 || isempty(base_dir), base_dir = 'results_ablation'; end
%     if ischar(func_names), func_names = {func_names}; end
% 
%     valid_fs = {'base','keff','eval'};
%     assert(ismember(lower(fs_variant), valid_fs), 'fs_variant must be base|keff|eval');
%     valid_which = {'full','interior','edge_gain'};
%     assert(ismember(lower(which_series), valid_which), 'which_series must be full|interior|edge_gain');
% 
%     series = struct([]); all_sNs = [];
% 
%     for i = 1:numel(func_names)
%         fname = func_names{i};
%         abldir = fullfile(base_dir, fname, 'edge_vs_interior');
%         if ~exist(abldir,'dir')
%             error('Directory not found: %s', abldir);
%         end
%         matfile = getLatestMat(abldir, sprintf('abl2d_%s_', fname));
%         fprintf('[compare] loading: %s\n', matfile);
%         S = load(matfile);
%         assert(isfield(S,'dim') && S.dim==2, 'This plotter expects dim=2 ablation files.');
% 
%         % Build sN for exactly this file’s levels
%         st = load('DiskPoissonNodesClustered.mat');
%         levs = S.start_nodes:S.end_nodes;
%         sNs  = nan(numel(levs),1);
%         for t = 1:numel(levs)
%             k = levs(t);
%             x = [st.fullintnodes{k}; st.bdrynodes{k}];
%             sNs(t) = nthroot(size(x,1), S.dim);
%         end
% 
%         % Pick full/interior series for this FS variant
%         switch lower(fs_variant)
%             case 'base'
%                 fullv = getfield_safe(S,'el2_full_fs_base');
%                 intv  = getfield_safe(S,'el2_int_fs_base');
%             case 'keff'
%                 fullv = getfield_safe(S,'el2_full_fs_keff');
%                 intv  = getfield_safe(S,'el2_int_fs_keff');
%             case 'eval'
%                 fullv = getfield_safe(S,'el2_full_fs_eval');
%                 intv  = getfield_safe(S,'el2_int_fs_eval');
%         end
% 
%         % Align lengths: truncate to common min length for this function
%         m = min([numel(sNs), numel(fullv), numel(intv)]);
%         sNs  = sNs(1:m);
%         fullv = fullv(1:m);
%         intv  = intv(1:m);
% 
%         fprintf('[sizes] %s: sNs=%d, full=%d, int=%d -> using %d\n', ...
%             fname, numel(sNs), numel(fullv), numel(intv), m);
% 
%         % Select the plotting series
%         switch lower(which_series)
%             case 'full', y = fullv(:);
%             case 'interior', y = intv(:);
%             case 'edge_gain'
%                 y = nan(m,1);
%                 ok = isfinite(fullv) & isfinite(intv) & (intv > 0);
%                 y(ok) = fullv(ok)./intv(ok);
%         end
% 
%         series(i).name  = fname;
%         series(i).label = label_for_func(fname); % <-- pretty math label
%         series(i).x = sNs(:);
%         series(i).y = y(:);
% 
%         all_sNs = [all_sNs; sNs(:)];
%     end
% 
%     xmin = min(all_sNs); xmax = max(all_sNs);
% 
%     % Build y-lims from what we will actually plot
%     vals = [];
%     for i = 1:numel(series)
%         yi = series(i).y;
%         vals = [vals; yi(isfinite(yi) & yi>0)];
%     end
%     if isempty(vals)
%         ymin = 1e-16; ymax = 1;
%     else
%         ymin = 10^floor(log10(min(vals)));
%         ymax = 10^ceil (log10(max(vals)));
%         if strcmpi(which_series,'edge_gain')
%             ymin = max(ymin, 1e-2);
%             ymax = min(max(ymax, 1), 1e4);
%         end
%     end
% 
%     % Plot 
%     h = figure; set(h,'Color','none'); ax = gca;
%     set(ax,'Color','none','FontSize',12,'LineWidth',1.2);
%     marks = {'-o','-s','-^','-x','-d','-v','-+'};
%     H = gobjects(0); L = {};
% 
%     for i = 1:numel(series)
%         x = series(i).x; y = series(i).y;
%         good = isfinite(x) & isfinite(y) & (y > 0);
%         xg = x(good); yg = y(good);
%         if numel(xg) >= 1 && numel(yg) == numel(xg)
%             mk = marks{1 + mod(i-1, numel(marks))};
%             hl = semilogy(xg, yg, mk, 'LineWidth',2,'MarkerFaceColor','w','MarkerSize',6); hold on;
%             H(end+1) = hl;
%             L{end+1} = series(i).label; % <-- pretty label
%         else
%             fprintf('[warn] %s has no valid points after masking; skipping.\n', series(i).name);
%         end
%     end
% 
%     ax.LineWidth = 2; ax.FontSize = 16;
%     ax.XAxis.FontWeight = 'bold'; ax.YAxis.FontWeight = 'bold';
%     xlabel('N^{1/2}','Interpreter','tex','FontSize',16,'FontWeight','bold');
% 
%     switch lower(which_series)
%         case 'full', yl = 'Relative l_2 error (full domain)';
%         case 'interior', yl = 'Relative l_2 error (interior only)';
%         case 'edge_gain',yl = 'Edge gain = E_{full}/E_{int}';
%     end
%     ylabel(yl,'Interpreter','tex','FontSize',16,'FontWeight','bold');
% 
%     if ~isempty(H)
%         legend(H,L,'Location','best','Interpreter','tex','FontWeight','bold','FontSize',14);
%     end
% 
%     switch lower(fs_variant)
%         case 'base', ttl = 'FS (K_t = 1e12)';
%         case 'keff', ttl = 'FS-eff';
%         case 'eval', ttl = 'FS-eval';
%     end
%     title(sprintf('%s — %s', ttl, series_title(which_series)), ...
%         'Interpreter','tex','FontWeight','bold','FontSize',16);
% 
%     xlim([xmin xmax]); ylim([ymin ymax]);
% 
%     outdir = fullfile(base_dir, func_names{1}, 'edge_vs_interior');
%     if ~exist(outdir,'dir'), outdir = pwd; end
%     outfile = fullfile(outdir, sprintf('compare_%s_%s.png', lower(fs_variant), lower(which_series)));
%     export_fig(gcf, outfile, '-png','-r300','-transparent');
%     close(h);
%     fprintf('[compare] saved: %s\n', outfile);
% end
% 
% % helpers
% 
% function v = getfield_safe(S, fname)
%     if isfield(S, fname), v = S.(fname);
%     else, error('Missing field "%s" in %s', fname, inputname(1)); end
% end
% 
% function matfile = getLatestMat(dirpath, prefix)
%     d = dir(fullfile(dirpath, [prefix, '*.mat']));
%     if isempty(d), error('No ablation .mat files found in %s', dirpath); end
%     [~,idx] = max([d.datenum]);
%     matfile = fullfile(dirpath, d(idx).name);
% end
% 
% function s = series_title(which_series)
%     switch lower(which_series)
%         case 'full', s = 'Full-domain';
%         case 'interior', s = 'Interior-only';
%         case 'edge_gain',s = 'Edge gain';
%     end
% end
% 
% function s = label_for_func(nm)
% % Return a math label (TeX) for the function name
% % r = \sqrt{x^2+y^2}, (a)_+ = max(a,0)
%     switch lower(nm)
%         case 'edge_bl_exp'
%             s = 'e^{-(1-r)/\tau}';
%         case 'edge_pow'
%             s = '(1-r)_{+}^{0.25}';
%         case 'mild_exp_xy'
%             s = 'e^{x+y}';
%         case 'mild_runge_r'
%             s = '(1+9r^2)^{-1}';
%         case 'interior_cusp'
%             s = '(x^2+y^2)^{3/2}';
%         otherwise
%             % fallback: clean underscores
%             s = strrep(nm,'_','\_');
%     end
% end