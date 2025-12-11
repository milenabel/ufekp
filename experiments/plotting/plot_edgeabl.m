function plot_edgeabl(func_names, which_series, base_dir)
% FS-only ablation plotter: one figure, functions in the legend.
% Usage:
%   plot_edgeabl({'edge_bl_exp','edge_pow','mild_exp_xy','mild_runge_r','interior_cusp'}, 'full');
%   plot_edgeabl('edge_bl_exp', 'edge_gain');
%
% Inputs
%   func_names   : string or cellstr of function IDs (folder names under results_ablation)
%   which_series : 'full' | 'interior' | 'edge_gain'
%   base_dir     : (optional) root results folder, default 'results_ablation'
%
% Expects MAT files produced by DriverRun2DAbl.m (FS-only).

    if nargin < 3 || isempty(base_dir), base_dir = 'results_ablation'; end
    if ischar(func_names), func_names = {func_names}; end

    valid_which = {'full','interior','edge_gain'};
    assert(ismember(lower(which_series), valid_which), 'which_series must be full|interior|edge_gain');

    series = struct([]); all_sNs = [];

    for i = 1:numel(func_names)
        fname = func_names{i};
        abldir = fullfile(base_dir, fname, 'edge_vs_interior');
        if ~exist(abldir,'dir')
            error('Directory not found: %s', abldir);
        end
        matfile = getLatestMat(abldir, sprintf('abl2d_%s_', fname));
        fprintf('[compare] loading: %s\n', matfile);
        S = load(matfile);
        assert(isfield(S,'dim') && S.dim==2, 'This plotter expects dim=2 ablation files.');

        % Build sN for exactly this file’s levels
        st = load('DiskPoissonNodesClustered.mat');
        levs = S.start_nodes:S.end_nodes;
        sNs  = nan(numel(levs),1);
        for t = 1:numel(levs)
            k = levs(t);
            x = [st.fullintnodes{k}; st.bdrynodes{k}];
            sNs(t) = nthroot(size(x,1), S.dim);
        end

        % FS baseline series from the file
        fullv = getfield_safe(S,'el2_full_fs_base');
        intv = getfield_safe(S,'el2_int_fs_base');

        % Align lengths: truncate to common min length for this function
        m = min([numel(sNs), numel(fullv), numel(intv)]);
        sNs = sNs(1:m);
        fullv = fullv(1:m);
        intv = intv(1:m);

        fprintf('[sizes] %s: sNs=%d, full=%d, int=%d -> using %d\n', ...
            fname, numel(sNs), numel(fullv), numel(intv), m);

        % Select the plotting series
        switch lower(which_series)
            case 'full'
                y = fullv(:);
            case 'interior'
                y = intv(:);
            case 'edge_gain'
                y = nan(m,1);
                ok = isfinite(fullv) & isfinite(intv) & (intv > 0);
                y(ok) = fullv(ok)./intv(ok);
        end

        series(i).name = fname;
        series(i).label = label_for_func(fname); % TeX pretty label
        series(i).x = sNs(:);
        series(i).y = y(:);

        all_sNs = [all_sNs; sNs(:)];
    end

    xmin = min(all_sNs); xmax = max(all_sNs);

    % Build y-lims from what we will actually plot
    vals = [];
    for i = 1:numel(series)
        yi = series(i).y;
        vals = [vals; yi(isfinite(yi) & yi>0)];
    end
    if isempty(vals)
        ymin = 1e-16; ymax = 1;
    else
        ymin = 10^floor(log10(min(vals)));
        ymax = 10^ceil (log10(max(vals)));
        if strcmpi(which_series,'edge_gain')
            ymin = max(ymin, 1e-2);
            ymax = min(max(ymax, 1), 1e4);
        end
    end

    % Plot
    h = figure; set(h,'Color','none'); ax = gca;
    set(ax,'Color','none','FontSize',12,'LineWidth',1.2);
    marks = {'-o','-s','-^','-x','-d','-v','-+'};
    H = gobjects(0); L = {};

    for i = 1:numel(series)
        x = series(i).x; y = series(i).y;
        good = isfinite(x) & isfinite(y) & (y > 0);
        xg = x(good); yg = y(good);
        if numel(xg) >= 1 && numel(yg) == numel(xg)
            mk = marks{1 + mod(i-1, numel(marks))};
            hl = semilogy(xg, yg, mk, 'LineWidth',2,'MarkerFaceColor','w','MarkerSize',6); hold on;
            H(end+1) = hl;
            L{end+1} = series(i).label; % pretty label
        else
            fprintf('[warn] %s has no valid points after masking; skipping.\n', series(i).name);
        end
    end

    ax.LineWidth = 2; 
    ax.FontSize = 16;
    ax.XAxis.FontWeight = 'bold'; 
    ax.YAxis.FontWeight = 'bold';
    xlabel('N^{1/2}','Interpreter','tex','FontSize',16,'FontWeight','bold');

    switch lower(which_series)
        case 'full', yl = 'Relative l_2 error (full domain)';
        case 'interior', yl = 'Relative l_2 error (interior only)';
        case 'edge_gain', yl = 'Edge gain = E_{full}/E_{int}';
    end
    ylabel(yl,'Interpreter','tex','FontSize',16,'FontWeight','bold');

    if ~isempty(H)
        legend(H,L,'Location','best','Interpreter','tex','FontWeight','bold','FontSize',14);
    end

    title(sprintf('FS (fixed \\epsilon) — %s', series_title(which_series)), 'Interpreter','tex','FontWeight','bold','FontSize',16);

    xlim([xmin xmax]); ylim([ymin ymax]);

    % Save beside the first function’s results folder for convenience
    outdir = fullfile(base_dir, func_names{1}, 'edge_vs_interior');
    if ~exist(outdir,'dir'), outdir = pwd; end
    outfile = fullfile(outdir, sprintf('compare_fs_%s.png', lower(which_series)));
    export_fig(gcf, outfile, '-png','-r300','-transparent');
    close(h);
    fprintf('[compare] saved: %s\n', outfile);
end


% helpers

function v = getfield_safe(S, fname)
    if isfield(S, fname), v = S.(fname);
    else, error('Missing field "%s" in %s', fname, inputname(1)); end
end

function matfile = getLatestMat(dirpath, prefix)
    d = dir(fullfile(dirpath, [prefix, '*.mat']));
    if isempty(d), error('No ablation .mat files found in %s', dirpath); end
    [~,idx] = max([d.datenum]);
    matfile = fullfile(dirpath, d(idx).name);
end

function s = series_title(which_series)
    switch lower(which_series)
        case 'full', s = 'Full-domain';
        case 'interior', s = 'Interior-only';
        case 'edge_gain',s = 'Edge gain';
    end
end

function s = label_for_func(nm)
% Return a math label (TeX) for the function name
% r = \sqrt{x^2+y^2}, (a)_+ = max(a,0)
    switch lower(nm)
        case 'edge_bl_exp', s = 'e^{-(1-r)/\tau}';
        case 'edge_pow', s = '(1-r)_{+}^{0.25}';
        case 'mild_exp_xy', s = 'e^{x+y}';
        case 'mild_runge_r', s = '(1+9r^2)^{-1}';
        case 'interior_cusp', s = '(x^2+y^2)^{3/2}';
        otherwise, s = strrep(nm,'_','\_');
    end
end
