function plot_edgeabl_all(func_name, base_dir)
% One plot, one function: FS (K_t = 1e12) only.
% Curves: full-domain error, interior-only error, and edge gain.
%
% Usage:
%   plot_edgeabl_all('edge_bl_exp');
%   plot_edgeabl_all('mild_exp_xy','results_ablation');

    if nargin < 2 || isempty(base_dir), base_dir = 'results_ablation'; end

    % Load latest ablation for this function
    abldir = fullfile(base_dir, func_name, 'edge_vs_interior');
    if ~exist(abldir,'dir')
        error('Directory not found: %s', abldir);
    end
    matfile = getLatestMat(abldir, sprintf('abl2d_%s_', func_name));
    fprintf('[FS-only] loading: %s\n', matfile);
    S = load(matfile);
    assert(isfield(S,'dim') && S.dim==2, 'expects dim=2 ablation file');

    % Build sN (N^{1/2}) for exactly these levels
    st = load('DiskPoissonNodesClustered.mat');
    levs = S.start_nodes:S.end_nodes;
    sNs  = nan(numel(levs),1);
    for t = 1:numel(levs)
        k = levs(t);
        x = [st.fullintnodes{k}; st.bdrynodes{k}];
        sNs(t) = nthroot(size(x,1), S.dim);
    end

    % Pull FS(base) series and align
    full_b = getfield_safe(S,'el2_full_fs_base');
    int_b = getfield_safe(S,'el2_int_fs_base');
    m = min([numel(sNs), numel(full_b), numel(int_b)]);
    sNs = sNs(1:m);
    full_b = full_b(1:m);
    int_b = int_b(1:m);
    edge_g = safe_ratio(full_b, int_b);

    % y-lims from all three curves (treat as mixed: errors + ratio)
    [ymin,ymax] = lims_from_mixed([full_b(:); int_b(:); edge_g(:)]);

    % Plot (same style)
    h = figure; set(h,'Color','none'); ax = gca;
    set(ax,'Color','none','FontSize',12,'LineWidth',1.2);
    marks = {'-o','-s','-^'};
    labs = {'Full-domain error','Interior-only error','Edge gain  E_{full}/E_{int}'};

    plot_one(sNs, full_b, marks{1});
    plot_one(sNs, int_b,  marks{2});
    plot_one(sNs, edge_g, marks{3});

    ax.LineWidth = 2; ax.FontSize = 16;
    ax.XAxis.FontWeight = 'bold'; ax.YAxis.FontWeight = 'bold';
    xlabel('N^{1/2}','Interpreter','tex','FontSize',16,'FontWeight','bold');
    ylabel('Relative l_2 error / Edge gain','Interpreter','tex','FontSize',16,'FontWeight','bold');
    legend(labs,'Location','best','Interpreter','tex','FontWeight','bold','FontSize',14);

    ttl = sprintf('Relative l_2 error vs. N^{1/d}');
    title(ttl,'Interpreter','tex','FontWeight','bold','FontSize',16);

    xlim([min(sNs) max(sNs)]); ylim([ymin ymax]);

    outfile = fullfile(abldir, sprintf('fs_only_%s_full_interior_edge.png', lower(func_name)));
    export_fig(gcf, outfile, '-png','-r300','-transparent');
    close(h);
    fprintf('[FS-only] saved: %s\n', outfile);
end

% helpers

function plot_one(x, y, mk)
    % Make both column vectors and length-match BEFORE masking
    x = x(:); y = y(:);
    m = min(numel(x), numel(y));
    x = x(1:m); y = y(1:m);
    good = isfinite(x) & isfinite(y) & (y > 0);
    xg = x(good); yg = y(good);
    if ~isempty(xg)
        semilogy(xg, yg, mk, 'LineWidth',2,'MarkerFaceColor','w','MarkerSize',6); hold on;
    end
end

function r = safe_ratio(a,b)
    a = a(:); b = b(:);
    m = min(numel(a), numel(b));
    a = a(1:m); b = b(1:m);
    r = nan(size(a));
    ok = isfinite(a) & isfinite(b) & b > 0;
    r(ok) = a(ok)./b(ok);
end

function [ymin,ymax] = lims_from_mixed(vals)
    v = vals(isfinite(vals) & vals > 0);
    if isempty(v)
        ymin = 1e-16; ymax = 1;
    else
        ymin = 10^floor(log10(min(v)));
        ymax = 10^ceil (log10(max(v)));
        % keep ratio range sane if it dominates
        ymax = min(max(ymax, 1e-2), 1e6);
    end
end

function v = getfield_safe(S, fname)
    if isfield(S, fname), v = S.(fname);
    else, error('Missing field "%s"', fname);
    end
end

function matfile = getLatestMat(dirpath, prefix)
    d = dir(fullfile(dirpath, [prefix, '*.mat']));
    if isempty(d), error('No ablation .mat files found in %s', dirpath); end
    [~,idx] = max([d.datenum]);
    matfile = fullfile(dirpath, d(idx).name);
end

function s = label_for_func(nm)
    switch lower(nm)
        case 'edge_bl_exp'
            s = 'e^{-(1-r)/\tau}';
        case 'edge_pow'
            s = '(1-r)_{+}^{0.25}';
        case 'mild_exp_xy'
            s = 'e^{x+y}';
        case 'mild_runge_r'
            s = '(1+9r^2)^{-1}';
        case 'interior_cusp'
            s = '(x^2+y^2)^{3/2}';
        otherwise
            s = strrep(nm,'_','\_');
    end
end
