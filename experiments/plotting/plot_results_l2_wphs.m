%% Pick a function to run
function plot_results_l2_wphs(function_name, subdir, base_results_dir, smoothness_idx, variety)
    % Usage:
    %   plot_results_l2_wphs('abs_1d','high','results',1,true);  
    %
    % variety (logical, default = true):
    %   true  -> plot PLS, Diag, FS with K_t = 1e12, 1e8, 1e4 (like 1D
    %   panels), PHS for FS only
    %   false -> plot PLS, Diag, FS(fs1) only (for 2D/3D), PHS
    % make_fs_family_plot = true;  % <-- set to false if you don't want the extra FS-only figure
    if nargin < 3 || isempty(base_results_dir), base_results_dir = fullfile('results'); end
    if nargin < 2 || isempty(subdir), subdir = 'high'; end
    if nargin < 4 || isempty(smoothness_idx), smoothness_idx = 1; end
    if nargin < 5 || isempty(variety), variety = true; end
    sm = max(1,min(3,smoothness_idx));

    % Load unified results
    res = loadResults_with_PHS(function_name, subdir, base_results_dir);
    res.sNs = res.sN;
    dim = res.dim;

    % Global x-limits
    xmin = min(res.sNs(:));
    xmax = max(res.sNs(:));

    % Which FS series to include
    include_all_Kt = variety || (dim==1);   % 1D always “all K_t”

    % Y-lims from the series we’ll actually plot
    err_vals = [];
    if isfield(res,'el2_poly'), err_vals = [err_vals; res.el2_poly(:)]; end
    if isfield(res,'el2_diag'), err_vals = [err_vals; res.el2_diag(:,sm)]; end
    if isfield(res,'el2_fs1'),  err_vals = [err_vals; res.el2_fs1(:,sm)]; end
    if include_all_Kt
        if isfield(res,'el2_fs2'), err_vals = [err_vals; res.el2_fs2(:,sm)]; end
        if isfield(res,'el2_fs3'), err_vals = [err_vals; res.el2_fs3(:,sm)]; end
        if isfield(res,'el2_fc1'), err_vals = [err_vals; res.el2_fc1(:,sm)]; end
        if isfield(res,'el2_fc2'), err_vals = [err_vals; res.el2_fc2(:,sm)]; end
        if isfield(res,'el2_fc3'), err_vals = [err_vals; res.el2_fc3(:,sm)]; end
    end
    if isfield(res,'el2_fs_eff'), err_vals = [err_vals; res.el2_fs_eff(:)]; end
    if isfield(res,'el2_fs_eval'), err_vals = [err_vals; res.el2_fs_eval(:)]; end
    % include PHS (any dim that has it)
    if isfield(res,'el2_phs'), err_vals = [err_vals; res.el2_phs(:)]; end

    err_vals = err_vals(isfinite(err_vals) & err_vals>0);
    if isempty(err_vals), err_min = 1e-16; err_max = 1;
    else
        err_min = 10^floor(log10(min(err_vals)));
        err_max = 10^ceil (log10(max(err_vals)));
    end

    % Main plot
    plotErrorVsN(res, sm, xmin, xmax, err_min, err_max, include_all_Kt);
end


%% Load results (+ optional PHS p=9)
function res = loadResults_with_PHS(function_name, subdir, base_results_dir)
    valid_subdirs = {'high','low','fixed','mixed'};
    if ~ismember(subdir, valid_subdirs)
        error('Invalid subdirectory: %s. Must be one of: high, low, fixed, mixed', subdir);
    end
    results_dir = fullfile(base_results_dir, function_name, subdir);
    matfile = fullfile(results_dir, sprintf('results_%s.mat', function_name));
    if ~exist(matfile,'file'), error('Results file not found: %s', matfile); end
    fprintf('[plot] loading: %s\n', matfile);
    data = load(matfile);
    data.results_dir = results_dir;

    % Attach PHS p=9 if present (el2_phs, a_time_phs, e_time_phs)
    phsfile = fullfile(results_dir, sprintf('results_%s_phsp_9.mat', function_name));
    if exist(phsfile,'file')
        Phs = load(phsfile);
        n = numel(data.sN(:));
        data.el2_phs = padToN(getVecOrNaN(Phs,'el2_phs'), n);
        data.a_time_phs = padToN(getVecOrNaN(Phs,'a_time_phs'), n);
        data.e_time_phs = padToN(getVecOrNaN(Phs,'e_time_phs'), n);
        fprintf('[plot] PHS p=9 attached: %s\n', phsfile);
    end
    res = data;
end

function v = getVecOrNaN(S, fld)
    if isfield(S,fld), v = S.(fld)(:); else, v = []; end
end
function v = padToN(v, n)
    v = v(:);
    if isempty(v), v = NaN(n,1); return; end
    m = numel(v);
    if m >= n, v = v(1:n); else, v = [v; NaN(n-m,1)]; end
end


%% Plot relative L2 error vs N^{1/d} (main panel) — with unified styles
function plotErrorVsN(res, sm, xmin, xmax, err_min, err_max, include_all_Kt)
    dim = res.dim;
    sNs = res.sNs(:);
    n = numel(sNs);

    styles = getPlotStyles(dim);   

    % Series
    el2_poly_s = NaN(n,1); if isfield(res,'el2_poly'), el2_poly_s = res.el2_poly; end
    el2_diag_s = NaN(n,1); if isfield(res,'el2_diag'), el2_diag_s = res.el2_diag(:,sm); end

    el2_fs1_s = NaN(n,1); if isfield(res,'el2_fs1'), el2_fs1_s = res.el2_fs1(:,sm); end
    el2_fs2_s = NaN(n,1); if isfield(res,'el2_fs2'), el2_fs2_s = res.el2_fs2(:,sm); end
    el2_fs3_s = NaN(n,1); if isfield(res,'el2_fs3'), el2_fs3_s = res.el2_fs3(:,sm); end

    el2_fc1_s = NaN(n,1); if isfield(res,'el2_fc1'), el2_fc1_s = res.el2_fc1(:,sm); end
    el2_fc2_s = NaN(n,1); if isfield(res,'el2_fc2'), el2_fc2_s = res.el2_fc2(:,sm); end
    el2_fc3_s = NaN(n,1); if isfield(res,'el2_fc3'), el2_fc3_s = res.el2_fc3(:,sm); end


    el2_fs_eff_s = NaN(n,1); if isfield(res,'el2_fs_eff'), el2_fs_eff_s = res.el2_fs_eff(:); end
    el2_fs_eval_s = NaN(n,1); if isfield(res,'el2_fs_eval'), el2_fs_eval_s = res.el2_fs_eval(:); end

    el2_phs_s = NaN(n,1); if isfield(res,'el2_phs'), el2_phs_s = res.el2_phs; end

    % 3D: keep only 5 points total, trimmed to last valid PHS
    if dim == 3
        % last index where PHS exists (if present)
        if any(isfinite(el2_phs_s) & el2_phs_s > 0)
            lastPHS = find(isfinite(el2_phs_s) & el2_phs_s > 0, 1, 'last');
        else
            lastPHS = n;  % no PHS -> just use available length
        end
        keep = 1:min(5, lastPHS);

        % apply to x and every series so all curves end at the same x
        sNs = sNs(keep);
        el2_poly_s = el2_poly_s(keep);
        el2_diag_s = el2_diag_s(keep);
        el2_fs1_s = el2_fs1_s(keep);
        el2_fs2_s = el2_fs2_s(keep);
        el2_fs3_s = el2_fs3_s(keep);
        el2_fs_eff_s = el2_fs_eff_s(keep);
        el2_fs_eval_s = el2_fs_eval_s(keep);
        el2_phs_s = el2_phs_s(keep);
        n = numel(sNs); % keep n consistent

    end

    % Exponential trendline (FS1)
    ok = ~(isnan(el2_fs1_s) | el2_fs1_s<=0);
    if any(ok)
        xx = sNs(ok); yy = el2_fs1_s(ok);
        coeff_exp = polyfit(xx, log(yy), 1);
        yfit_exp = exp(coeff_exp(2)) * exp( coeff_exp(1) * xx );
        fprintf('Exponential fit (FS1): E(N) = %.2f · exp(%.2f·N^{1/d})\n', exp(coeff_exp(2)), coeff_exp(1));
    else
        xx = []; yfit_exp = [];
    end

    % Plot
    h = figure; set(h,'Color','none'); ax = gca;
    set(ax,'Color','none','FontSize',12,'LineWidth',1.2);
    H = gobjects(0); L = {};

    % helper using the style map
    function add_series(y, key)
        st = styles.(key);
        good = isfinite(y) & y>0;
        if any(good)
            hl = semilogy(sNs, y, 'LineWidth',2, 'LineStyle',st.LineStyle, 'Marker',st.Marker, 'MarkerSize',st.MarkerSize, 'MarkerFaceColor',st.MarkerFaceColor, 'Color',st.Color);
            hold on;
            H(end+1) = hl; 
            L{end+1} = st.Legend; 
        end
    end

    % PLS, Diag
    add_series(el2_poly_s, 'PLS');
    add_series(el2_diag_s, 'Diag');

    % FS K_t series
    if include_all_Kt
        add_series(el2_fs1_s, 'FS1');
        add_series(el2_fs2_s, 'FS2');
        add_series(el2_fs3_s, 'FS3');
    else
        add_series(el2_fs1_s, 'FS1'); % shown as “FS” in legend text below if you prefer
    end

    % Overlays
    add_series(el2_fs_eff_s, 'FSeff');
    add_series(el2_fs_eval_s, 'FSeval');

    % PHS+poly (same style everywhere)
    add_series(el2_phs_s, 'PHS');

    % Trendline
    if ~isempty(yfit_exp)
        st = styles.Trend;
        hfit = semilogy(xx, yfit_exp, 'LineWidth',2, 'LineStyle',st.LineStyle, 'Color',st.Color);
        H(end+1) = hfit; L{end+1} = st.Legend;
    end

    ax.LineWidth = 2; ax.FontSize = 16;
    ax.XAxis.FontWeight = 'bold'; ax.YAxis.FontWeight = 'bold';
    if dim==1, xlabel('N','Interpreter','tex','FontSize',16,'FontWeight','bold');
    else,      xlabel(sprintf('N^{1/%d}', dim),'Interpreter','tex','FontSize',16,'FontWeight','bold'); end
    ylabel('Relative l_2 error','Interpreter','tex','FontSize',16,'FontWeight','bold');

    if ~isempty(H), legend(H,L,'Location','best','Interpreter','tex','FontWeight','bold','FontSize',14); end

    if sm==1, ttl = 'Relative l_2 error vs. N^{1/d}, C^2(R^3) Wendland Kernel';
    elseif sm==2, ttl = 'Relative l_2 error vs. N^{1/d}, C^4(R^3) Wendland Kernel';
    else ttl = 'Relative l_2 error vs. N^{1/d}, C^6(R^3) Wendland Kernel';
    end
    title(ttl,'Interpreter','tex','FontWeight','bold','FontSize',16);
    xlim([min(sNs) max(sNs)]); 
    ylim([err_min err_max]);


    if include_all_Kt, outname = sprintf('error_vs_N_fs_all_s%d.png', sm);
    else, outname = sprintf('error_vs_N_all_s%d.png', sm);
    end
    export_fig(gcf, fullfile(res.results_dir, outname), '-png','-r300','-transparent');
    close(h);

        % 1D FC panel that shows the same PHS curve as FS 
    if dim == 1
        h = figure; set(h,'Color','none'); ax = gca;
        set(ax,'Color','none','FontSize',12,'LineWidth',1.2);
        H = gobjects(0); L = {};

        % helper (re-use nested add_series)
        % PLS, Diag
        add_series(el2_poly_s, 'PLS');
        add_series(el2_diag_s, 'Diag');

        % FC family
        add_series(el2_fc1_s, 'FC1');
        if include_all_Kt
            add_series(el2_fc2_s, 'FC2');
            add_series(el2_fc3_s, 'FC3');
        end

        % PHS+poly — EXACT SAME CURVE as on the FS plot
        add_series(el2_phs_s, 'PHS');

        if ~isempty(yfit_exp)
            st = styles.Trend;
            hfit = semilogy(xx, yfit_exp, 'LineWidth',2, 'LineStyle',st.LineStyle, 'Color',st.Color);
            H(end+1) = hfit; L{end+1} = st.Legend;
        end

        ax.LineWidth = 2; ax.FontSize = 16;
        ax.XAxis.FontWeight = 'bold'; ax.YAxis.FontWeight = 'bold';
        xlabel('N','Interpreter','tex','FontSize',16,'FontWeight','bold');
        ylabel('Relative l_2 error','Interpreter','tex','FontSize',16,'FontWeight','bold');
        if ~isempty(H), legend(H,L,'Location','best','Interpreter','tex','FontWeight','bold','FontSize',14); end

        if sm==1, ttl = 'Relative l_2 error vs. N, C^2(R^3) Wendland Kernel';
        elseif sm==2, ttl = 'Relative l_2 error vs. N, C^4(R^3) Wendland Kernel';
        else ttl = 'Relative l_2 error vs. N, C^6(R^3) Wendland Kernel';
        end
        title(ttl,'Interpreter','tex','FontWeight','bold','FontSize',16);
        xlim([min(sNs) max(sNs)]);
        ylim([err_min err_max]);

        if include_all_Kt, outname_fc = sprintf('error_vs_N_fc_all_s%d.png', sm);
        else, outname_fc = sprintf('error_vs_N_fc_s%d.png', sm);
        end
        export_fig(gcf, fullfile(res.results_dir, outname_fc), '-png','-r300','-transparent');
        close(h);
    end
end


%% Style map (consistent across all figures)
function S = getPlotStyles(dim)
    % MATLAB default lines (for reference):
    % blue, orange, yellow, purple, green, light blue, etc.
    S.PLS = struct('Legend','PLS', 'Color',[0.000 0.447 0.741],'Marker','o','MarkerSize',6,'MarkerFaceColor','w','LineStyle','-');
    S.Diag = struct('Legend','Diag', 'Color',[0.850 0.325 0.098],'Marker','s','MarkerSize',6,'MarkerFaceColor','w','LineStyle','-');
    if dim==1
        S.FS1 = struct('Legend','K_t = 1e12','Color',[0.929 0.694 0.125],'Marker','^','MarkerSize',6,'MarkerFaceColor','w','LineStyle','-');
    else
        S.FS1 = struct('Legend','FS','Color',[0.929 0.694 0.125],'Marker','^','MarkerSize',6,'MarkerFaceColor','w','LineStyle','-');
    end
    S.FS2 = struct('Legend','K_t = 1e8','Color',[0.494 0.184 0.556],'Marker','x','MarkerSize',8,'MarkerFaceColor','none','LineStyle','-');
    S.FS3 = struct('Legend','K_t = 1e4','Color',[0.466 0.674 0.188],'Marker','+','MarkerSize',8,'MarkerFaceColor','none', 'LineStyle','-');
    S.FC1 = struct('Legend','K_t = 1e12','Color',[0.929 0.694 0.125],'Marker','^','MarkerSize',6,'MarkerFaceColor','w','LineStyle','-');
    S.FC2 = struct('Legend','K_t = 1e8','Color',[0.494 0.184 0.556],'Marker','x','MarkerSize',8,'MarkerFaceColor','none','LineStyle','-');
    S.FC3 = struct('Legend','K_t = 1e4','Color',[0.466 0.674 0.188],'Marker','+','MarkerSize',8,'MarkerFaceColor','none', 'LineStyle','-');

    % PHS+poly: make it the SAME everywhere (match your timing fig: light blue circle)
    S.PHS = struct('Legend','PHS+poly','Color',[0.000 0.75 1.000],'Marker','>', 'MarkerSize',6,'MarkerFaceColor','w','LineStyle','-');
    % Optional overlays
    S.FSeff = struct('Legend','FS-eff','Color',[0.300 0.300 0.300],'Marker','d','MarkerSize',6,'MarkerFaceColor','w','LineStyle','-');
    S.FSeval = struct('Legend','FS-eval','Color',[0.150 0.150 0.150],'Marker','v','MarkerSize',6,'MarkerFaceColor','w','LineStyle','-');
    % Trendline
    S.Trend = struct('Legend','exp trendline','Color',[.5 .5 .5],'Marker','none','MarkerSize',6,'MarkerFaceColor','none','LineStyle','--');
end