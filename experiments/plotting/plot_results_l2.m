%% Pick a function to run
function plot_results_l2(function_name, subdir, base_results_dir, smoothness_idx, variety)
    % Usage:
    %   plot_results_l2('abs_1d','high')
    %   plot_results_l2('xy_p_2d','high','results',1,false)  % 2D/3D: PLS/Diag/FS-only
    %
    % variety (logical, default = true):
    %   true  -> plot PLS, Diag, FS with K_t = 1e12, 1e8, 1e4 (like 1D panels)
    %   false -> plot PLS, Diag, FS(fs1) only (for 2D/3D)
    %
    % We also overlay FS-κeff and FS-eval by default in the main plot.
    % To show ONLY PLS/Diag/FS(fs1), comment the two overlay lines (marked).

    if nargin < 3 || isempty(base_results_dir), base_results_dir = fullfile('results'); end
    if nargin < 2 || isempty(subdir), subdir = 'high'; end
    if nargin < 4 || isempty(smoothness_idx), smoothness_idx = 1; end
    if nargin < 5 || isempty(variety), variety = true; end
    sm = max(1,min(3,smoothness_idx));

    % Load data
    res = loadResults(function_name, subdir, base_results_dir);
    res.sNs = res.sN;
    dim = res.dim;

    % Compute global x‐limits 
    xmin = min(res.sNs(:));
    xmax = max(res.sNs(:));

    % Decide which series will be included in the main plot
    include_all_Kt = variety || (dim==1); % 1D always “all K_t”

    % Build y-limits from only the series we’ll actually plot
    err_vals = [];
    if isfield(res,'el2_poly'), err_vals = [err_vals; res.el2_poly(:)]; end
    if isfield(res,'el2_diag'), err_vals = [err_vals; res.el2_diag(:,sm)]; end
        if isfield(res,'el2_fs1'), err_vals = [err_vals; res.el2_fs1(:,sm)]; end
    if include_all_Kt
        if isfield(res,'el2_fs2'), err_vals = [err_vals; res.el2_fs2(:,sm)]; end
        if isfield(res,'el2_fs3'), err_vals = [err_vals; res.el2_fs3(:,sm)]; end
    end

    % 1D FC variants (only dimension 1 has FC data)
    if dim == 1
        if isfield(res,'el2_fc1'), err_vals = [err_vals; res.el2_fc1(:,sm)]; end
        if include_all_Kt
            if isfield(res,'el2_fc2'), err_vals = [err_vals; res.el2_fc2(:,sm)]; end
            if isfield(res,'el2_fc3'), err_vals = [err_vals; res.el2_fc3(:,sm)]; end
        end
    end

    % FS-κeff / FS-eval overlays
    if isfield(res,'el2_fs_eff'),  err_vals = [err_vals; res.el2_fs_eff(:)];  end  % (comment to hide FS-κeff / FS-eval)
    if isfield(res,'el2_fs_eval'), err_vals = [err_vals; res.el2_fs_eval(:)]; end  % (comment to hide FS-κeff / FS-eval)

    err_vals = err_vals(isfinite(err_vals) & err_vals>0);
    if isempty(err_vals), err_min = 1e-16; err_max = 1; else
        err_min = 10^floor(log10(min(err_vals)));
        err_max = 10^ceil (log10(max(err_vals)));
    end

    % Main plot (PLS/Diag + FS, with optional overlays)
    plotErrorVsN(res, sm, xmin, xmax, err_min, err_max, include_all_Kt);
end


%% Load results files
function res = loadResults(function_name, subdir, base_results_dir)
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
    res = data;
end


%% Plot relative L2 error vs N^{1/d} (main panel)
function plotErrorVsN(res, sm, xmin, xmax, err_min, err_max, include_all_Kt)
    dim = res.dim;
    sNs = res.sNs(:);
    n = numel(sNs);

    % Series (fill with NaNs if missing so plotting/limits stay consistent)
    el2_poly_s = NaN(n,1); if isfield(res,'el2_poly'), el2_poly_s = res.el2_poly; end
    el2_diag_s = NaN(n,1); if isfield(res,'el2_diag'), el2_diag_s = res.el2_diag(:,sm); end

    el2_fs1_s = NaN(n,1); if isfield(res,'el2_fs1'), el2_fs1_s = res.el2_fs1(:,sm); end
    el2_fs2_s = NaN(n,1); if isfield(res,'el2_fs2'), el2_fs2_s = res.el2_fs2(:,sm); end
    el2_fs3_s = NaN(n,1); if isfield(res,'el2_fs3'), el2_fs3_s = res.el2_fs3(:,sm); end

    % 1D FC variants (leave as NaN for higher dims)
    el2_fc1_s = NaN(n,1); if isfield(res,'el2_fc1'), el2_fc1_s = res.el2_fc1(:,sm); end
    el2_fc2_s = NaN(n,1); if isfield(res,'el2_fc2'), el2_fc2_s = res.el2_fc2(:,sm); end
    el2_fc3_s = NaN(n,1); if isfield(res,'el2_fc3'), el2_fc3_s = res.el2_fc3(:,sm); end

    % Diagnostic-tuned FS variants
    el2_fs_eff_s  = NaN(n,1); if isfield(res,'el2_fs_eff'),  el2_fs_eff_s  = res.el2_fs_eff(:);  end  % (comment to hide overlay)
    el2_fs_eval_s = NaN(n,1); if isfield(res,'el2_fs_eval'), el2_fs_eval_s = res.el2_fs_eval(:); end  % (comment to hide overlay)


    % Exponential trendline: fit to FS1
    ok = ~(isnan(el2_fs1_s) | el2_fs1_s<=0);
    if any(ok)
        xx = sNs(ok); yy = el2_fs1_s(ok);
        coeff_exp = polyfit(xx, log(yy), 1);
        yfit_exp = exp(coeff_exp(2)) * exp( coeff_exp(1) * xx );
        fprintf('Exponential fit (FS1): E(N) = %.2f · exp(%.2f·N^{1/d})\n', exp(coeff_exp(2)), coeff_exp(1));
    else
        xx = []; yfit_exp = [];
        fprintf('Exponential fit (FS1): insufficient data.\n');
    end

    % Plot
    h = figure; set(h,'Color','none'); ax = gca;
    set(ax,'Color','none','FontSize',12,'LineWidth',1.2);
    marks = {'-o','-s','-^','-x','-+','-d','-v'}; 
    if dim==1
        marks = {'-o','-s','-^','-x','-+','-d','-v','->','-<','-p'};
    end
    H = gobjects(0); L = {};

    % helper
    function add_if_valid(y, mk, lab)
        good = isfinite(y) & y>0;
        if any(good)
            hl = semilogy(sNs, y, mk, 'LineWidth',2,'MarkerFaceColor','w','MarkerSize',6); hold on;
            H(end+1) = hl; L{end+1} = lab;
        end
    end

    % PLS, Diag
    add_if_valid(el2_poly_s, marks{1}, 'PLS');
    add_if_valid(el2_diag_s, marks{2}, 'Diag');

        % FS K_t series
    if include_all_Kt
        add_if_valid(el2_fs1_s, marks{3}, 'K_t = 1e12');
        add_if_valid(el2_fs2_s, marks{4}, 'K_t = 1e8');
        add_if_valid(el2_fs3_s, marks{5}, 'K_t = 1e4');
    else
        add_if_valid(el2_fs1_s, marks{3}, 'FS');
    end

    % 1D FC family (only meaningful when dim==1)
    if dim == 1
        if include_all_Kt
            add_if_valid(el2_fc1_s, marks{8}, 'FC, K_t = 1e12');
            add_if_valid(el2_fc2_s, marks{9}, 'FC, K_t = 1e8');
            add_if_valid(el2_fc3_s, marks{10}, 'FC, K_t = 1e4');
        else
            add_if_valid(el2_fc1_s, marks{8}, 'FC');
        end
    end

    % Overlays (comment to hide)
    add_if_valid(el2_fs_eff_s,  marks{6}, 'FS-eff');
    add_if_valid(el2_fs_eval_s, marks{7}, 'FS-eval');

    % Trendline
    if ~isempty(yfit_exp)
        hfit = semilogy(xx, yfit_exp, '--','LineWidth',2,'Color',[.5 .5 .5]);
        H(end+1) = hfit; L{end+1} = 'exp trendline';
    end

    ax.LineWidth = 2; ax.FontSize = 16;
    ax.XAxis.FontWeight = 'bold'; ax.YAxis.FontWeight = 'bold';

    if dim==1
        xlabel('N','Interpreter','tex','FontSize',16,'FontWeight','bold');
    else
        xlabel(sprintf('N^{1/%d}', dim),'Interpreter','tex','FontSize',16,'FontWeight','bold');
    end
    ylabel('Relative l_2 error','Interpreter','tex','FontSize',16,'FontWeight','bold');

    if ~isempty(H), legend(H,L,'Location','best','Interpreter','tex','FontWeight','bold','FontSize',14); end

    % Title + axes
    if sm==1, title('Relative l_2 error vs. N^{1/d}, C^2(R^3) Wendland Kernel','Interpreter','tex','FontWeight','bold','FontSize',16);
    elseif sm==2, title('Relative l_2 error vs. N^{1/d}, C^4(R^3) Wendland Kernel','Interpreter','tex','FontWeight','bold','FontSize',16);
    else, title('Relative l_2 error vs. N^{1/d}, C^6(R^3) Wendland Kernel','Interpreter','tex','FontWeight','bold','FontSize',16);
    end
    xlim([xmin xmax]); ylim([err_min err_max]);

    % Filename unchanged
    if include_all_Kt, outname = sprintf('error_vs_N_fs_s%d.png', sm);
    else, outname = sprintf('error_vs_N_combo_s%d.png', sm);
    end
    export_fig(gcf, fullfile(res.results_dir, outname), '-png','-r300','-transparent');
    close(h);
end