%% plot_efficiency_wphs.m — timings with optional PHS+poly (all dims)
function plot_efficiency_wphs(function_name, subdir, base_results_dir, smoothness_idx, variety)
% Usage:
%   plot_efficiency_wphs('abs_1d','high')
%   plot_efficiency_wphs('xy_p_2d','high','results',1,false)
%
% If a PHS p=9 file exists, include “PHS+poly” in FS/COMBO timing panels
% in ANY dimension (1D/2D/3D).

if nargin < 3 || isempty(base_results_dir), base_results_dir = 'results'; end
if nargin < 2 || isempty(subdir), subdir = 'high'; end
if nargin < 4 || isempty(smoothness_idx), smoothness_idx = 1; end
if nargin < 5 || isempty(variety), variety = true; end
sm = max(1, min(3, smoothness_idx));  % 1=C^2, 2=C^4, 3=C^6

% load unified (+ optional PHS p=9)
res = loadResults_with_PHS(function_name, subdir, base_results_dir);
res.sNs = res.sN;
results_dir = res.results_dir;

% limits (consider PHS if present)
xmin = min(res.sNs(:));
xmax = max(res.sNs(:));
[at_min, at_max] = globalYLim_times_with_PHS(res, sm, 'assembly');
[et_min, et_max] = globalYLim_times_with_PHS(res, sm, 'eval');

% plotting
if variety
    % FS panels (PLS + Diag + FS) (+ PHS if available)
    plotAssemblyOrEval_FSFC_with_PHS(res, sm, xmin, xmax, at_min, at_max, 'assembly','fs');
    plotAssemblyOrEval_FSFC_with_PHS(res, sm, xmin, xmax, et_min, et_max, 'eval','fs');

    % FC panels (unchanged, no PHS)
    plotAssemblyOrEval_FSFC_with_PHS(res, sm, xmin, xmax, at_min, at_max, 'assembly','fc');
    plotAssemblyOrEval_FSFC_with_PHS(res, sm, xmin, xmax, et_min, et_max, 'eval','fc');
else
    % one combined panel: PLS, Diag, FS(fs1) (+ PHS if available), FC(vs1)
    plotAssemblyOrEval_FSFC_with_PHS(res, sm, xmin, xmax, at_min, at_max, 'assembly','combo');
    plotAssemblyOrEval_FSFC_with_PHS(res, sm, xmin, xmax, et_min, et_max, 'eval','combo');
end

% CSV (unchanged)
T = build_efficiency_table_long(res, sm);
writetable(T, fullfile(results_dir,sprintf('efficiency_table_%d.csv', sm)));
fprintf('Wrote: %s\n', fullfile(results_dir,sprintf('efficiency_table_%d.csv', sm)));
end


%% unified loader (+ optional PHS p=9)
function res = loadResults_with_PHS(function_name, subdir, base_results_dir)
valid_subdirs = {'high','low','fixed','mixed'};
if ~ismember(subdir, valid_subdirs)
    error('Invalid subdirectory: %s. Must be one of: high, low, fixed, mixed', subdir);
end
results_dir = fullfile(base_results_dir, function_name, subdir);
matfile     = fullfile(results_dir, sprintf('results_%s.mat', function_name));
if ~exist(matfile,'file'), error('Results file not found: %s', matfile); end
data = load(matfile);
data.results_dir = results_dir;

% optional PHS p=9
phsfile = fullfile(results_dir, sprintf('results_%s_phsp_9.mat', function_name));
if exist(phsfile,'file')
    Phs = load(phsfile);
    n  = numel(data.sN(:));
    data.el2_phs    = padToN(getVecOrNaN(Phs,'el2_phs'),    n);
    data.a_time_phs = padToN(getVecOrNaN(Phs,'a_time_phs'), n);
    data.e_time_phs = padToN(getVecOrNaN(Phs,'e_time_phs'), n);
    fprintf('[eff] PHS p=9 attached: %s\n', phsfile);
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


%% Plotting (with optional PHS)
function plotAssemblyOrEval_FSFC_with_PHS(res, sm, xmin, xmax, y_min, y_max, kind, which)
    dim  = res.dim; 
    sNs = res.sNs(:);
    n = numel(sNs);

    styles = getPlotStyles(dim);  

    % poly / diag
    if strcmpi(kind,'assembly')
        polyv = isfield(res,'a_time_poly') * res.a_time_poly(:); if isempty(polyv), polyv = NaN(n,1); end
        diagv = isfield(res,'a_time_diag') * res.a_time_diag(:,sm); if isempty(diagv), diagv = NaN(n,1); end
    else
        polyv = isfield(res,'e_time_poly') * res.e_time_poly(:); if isempty(polyv), polyv = NaN(n,1); end
        diagv = isfield(res,'e_time_diag') * res.e_time_diag(:,sm); if isempty(diagv), diagv = NaN(n,1); end
    end

    % FS/FC/COMBO series
    if strcmpi(which,'fs') || strcmpi(which,'fc')
        labtag = iff(strcmpi(which,'fs'),'fs','vs');
        S = cell(1,3);
        for j=1:3
            tag = sprintf('%s%d', labtag, j);
            fld = sprintf('%s_time_%s', lower(kind(1)), tag);  % a_time_* or e_time_*
            if isfield(res, fld), S{j} = res.(fld)(:,sm); else, S{j} = NaN(n,1); end
        end
    else
        % combo: FS := fs1, FC := vs1
        if strcmpi(kind,'assembly')
            FS1 = isfield(res,'a_time_fs1') * res.a_time_fs1(:,sm); if isempty(FS1), FS1 = NaN(n,1); end
            VS1 = isfield(res,'a_time_vs1') * res.a_time_vs1(:,sm); if isempty(VS1), VS1 = NaN(n,1); end
        else
            FS1 = isfield(res,'e_time_fs1') * res.e_time_fs1(:,sm); if isempty(FS1), FS1 = NaN(n,1); end
            VS1 = isfield(res,'e_time_vs1') * res.e_time_vs1(:,sm); if isempty(VS1), VS1 = NaN(n,1); end
        end
    end

    % PHS timings (any dim)
    if strcmpi(kind,'assembly')
        y_phs = isfield(res,'a_time_phs') * res.a_time_phs(:); if isempty(y_phs), y_phs = NaN(n,1); end
    else
        y_phs = isfield(res,'e_time_phs') * res.e_time_phs(:); if isempty(y_phs), y_phs = NaN(n,1); end
    end

    % 3D trim: cut all series to the last valid PHS point and to 5 pts ---
    if dim == 3
        if any(isfinite(y_phs) & y_phs > 0)
            lastPHS = find(isfinite(y_phs) & y_phs > 0, 1, 'last');
        else
            lastPHS = n;  % no PHS timings: fall back to what we have
        end
        keep = 1:min(5, lastPHS);

        % apply trimming to x and every series we might plot
        sNs   = sNs(keep);
        polyv = polyv(keep);
        diagv = diagv(keep);
        if exist('S','var')     % fs/fc panels
            for j = 1:3, S{j} = S{j}(keep); end
        else                     % combo panel
            FS1 = FS1(keep);
            VS1 = VS1(keep);
        end
        y_phs = y_phs(keep);
        n = numel(sNs);  % keep n consistent after trimming
    end

    % figure
    h = figure; set(h,'Color','none'); ax = gca;
    set(ax,'Color','none','FontSize',12,'LineWidth',1.2);
    H = gobjects(0); L = {};

    % helper using the shared style map
    function add_series(y, key)
        st = styles.(key);
        g = isfinite(y) & y>0;
        if any(g)
            hl = semilogy(sNs, y, ...
                'LineWidth',2, ...
                'LineStyle',st.LineStyle, ...
                'Marker',st.Marker, ...
                'MarkerSize',st.MarkerSize, ...
                'MarkerFaceColor',st.MarkerFaceColor, ...
                'Color',st.Color); hold on;
            H(end+1) = hl; L{end+1} = st.Legend;
        end
    end

    % draw
    add_series(polyv, 'PLS');
    add_series(diagv, 'Diag');

    if strcmpi(which,'combo')
        add_series(FS1, 'FS1');   
    else
        add_series(S{1}, 'FS1');  % K_t = 1e12
        add_series(S{2}, 'FS2');  % K_t = 1e8
        add_series(S{3}, 'FS3');  % K_t = 1e4
    end

    % PHS+poly — SAME style in all plots
    add_series(y_phs, 'PHS');

    % axes/labels
    ax.LineWidth = 2; ax.FontSize = 16;
    ax.XAxis.FontWeight = 'bold'; ax.YAxis.FontWeight = 'bold';
    if dim==1, xlabel('N','Interpreter','tex','FontSize',16,'FontWeight','bold');
    else,      xlabel(sprintf('N^{1/%d}', dim),'Interpreter','tex','FontSize',16,'FontWeight','bold'); end
    if strcmpi(kind,'assembly')
        ylabel('Assembly and Solve time (s)','Interpreter','tex','FontSize',16,'FontWeight','bold');
        ttl = 'Assembly and Solve time';
    else
        ylabel('Evaluation time (s)','Interpreter','tex','FontSize',16,'FontWeight','bold');
        ttl = 'Evaluation time';
    end
    title(texTitle(ttl, sm), 'Interpreter','tex','FontWeight','bold','FontSize',16);

    if ~isempty(H), legend(H,L,'Location','best','Interpreter','tex','FontWeight','bold','FontSize',14); end
    xlim([min(sNs) max(sNs)]); 
    ylim([y_min y_max]);

    if strcmpi(which,'combo')
        if strcmpi(kind,'assembly'), outname = sprintf('assembly_time_vs_N_all_s%d.png', sm);
        else,                        outname = sprintf('eval_time_vs_N_all_s%d.png', sm);
        end
    else
        if strcmpi(kind,'assembly'), outname = sprintf('assembly_time_vs_N_%s_all_s%d.png', which, sm);
        else,                        outname = sprintf('eval_time_vs_N_%s_all_s%d.png', which, sm);
        end
    end
    export_fig(gcf, fullfile(res.results_dir, outname), '-png','-r300','-transparent');
    close(h);
end

% tiny inline helper
function out = iff(cond,a,b), if cond, out=a; else, out=b; end, end


%% Y-lims (consider PHS if present)
function [ymin, ymax] = globalYLim_times_with_PHS(res, sm, kind)
vals = [];
switch lower(kind)
    case 'assembly'
        if isfield(res,'a_time_poly'), vals = [vals; res.a_time_poly(:)]; end
        if isfield(res,'a_time_diag'), vals = [vals; res.a_time_diag(:,sm)]; end
        for tag = ["fs1","fs2","fs3","vs1","vs2","vs3"]
            f = sprintf('a_time_%s', tag);
            if isfield(res,f), vals = [vals; res.(f)(:,sm)]; end
        end
        % NEW
        if isfield(res,'a_time_phs'), vals = [vals; res.a_time_phs(:)]; end
    case 'eval'
        if isfield(res,'e_time_poly'), vals = [vals; res.e_time_poly(:)]; end
        if isfield(res,'e_time_diag'), vals = [vals; res.e_time_diag(:,sm)]; end
        for tag = ["fs1","fs2","fs3","vs1","vs2","vs3"]
            f = sprintf('e_time_%s', tag);
            if isfield(res,f), vals = [vals; res.(f)(:,sm)]; end
        end
        % NEW
        if isfield(res,'e_time_phs'), vals = [vals; res.e_time_phs(:)]; end
end
vals = vals(isfinite(vals) & vals>0);
if isempty(vals), ymin=1e-3; ymax=1; return; end
ymin = 10^floor(log10(min(vals)));
ymax = 10^ceil (log10(max(vals)));
end

function t = texTitle(base, sm)
switch sm
    case 1, kern = 'C^2(R^3)';
    case 2, kern = 'C^4(R^3)';
    case 3, kern = 'C^6(R^3)';
end
t = sprintf('%s, %s Wendland Kernel', base, kern);
end


%% CSV Table (unchanged from your version)
function T = build_efficiency_table_long(res, sm)
N1 = res.sNs(:);
if isfield(res,'sN') && isfield(res,'dim')
    N = round(res.sN(:).^res.dim);
else
    N = (1:numel(res.sNs(:))).';
end

mth = {}; Kt = []; Asm = []; Evl = []; Sp = []; Cd = []; Ncol = []; N1col = [];

    function add(methodName, Ktval, a, e, sp, cd)
        n = numel(N1);
        if numel(a) < n, a  = [a;  nan(n-numel(a),1)]; end
        if numel(e) < n, e  = [e;  nan(n-numel(e),1)]; end
        if numel(sp)< n, sp = [sp; nan(n-numel(sp),1)]; end
        if numel(cd)< n, cd = [cd; nan(n-numel(cd),1)]; end
        mth = [mth; repmat({methodName}, n, 1)];
        Kt  = [Kt;  repmat(Ktval, n, 1)];
        Asm = [Asm; a];
        Evl = [Evl; e];
        Sp  = [Sp;  sp];
        Cd  = [Cd;  cd];
        Ncol= [Ncol; N];
        N1col=[N1col;N1];
    end

% PLS
if isfield(res,'a_time_poly') && isfield(res,'e_time_poly')
    add('PLS', NaN, res.a_time_poly(:), res.e_time_poly(:), NaN(size(N1)), NaN(size(N1)));
end

% Diag
if isfield(res,'a_time_diag') && isfield(res,'e_time_diag')
    add('Diag', NaN, res.a_time_diag(:,sm), res.e_time_diag(:,sm), NaN(size(N1)), NaN(size(N1)));
end

% FS (Kt = 1e12, 1e8, 1e4)
KtVals = [1e12, 1e8, 1e4];
for j=1:3
    tag = sprintf('fs%d', j);
    a_fld = sprintf('a_time_%s', tag);
    e_fld = sprintf('e_time_%s', tag);
    s_fld = sprintf('sparsity_%s', tag);
    c_fld = sprintf('cond_%s', tag);
    if isfield(res, a_fld), a = res.(a_fld)(:,sm); else, a = NaN(numel(N1),1); end
    if isfield(res, e_fld), e = res.(e_fld)(:,sm); else, e = NaN(numel(N1),1); end
    if isfield(res, s_fld), sp = res.(s_fld)(:,sm); else, sp = NaN(numel(N1),1); end
    if isfield(res, c_fld), cd = res.(c_fld)(:,sm); else, cd = NaN(numel(N1),1); end
    if ~all(isnan(a)) || ~all(isnan(e))
        add('FS', KtVals(j), a, e, sp, cd);
    end
end

% FC (Kt = 1e12, 1e8, 1e4)
for j=1:3
    tag = sprintf('vs%d', j);
    a_fld = sprintf('a_time_%s', tag);
    e_fld = sprintf('e_time_%s', tag);
    s_fld = sprintf('sparsity_%s', tag);
    c_fld = sprintf('cond_%s', tag);
    hasA = isfield(res, a_fld);
    hasE = isfield(res, e_fld);
    if ~(hasA || hasE), continue; end
    if hasA, a = res.(a_fld)(:,sm); else, a = NaN(numel(N1),1); end
    if hasE, e = res.(e_fld)(:,sm); else, e = NaN(numel(N1),1); end
    if isfield(res, s_fld), sp = res.(s_fld)(:,sm); else, sp = NaN(numel(N1),1); end
    if isfield(res, c_fld), cd = res.(c_fld)(:,sm); else, cd = NaN(numel(N1),1); end
    if ~all(isnan(a)) || ~all(isnan(e))
        add('FC', KtVals(j), a, e, sp, cd);
    end
end

T = table(Ncol, N1col, string(mth), Kt, Asm, Evl, Sp, Cd, ...
    'VariableNames',{'N','N1overd','Method','Kt','AsmSolve_s','Eval_s','Sparsity','Cond'});

% sort nicely: by N, then Method (PLS,Diag,FS,FC), then Kt (NaN last)
[~, ordMethod] = ismember(T.Method, ["PLS","Diag","FS","FC"]);
ordMethod(ordMethod==0) = 99;
[~, order] = sortrows([T.N, ordMethod, isnan(T.Kt), T.Kt], [1 2 -3 4]);
T = T(order,:);
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
    S.FS2 = struct('Legend','K_t = 1e8', 'Color',[0.494 0.184 0.556],'Marker','x','MarkerSize',8,'MarkerFaceColor','none','LineStyle','-');
    S.FS3 = struct('Legend','K_t = 1e4', 'Color',[0.466 0.674 0.188],'Marker','+','MarkerSize',8,'MarkerFaceColor','none', 'LineStyle','-');
    % PHS+poly: make it the SAME everywhere (match your timing fig: light blue circle)
    S.PHS = struct('Legend','PHS+poly',  'Color',[0.000 0.75 1.000],'Marker','>', 'MarkerSize',6,'MarkerFaceColor','w','LineStyle','-');
    % Optional overlays
    S.FSeff = struct('Legend','FS-eff', 'Color',[0.300 0.300 0.300],'Marker','d','MarkerSize',6,'MarkerFaceColor','w','LineStyle','-');
    S.FSeval = struct('Legend','FS-eval', 'Color',[0.150 0.150 0.150],'Marker','v','MarkerSize',6,'MarkerFaceColor','w','LineStyle','-');
    % Trendline
    S.Trend = struct('Legend','exp trendline','Color',[.5 .5 .5],'Marker','none','MarkerSize',6,'MarkerFaceColor','none','LineStyle','--');
end