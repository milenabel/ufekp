function plot_compare_phs_vs_unified(function_name, subdir, base_results_dir, phs_p_list)
% Plot PLS / Diag / FS(mid) vs PHS+poly (p in phs_p_list).
%
% Usage:
%   plot_compare_phs_style('xy_p_2d'); % auto-detect PHS p in folder
%   plot_compare_phs_style('xy_p_2d','high','results',[5 7 9]);
%
% Writes to <results>/<function>/<subdir>/figs/:
%   - error_vs_N_compare.png
%   - assembly_time_vs_N_compare.png
%   - eval_time_vs_N_compare.png

if nargin < 3 || isempty(base_results_dir), base_results_dir = 'results'; end
if nargin < 2 || isempty(subdir), subdir = 'high'; end

% Accept a single name or a list; if a list, loop and return
if iscell(function_name)
    for k = 1:numel(function_name)
        plot_compare_phs_vs_unified(function_name{k}, subdir, base_results_dir, phs_p_list);
    end
    return;
end
% If string, make it char; if numeric, error
if isstring(function_name), function_name = char(function_name); end
if ~ischar(function_name) || isempty(function_name)
    error('function_name must be a nonempty char or string.');
end

% Load baseline unified results (has PLS/Diag/FS + sN/dim)
res = loadResults(function_name, subdir, base_results_dir);
results_dir = res.results_dir;
figs_dir = fullfile(results_dir,'figs');
if ~exist(figs_dir,'dir'), mkdir(figs_dir); end

% Find PHS result files (results_<fun>_phsp_<p>.mat)
if nargin < 4 || isempty(phs_p_list)
    listing = dir(fullfile(results_dir, sprintf('results_%s_phsp_*.mat', function_name)));
    phs_p_list = [];
    for j=1:numel(listing)
        t = regexp(listing(j).name,'phsp_(\d+)\.mat$','tokens','once');
        if ~isempty(t), phs_p_list(end+1) = str2double(t{1}); end %#ok<AGROW>
    end
    phs_p_list = unique(phs_p_list);
end

PHS = cell(1,numel(phs_p_list));
for j=1:numel(phs_p_list)
    fP = fullfile(results_dir, sprintf('results_%s_phsp_%d.mat', function_name, phs_p_list(j)));
    if exist(fP,'file'), PHS{j} = load(fP);
    else, warning('Missing PHS file: %s', fP); PHS{j} = []; end
end

% SERIES 
sNs = res.sN(:); n = numel(sNs); dim = res.dim;

% PLS
E_pls = getOrNaN(res,'el2_poly', n);
TA_pls = getOrNaN(res,'a_time_poly', n);
TE_pls = getOrNaN(res,'e_time_poly', n);

% Diag
E_diag  = getColOrNaN(res,'el2_diag',  1, n);
TA_diag = getColOrNaN(res,'a_time_diag',1, n);
TE_diag = getColOrNaN(res,'e_time_diag',1, n);

% FS(mid): prefer fs1 (Kt=1e12); else fs2 (1e8); else fs3 (1e4)
[fs_label, E_fs, TA_fs, TE_fs] = pickFSmid(res, n);

% PHS p-list (each aligned to baseline sN)
PHSlab = cell(1,numel(phs_p_list));
E_phs  = cell(1,numel(phs_p_list));
TA_phs = cell(1,numel(phs_p_list));
TE_phs = cell(1,numel(phs_p_list));
for j=1:numel(phs_p_list)
    PHSlab{j} = sprintf('PHS+poly (p=%d)', phs_p_list(j));
    if isempty(PHS{j})
        E_phs{j}  = NaN(n,1); TA_phs{j} = NaN(n,1); TE_phs{j} = NaN(n,1);
        continue;
    end
    % Map PHS (which may have shorter sN) onto baseline length n
    E_phs{j}  = padToN(getVecOrNaN(PHS{j},'el2_phs'),  n);
    TA_phs{j} = padToN(getVecOrNaN(PHS{j},'a_time_phs'),n);
    TE_phs{j} = padToN(getVecOrNaN(PHS{j},'e_time_phs'),n);
end

% 3D: keep only 5 points across ALL series
if dim == 3
    keep = 1:min(5, numel(sNs));
    sNs    = sNs(keep);

    E_pls = E_pls(keep); TA_pls = TA_pls(keep); TE_pls = TE_pls(keep);
    E_diag = E_diag(keep); TA_diag= TA_diag(keep); TE_diag= TE_diag(keep);
    E_fs = E_fs(keep); TA_fs  = TA_fs(keep); TE_fs = TE_fs(keep);

    for j = 1:numel(phs_p_list)
        if ~isempty(E_phs{j}), E_phs{j} = E_phs{j}(keep); end
        if ~isempty(TA_phs{j}), TA_phs{j} = TA_phs{j}(keep); end
        if ~isempty(TE_phs{j}), TE_phs{j} = TE_phs{j}(keep); end
    end
end

% Global y-lims (log decades) *after* trimming
[e_ymin,e_ymax] = decade_limits_smart([E_pls;E_diag;E_fs;cat(1,E_phs{:})]);
[ta_ymin,ta_ymax] = decade_limits_smart([TA_pls;TA_diag;TA_fs;cat(1,TA_phs{:})]);
[te_ymin,te_ymax] = decade_limits_smart([TE_pls;TE_diag;TE_fs;cat(1,TE_phs{:})]);

xmin = min(sNs); xmax = max(sNs);


% Marks/style (match plot_efficiency)
marks = {'-o','-s','-x','-+','-d','-^'};  % PLS, Diag, FS, PHS p=5,7,9...
lw = 2;

% Plot 1: Relative l2 error
h = figure; set(h,'Color','none'); ax = gca;
set(ax,'Color','none','FontSize',12,'LineWidth',1.2);
semilogy(sNs, E_pls,  marks{1}, 'LineWidth',lw); hold on;
semilogy(sNs, E_diag, marks{2}, 'LineWidth',lw);
semilogy(sNs, E_fs,   marks{3}, 'LineWidth',lw);
for j=1:numel(phs_p_list)
    mIdx = min(3+j, numel(marks));
    semilogy(sNs, E_phs{j}, marks{mIdx}, 'LineWidth',lw);
end
style_axes(ax, dim, 'Relative l_2 error');
xlim([xmin xmax]); ylim([e_ymin e_ymax]);
leg = [{'PLS','Diag','FS'}, PHSlab];
legend(leg, 'Location','best','Interpreter','tex','FontWeight','bold','FontSize',14);
if dim==1
    title(sprintf('Relative l_2 error vs. N'), 'Interpreter','tex','FontWeight','bold','FontSize',16);
else
    title(sprintf('Relative l_2 error vs. N^{1/%d}', dim), 'Interpreter','tex','FontWeight','bold','FontSize',16);
end
export_fig(gcf, fullfile(figs_dir, 'error_vs_N_compare.png'), '-png','-r300','-transparent');
close(h);

% Plot 2: Assembly+Solve
h = figure; set(h,'Color','none'); ax = gca;
set(ax,'Color','none','FontSize',12,'LineWidth',1.2);
semilogy(sNs, TA_pls,  marks{1}, 'LineWidth',lw); hold on;
semilogy(sNs, TA_diag, marks{2}, 'LineWidth',lw);
semilogy(sNs, TA_fs,   marks{3}, 'LineWidth',lw);
for j=1:numel(phs_p_list)
    mIdx = min(3+j, numel(marks));
    semilogy(sNs, TA_phs{j}, marks{mIdx}, 'LineWidth',lw);
end
style_axes(ax, dim, 'Assembly and Solve time (s)');
xlim([xmin xmax]); ylim([ta_ymin ta_ymax]);
legend(leg, 'Location','best','Interpreter','tex','FontWeight','bold','FontSize',14);
title('Assembly and Solve time', 'Interpreter','tex','FontWeight','bold','FontSize',16);
export_fig(gcf, fullfile(figs_dir, 'assembly_time_vs_N_compare.png'), '-png','-r300','-transparent');
close(h);

% Plot 3: Evaluation time
h = figure; set(h,'Color','none'); ax = gca;
set(ax,'Color','none','FontSize',12,'LineWidth',1.2);
semilogy(sNs, TE_pls,  marks{1}, 'LineWidth',lw); hold on;
semilogy(sNs, TE_diag, marks{2}, 'LineWidth',lw);
semilogy(sNs, TE_fs,   marks{3}, 'LineWidth',lw);
for j=1:numel(phs_p_list)
    mIdx = min(3+j, numel(marks));
    semilogy(sNs, TE_phs{j}, marks{mIdx}, 'LineWidth',lw);
end
style_axes(ax, dim, 'Evaluation time (s)');
xlim([xmin xmax]); ylim([te_ymin te_ymax]);
legend(leg, 'Location','best','Interpreter','tex','FontWeight','bold','FontSize',14);
title('Evaluation time', 'Interpreter','tex','FontWeight','bold','FontSize',16);
export_fig(gcf, fullfile(figs_dir, 'eval_time_vs_N_compare.png'), '-png','-r300','-transparent');
close(h);

fprintf('Saved figs to %s\n', figs_dir);
end

% helpers
function v = getOrNaN(S, fld, n)
if isfield(S,fld), v = S.(fld)(:,1); else, v = NaN(n,1); end
end

function v = getVecOrNaN(S, fld)
if isfield(S,fld), v = S.(fld)(:); else, v = []; end
end

function v = padToN(v, n)
v = v(:);
if isempty(v), v = NaN(n,1); return; end
m = numel(v);
if m >= n, v = v(1:n);
else, v = [v; NaN(n-m,1)];
end
end

function v = getColOrNaN(S, fld, j, n)
if isfield(S,fld)
    V = S.(fld);
    if size(V,2) >= j, v = V(:,j);
    else, v = NaN(n,1);
    end
else
    v = NaN(n,1);
end
end

function [lab, E, TA, TE] = pickFSmid(S, n)
% Prefer fs1; fallback fs2; fallback fs3; else NaNs.
if all(isfield(S,{'el2_fs1','a_time_fs1','e_time_fs1'}))
    lab = 'FS (K_t = 1e12)';
    E  = S.el2_fs1(:,1); TA = S.a_time_fs1(:,1); TE = S.e_time_fs1(:,1);
elseif all(isfield(S,{'el2_fs2','a_time_fs2','e_time_fs2'}))
    lab = 'FS (K_t = 1e8)';
    E  = S.el2_fs2(:,1); TA = S.a_time_fs2(:,1); TE = S.e_time_fs2(:,1);
elseif all(isfield(S,{'el2_fs3','a_time_fs3','e_time_fs3'}))
    lab = 'FS (K_t = 1e4)';
    E  = S.el2_fs3(:,1); TA = S.a_time_fs3(:,1); TE = S.e_time_fs3(:,1);
else
    lab = 'FS';
    E = NaN(n,1); TA = NaN(n,1); TE = NaN(n,1);
end
end

function style_axes(ax, d, ylab)
ax.LineWidth = 2; ax.FontSize = 16;
ax.XAxis.FontWeight = 'bold'; ax.YAxis.FontWeight = 'bold';
if d==1
    xlabel('N','Interpreter','tex','FontSize',16,'FontWeight','bold');
else
    xlabel(sprintf('N^{1/%d}', d),'Interpreter','tex','FontSize',16,'FontWeight','bold');
end
ylabel(ylab,'Interpreter','tex','FontSize',16,'FontWeight','bold');
end

function [ymin, ymax] = decade_limits_smart(v)
v = v(isfinite(v) & v>0);
if isempty(v), ymin = 1e-6; ymax = 1; return; end
ymin = 10^floor(log10(min(v)));
ymax = 10^ceil (log10(max(v)));
if ymin==ymax, ymin = ymin/10; ymax = ymax*10; end
end

function res = loadResults(function_name, subdir, base_results_dir)
valid_subdirs = {'high','low','fixed','mixed'};
if ~ismember(subdir, valid_subdirs)
    error('Invalid subdirectory: %s. Must be one of: high, low, fixed, mixed', subdir);
end
results_dir = fullfile(base_results_dir, function_name, subdir);
matfile = fullfile(results_dir, sprintf('results_%s.mat', function_name));
if ~exist(matfile,'file'), error('Results file not found: %s', matfile); end
data = load(matfile);
data.results_dir = results_dir;
res = data;
end

