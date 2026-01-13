function plot_l2_manifold(subdir, base_results_dir)
% plot_l2_manifolds('high','results')
% For each manifold in {'boundary_test','sphere','torus'} plot:
%   PLS (el2_poly), FS (el2), PHS+poly (el2_phs) vs N^{1/2} + FS exp trend.
% Saves error_vs_N_<manifold>.png into that manifold's <results_dir>.

if nargin < 1 || isempty(subdir), subdir = 'high';  end
if nargin < 2 || isempty(base_results_dir), base_results_dir = 'results'; end

manifolds = {'boundary_test','sphere','torus'};
for k = 1:numel(manifolds)
    fn = manifolds{k};
    try
        R = loadResults_with_PHS_manifold(fn, subdir, base_results_dir);
        plotErrorVsN_manifold(R);
    catch ME
        warning('Skipping %s: %s', fn, ME.message);
    end
end
end


%% Load results (FS+PLS from any base file; attach PHS from same/sibling file)
function res = loadResults_with_PHS_manifold(manifold, subdir, base_results_dir)
results_dir = fullfile(base_results_dir, manifold, subdir);
pat = fullfile(results_dir, sprintf('results_%s*.mat', manifold));
D = dir(pat);
if isempty(D)
    error('No results files matching %s', pat);
end

% Pick a base file that has sN/Ns + el2 + el2_poly
baseIdx = NaN;
for j = 1:numel(D)
    S = load(fullfile(results_dir, D(j).name));
    if (isfield(S,'sN') || isfield(S,'Ns')) && isfield(S,'el2') && isfield(S,'el2_poly')
        baseIdx = j; break
    end
end
if isnan(baseIdx)
    error('No MAT in %s has required fields sN/Ns + el2 + el2_poly.', results_dir);
end
Sbase = load(fullfile(results_dir, D(baseIdx).name));

% sN from sN or sqrt(Ns)
if isfield(Sbase,'sN') && ~isempty(Sbase.sN)
    sN = Sbase.sN(:);
elseif isfield(Sbase,'Ns') && ~isempty(Sbase.Ns)
    sN = sqrt(Sbase.Ns(:));
else
    error('Base file lacks sN and Ns: %s', D(baseIdx).name);
end

% FS + PLS from base
res.results_dir = results_dir;
res.function_name = manifold;
res.sN = sN(:);
res.el2 = getv(Sbase,'el2');
res.el2_poly = getv(Sbase,'el2_poly');
srcFS = D(baseIdx).name; 
srcPLS = D(baseIdx).name;

% PHS: prefer base; otherwise search other files
res.el2_phs = getv(Sbase,'el2_phs'); 
srcPHS = D(baseIdx).name;
if all(~isfinite(res.el2_phs))
    found = false;
    for j = 1:numel(D)
        if j == baseIdx, continue; end
        S = load(fullfile(results_dir, D(j).name));
        if isfield(S,'el2_phs') && ~isempty(S.el2_phs)
            res.el2_phs = S.el2_phs(:);
            srcPHS = D(j).name;
            found = true;
            break
        end
    end
    if ~found, res.el2_phs = NaN(size(res.sN)); end
end

% pad/trim all to |sN|
n = numel(res.sN);
res.el2 = padToN(res.el2, n);
res.el2_poly = padToN(res.el2_poly, n);
res.el2_phs = padToN(res.el2_phs, n);

fprintf('[plot] %s\n   FS : %s\n   PLS: %s\n   PHS: %s\n', manifold, srcFS, srcPLS, srcPHS);
fprintf('finite pts — FS:%d  PLS:%d  PHS:%d\n', nnz(isfinite(res.el2)), nnz(isfinite(res.el2_poly)), nnz(isfinite(res.el2_phs)));
end

function v = getv(S, fld)
if isfield(S,fld), v = S.(fld)(:); else, v = []; end
end

function v = padToN(v, n)
v = v(:); 
m = numel(v);
if m==0, v = NaN(n,1); return; end
if m<n, v = [v; NaN(n-m,1)]; else, v = v(1:n); end
end


%% Plot relative L2 error vs N^{1/d} (manifold panel): FS, PLS, PHS+poly + trend
function plotErrorVsN_manifold(res)
    % styles
    cPLS = [0.000 0.447 0.741]; mPLS='o';
    cFS = [0.850 0.325 0.098]; mFS='^';
    cPHS = [0.000 0.75  1.000]; mPHS='>';
    cTr = [0.5 0.5 0.5];
    ms = 7;                
    lw = 2.2;              

    x = res.sN(:);
    yFS = res.el2(:);
    yPLS = res.el2_poly(:);
    yPHS = res.el2_phs(:);

    % log-safe (drop non-positive/inf)
    yFS (yFS <=0 | ~isfinite(yFS )) = NaN;
    yPLS(yPLS<=0 | ~isfinite(yPLS)) = NaN;
    yPHS(yPHS<=0 | ~isfinite(yPHS)) = NaN;

    % FS exponential trendline
    ok = isfinite(x) & isfinite(yFS);
    if nnz(ok) >= 2
        xx = x(ok); yy = yFS(ok);
        p  = polyfit(xx, log(yy), 1);
        yfit = exp(p(2)) * exp(p(1) * xx);
    else
        p = [NaN,NaN]; yfit = [];
    end

    % figure/axes
    h = figure; set(h,'Color','none','InvertHardcopy','off');
    ax = gca;
    set(ax,'Color','none','LineWidth',2,'FontSize',16,'FontName','Helvetica');
    ax.FontWeight = 'bold';        
    ax.XAxis.FontWeight = 'bold';
    ax.YAxis.FontWeight = 'bold';
    set(ax,'XScale','linear','YScale','log'); % x = N^{1/2} (linear), y = error (log)
    ax.TickLength = [0.02 0.02];
    ax.Box = 'on';

    hold on;
    H = gobjects(0); L = {};

    % series 
    if any(isfinite(yFS))
        H(end+1) = semilogy(x, yFS, '-', 'Color', cFS,  'Marker', mFS, 'MarkerSize', ms, 'MarkerFaceColor', 'w', 'LineWidth', lw);  L{end+1}='FS';
    end
    if any(isfinite(yPLS))
        H(end+1) = semilogy(x, yPLS,'-', 'Color', cPLS, 'Marker', mPLS, 'MarkerSize', ms, 'MarkerFaceColor', 'w', 'LineWidth', lw);  L{end+1}='PLS';
    end
    if any(isfinite(yPHS))
        H(end+1) = semilogy(x, yPHS,'-', 'Color', cPHS, 'Marker', mPHS, 'MarkerSize', ms, 'MarkerFaceColor', 'w', 'LineWidth', lw);  L{end+1}='PHS+poly';
    end
    if ~isempty(yfit)
        H(end+1) = semilogy(xx, yfit, '--', 'Color', cTr, 'LineWidth', 2.0); L{end+1}='exp trendline';
    end

    % labels/title
    xlabel('N^{1/2}','Interpreter','tex','FontSize',16,'FontWeight','bold');
    ylabel('Relative l_2 error','Interpreter','tex','FontSize',16,'FontWeight','bold');
    title('Relative l_2 error vs. N^{1/d}, C^2(R^3) Wendland Kernel', 'Interpreter','tex','FontSize',16,'FontWeight','bold');

    % legend
    if ~isempty(H)
        lgd = legend(H,L,'Location','best','Interpreter','tex');
        set(lgd,'FontSize',14,'FontWeight','bold','Box','on','LineWidth',1.5,'Color','w');
    end

    % limits
    xgood = x(isfinite(x)); if isempty(xgood), xgood = [1;1]; end
    vals = [yFS; yPLS; yPHS]; vals = vals(isfinite(vals));
    if isempty(vals), ymin = 1e-16; ymax = 1;
    else, ymin = 10^floor(log10(min(vals))); ymax = 10^ceil(log10(max(vals)));
    end
    xlim([min(xgood) max(xgood)]); ylim([ymin ymax]);

    % export
    outpng = fullfile(res.results_dir, sprintf('error_vs_N_%s.png', res.function_name));
    export_fig(gcf, outpng, '-png','-r300','-transparent'); close(h);

    if all(isfinite(p))
        fprintf('[%s] FS exp fit: E(N) = %.2f * exp(%.2f * N^{1/d}) -> %s\n', res.function_name, exp(p(2)), p(1), outpng);
    else
        fprintf('[%s] FS exp fit skipped -> %s\n', res.function_name, outpng);
    end
end