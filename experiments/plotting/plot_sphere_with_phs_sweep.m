function plot_sphere_with_phs_sweep(subdir, base_results_dir)
% Plot FS, PLS, and ALL available PHS+poly (different degrees p) for SPHERE.
% Scans the folder for results_sphere_phsp_*.mat and overlays their curves.
%
% Usage:
%   plot_sphere_with_phs_sweep('high','results')
%
% Output:
%   results/sphere/<subdir>/error_vs_N_sphere_phs_sweep.png

    if nargin < 1 || isempty(subdir), subdir = 'high'; end
    if nargin < 2 || isempty(base_results_dir), base_results_dir = 'results'; end

    results_dir = fullfile(base_results_dir, 'sphere', subdir);

    % load base file (FS + PLS + maybe PHS)
    pat_base = fullfile(results_dir, 'results_sphere*.mat');
    D = dir(pat_base);
    assert(~isempty(D), 'No results_sphere*.mat files in %s', results_dir);

    % choose a base file that has sN/Ns + el2 + el2_poly
    baseIdx = NaN;
    for j = 1:numel(D)
        S = load(fullfile(results_dir, D(j).name));
        if (isfield(S,'sN') || isfield(S,'Ns')) && isfield(S,'el2') && isfield(S,'el2_poly')
            baseIdx = j; break
        end
    end
    assert(~isnan(baseIdx), 'No base file with (sN or Ns) + el2 + el2_poly found.');

    Sbase = load(fullfile(results_dir, D(baseIdx).name));

    % x-axis: sN
    if isfield(Sbase,'sN') && ~isempty(Sbase.sN)
        sN = Sbase.sN(:);
    elseif isfield(Sbase,'Ns') && ~isempty(Sbase.Ns)
        sN = sqrt(Sbase.Ns(:));
    else
        error('Base file lacks sN and Ns: %s', D(baseIdx).name);
    end
    n = numel(sN);

    % FS/PLS from base
    el2_fs = padToN(getv(Sbase,'el2'), n);
    el2_pls = padToN(getv(Sbase,'el2_poly'), n);

    % gather ALL PHS series by degree
    Dp = dir(fullfile(results_dir, 'results_sphere_phsp_*.mat'));
    % also include PHS in base (if present & non-NaN)
    PHScurves = struct('p', {}, 'y', {});
    if isfield(Sbase,'el2_phs') && any(isfinite(Sbase.el2_phs(:)))
        PHScurves(end+1).p = getPfromName(D(baseIdx).name); % may be NaN; will rename to 'base' if needed
        PHScurves(end).y = padToN(Sbase.el2_phs(:), n);
    end
    % sweep files
    for j = 1:numel(Dp)
        Sp = load(fullfile(results_dir, Dp(j).name));
        if isfield(Sp,'el2_phs') && ~isempty(Sp.el2_phs)
            pdeg = getPfromName(Dp(j).name);
            PHScurves(end+1).p = pdeg;
            PHScurves(end).y = padToN(Sp.el2_phs(:), n);
        end
    end

    % sort PHS series by degree (NaN last)
    if ~isempty(PHScurves)
        [~, ord] = sort( arrayfun(@(t) t.p, PHScurves) );
        PHScurves = PHScurves(ord);
    end

    % y-limits (log-safe)
    collect = [el2_fs; el2_pls];
    for j = 1:numel(PHScurves), collect = [collect; PHScurves(j).y(:)]; %#ok<AGROW>
    end
    collect = collect(isfinite(collect) & collect>0);
    if isempty(collect), ymin = 1e-16; ymax = 1;
    else
        ymin = 10^floor(log10(min(collect)));
        ymax = 10^ceil (log10(max(collect)));
    end

    % exponential trendline on FS (finite points)
    ok = isfinite(sN) & isfinite(el2_fs) & el2_fs>0;
    if nnz(ok) >= 2
        xx = sN(ok); 
        yy = el2_fs(ok);
        p  = polyfit(xx, log(yy), 1);
        yfit = exp(p(2)) * exp(p(1) * xx);
        fit_txt = sprintf('exp trendline');
    else
        xx = []; 
        yfit = []; 
        fit_txt = '';
    end

    % styling
    cPLS = [0.000 0.447 0.741]; 
    mPLS='o';
    cFS = [0.850 0.325 0.098]; 
    mFS='^';
    cTr = [0.5 0.5 0.5];
    ms = 7;  
    lw = 2.2;

    % distinct colormap for PHS sweep
    phs_colors = lines(max(1, numel(PHScurves))); 

    %  figure
    h = figure; set(h,'Color','none','InvertHardcopy','off');
    ax = gca;
    set(ax,'Color','none','LineWidth',2,'FontSize',16,'FontName','Helvetica');
    ax.FontWeight = 'bold';
    ax.XAxis.FontWeight = 'bold';
    ax.YAxis.FontWeight = 'bold';
    set(ax,'XScale','linear','YScale','log');
    ax.TickLength = [0.02 0.02];
    ax.Box = 'on'; 
    hold on;

    H = gobjects(0); 
    L = {};

    % FS & PLS
    if any(isfinite(el2_fs))
        H(end+1) = semilogy(sN, el2_fs, '-', 'Color', cFS,  'Marker', mFS, 'MarkerSize', ms, 'MarkerFaceColor', 'w', 'LineWidth', lw);  L{end+1}='FS';
    end
    if any(isfinite(el2_pls))
        H(end+1) = semilogy(sN, el2_pls,'-', 'Color', cPLS, 'Marker', mPLS, 'MarkerSize', ms, 'MarkerFaceColor', 'w', 'LineWidth', lw);  L{end+1}='PLS';
    end

    % PHS curves
    for j = 1:numel(PHScurves)
        y = PHScurves(j).y(:);
        if any(isfinite(y) & y>0)
            name = PHScurves(j).p;
            if isnan(name), leg = 'PHS+poly (base)'; else, leg = sprintf('PHS+poly p=%d', name); end
            H(end+1) = semilogy(sN, y, '-', 'Color', phs_colors(j,:), 'Marker', '>', 'MarkerSize', ms, 'MarkerFaceColor','w', 'LineWidth', lw);
            L{end+1} = leg;
        end
    end

    % trendline
    if ~isempty(yfit)
        H(end+1) = semilogy(xx, yfit, '--', 'Color', cTr, 'LineWidth', 2.0); L{end+1}=fit_txt;
    end

    % labels and title
    xlabel('N^{1/2}','Interpreter','tex','FontSize',16,'FontWeight','bold');
    ylabel('Relative l_2 error','Interpreter','tex','FontSize',16,'FontWeight','bold');
    title('Sphere: Relative l_2 error vs. N^{1/d}, C^2(R^3) Wendland Kernel','Interpreter','tex','FontSize',18,'FontWeight','bold');

    % legend
    if ~isempty(H)
        lgd = legend(H,L,'Location','best','Interpreter','tex');
        set(lgd,'FontSize',14,'FontWeight','bold','Box','on','LineWidth',1.5,'Color','w');
    end

    % limits
    xgood = sN(isfinite(sN)); if isempty(xgood), xgood = [1;1]; end
    xlim([min(xgood) max(xgood)]); ylim([ymin ymax]);

    % export
    outpng = fullfile(results_dir, 'error_vs_N_sphere_phs_sweep.png');
    export_fig(gcf, outpng, '-png','-r300','-transparent');
    close(h);

    if ~isempty(yfit)
        fprintf('[sphere] FS exp fit: E(N) = %.2f * exp(%.2f * N^{1/d}) -> %s\n', exp(p(2)), p(1), outpng);
    else
        fprintf('[sphere] FS exp fit skipped -> %s\n', outpng);
    end

    % helpers
    function v = getv(S, fld)
        if isfield(S,fld), v = S.(fld)(:); else, v = []; end
    end
    function v = padToN(v, n0)
        v = v(:); m = numel(v);
        if m==0, v = NaN(n0,1); return; end
        if m<n0,  v = [v; NaN(n0-m,1)]; else, v = v(1:n0); end
    end
    function p = getPfromName(fname)
        % parse ..._phsp_<p>.mat
        p = NaN;
        tok = regexp(fname, 'phsp_(\d+)\.mat$', 'tokens','once');
        if ~isempty(tok), p = str2double(tok{1}); end
    end
end
