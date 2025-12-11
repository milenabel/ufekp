function sweep_phs_only_sphere(p_list, subdir, base_results_dir)
% Run ONLY PHS+poly on the sphere for multiple degrees p and save results.
% Making sure the picked PHS degree is the best one
% One MAT file per degree:
%   results_sphere_phsp_<p>.mat
%
% Usage:
%   sweep_phs_only_sphere(1:9);  % defaults: subdir='high', base='results'
%   sweep_phs_only_sphere([1 3 5], 'high', 'results');
%
% Notes:
% - Classic PHS uses odd degrees (1,3,5,7,9). You can pass even p if desired.
% - This does NOT recompute PLS/CSRBF; it only writes el2_phs* series.

    if nargin < 1 || isempty(p_list),        p_list = [1, 3, 5, 7]; end
    if nargin < 2 || isempty(subdir),        subdir = 'high'; end
    if nargin < 3 || isempty(base_results_dir), base_results_dir = 'results'; end

    rng(0); % reproducible node sets 

    % grid sizes 
    Ns = [100, 200, 400, 800, 1600, 3200, 6400, 12800];
    Ns = max(1, round(Ns));        % ensure valid integers
    sN = sqrt(Ns);

    % evaluation set on the sphere (spiral_points on S^2)
    Ne = 15000;
    xe = spiral_points(1, [0,0,0], Ne);

    % target function 
    f = @(x) (x(:,1).^2).*abs(x(:,1)) + ...
             (x(:,2).^2).*abs(x(:,2)) + ...
             (x(:,3).^2).*abs(x(:,3));
    f_xe = f(xe); f_xe(isnan(f_xe)) = 0;

    % numeric control for the KKT pseudoinverse
    pinv_tol = 1e-11;

    % io
    function_name = 'sphere';
    results_dir = fullfile(base_results_dir, function_name, subdir);
    if ~exist(results_dir, 'dir'), mkdir(results_dir); end

    % sweep degrees
    for p_phs = p_list(:).'
        fprintf('[sphere/PHS] degree p = %d\n', p_phs);

        K = numel(Ns);
        el2_phs    = nan(K,1);
        elinf_phs  = nan(K,1);
        a_time_phs = nan(K,1);
        e_time_phs = nan(K,1);
        cond_phs   = nan(K,1);
        spars_phs  = nan(K,1);
        ell_store  = nan(K,1);

        for k = 1:K
            Nk = Ns(k);
            % nodes on the sphere
            x = spiral_points(1, [0,0,0], Nk);

            % values at nodes
            f_x = f(x); f_x(isnan(f_x)) = 0;

            % ambient poly degree in R^3 (your rule)
            ell = floor(nthroot(size(x,1), 3));
            ell_store(k) = ell;

            % PHS+poly solve (dense KKT + pseudoinverse)
            [el2_phs(k), elinf_phs(k), a_time_phs(k), e_time_phs(k), ...
                ~, cond_phs(k), ~, spars_phs(k)] = ...
                PHS_poly_man(x, f_x, ell, xe, p_phs, f_xe, pinv_tol);
        end

        % save one file per degree
        timestamp = datestr(datetime('now'), 'yyyyMMdd_HHmmss');
        out = fullfile(results_dir, sprintf('results_%s_phsp_%d.mat', function_name, p_phs));
        save(out, 'function_name','Ns','sN','p_phs','pinv_tol','timestamp', ...
                  'el2_phs','elinf_phs','a_time_phs','e_time_phs','cond_phs','spars_phs','ell_store');
        fprintf('  -> wrote %s\n', out);
    end

    fprintf('[sphere/PHS] done. Files in: %s\n', results_dir);
end