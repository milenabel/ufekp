%% Pick a function to run
function plotnb_results(function_name, subdir, base_results_dir)
    if nargin < 3 || isempty(base_results_dir)
        base_results_dir = fullfile('code/Approximation/CleanLagrangeApprox/results');
    end
    if nargin < 2 || isempty(subdir)
        subdir = 'high'; % Default to 'high' if not specified
    end
    % Load data
    res = loadResults(function_name, subdir, base_results_dir);
    res.sNs = res.sN;

    % Compute global x‐limits 
    xmin = min(res.sNs);
    xmax = max(res.sNs);

    % Compute global y‐limits for each plot‐type:
    % relative L2 error
    all_errs = [res.el2_poly; res.el2_diag(:); res.el2_fs1(:); res.el2_fs2(:); res.el2_fs3(:); res.el2_vs1(:); res.el2_vs2(:); res.el2_vs3(:)];
    err_min = min(all_errs(:));
    err_max = max(all_errs(:));

    % assembly time
    all_atime = [res.a_time_fs1(:); res.a_time_fs2(:); res.a_time_fs3(:); res.a_time_vs1(:); res.a_time_vs2(:); res.a_time_vs3(:)];
    at_min = min(all_atime(:));
    at_max = max(all_atime(:));

    % evaluation time
    all_etime = [res.e_time_fs1(:); res.e_time_fs2(:); res.e_time_fs3(:); res.e_time_vs1(:); res.e_time_vs2(:); res.e_time_vs3(:)];
    et_min = min(all_etime(:));
    et_max = max(all_etime(:));

    % sparsity
    all_sp = [res.sparsity_fs1(:); res.sparsity_fs2(:); res.sparsity_fs3(:); res.sparsity_vs1(:); res.sparsity_vs2(:); res.sparsity_vs3(:)];
    sp_min = min(all_sp(:));
    sp_max = max(all_sp(:));

    % evaluation‐time y‐limits (never zero so no need to drop zeros)
    et_min = 10^floor(log10(et_min));
    et_max = 10^ceil (log10(et_max));

    % Create all plots (sor all smoothenesses)
    for sm = 1:3
        plotErrorVsN(res, sm, xmin, xmax, err_min, err_max);
        plotAssemblyTimeVsN(res, sm, xmin, xmax, at_min, at_max);
        plotEvalTimeVsN(res, sm, xmin, xmax, et_min, et_max);
        plotSparsityVsN(res, sm, xmin, xmax, sp_min, sp_max);
    end
end

%% Load results files
function res = loadResults(function_name, subdir, base_results_dir)
    valid_subdirs = {'high', 'low', 'fixed', 'mixed'};
    if ~ismember(subdir, valid_subdirs)
        error('Invalid subdirectory: %s. Must be one of: high, low, fixed', subdir);
    end
    results_dir = fullfile(base_results_dir, function_name, subdir);
    matfile = fullfile(results_dir, sprintf('results_%s.mat', function_name));

    if ~exist(matfile,'file')
        error('Results file not found: %s', matfile);
    end
    data = load(matfile);
    data.results_dir = results_dir;
    res = data;
end

%% Plot relative L2 error vs N^{1/d}
function plotErrorVsN(res, sm, xmin, xmax, err_min, err_max)
    dim = res.dim;
    sNs = res.sNs;

    el2_poly_s = res.el2_poly;
    el2_diag_s = res.el2_diag(:,sm);
    el2_fs1_s = res.el2_fs1(:,sm);
    el2_fs2_s = res.el2_fs2(:,sm);
    el2_fs3_s = res.el2_fs3(:,sm);
    el2_vs1_s = res.el2_vs1(:,sm);
    el2_vs2_s = res.el2_vs2(:,sm);
    el2_vs3_s = res.el2_vs3(:,sm);

    h = figure;  %grid on;
    set(h, 'Color', 'none');          
    ax = gca;
    set(ax, 'Color', 'none','LineWidth', 1.2);     
    marks = {'-o','-s','-^','-x','-+'};
    semilogy(sNs, el2_poly_s, marks{1}, 'LineWidth',1.2);
    hold on;
    semilogy(sNs, el2_diag_s, marks{2}, 'LineWidth',1.2);
    semilogy(sNs, el2_fs1_s,  marks{3}, 'LineWidth',1.2);
    semilogy(sNs, el2_fs2_s,  marks{4}, 'LineWidth',1.2);
    semilogy(sNs, el2_fs3_s,  marks{5}, 'LineWidth',1.2);

    if dim==1
        xlabel(sprintf('$N$', dim),'Interpreter','latex','FontSize',16);
    else
        xlabel(sprintf('$N^{1/%d}$', dim),'Interpreter','latex','FontSize',16);
    end
    ylabel('Relative $\ell_2$ error');
    legend({'PLS','Diag','$K_t$ = $1e12$','$K_t$ = $1e8$','$K_t$ = $1e4$'}, 'Location','best','Interpreter','latex','FontSize',14);

    if sm == 1
        title(sprintf('$\\ell_2$ error vs. $N^{1/d}$, $C^2(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18);
    elseif sm==2
        title(sprintf('$\\ell_2$ error vs. $N^{1/d}$, $C^4(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18);
    elseif sm==3
        title(sprintf('$\\ell_2$ error vs. $N^{1/d}$, $C^6(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18); 
    end

    % force all plots to use the same x‐ and y-ranges
    xlim([xmin xmax]);
    ylim([err_min err_max]);

    %set(gca,'FontSize',12,'FontWeight', 'bold');
    export_fig(gcf, fullfile(res.results_dir,sprintf('error_vs_N_fs_s%d.png',sm)),'-png','-r300','-transparent');
    close(h);


    h = figure;  %grid on;
    set(h, 'Color', 'none');          
    ax = gca;
    set(ax, 'Color', 'none','FontSize', 12, 'LineWidth', 1.2, 'TickLabelInterpreter', 'latex');
    marks = {'-o','-s','-^','-x','-+'};
    semilogy(sNs, el2_poly_s, marks{1}, 'LineWidth',1.2);
    hold on;
    semilogy(sNs, el2_diag_s, marks{2}, 'LineWidth',1.2);
    semilogy(sNs, el2_vs1_s,  marks{3}, 'LineWidth',1.2);
    semilogy(sNs, el2_vs2_s,  marks{4}, 'LineWidth',1.2);
    semilogy(sNs, el2_vs3_s,  marks{5}, 'LineWidth',1.2);

    if dim==1
        xlabel(sprintf('$N$', dim),'Interpreter','latex','FontSize',16);
    else
        xlabel(sprintf('$N^{1/%d}$', dim),'Interpreter','latex','FontSize',16);
    end
    ylabel('Relative $\ell_2$ error','Interpreter','latex','FontSize',16);
    legend({'PLS','Diag','$K_t$ = $1e12$','$K_t$ = $1e8$','$K_t$ = $1e4$'}, 'Location','best','Interpreter','latex','FontSize',14);

    if sm == 1
        title(sprintf('$\\ell_2$ error vs. $N^{1/d}$, $C^2(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18);
    elseif sm==2
        title(sprintf('$\\ell_2$ error vs. $N^{1/d}$, $C^4(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18);
    elseif sm==3
        title(sprintf('$\\ell_2$ error vs. $N^{1/d}$, $C^6(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18); 
    end

    % force all plots to use the same x‐ and y-ranges
    xlim([xmin xmax]);
    ylim([err_min err_max]);

    % set(gca,'FontSize',12,'FontWeight', 'bold');
    export_fig(gcf, fullfile(res.results_dir,sprintf('error_vs_N_vs_s%d.png',sm)),'-png','-r300','-transparent');
    close(h);
end

%% Plot assembly time vs N^{1/d}
function plotAssemblyTimeVsN(res, sm, xmin, xmax, at_min, at_max)
    dim = res.dim;
    sNs = res.sNs;

    a_fs1 = res.a_time_fs1(:,sm);
    a_fs2 = res.a_time_fs2(:,sm);
    a_fs3 = res.a_time_fs3(:,sm);
    a_vs1 = res.a_time_vs1(:,sm);
    a_vs2 = res.a_time_vs2(:,sm);
    a_vs3 = res.a_time_vs3(:,sm);

    h = figure; grid on;
    set(h, 'Color', 'none');          
    ax = gca;
    set(ax, 'Color', 'none','FontSize', 12, 'LineWidth', 1.2);
    marks = {'-o','-s','-^','-x','-+'};
    semilogy(sNs, a_fs1, marks{1}, 'LineWidth',1.2);
    hold on; 
    semilogy(sNs, a_fs2, marks{2}, 'LineWidth',1.2);
    semilogy(sNs, a_fs3, marks{3}, 'LineWidth',1.2);
    if dim==1
        xlabel(sprintf('$N$', dim), 'FontSize',16,'Interpreter','latex');
    else
        xlabel(sprintf('$N^{1/%d}$', dim), 'FontSize',16,'Interpreter','latex');
    end
    ylabel('Assembly and Solve time (s)','FontSize',16,'Interpreter','latex');
    legend({'$K_t$ = $1e12$','$K_t$ = $1e8$','$K_t$ = $1e4$'}, 'Location','northwest','FontSize',14,'Interpreter','latex');
    if sm==1
        title(sprintf('Assembly and Solve time, $C^2(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18);
    elseif sm==2
        title(sprintf('Assembly and Solve time, $C^4(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18);
    elseif sm==3
        title(sprintf('Assembly and Solve time, $C^6(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18);
    end

    % force all plots to use the same x‐ and y-ranges
    xlim([xmin xmax]);
    ylim([at_min at_max]);

    % set(gca,'FontSize',12,'FontWeight', 'bold');
    export_fig(gcf, fullfile(res.results_dir,sprintf('assembly_time_vs_N_fs_s%d.png',sm)), '-png','-r300','-transparent');
    close(h);


    h = figure; grid on;
    set(h, 'Color', 'none');          
    ax = gca;
    set(ax, 'Color', 'none','FontSize', 12, 'LineWidth', 1.2);
    marks = {'-o','-s','-^','-x','-+'};
    semilogy(sNs, a_vs1, marks{1}, 'LineWidth',1.2);
    hold on; 
    semilogy(sNs, a_vs2, marks{2}, 'LineWidth',1.2);
    semilogy(sNs, a_vs3, marks{3}, 'LineWidth',1.2);
    if dim==1
        xlabel(sprintf('$N$', dim), 'FontSize',16,'Interpreter','latex');
    else
        xlabel(sprintf('$N^{1/%d}$', dim), 'FontSize',16,'Interpreter','latex');
    end
    ylabel('Assembly and Solve time (s)','FontSize',16,'Interpreter','latex');
    legend({'$K_t$ = $1e12$','$K_t$ = $1e8$','$K_t$ = $1e4$'}, 'Location','northwest','FontSize',14,'Interpreter','latex');
    if sm==1
        title(sprintf('Assembly and Solve time, $C^2(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18);
    elseif sm==2
        title(sprintf('Assembly and Solve time, $C^4(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18);
    elseif sm==3
        title(sprintf('Assembly and Solve time, $C^6(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18);
    end

    % force all plots to use the same x‐ and y-ranges
    xlim([xmin xmax]);
    ylim([at_min at_max]);

    % set(gca,'FontSize',12,'FontWeight', 'bold');
    export_fig(gcf, fullfile(res.results_dir,sprintf('assembly_time_vs_N_vs_s%d.png',sm)), '-png','-r300','-transparent');
    close(h);
end

%% Plot evaluation time vs N^{1/d}
function plotEvalTimeVsN(res, sm, xmin, xmax, et_min, et_max)
    dim = res.dim;
    sNs = res.sNs;

    e_fs1 = res.e_time_fs1(:,sm);
    e_fs2 = res.e_time_fs2(:,sm);
    e_fs3 = res.e_time_fs3(:,sm);
    e_vs1 = res.e_time_vs1(:,sm);
    e_vs2 = res.e_time_vs2(:,sm);
    e_vs3 = res.e_time_vs3(:,sm);

    h = figure; grid on;
    set(h, 'Color', 'none');          
    ax = gca;
    set(ax, 'Color', 'none','FontSize', 12, 'LineWidth', 1.2);
    marks = {'-o','-s','-^','-x','-+'};
    semilogy(sNs, e_fs1, marks{1}, 'LineWidth',1.2);
    hold on;
    semilogy(sNs, e_fs2, marks{2}, 'LineWidth',1.2);
    semilogy(sNs, e_fs3, marks{3}, 'LineWidth',1.2);
    if dim==1
        xlabel(sprintf('$N$', dim), 'FontSize',16,'Interpreter','latex');
    else
        xlabel(sprintf('$N^{1/%d}$', dim), 'FontSize',16,'Interpreter','latex');
    end
    ylabel('Evaluation time (s)','FontSize',16,'Interpreter','latex');
    legend({'$K_t$ = $1e12$','$K_t$ = $1e8$','$K_t$ = $1e4$'}, 'Location','best','FontSize',14,'Interpreter','latex');
    if sm==1
        title(sprintf('Evaluation time, $C^2(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18);
    elseif sm==2
        title(sprintf('Evaluation time, $C^4(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18);
    elseif sm==3
        title(sprintf('Evaluation time, $C^6(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18);
    end

    % force all plots to use the same x‐ and y-ranges
    xlim([xmin xmax]);
    ylim([et_min et_max]);

    % set(gca,'FontSize',12,'FontWeight', 'bold');
    export_fig(gcf, fullfile(res.results_dir,sprintf('eval_time_vs_N_fs_s%d.png',sm)), '-png','-r300','-transparent');
    close(h);

    h = figure; grid on;
    set(h, 'Color', 'none');          
    ax = gca;
    set(ax, 'Color', 'none','FontSize', 12, 'LineWidth', 1.2);
    marks = {'-o','-s','-^','-x','-+'};
    semilogy(sNs, e_vs1, marks{1}, 'LineWidth',1.2);
    hold on; 
    semilogy(sNs, e_vs2, marks{2}, 'LineWidth',1.2);
    semilogy(sNs, e_vs3, marks{3}, 'LineWidth',1.2);
    if dim==1
        xlabel(sprintf('$N$', dim), 'FontSize',16,'Interpreter','latex');
    else
        xlabel(sprintf('$N^{1/%d}$', dim), 'FontSize',16,'Interpreter','latex');
    end
    ylabel('Evaluation time (s)','FontSize',16,'Interpreter','latex');
    legend({'$K_t$ = $1e12$','$K_t$ = $1e8$','$K_t$ = $1e4$'}, 'Location','best','FontSize',14,'Interpreter','latex');
    if sm==1
        title(sprintf('Evaluation time, $C^2(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18);
    elseif sm==2
        title(sprintf('Evaluation time, $C^4(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18);
    elseif sm==3
        title(sprintf('Evaluation time, $C^6(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18);
    end

    % force all plots to use the same x‐ and y-ranges
    xlim([xmin xmax]);
    ylim([et_min et_max]);

    % set(gca,'FontSize',12,'FontWeight', 'bold');
    export_fig(gcf, fullfile(res.results_dir,sprintf('eval_time_vs_N_vs_s%d.png',sm)), '-png','-r300','-transparent');
    close(h);
end

%% Plot error vs sparsity (tiled, one panel per smoothness)
function plotSparsityVsN(res, sm, xmin, xmax, sp_min, sp_max)
    dim = res.dim;
    sNs = res.sNs;

    sp_fs1 = res.sparsity_fs1(:,sm);
    sp_fs2 = res.sparsity_fs2(:,sm);
    sp_fs3 = res.sparsity_fs3(:,sm);
    sp_vs1 = res.sparsity_vs1(:,sm);
    sp_vs2 = res.sparsity_vs2(:,sm);
    sp_vs3 = res.sparsity_vs3(:,sm);

    h = figure; grid on;
    set(h, 'Color', 'none');          
    ax = gca;
    set(ax, 'Color', 'none','FontSize', 12, 'LineWidth', 1.2);
    marks = {'-o','-s','-^','-x','-+'};
    plot(sNs, sp_fs1, marks{1}, 'LineWidth',1.2);
    hold on; 
    plot(sNs, sp_fs2, marks{2}, 'LineWidth',1.2);
    plot(sNs, sp_fs3, marks{3}, 'LineWidth',1.2);

    if dim==1
        xlabel('$N$','Interpreter','latex','FontSize',16);
    else
        xlabel(sprintf('$N^{1/%d}$',dim),'Interpreter','latex','FontSize',16);
    end
    ylabel('Achieved sparsity','FontSize',16,'Interpreter','latex');
    legend({'$K_t$ = $1e12$','$K_t$ = $1e8$','$K_t$ = $1e4$'}, 'Location','southeast','FontSize',14,'Interpreter','latex');
    if sm==1
        title(sprintf('Sparsity vs. $N^{1/d}$, $C^2(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18);
    elseif sm==2
        title(sprintf('Sparsity vs. $N^{1/d}$, $C^4(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18);
    elseif sm==3
        title(sprintf('Sparsity vs. $N^{1/d}$, $C^6(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18);
    end

    % force all plots to use the same x‐ and y-ranges:
    xlim([xmin xmax]);
    ylim([sp_min sp_max]);

    % set(gca,'FontSize',12,'FontWeight', 'bold');

    export_fig(gcf, fullfile(res.results_dir,sprintf('sparsity_vs_N_fs_s%d.png', sm)), '-png','-r300','-transparent');
    close(h);

    h = figure; grid on;
    set(h, 'Color', 'none');          
    ax = gca;
    set(ax, 'Color', 'none','FontSize', 12, 'LineWidth', 1.2);
    marks = {'-o','-s','-^','-x','-+'};
    plot(sNs, sp_fs1, marks{1}, 'LineWidth',1.2);
    hold on; 
    plot(sNs, sp_vs2, marks{2}, 'LineWidth',1.2);
    plot(sNs, sp_vs3, marks{3}, 'LineWidth',1.2);

    if dim==1
        xlabel('$N$','Interpreter','latex','FontSize',16);
    else
        xlabel(sprintf('$N^{1/%d}$',dim),'Interpreter','latex','FontSize',16);
    end
    ylabel('Achieved sparsity','FontSize',16,'Interpreter','latex');
    legend({'$K_t$ = $1e12$','$K_t$ = $1e8$','$K_t$ = $1e4$'}, 'Location','southeast','FontSize',14,'Interpreter','latex');
    if sm==1
        title(sprintf('Sparsity vs. $N^{1/d}$, $C^2(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18);
    elseif sm==2
        title(sprintf('Sparsity vs. $N^{1/d}$, $C^4(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18);
    elseif sm==3
        title(sprintf('Sparsity vs. $N^{1/d}$, $C^6(\\mathrm{R}^3)$ Wendland Kernel'),'Interpreter','latex','FontSize',18);
    end

    % force all plots to use the same x‐ and y-ranges:
    xlim([xmin xmax]);
    ylim([sp_min sp_max]);

    export_fig(gcf, fullfile(res.results_dir,sprintf('sparsity_vs_N_vs_s%d.png', sm)), '-png','-r300','-transparent');
    close(h);
end
