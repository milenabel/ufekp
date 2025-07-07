%% Pick a function to run
function plot_results(function_name, base_results_dir)
    if nargin < 2 || isempty(base_results_dir)
        base_results_dir = fullfile('csrbfjacobi/code/Approximation/CleanLagrangeApprox/results/');
    end
    % Load data
    res = loadResults(function_name, base_results_dir);

    % Sort once by N
    [res.sNs, res.I] = sort(res.sN);

    % Reorder all values vs N
    res.el2_poly_s = res.el2_poly(res.I);
    res.el2_diag_s = res.el2_diag(res.I,:);    
    res.el2_fs1_s = res.el2_fs1(res.I,:);
    res.el2_fs2_s = res.el2_fs2(res.I,:);
    res.el2_fs3_s = res.el2_fs3(res.I,:);
    res.el2_vs1_s = res.el2_vs1(res.I,:);
    res.el2_vs2_s = res.el2_vs2(res.I,:);
    res.el2_vs3_s = res.el2_vs3(res.I,:);
    % assembly times
    res.a_time_fs1_s = res.a_time_fs1(res.I,:);
    res.a_time_fs2_s = res.a_time_fs2(res.I,:);
    res.a_time_fs3_s = res.a_time_fs3(res.I,:);
    res.a_time_vs1_s = res.a_time_vs1(res.I,:);
    res.a_time_vs2_s = res.a_time_vs2(res.I,:);
    res.a_time_vs3_s = res.a_time_vs3(res.I,:);
    % eval times
    res.e_time_fs1_s = res.e_time_fs1(res.I,:);
    res.e_time_fs2_s = res.e_time_fs2(res.I,:);
    res.e_time_fs3_s = res.e_time_fs3(res.I,:);
    res.e_time_vs1_s = res.e_time_vs1(res.I,:);
    res.e_time_vs2_s = res.e_time_vs2(res.I,:);
    res.e_time_vs3_s = res.e_time_vs3(res.I,:);
    % sparsities
    res.sp_fs1_s = res.sparsity_fs1(res.I,:);
    res.sp_fs2_s = res.sparsity_fs2(res.I,:);
    res.sp_fs3_s = res.sparsity_fs3(res.I,:);
    res.sp_vs1_s = res.sparsity_vs1(res.I,:);
    res.sp_vs2_s = res.sparsity_vs2(res.I,:);
    res.sp_vs3_s = res.sparsity_vs3(res.I,:);
    % Create all plots (sor all smoothenesses)
    for sm = 1:3
        plotErrorVsN(res, sm);
        plotAssemblyTimeVsN(res, sm);
        plotEvalTimeVsN(res, sm);
        plotSparsityVsN(res, sm);
    end
end

%% Load results files
function res = loadResults(function_name, base_results_dir)
    results_dir = fullfile(base_results_dir, function_name);
    matfile = fullfile(results_dir, sprintf('results_%s.mat', function_name));
    if ~exist(matfile,'file')
        error('Results file not found: %s', matfile);
    end
    data = load(matfile);
    data.results_dir = results_dir;
    res = data;
end

%% Plot relative L2 error vs N^{1/d}
function plotErrorVsN(res, sm)
    dim = res.dim;
    sNs = res.sNs;

    el2_poly_s = res.el2_poly_s;
    el2_diag_s = res.el2_diag_s(:,sm);
    el2_fs1_s = res.el2_fs1_s(:,sm);
    el2_fs2_s = res.el2_fs2_s(:,sm);
    el2_fs3_s = res.el2_fs3_s(:,sm);
    el2_vs1_s = res.el2_vs1_s(:,sm);
    el2_vs2_s = res.el2_vs2_s(:,sm);
    el2_vs3_s = res.el2_vs3_s(:,sm);

    h = figure;  %grid on;
    % marks = {'-o','-s','-^','--o','--s','--^','-.x','-.+'};
    marks = {'-o','-s','-^','-x','-+'};
    semilogy(sNs, el2_poly_s, marks{1}, 'LineWidth',1.2);
    hold on;
    semilogy(sNs, el2_diag_s, marks{2}, 'LineWidth',1.2);
    semilogy(sNs, el2_fs1_s,  marks{3}, 'LineWidth',1.2);
    semilogy(sNs, el2_fs2_s,  marks{4}, 'LineWidth',1.2);
    semilogy(sNs, el2_fs3_s,  marks{5}, 'LineWidth',1.2);

    if dim==1
        xlabel(sprintf('$N$', dim), 'FontSize',14,'Interpreter','latex');
    else
        xlabel(sprintf('$N^{1/%d}$', dim), 'FontSize',14,'Interpreter','latex');
    end
    ylabel('Relative $\ell_2$ error','FontSize',14,'Interpreter','latex');
    legend({'PLS','Diag','$K_t$ = $1e12$','$K_t$ = $1e8$','$K_t$ = $1e4$'}, 'Location','best','FontSize',10,'Interpreter','latex');

    title(sprintf('$\\ell_2$ error vs. $N^{1/d}$, smoothness=%d',sm),'Interpreter','latex','FontSize',16);
    set(gca,'FontSize',12);
    export_fig(gcf, fullfile(res.results_dir,sprintf('error_vs_N_fs_s%d.png',sm)),'-png','-r300');
    close(h);


    h = figure;  %grid on;
    marks = {'-o','-s','-^','-x','-+'};
    semilogy(sNs, el2_poly_s, marks{1}, 'LineWidth',1.2);
    hold on;
    semilogy(sNs, el2_diag_s, marks{2}, 'LineWidth',1.2);
    semilogy(sNs, el2_vs1_s,  marks{3}, 'LineWidth',1.2);
    semilogy(sNs, el2_vs2_s,  marks{4}, 'LineWidth',1.2);
    semilogy(sNs, el2_vs3_s,  marks{5}, 'LineWidth',1.2);

    if dim==1
        xlabel(sprintf('$N$', dim), 'FontSize',14,'Interpreter','latex');
    else
        xlabel(sprintf('$N^{1/%d}$', dim), 'FontSize',14,'Interpreter','latex');
    end
    ylabel('Relative $\\ell_2$ error','FontSize',14,'Interpreter','latex');
    legend({'PLS','Diag','$K_t$ = $1e12$','$K_t$ = $1e8$','$K_t$ = $1e4$'}, 'Location','best','FontSize',10,'Interpreter','latex');
    title(sprintf('$\\ell_2$ error vs. $N^{1/d}$, smoothness=%d',sm),'Interpreter','latex','FontSize',16);
    set(gca,'FontSize',12);
    export_fig(gcf, fullfile(res.results_dir,sprintf('error_vs_N_vs_s%d.png',sm)),'-png','-r300');
    close(h);
end

%% Plot assembly time vs N^{1/d}
function plotAssemblyTimeVsN(res, sm)
    dim = res.dim;
    sN = res.sN;

    a_fs1 = res.a_time_fs1_s(:,sm);
    a_fs2 = res.a_time_fs2_s(:,sm);
    a_fs3 = res.a_time_fs3_s(:,sm);
    a_vs1 = res.a_time_vs1_s(:,sm);
    a_vs2 = res.a_time_vs2_s(:,sm);
    a_vs3 = res.a_time_vs3_s(:,sm);

    h = figure; grid on;
    marks = {'-o','-s','-^','-x','-+'};
    semilogy(sN, a_fs1, marks{1}, 'LineWidth',1.2);
    hold on; 
    semilogy(sN, a_fs2, marks{2}, 'LineWidth',1.2);
    semilogy(sN, a_fs3, marks{3}, 'LineWidth',1.2);
    if dim==1
        xlabel(sprintf('$N$', dim), 'FontSize',14,'Interpreter','latex');
    else
        xlabel(sprintf('$N^{1/%d}$', dim), 'FontSize',14,'Interpreter','latex');
    end
    ylabel('Assembly time (s)','FontSize',14);
    legend({'$K_t$ = $1e12$','$K_t$ = $1e8$','$K_t$ = $1e4$'}, 'Location','northwest','Interpreter','latex');
    title(sprintf('Assembly time, smoothness=$%d$',sm),'Interpreter','latex','FontSize',16);
    set(gca,'FontSize',12);
    export_fig(gcf, fullfile(res.results_dir,sprintf('assembly_time_vs_N_fs_s%d.png',sm)), '-png','-r300');
    close(h);

    h = figure; grid on;
    marks = {'-o','-s','-^','-x','-+'};
    semilogy(sN, a_vs1, marks{1}, 'LineWidth',1.2);
    hold on; 
    semilogy(sN, a_vs2, marks{2}, 'LineWidth',1.2);
    semilogy(sN, a_vs3, marks{3}, 'LineWidth',1.2);
    if dim==1
        xlabel(sprintf('$N$', dim), 'FontSize',14,'Interpreter','latex');
    else
        xlabel(sprintf('$N^{1/%d}$', dim), 'FontSize',14,'Interpreter','latex');
    end
    ylabel('Assembly time (s)','FontSize',14);
    legend({'$K_t$ = $1e12$','$K_t$ = $1e8$','$K_t$ = $1e4$'}, 'Location','northwest','Interpreter','latex');
    title(sprintf('Assembly time, smoothness=$%d$',sm),'Interpreter','latex','FontSize',16);
    set(gca,'FontSize',12);
    export_fig(gcf, fullfile(res.results_dir,sprintf('assembly_time_vs_N_vs_s%d.png',sm)), '-png','-r300');
    close(h);
end

%% Plot evaluation time vs N^{1/d}
function plotEvalTimeVsN(res, sm)
    dim = res.dim;
    sN = res.sN;

    e_fs1 = res.e_time_fs1_s(:,sm);
    e_fs2 = res.e_time_fs2_s(:,sm);
    e_fs3 = res.e_time_fs3_s(:,sm);
    e_vs1 = res.e_time_vs1_s(:,sm);
    e_vs2 = res.e_time_vs2_s(:,sm);
    e_vs3 = res.e_time_vs3_s(:,sm);

    h = figure; grid on;
    marks = {'-o','-s','-^','-x','-+'};
    semilogy(sN, e_fs1, marks{1}, 'LineWidth',1.2);
    hold on;
    semilogy(sN, e_fs2, marks{2}, 'LineWidth',1.2);
    semilogy(sN, e_fs3, marks{3}, 'LineWidth',1.2);
    if dim==1
        xlabel(sprintf('$N$', dim), 'FontSize',14,'Interpreter','latex');
    else
        xlabel(sprintf('$N^{1/%d}$', dim), 'FontSize',14,'Interpreter','latex');
    end
    ylabel('Evaluation time (s)','FontSize',14);
    legend({'$K_t$ = $1e12$','$K_t$ = $1e8$','$K_t$ = $1e4$'}, 'Location','northwest','Interpreter','latex');
    title(sprintf('Evaluation time, smoothness=$%d$',sm),'Interpreter','latex','FontSize',16);
    set(gca,'FontSize',12);
    export_fig(gcf, fullfile(res.results_dir,sprintf('eval_time_vs_N_fs_s%d.png',sm)), '-png','-r300');
    close(h);

    h = figure; grid on;
    marks = {'-o','-s','-^','-x','-+'};
    semilogy(sN, e_vs1, marks{1}, 'LineWidth',1.2);
    hold on; 
    semilogy(sN, e_vs2, marks{2}, 'LineWidth',1.2);
    semilogy(sN, e_vs3, marks{3}, 'LineWidth',1.2);
    if dim==1
        xlabel(sprintf('$N$', dim), 'FontSize',14,'Interpreter','latex');
    else
        xlabel(sprintf('$N^{1/%d}$', dim), 'FontSize',14,'Interpreter','latex');
    end
    ylabel('Evaluation time (s)','FontSize',14);
    legend({'$K_t$ = $1e12$','$K_t$ = $1e8$','$K_t$ = $1e4$'}, 'Location','northwest','Interpreter','latex');
    title(sprintf('Evaluation time, smoothness=$%d$',sm),'Interpreter','latex','FontSize',16);
    set(gca,'FontSize',12);
    export_fig(gcf, fullfile(res.results_dir,sprintf('eval_time_vs_N_vs_s%d.png',sm)), '-png','-r300');
    close(h);
end

%% Plot error vs sparsity (tiled, one panel per smoothness)
function plotSparsityVsN(res, sm)
    dim = res.dim;
    sN = res.sN;

    sp_fs1 = res.sp_fs1_s(:,sm);
    sp_fs2 = res.sp_fs2_s(:,sm);
    sp_fs3 = res.sp_fs3_s(:,sm);
    sp_vs1 = res.sp_vs1_s(:,sm);
    sp_vs2 = res.sp_vs2_s(:,sm);
    sp_vs3 = res.sp_vs3_s(:,sm);

    h = figure; grid on;
    marks = {'-o','-s','-^','-x','-+'};
    semilogy(sN, sp_fs1, marks{1}, 'LineWidth',1.2);
    hold on; 
    semilogy(sN, sp_fs2, marks{2}, 'LineWidth',1.2);
    semilogy(sN, sp_fs3, marks{3}, 'LineWidth',1.2);

    if dim==1
        xlabel('$N$','Interpreter','latex','FontSize',14);
    else
        xlabel(sprintf('$N^{1/%d}$',dim),'Interpreter','latex','FontSize',14);
    end
    ylabel('Achieved sparsity','FontSize',14);
    legend({'$K_t$ = $1e12$','$K_t$ = $1e8$','$K_t$ = $1e4$'}, 'Location','southeast','Interpreter','latex');
    title(sprintf('Sparsity vs. $N^{1/d}$, smoothness=%d',sm),'Interpreter','latex','FontSize',16);
    set(gca,'FontSize',12);

    export_fig(gcf, fullfile(res.results_dir,sprintf('sparsity_vs_N_fs_s%d.png', sm)), '-png','-r300');
    close(h);

    h = figure; grid on;
    marks = {'-o','-s','-^','-x','-+'};
    semilogy(sN, sp_fs1, marks{1}, 'LineWidth',1.2);
    hold on; 
    semilogy(sN, sp_vs2, marks{2}, 'LineWidth',1.2);
    semilogy(sN, sp_vs3, marks{3}, 'LineWidth',1.2);

    if dim==1
        xlabel('$N$','Interpreter','latex','FontSize',14);
    else
        xlabel(sprintf('$N^{1/%d}$',dim),'Interpreter','latex','FontSize',14);
    end
    ylabel('Achieved sparsity','FontSize',14);
    legend({'$K_t$ = $1e12$','$K_t$ = $1e8$','$K_t$ = $1e4$'}, 'Location','southeast','Interpreter','latex');
    title(sprintf('Sparsity vs. $N^{1/d}$, smoothness=%d',sm),'Interpreter','latex','FontSize',16);
    set(gca,'FontSize',12);

    export_fig(gcf, fullfile(res.results_dir,sprintf('sparsity_vs_N_vs_s%d.png', sm)), '-png','-r300');
    close(h);
end