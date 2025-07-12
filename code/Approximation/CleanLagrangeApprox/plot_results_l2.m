%% Pick a function to run
function plot_results_l2(function_name, subdir, base_results_dir)
    if nargin < 3 || isempty(base_results_dir)
        base_results_dir = fullfile('results');
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
    all_errs = [res.el2_poly; res.el2(:)];
    % all_errs = [res.el2_poly; res.el2_diag(:); res.el2_fs1(:)];
    err_min = min(all_errs(:));
    err_max = max(all_errs(:));
    % 
    % % assembly time
    % all_atime = [res.a_time_fs1(:); res.a_time_fs2(:); res.a_time_fs3(:); res.a_time_vs1(:); res.a_time_vs2(:); res.a_time_vs3(:)];
    % % all_atime = [res.a_time_fs1(:),res.a_time_poly(:),res.a_time_diag(:)];
    % at_min = min(all_atime(:));
    % at_max = max(all_atime(:));
    % 
    % % evaluation time
    % all_etime = [res.e_time_fs1(:); res.e_time_fs2(:); res.e_time_fs3(:); res.e_time_vs1(:); res.e_time_vs2(:); res.e_time_vs3(:)];
    % all_etime = [res.e_time_fs1(:),res.e_time_poly(:),res.e_time_diag(:)];
    % et_min = min(all_etime(:));
    % et_max = max(all_etime(:));
    % 
    % % sparsity
    % all_sp = [res.sparsity_fs1(:); res.sparsity_fs2(:); res.sparsity_fs3(:); res.sparsity_vs1(:); res.sparsity_vs2(:); res.sparsity_vs3(:)];
    % all_sp = [res.sparsity_fs1(:)];
    % sp_min = min(all_sp(:));
    % sp_max = max(all_sp(:));
    % 
    % % evaluation‐time y‐limits (never zero so no need to drop zeros)
    % et_min = 10^floor(log10(et_min));
    % et_max = 10^ceil (log10(et_max));

    sm=1;

    % Create all plots (sor all smoothenesses)
    plotErrorVsN(res, sm, xmin, xmax, err_min, err_max);

    % for sm = 1:3
    %     plotErrorVsN(res, sm, xmin, xmax, err_min, err_max);
    %     % plotAssemblyTimeVsN(res, sm, xmin, xmax, at_min, at_max);
    %     % plotEvalTimeVsN(res, sm, xmin, xmax, et_min, et_max);
    %     % plotSparsityVsN(res, sm, xmin, xmax, sp_min, sp_max);
    % end
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
    % dim = res.dim;
    dim = 2;
    sNs = res.sNs;

    el2_poly_s = res.el2_poly;
    % el2_diag_s = res.el2_diag(:,sm);
    el2_s = res.el2(:,sm);
    % el2_fs1_s = res.el2_fs1(:,sm);
    % el2_fs2_s = res.el2_fs2(:,sm);
    % el2_fs3_s = res.el2_fs3(:,sm);
    % el2_vs1_s = res.el2_vs1(:,sm);
    % el2_vs2_s = res.el2_vs2(:,sm);
    % el2_vs3_s = res.el2_vs3(:,sm);

    xx = sNs(:);
    yy = el2_s(:);    
    % yy = el2_fs1_s(:);
    ok = isnan(yy)|yy<=0;     
    xx(ok)=[];  yy(ok)=[];
    
    coeff_exp = polyfit(xx, log(yy), 1);
    a = coeff_exp(1);
    b = coeff_exp(2);
    fprintf('Exponential fit:   E(N) ≃ %.2f · exp(%.2f·N^{1/d})\n',  exp(b), a);
    yfit_exp = exp(b) * exp( a * xx );

    h = figure;  %grid on;
    set(h, 'Color', 'none');          
    ax = gca;
    set(ax, 'Color', 'none','FontSize', 12,'LineWidth', 1.2);     
    marks = {'-o','-s','-^','-x','-+'};
    semilogy(sNs, el2_poly_s, marks{1}, 'LineWidth',2);
    hold on;
    % semilogy(sNs, el2_diag_s, marks{2}, 'LineWidth',2);
    semilogy(sNs, el2_s,  marks{3}, 'LineWidth',2);
    % semilogy(sNs, el2_fs1_s,  marks{4}, 'LineWidth',2);
    % semilogy(sNs, el2_fs3_s,  marks{5}, 'LineWidth',2);
    semilogy(xx, yfit_exp, '--','LineWidth',2, 'Color',[.5 .5 .5]);

    ax.LineWidth = 2;
    ax.FontSize = 16;
    ax.XAxis.FontWeight = 'bold';
    ax.YAxis.FontWeight = 'bold';


    if dim==1
        xlabel('N','Interpreter','tex','FontSize',16,'FontWeight','bold');
    else
        xlabel(sprintf('N^{1/%d}', dim),'Interpreter','tex','FontSize',16,'FontWeight','bold');
    end
    ylabel('Relative l_2 error','Interpreter','tex','FontSize',16,'FontWeight','bold');
    % legend({'PLS','Diag','K_t = 1e12','K_t = 1e8','K_t = 1e4', 'exp trendline'}, 'Location','northeast','Interpreter','tex','FontWeight','bold','FontSize',16);
    legend({'PLS','FS','exp trendline'}, 'Location','best','Interpreter','tex','FontWeight','bold','FontSize',14);
    % legend({'PLS','Diag','FS','exp trendline'}, 'Location','best','Interpreter','tex','FontWeight','bold','FontSize',14);

    if sm == 1
        title(sprintf('Relative l_2 error vs. N^{1/d}, C^2(R^3) Wendland Kernel'),'Interpreter','tex','FontWeight','bold','FontSize',16);
    elseif sm==2
        title(sprintf('Relative l_2 error vs. N^{1/d}, C^4(R^3) Wendland Kernel'),'Interpreter','tex','FontWeight','bold','FontSize',16);
    elseif sm==3
        title(sprintf('Relative l_2 error vs. N^{1/d}, C^6(R^3) Wendland Kernel'),'Interpreter','tex','FontWeight','bold','FontSize',16); 
    end

    fprintf('Exponential fit:   E(N) = %.2f · exp(%.2f·N^{1/d})\n',  exp(b), a);


    % force all plots to use the same x‐ and y-ranges
    xlim([xmin xmax]);
    ylim([err_min err_max]);

    %set(gca,'FontSize',12,'FontWeight', 'bold');
    export_fig(gcf, fullfile(res.results_dir,sprintf('error_vs_N_fs_s%d.png',sm)),'-png','-r300','-transparent');
    close(h);

    % xx = sNs(:);
    % yy = el2_fs1_s(:);        
    % ok = isnan(yy)|yy<=0;     
    % xx(ok)=[];  yy(ok)=[];
    % 
    % coeff_exp = polyfit(xx, log(yy), 1);
    % a = coeff_exp(1);
    % b = coeff_exp(2);
    % fprintf('Exponential fit:   E(N) ≃ %.2f · exp(%.2f·N^{1/d})\n',  exp(b), a);
    % yfit_exp = exp(b) * exp( a * xx );
    % 
    % h = figure;  %grid on;
    % set(h, 'Color', 'none');          
    % ax = gca;
    % set(ax, 'Color', 'none','FontSize', 12, 'LineWidth', 1.2, 'TickLabelInterpreter', 'latex');
    % marks = {'-o','-s','-^','-x','-+'};
    % semilogy(sNs, el2_poly_s, marks{1}, 'LineWidth',1.2);
    % hold on;
    % semilogy(sNs, el2_diag_s, marks{2}, 'LineWidth',1.2);
    % semilogy(sNs, el2_vs1_s,  marks{3}, 'LineWidth',1.2);
    % semilogy(sNs, el2_vs2_s,  marks{4}, 'LineWidth',1.2);
    % semilogy(sNs, el2_vs3_s,  marks{5}, 'LineWidth',1.2);
    % semilogy(xx, yfit_exp, '--','LineWidth',1.5, 'Color',[.5 .5 .5]);
    % 
    % 
    % if dim==1
    %     xlabel('N','Interpreter','tex','FontSize',14,'FontWeight','bold');
    % else
    %     xlabel(sprintf('N^{1/%d}', dim),'Interpreter','tex','FontSize',14,'FontWeight','bold');
    % end
    % ylabel('Relative l_2 error','Interpreter','tex','FontSize',14,'FontWeight','bold');
    % legend({'PLS','Diag','K_t = 1e12','K_t = 1e8','K_t = 1e4','exp trendline'}, 'Location','northeast','Interpreter','tex','FontWeight','bold','FontSize',14);
    % % legend({'PLS','Diag','FS', 'exp trendline'}, 'Location','best','Interpreter','tex','FontWeight','bold','FontSize',14);
    % 
    % if sm == 1
    %     title(sprintf('Relative l_2 error vs. N^{1/d}, C^2(R^3) Wendland Kernel'),'Interpreter','tex','FontWeight','bold','FontSize',14);
    % elseif sm==2
    %     title(sprintf('Relative l_2 error vs. N^{1/d}, C^4(R^3) Wendland Kernel'),'Interpreter','tex','FontWeight','bold','FontSize',14);
    % elseif sm==3
    %     title(sprintf('Relative l_2 error vs. N^{1/d}, C^6(R^3) Wendland Kernel'),'Interpreter','tex','FontWeight','bold','FontSize',14); 
    % end
    % 
    % % force all plots to use the same x‐ and y-ranges
    % xlim([xmin xmax]);
    % ylim([err_min err_max]);
    % 
    % % set(gca,'FontSize',12,'FontWeight', 'bold');
    % export_fig(gcf, fullfile(res.results_dir,sprintf('error_vs_N_vs_s%d.png',sm)),'-png','-r300','-transparent');
    % close(h);
end

