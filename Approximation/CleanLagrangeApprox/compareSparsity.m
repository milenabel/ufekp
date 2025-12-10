function compareSparsity(function_list, strategy, smoothness, base_results_dir)
    if nargin<4 || isempty(base_results_dir)
        base_results_dir = fullfile('csrbfjacobi/code/Approximation/CleanLagrangeApprox/results/');
    end

    % a few nice markers/colors
    markers = {'-o','-s','-^','--d','-.x','-.+'};
    colors  = lines(numel(function_list));

    figure; hold on; grid on;
    for k = 1:numel(function_list)
        fn = function_list{k};
        res = loadResults(fn, base_results_dir);  
        S = res.(['sparsity_' strategy])(:, smoothness);
        E = res.(['el2_'       strategy])(:, smoothness);
        plot(S, E, markers{mod(k-1,numel(markers))+1}, ...
             'Color', colors(k,:), ...
             'LineWidth',1.2, ...
             'DisplayName', fn);
    end
    xlabel('Achieved sparsity','FontSize',14);
    ylabel('Relative \ell_2 error','FontSize',14);
    title( sprintf('Strategy %s at Wendland C%d', upper(strategy), smoothness*2), ...
           'FontSize',16 );
    legend('Location','best','Interpreter','none');
    set(gca,'FontSize',12);

    % export if you like
    export_fig(gcf, fullfile(base_results_dir,...
                 sprintf('compare_sparsity_%s_C%d.png',strategy,smoothness)), ...
               '-png','-r300');
end
