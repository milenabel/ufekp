h = figure;  %grid on;    
set(h, 'Color', 'none');          
ax = gca;
set(ax, 'Color', 'none','FontSize', 12, 'LineWidth', 1.2, 'TickLabelInterpreter', 'latex');
semilogy(sqrt(N), el2, '-o', 'LineWidth',1.2);
xlabel('N^{1/2}','FontSize',14,'FontWeight','bold');    
ylabel('Relative l_2 error','FontSize',14,'FontWeight','bold'); 
title('Spatial Convergence, C^2(R^3) Wendland kernel','FontWeight','bold','FontSize',14);
export_fig torus_error_vs_N.png -png -r300 -transparent
   % close(h);