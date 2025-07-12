%% Test our method on manifolds. The polynomial matrix should be
%% rank-deficient on any algebraic Riemannian manifold.
%% We will tackle this with a rank-revealing QR factorization of
%% matrices that involve P.

%% Numbers of points on the torus for a convergence study
n = [5,10,15,20,25,30,35,40];

%% Number of eval points
ne = 50;
xe = computeHexagonTorusNodes(ne);

%% Target function
rb = 1/3;
%f = @(x) exp( (x(:,1) + x(:,2) + x(:,3)));
%f = @(x) exp( (x(:,1) + x(:,2) + x(:,3)).^2);
%f = @(x) (acos(x(:,1))<rb).* (1 + cos(pi*acos(x(:,1))/rb))/2; %C^1
f = @(x) (x(:,1).^2).*abs(x(:,1)) + x(:,2).^2.*abs(x(:,2)) + x(:,3).^2.*abs(x(:,3)); 

%% Evaluate this target at the eval points
f_xe = f(xe);
f_xe(isnan(f_xe)) = 0;

%% Define a Wendland RBF
rbf = @(e,r) -r.^4.*(4.*r - 5*spones(r)); %C2
%rbf = @(e,r) (r.^6.*(35.*r.^2 - 88.*r + 56*spones(r)))./3; %C4
%rbf = @(e,r) r.^8.*(66*spones(r) - 154*r + 121*r.^2 - 32*r.^3); %C6

%% Convergence study
for k=1:length(n)
    %% Generate the interpolation nodes
    x = computeHexagonTorusNodes(n(k));
    Ns(k,1) = length(x);

    %% Get a tree
    tree = KDTreeSearcher(x);

    %% Evaluate the target function there
    f_x = f(x);  
    f_x(isnan(f_x)) = 0;

    %% Get the polynomial degree. We'll use ambient polynomials in R^3
    ell = floor(nthroot(length(x),3));   

    %% Do poly least squares with rank-revealing QR
    [el2_poly(k,1),elinf_poly(k),a_time_poly(k),e_time_poly(k),~] = PLSManifold(x,f_x,ell,xe,0,f_xe);
    

    %% Do a CSRBF + poly fit that uses a rank-revealing QR
   [el2(k,1),elinf(k,1),a_time(k,1),e_time(k,1),~,cond_num(k,1),~, sparsity(k,1)] = CSRBFGenManifold(x,f_x,ell,xe,0,rbf,7,tree,f_xe);   
end

%% Save results 
% Timestamp for uniqueness

sN = sqrt(n);
function_name = 'torus';
timestamp = datestr(datetime('now'), 'yyyyMMdd_HHmmss');
results_dir = fullfile('results/', sprintf('%s', function_name),'/high');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end
results_filename = fullfile(results_dir, sprintf('results_%s.mat', function_name));


% Save everything
save(results_filename, ...
    'el2_poly', 'elinf_poly', 'a_time_poly', 'e_time_poly', 'coeff_exp', 'ell', ...
    'el2', 'elinf', 'a_time', 'e_time', 'cond_num', 'sparsity', ...
    'sN', 'function_name', 'timestamp');

%% Plot
xx = sN(:);
yy = el2(:); 
ok = isnan(yy)|yy<=0;
xx(ok)=[];  yy(ok)=[];

coeff_exp = polyfit(xx, log(yy), 1);
a = coeff_exp(1);
b = coeff_exp(2);
fprintf('Exponential fit:   E(N) ≃ %.2f · exp(%.2f·N^{1/d})\n',  exp(b), a);
yfit_exp = exp(b) * exp( a * xx );

h=figure;
set(h, 'Color', 'none');          
ax = gca;
set(ax, 'Color', 'none','FontSize', 12,'LineWidth', 1.2,'TickLabelInterpreter', 'latex');     
marks = {'-o','-s','-^','-x','-+'};
semilogy(sqrt(n), el2(:),marks{1}, 'LineWidth',1.2);
hold on;
semilogy(sqrt(n), el2_poly(:),marks{2}, 'LineWidth',1.2);
semilogy(xx, yfit_exp, '--','LineWidth',1.5, 'Color',[.5 .5 .5]);
xlabel('N^{1/2}','Interpreter','tex','FontSize',14,'FontWeight','bold');
ylabel('Relative l_2 error','Interpreter','tex','FontSize',14,'FontWeight','bold');
legend({'FS','PLS', 'exp trendline'}, 'Location','best','Interpreter','tex','FontWeight','bold','FontSize',14);
title(sprintf('Relative l_2 error vs. N^{1/d}, C^2(R^3) Wendland Kernel'),'Interpreter','tex','FontWeight','bold','FontSize',14);
export_fig(gcf, fullfile(results_dir,'torus_error_vs_N.png'),'-png','-r300','-transparent');
close(h);