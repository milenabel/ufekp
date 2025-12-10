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


% PHS+poly controls 
ep = 7;        
p_phs = 3;        
pinv_tol = 1e-11;   
Ns = zeros(numel(n), 1);   


%% Convergence study
for k=1:length(n)
    %% Generate the interpolation nodes
    x = computeHexagonTorusNodes(n(k));
    Ns(k) = size(x,1);

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
    [el2(k,1),elinf(k,1),a_time(k,1),e_time(k,1),~,cond_num(k,1),~, sparsity(k,1)] = CSRBFGenManifold(x,f_x,ell,xe,0,rbf,ep,tree,f_xe);   

    %% Do a PHS+poly fit that uses a dense KKT with pseudoinverse
    [el2_phs(k,1),elinf_phs(k,1),a_time_phs(k,1),e_time_phs(k,1),~,cond_phs(k,1),~,spars_phs(k,1)] = PHS_poly_man(x, f_x, ell, xe, p_phs, f_xe, pinv_tol);

end

%% Save results 
% Timestamp for uniqueness

sN = sqrt(Ns);
function_name = 'torus';
timestamp = datestr(datetime('now'), 'yyyyMMdd_HHmmss');
results_dir = fullfile('results/', sprintf('%s', function_name),'/high');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end
results_filename = fullfile(results_dir, sprintf('results_%s.mat', function_name));

%% Plot
xx = sN(:);
yy = el2(:); 
% ok = isnan(yy)|yy<=0;
% xx(ok)=[];  yy(ok)=[];


valid = isfinite(xx) & isfinite(yy) & (yy > 0);
xx = xx(valid);
yy = yy(valid);

if numel(xx) >= 2
    coeff_exp = polyfit(xx, log(yy), 1);
    a = coeff_exp(1); b = coeff_exp(2);
    yfit_exp = exp(b) * exp(a * xx);
    fprintf('Exponential fit:   E(N) ≃ %.2f · exp(%.2f·N^{1/d})\n', exp(b), a);
else
    a = NaN; b = NaN; yfit_exp = [];
    fprintf('Exponential fit skipped (insufficient valid points)\n');
end

% coeff_exp = polyfit(xx, log(yy), 1);
% a = coeff_exp(1);
% b = coeff_exp(2);
fprintf('Exponential fit:   E(N) ≃ %.2f · exp(%.2f·N^{1/d})\n',  exp(b), a);
% yfit_exp = exp(b) * exp( a * xx );

% Save everything
save(results_filename, ...
    'el2_poly', 'elinf_poly', 'a_time_poly', 'e_time_poly', 'coeff_exp', 'ell', ...
    'el2', 'elinf', 'a_time', 'e_time', 'cond_num', 'sparsity', ...
    'el2_phs','elinf_phs','a_time_phs','e_time_phs','cond_phs','spars_phs', ...
    'a','b','ep','p_phs','pinv_tol', ...
    'Ns','sN', 'function_name', 'timestamp');

h=figure;
set(h, 'Color', 'none');  
set(groot,'defaultAxesTickLabelInterpreter','tex');
set(groot,'defaultTextInterpreter','tex');
set(groot,'defaultLegendInterpreter','tex');
ax = gca;
ax.TickLabelInterpreter = 'tex'; 
ax.Color = 'none';
ax.FontWeight = 'bold';
ax.FontSize = 16;
ax.LineWidth = 2;
% ax.XAxis.FontWeight = 'bold'; 
% ax.YAxis.FontWeight = 'bold';
ax.XTickLabelMode = 'auto';
ax.YTickLabelMode = 'auto';
drawnow;
% Styles 
% PLS
cPLS = [0.000 0.447 0.741]; mPLS = 'o'; msPLS = 6; mfPLS = 'w'; lsPLS = '-';
% FS label for manifolds = 'FS'
cFS  = [0.850 0.325 0.098]; mFS  = '^'; msFS  = 6; mfFS  = 'w'; lsFS  = '-';
% PHS+poly 
cPHS = [0.000 0.75  1.000]; mPHS = '>'; msPHS = 6; mfPHS = 'w'; lsPHS = '-';
% Trend
cTr  = [0.5 0.5 0.5]; lsTr = '--';
semilogy(sqrt(Ns), el2(:),lsFS,'Color',cFS, 'Marker',mFS, 'MarkerSize',msFS, 'MarkerFaceColor',mfFS, 'LineWidth',1.2); 
hold on;
semilogy(sqrt(Ns), el2_poly(:),lsPLS,'Color',cPLS,'Marker',mPLS,'MarkerSize',msPLS,'MarkerFaceColor',mfPLS,'LineWidth',1.2);
% PHS+poly curve
semilogy(sqrt(Ns), el2_phs(:), lsPHS,'Color',cPHS,'Marker',mPHS,'MarkerSize',msPHS,'MarkerFaceColor',mfPHS,'LineWidth',1.2);
semilogy(xx, yfit_exp, '--','LineWidth',1.5, 'Color',[.5 .5 .5]);
xlabel('N^{1/2}','Interpreter','tex','FontSize',16,'FontWeight','bold');
ylabel('Relative l_2 error','Interpreter','tex','FontSize',16,'FontWeight','bold');
legend({'FS','PLS','PHS+poly','exp trendline'}, 'Location','best','Interpreter','tex','FontWeight','bold','FontSize',14);
title(sprintf('Relative l_2 error vs. N^{1/d}, C^2(R^3) Wendland Kernel'),'Interpreter','tex','FontWeight','bold','FontSize',16);
export_fig(gcf, fullfile(results_dir,'torus_error_vs_N.png'),'-png','-r300','-transparent');
close(h);