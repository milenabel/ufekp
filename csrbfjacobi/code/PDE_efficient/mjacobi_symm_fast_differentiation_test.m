% This script 
% (1) generates an interpolatory approximation via an expansion in
% multivariate polynomials that are tensor products of univariate Jacobi
% polynomials with symmetric parameters (alph, bet) = (alph, alph).
%
% (2) uses a "fast" algorithm to transform the expansion coefficients to
% expansion coefficients of a derivative.
%
% (3) computes errors for this coefficient transform versus a pointwise derivative.
%
% Note that (3) is not a computation of error against the actual function
% derivative, but instead against the derivative of the interpolatory
% approximation. I.e., this should essentially be machine precision.

clear
close all

%% Problem setup
dim = 2;
region_ind = @(xx) sum(xx.^2, 2) <= 1;
f = @(xx) sin(xx(:,1) + 2*xx(:,2));

% Jacobi polynomial parameter
alph = 0;
recurrence = @(N) jacobi_recurrence(N, alph, alph);

%% Selecting a grid for interpolation: this selection process is pretty good, but not optimal.
M = 1e5;
a = total_degree_indices(dim, 100);
%a = hyperbolic_cross_indices(dim, 15); % Hyperbolic cross works, too
N = size(a, 1);
% x_candidates = zeros([0 dim]);
% while size(x_candidates,1) < M
%   temp = rand([M dim])*2 - 1;
%   flags = region_ind(temp);
%   x_candidates = [x_candidates; temp(flags,:)];
% end

% st = load('DiskHolePoissonNodes.mat');
% xi = st.fullintnodes{5};
% xbo = st.obdrynodes{5};
% xbi = st.ibdrynodes{5};
% x_candidates = [xi;xbo;xbi];

st = load('DomainPoissonNodes2.mat');
xi = st.fullintnodes{6};
xb = st.bdrynodes{6};
x_candidates = [xi;xb];

V = mpoly_eval(x_candidates, a, recurrence);
weights = 1./sqrt(sum(V.^2, 2));
V = repmat(weights, [1 N]).*V;

[~,~,e] = qr(V', 0);

% Identify interpolation grid and update weights and V
x = x_candidates(e(1:N),:);
weights = weights(e(1:N));
V = V(e(1:N),:);

%% Generate an interpolatory approximation in coefficient space
c = V\(weights.*f(x));

%% Effect various differentiations in coefficient space
% One-time setup for differentiation
[a_structure, dimension_jumps] = mindex_derivative_analysis(a);
Cmat = jacobi_sparse_diffmat_inv(size(a_structure, 1), alph);
[lc,uc,pc,qc] = lu(Cmat);

% Utilize "fast" algorithm
tic
dx   = mjacobi_symm_faster_differentiation(c, a, alph, [1 0], a_structure, dimension_jumps, lc,uc,pc,qc);
dy   = mjacobi_symm_faster_differentiation(c, a, alph, [0 1], a_structure, dimension_jumps, lc,uc,pc,qc);
dxy  = mjacobi_symm_faster_differentiation(c, a, alph, [1 1], a_structure, dimension_jumps, lc,uc,pc,qc);
dxxy = mjacobi_symm_faster_differentiation(c, a, alph, [2 1], a_structure, dimension_jumps, lc,uc,pc,qc);
toc
%% Compare coefficient-space differentiation versus explicit pointwise differentiation:
errx   = norm( mpoly_eval(x, a, recurrence)*dx   - mpoly_eval(x, a, recurrence, [1 0])*c );
erry   = norm( mpoly_eval(x, a, recurrence)*dy   - mpoly_eval(x, a, recurrence, [0 1])*c );
errxy  = norm( mpoly_eval(x, a, recurrence)*dxy  - mpoly_eval(x, a, recurrence, [1 1])*c );
errxxy = norm( mpoly_eval(x, a, recurrence)*dxxy - mpoly_eval(x, a, recurrence, [2 1])*c );

fprintf('Errors are \n%1.4e\n%1.4e\n%1.4e\n%1.4e\n', errx, erry, errxy, errxy);
