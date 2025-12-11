function [el2, elinf, a_time, e_time, c_poly, cond_num, c_rbf, sparsity] = ...
    PHS_poly_man(x, y, ell, xe, p, ye_true, pinv_tol)
%PHS_poly_man  Global PHS + polynomial augmentation on manifolds (KKT + pseudoinverse).
%
% [el2, elinf, a_time, e_time, c_poly, cond_num, c_rbf, sparsity] =
%     PHS_poly_man(x, y, ell, xe, p, ye_true, pinv_tol)
%
% Inputs
%   x         : N-by-d manifold nodes (ambient coords in R^3)
%   y         : N-by-1 data values at x
%   ell       : polynomial total degree in ambient R^3
%   xe        : Ne-by-d evaluation points (ambient coords in R^3)
%   p         : odd PHS power (e.g., 5, 7, 9)  -> phi(r) = r.^p
%   ye_true   : Ne-by-1 ground truth at xe (for error reporting)
%   pinv_tol  : (optional) tolerance for pseudoinverse (default 1e-11)
%
% Outputs
%   el2, elinf: relative L2 / Linf errors on xe
%   a_time    : assembly+solve wall time
%   e_time    : evaluation wall time
%   c_poly    : polynomial coefficients (size M-by-1)
%   cond_num  : condest of KKT (approx; NaN if too costly)
%   c_rbf     : RBF coefficients (size N-by-1)
%   sparsity  : 1 - nnz(G)/numel(G)  (dense -> ~0)

if nargin < 7 || isempty(pinv_tol), pinv_tol = 1e-11; end

N = size(x,1);
dim = size(x,2); %#ok<NASGU>  % ambient dimension (expected 3)

% Build ambient polynomial Vandermonde V
alph = 0;  % Legendre
recurrence = @(Ndeg) jacobi_recurrence(Ndeg, alph, alph);
a = total_degree_indices(3, ell);            % multi-indices in R^3
[xc, c1, c2] = affine_scale(x);
V  = mpoly_eval(xc, a, recurrence);
M  = size(V,2);

% Dense PHS kernel blocks
tic
K  = phs_kernel(x, x, p);          % N x N
G  = [K, V; V.', zeros(M,M)];      % (N+M) x (N+M)
rhs = [y; zeros(M,1)];

% Pseudoinverse solve to handle rank deficiency in V on manifolds
sol = pinv(G, pinv_tol) * rhs;

c_rbf = sol(1:N);
c_poly = sol(N+1:end);

a_time = toc;

% crude sparsity / conditioning measures
sparsity = 1 - nnz(G)/numel(G);
try
    % condest is for sparse; for dense G this may be expensive.
    % Use rcond surrogate if condest errors.
    cond_num = 1 / rcond(G);
catch
    cond_num = NaN;
end

% Evaluate on xe
tic
[xe_c, ~, ~] = affine_scale(xe, c1, c2);
Ke = phs_kernel(xe, x, p);         % Ne x N
Ve = mpoly_eval(xe_c, a, recurrence);
ye = Ke * c_rbf + Ve * c_poly;
e_time = toc;

% Errors
el2 = norm(ye - ye_true) / max(1e-16, norm(ye_true));
elinf = norm(ye - ye_true, inf) / max(1e-16, norm(ye_true, inf));

end

% helpers

function [xc, c1, c2] = affine_scale(x, c1_in, c2_in)
% Affinely scale each coordinate of x to [-1,1] using xmin/xmax.
if nargin < 2
    xmax = max(x, [], 1);
    xmin = min(x, [], 1);
    c1 = 0.5 .* (xmax + xmin);
    c2 = 2 ./ max(xmax - xmin, eps);
else
    c1 = c1_in;
    c2 = c2_in;
end
xc = (x - c1) .* c2;
end

function K = phs_kernel(X, Y, p)
% Full PHS kernel matrix K_ij = ||X_i - Y_j||^p (odd p).
XX = sum(X.^2, 2);
YY = sum(Y.^2, 2);
D2 = max(XX + YY.' - 2*(X*Y.'), 0);  % squared distances, guarded
R = sqrt(D2);
% By convention r^p with r=0 gives 0; no logs needed for odd p.
K = R.^p;
end
