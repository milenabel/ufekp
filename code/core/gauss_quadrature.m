function[x,w] = gauss_quadrature(a,b,N)
% gauss_quadrature -- Computes Gauss quadrature from recurrence coefficients
%
% [x,w] = gauss_quadrature(a,b,N)
%
%     Computes N Gauss quadrature nodes (x) and weights (w) from length-N arrays
%     a and b containing the standard orthonormal recurrence coefficients of the
%     family. 

if N < 1
  x = [];
  w = [];
  return
end

a = a(:); b = b(:);
if nargin < 3
  N = min(length(a), length(b));
end

% Form Jacobi matrix
J = spdiags([sqrt([b(2:N);0]) a(1:N) sqrt(b(1:N))], -1:1, N, N);
x = eig(full(J));

w = poly_eval(a, b, x, N-1);
w = 1./sum(w.^2,2);
w(~isfinite(w)) = 0;
