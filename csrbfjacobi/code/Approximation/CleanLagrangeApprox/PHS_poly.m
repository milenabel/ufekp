function [el2, elinf, a_time, e_time, coeffs] = PHS_poly(x, y, ell, xe, alph, ye_true, m)
%PHS_poly  Global PHS(r^(2m+1)) + polynomial interpolant (dense KKT)
%
% Usage:
%   [el2, elinf, a_time, e_time, coeffs] = ...
%       PHS_poly(x, y, ell, xe, alph, ye_true, m)
%
% Inputs:
%   x        : N x d data sites
%   y        : N x 1 samples
%   ell      : polynomial total degree
%   xe       : Ne x d evaluation sites
%   alph     : Jacobi alpha (Legendre -> alph = 0), matches PLS
%   ye_true  : Ne x 1 true values at xe (for errors)
%   m        : (optional) PHS power parameter: phi(r)=r^(2m+1)
%              If omitted, defaults to m = max(1, floor(ell))
%
% Outputs:
%   el2, elinf : relative L2 / Linf errors on xe
%   a_time     : wall time (assembly+solve)
%   e_time     : wall time (evaluation)
%   coeffs     : struct with fields:
%                  c : N x 1 RBF coefficients
%                  d : M x 1 poly coefficients
%                  c1,c2 : scaling used for Vandermonde (like PLS)
%                  a : multi-index for poly basis
%                  recurrence : handle used

    if nargin < 7 || isempty(m)
        m = max(1, floor(ell));  % ensure CPD order (m+1) ≤ ℓ+1  → ℓ ≥ m
    end

    p = 2*m + 1;  % PHS degree (must be odd: 5,7,9,...)

    [N, dim] = size(x);
    Ne = size(xe,1);

    % Build poly Vandermonde on scaled coordinates (like PLS.m)
    % Scale x to [-1,1]^d for orthogonal polys
    xmax = max(x, [], 1);
    xmin = min(x, [], 1);
    c1 = 0.5 .* (xmax + xmin);
    c2 = 2 ./ (xmax - xmin);
    xc  = (x  - repmat(c1, [N  1])) .* repmat(c2, [N  1]);
    xce = (xe - repmat(c1, [Ne 1])) .* repmat(c2, [Ne 1]);

    recurrence = @(Nmax) jacobi_recurrence(Nmax, alph, alph); % Legendre if alph=0
    a = total_degree_indices(dim, ell); % M x dim multi-index
    P = mpoly_eval(xc, a, recurrence); % N x M
    Pe = mpoly_eval(xce, a, recurrence); % Ne x M
    M = size(P,2);

    % Enforce unisolvency: ℓ ≥ m for r^(2m+1)
    if ell < m
        warning('PHS_poly: raising ell from %d to %d to satisfy ℓ ≥ m.', ell, m);
        ell = m;
        a = total_degree_indices(dim, ell);
        P = mpoly_eval(xc, a, recurrence);
        Pe = mpoly_eval(xce, a, recurrence);
        M = size(P,2);
    end

    % Build dense PHS kernel blocks
    % phi(r) = r^(2m+1), where 2m+1 = p
    % Distances
    rd = max(bsxfun(@plus,sum(x.*x,2),sum(x.*x,2)') - 2*(x*x'),0);
    A = sqrt(rd) .^ p; % N x N

    % Assemble and solve KKT system
    % [ A  P ] [c] = [ y ]
    % [ P' 0 ] [d]   [ 0 ]
    K = [A, P; P.', zeros(M,M)];
    rhs = [y; zeros(M,1)];

    tic
    z = K \ rhs;
    a_time = toc;

    c = z(1:N);
    d = z(N+1:end);

    % Free big assembly-time arrays 
    clear K rhs A rd P

    % Evaluate on xe
    tic
    rde = max(bsxfun(@plus,sum(xe.*xe,2),sum(x.*x,2)') - 2*(xe*x'),0);

    Ae = sqrt(rde) .^ p;  % Ne x N
    ye = Ae * c + Pe * d;
    e_time = toc;

    % Relative errors
    el2 = norm(ye - ye_true)   / norm(ye_true);
    elinf = norm(ye - ye_true,inf) / norm(ye_true,inf);

    % Pack coeffs
    coeffs = struct('c', c, 'd', d, 'c1', c1, 'c2', c2, 'a', a, ...
                    'recurrence', recurrence, 'm', m,'phsdeg', p);
end