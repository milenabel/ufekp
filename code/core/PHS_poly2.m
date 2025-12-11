function [el2, elinf, a_time, e_time, coeffs] = PHS_poly2(x, y, ell, xe, alph, ye_true, m)
%PHS_poly2  Global PHS(r^(2m+1)) + polynomial via LU-Schur (no KKT assemble)
%
% Same interface as PHS_poly.m but uses:
%   [L,U,P] = lu(A), with P*A = L*U  (row pivoting, call P as Q in advisor note)
% and Schur complement:
%   S = V.' * A^{-1} * V,   rhs = V.' * A^{-1} * y
%   solve S*d = rhs;  then   c = A^{-1} * (y - V*d)
%
% Inputs/Outputs: identical to PHS_poly.m

    if nargin < 7 || isempty(m)
        m = max(1, floor(ell));  % ensure ℓ ≥ m for CPD unisolvency
    end
    p = 2*m + 1;

    [N, dim] = size(x);
    Ne = size(xe,1);

    % Polynomial (PLS-style) scaling and Vandermonde
    xmax = max(x, [], 1);
    xmin = min(x, [], 1);
    c1 = 0.5 .* (xmax + xmin);
    c2 = 2 ./ (xmax - xmin);
    xc  = (x  - repmat(c1, [N  1])) .* repmat(c2, [N  1]);
    xce = (xe - repmat(c1, [Ne 1])) .* repmat(c2, [Ne 1]);

    recurrence = @(Nmax) jacobi_recurrence(Nmax, alph, alph); % Legendre if alph=0
    a = total_degree_indices(dim, ell);
    V  = mpoly_eval(xc,  a, recurrence);   % N x M
    Ve = mpoly_eval(xce, a, recurrence);   % Ne x M
    M  = size(V,2);

    % Enforce ℓ ≥ m (unisolvency for r^(2m+1))
    if ell < m
        warning('PHS_poly2: raising ell from %d to %d to satisfy ℓ ≥ m.', ell, m);
        ell = m;
        a  = total_degree_indices(dim, ell);
        V  = mpoly_eval(xc,  a, recurrence);
        Ve = mpoly_eval(xce, a, recurrence);
        M  = size(V,2);
    end

   
    % Dense PHS kernel A (phi(r) = r^p, p odd)
    rd = max(bsxfun(@plus,sum(x.*x,2),sum(x.*x,2)') - 2*(x*x'),0);
    A  = sqrt(rd) .^ p;  % N x N

    % LU with row pivoting: P*A = L*U
    tic;
    [L,U,Pmat] = lu(A);  
    % A^{-1} * v   =   U \ ( L \ (Pmat * v) )

    % Build A^{-1} * y and A^{-1} * V without forming A^{-1}
    Ay_inv_y = U \ ( L \ (Pmat * y) );   % N x 1
    Ay_inv_V = U \ ( L \ (Pmat * V) );   % N x M

    % Schur system: (V.' * A^{-1} * V) d = V.' * A^{-1} * y
    S   = V.' * Ay_inv_V;               % M x M (symmetric)
    rhs = V.' * Ay_inv_y;               % M x 1

    % Solve for polynomial coeffs d (use backslash as requested)
    d = S \ rhs;

    % Then c = A^{-1} (y - V*d)
    c = U \ ( L \ (Pmat * (y - V*d)) );
    a_time = toc;

    % Evaluate on xe
    tic;
    rde = max(bsxfun(@plus,sum(xe.*xe,2),sum(x.*x,2)') - 2*(xe*x'),0);
    Ae  = sqrt(rde) .^ p;                % Ne x N
    ye  = Ae * c + Ve * d;
    e_time = toc;

    % Relative errors
    el2   = norm(ye - ye_true)   / norm(ye_true);
    elinf = norm(ye - ye_true,inf) / norm(ye_true,inf);

    % Pack coeffs (for parity with PHS_poly)
    coeffs = struct('c', c, 'd', d, 'c1', c1, 'c2', c2, 'a', a, ...
                    'recurrence', recurrence, 'm', m, 'phsdeg', p, ...
                    'solver','LU-Schur', 'note','P*A=L*U; Schur S=V''A^{-1}V');
end
