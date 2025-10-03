function[d, varargout] = mjacobi_symm_faster_differentiation(c, a, alph, dorder, varargin)
% [d, [a_structure, dimension_jumps,Cmat]] = 
%     mjacobi_symm_fast_differentiation(c, a, alph, dorder, [a_structure, dimension_jumps,Cmat])
%
% With c a vector of expansion coefficients in a multivariate Jacobi polynomial
% family with symmetric parameters (alph,bet) = (alph, alph), this computes a
% vector of coefficients d corresonding to the derivative of the polynomial
% defined by c.
%
% The expansion is defined by the multi-indices a.
%
% The derivative computed is given by dorder. E.g., in 3 dimensions, a fourth
% derivative characterized by one derivative in each of x and z and two
% derivatives in y has dorder = [1 2 1].
%
% The complexity of this procedure is formally linear in the maximum degree in
% any direction, and also has complexity linear in the number of multi-indices
% in the "deflated" index set.

if nargin < 6
  [a_structure, dimension_jumps] = mindex_derivative_analysis(a);
else
  [a_structure, dimension_jumps] = deal(varargin{1}, varargin{2});
end

if nargin < 7
%  N = size(a_structure, 1);
%
%  C2 = jacobi_iterated_connection_coeffs(N, alph);
%  Connmat = spdiags(C2, [0 2], N, N+2); 
%  Connmat(:,N+1:N+2) = [];
%
%  temp = jacobi_derivative_coefficients(alph, alph, N);
%  Cmat = spdiags(1./temp(2:N), 0, N-1,N-1)*Connmat(1:N-1,1:N-1);

  Cmat = jacobi_sparse_diffmat_inv(size(a_structure, 1), alph);
else
  lc = varargin{3};
  uc = varargin{4};
  pc = varargin{5};
  qc = varargin{6};
end

if all(dorder == 0)
  d = c;
  return
end


N = size(a_structure, 1);
d = [c; zeros([1 size(c,2)])];
for q = 1:length(dorder)

  if dorder(q) > 0

    colinds = dimension_jumps(q):(dimension_jumps(q+1)-1);

    D = d(a_structure(:,colinds));

    for qq = 1:dorder(q)
      rc = D(2:N,:);
      sc = qc*(uc\(lc\(pc*rc)));
      D = [sc; zeros([1 size(D,2)])];
    end

    d(a_structure(:,colinds)) = D;

  end

end

d(end) = [];
