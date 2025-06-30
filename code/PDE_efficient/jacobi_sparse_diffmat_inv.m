function[Cmat] = jacobi_sparse_diffmat_inv(N, alph);
% [Cmat] = jacobi_sparse_diffmat_inv(N, alph);
%
% Returns an (N-1) x (N-1) matrix Cmat, for which the following operation on a
% size-N vector of Jacobi(alph,alph) expansion coefficients effects
% differentiation:
%
% d = [Cmat\c(2:N); 0];

C2 = jacobi_iterated_connection_coeffs(N, alph);
Connmat = spdiags(C2, [0 2], N, N+2); 
Connmat(:,N+1:N+2) = [];

temp = jacobi_derivative_coefficients(alph, alph, N);
Cmat = spdiags(1./temp(2:N), 0, N-1,N-1)*Connmat(1:N-1,1:N-1);
