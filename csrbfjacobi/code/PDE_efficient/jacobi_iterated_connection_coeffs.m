function[C] = jacobi_iterated_connection_coeffs(N, alph)
% C = jacobi_iterated_connection_coeffs(N, alph)
%
% Returns an N x 2 matrix C such that
%
% (1 - x^2) p_n^{(alph+1,alph+1)}(x) =   C(n+1, 1) p_n^{(alph,alph)}(x) +
%                                        C(n+1, 2) p_{n+2}^{(alph,alph)}(x).

C = zeros([N 2]);

C1 = jacobi_connection_coeffs(N, alph, alph+1, 1);
C2 = jacobi_connection_coeffs(N+1, alph, alph, -1);

C = [C1(:,1).*C2(1:N,1) C1(:,2).*C2(2:N+1,2)];
