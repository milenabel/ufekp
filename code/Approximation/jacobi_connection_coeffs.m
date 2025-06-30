function[C] = jacobi_connection_coeffs(N, alph, bet, sgn)
% C = jacobi_connection_coeffs(N, alph, bet, sgn)
%
% If sgn is +1, returns an N x 2 matrix C such that 
%
% (1 - x) p_n^{(alph+1,bet)}(x) =   C(n+1, 1) p_n^{(alph,bet)}(x) +
%                                   C(n+1, 2) p_{n+1}^{(alph,bet)}(x)
%
% If sgn is -1, returns an N x 2 matrix C such that 
%
% (1 + x) p_n^{(alph,bet+1)}(x) =   C(n+1, 1) p_n^{(alph,bet)}(x) +
%                                   C(n+1, 2) p_{n+1}^{(alph,bet)}(x)
%
% The input sgn is either +1 or -1 (the input is rounded accordingly). The
% parameters (alph, bet) are larger than -1. N is a positive integer.

C = zeros([N 2]);
sgn = sign(sgn);
z = sgn;

[a,b] = jacobi_recurrence(N+1, alph, bet);
r = ratio_eval(a, b, z, N);
U = sqrt(sqrt(b(2:N+1)));

C(:,2) = -sgn*U./sqrt(sgn*r);
C(:,1) = U.*sqrt(sgn*r);
