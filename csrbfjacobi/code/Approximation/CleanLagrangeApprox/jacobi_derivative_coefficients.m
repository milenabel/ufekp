function[c] = jacobi_derivative_coefficients(alph, bet, N)
% c = jacobi_derivative_coefficients(alph, bet, N)
%
% Returns the first N coefficients c_0, ..., c_{N-1} for the identity
%
% d/dx p^{(alph,bet)}_j(x) = c_j p^{(alph+1,bet+1)}_{j-1}(x).
%
% Here, p^{(alph,bet)} are the Jacobi polynomials with parameters (alph,bet)
% that are orthonormal on the interval [-1,1] with respect to the
% *unnormalized* density
%
% w(x) = (1-x)^alph * (1+x)^bet

n = (0:(N-1)).';

%c = [0; (n + alph + bet + 1) .* sqrt( n./(n + alph + bet + 1))];
c = sqrt(n.*(n+alph+bet+1));
