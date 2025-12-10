function [U, S, V] = fast_rand_svd(A,r,e)
%FAST_RAND_SVD Fast Randomized computation of SVD.
% Given an m × n matrix A, a tolerance ?, and an integer r
%(e.g. r = 10), the following scheme computes an orthonormal
% matrix Q such that (4.2) holds with probability at least 1 ? min{m, n}10?r.

% r is probabilty 1-10^-r
% e is our error tolerance variable
% A is matrix we want to find svd of


[M,N] = size(A);
%Draw standard Gaussian vectors ?(1), . . . , ?(r) of length n.
W = randn(N, r);

%For i = 1,2,...,r, compute y(i) = A?(i).
Y = A*W;
%initialize j for loop
% j = 0;

%Q(0) =[ ], the m×0 empty matrix.
maxit = N;
Q = zeros(M,maxit);

% mx = vecnorm(Y);
% 
% %error threshold for loop
th = e/(10*sqrt(2/pi));
W = randn(N,maxit);

I = speye(M,M);
for j=1:maxit
    
    
    %Overwrite y(j) by I - Q(j-1)(Q(j-1))* y(j) 
    Y(:,j) = (I - Q*Q')*Y(:,j);
    
    %normalize and put into little q temporarily
    q = Y(:,j)/norm(Y(:,j));
    
    %add to big Q
    Q(:,j) = q;
    
    %make another random vector and compute something similar to above
    w = W(:,j);
    Y(:, j+r) = (I - Q*Q.')*(A*w);        
    
    Y(:,j+1:j+r-1) = (I - q*q')*Y(:,j+1:j+r-1);
    
 
    %recalculate max
    mx = max(vecnorm(Y(:,j+1:j+r)));
    if mx<=th
        break;
    end
end

%calculate SVD with truncated matrix
B = Q'*A;
[Ub, S, V] = svd(B,'econ');
S = diag(S);
U = Q*Ub;


end

