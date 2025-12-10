function [x] = fast_rand_svd_solve(A,b,r,e)
%FAST_RAND_SVD Fast Randomized computation of SVD + linear solve.
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
j = 0;

%Q(0) =[ ], the m×0 empty matrix.
Q = zeros(M,0);

mx = max(vecnorm(Y));

%error threshold for loop
th = e/(10*sqrt(2/pi));

while mx > th
    j = j+1;
    
    %Overwrite y(j) by I - Q(j-1)(Q(j-1))* y(j) 
    Y(:,j) = Y(:,j) - Q*(Q'*Y(:,j));
    
    %normalize and put into little q temporarily
    q = Y(:,j)/norm(Y(:,j));
    
    %add to big Q
    Q = [Q, q];
    
    %make another random vector and compute something similar to above
    w = randn(N, 1);
    Y(:, j+r) = A*w - Q*(Q.'*(A*w));
        
    for i = j+1:j+r-1
        Y(:,i) = Y(:,i) - q*(q'*Y(:,i));
    end
 
    %recalculate max
    mx = max(vecnorm(Y(:,j+1:j+r)));
end

%calculate SVD with truncated matrix
[Ub, S, V] = svd(Q'*A,'econ');
x = V*(1./diag(S).*(Ub.'*(Q.'*b)));


end

