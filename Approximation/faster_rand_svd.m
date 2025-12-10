function [U, S, V] = faster_rand_svd(A,k)
%FAST_RAND_SVD Fast Randomized computation of SVD.
[~,n] = size(A);
l = k+2;
R = 2*rand(n,l)-1 + ~isreal(A)*1i*(2*rand(n,l)-1);
%   Apply A to a random matrix, obtaining Q.
Q = A*R;
[Q,~] = lu(Q);
its = n;
for it = 1:its
    Q = A'*Q;
    [Q,~] = lu(Q);
    Q = A*Q;
    if(it < its)
        [Q,~] = lu(Q);
    else
        [Q,~,~] = qr(Q,0);
    end

end
%   SVD Q' applied to the centered A to obtain approximations
%   to the singular values and right singular vectors of the A;
R = Q';
[R,S,V] = svd((A'*R')','econ');
U = Q*R;
U = U(:,1:k);
V = V(:,1:k);
S = S(1:k,1:k);
S = diag(S);