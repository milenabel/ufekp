function kappa = rbf_effective_cond(K, P)
% Effective condition number of RBF block restricted to null(P^T).
% K: n×n SPD kernel Gram; P: n×m polynomial constraints.
    if isempty(P), kappa = cond(K); return; end
    [Q,~] = qr(P,0); % Q = [Q1 Q2], range(Q1)=range(P)
    r = rank(P);
    Q2 = Q(:, r+1:end); % orthonormal basis of null(P^T)
    Keff = Q2'*(K*Q2); % SPD on constrained subspace
    Keff = (Keff+Keff')/2;
    kappa = cond(Keff);
end
