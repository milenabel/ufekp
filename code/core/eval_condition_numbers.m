function [lam1_max, lam1_med, lam2_max, lam2_med] = ...
         eval_condition_numbers(solveTt, Kx_fun, Px_fun, Xeval)
% Evaluation "condition number": norms of the cardinal vector L(x*).
% solveTt : function handle that solves A' z = rhs for many rhs
% Kx_fun : function handle @(xstar) k(xstar,X) (n×1)
% Px_fun : function handle @(xstar) p(xstar) (m×1)
% Xeval : q×d evaluation points to sample
%
% Returns max/median of ||L||_1 and ||L||_2 over Xeval.

    q = size(Xeval,1);
    lam1 = zeros(q,1); lam2 = zeros(q,1);
    for i=1:q
        xstar = Xeval(i,:);
        rhs = [Kx_fun(xstar); Px_fun(xstar).'];
        z = solveTt(rhs); % z = [L; μ]
        n = numel(Kx_fun(xstar));
        L = z(1:n);
        lam1(i) = norm(L,1);
        lam2(i) = norm(L,2);
    end
    lam1_max = max(lam1); lam1_med = median(lam1);
    lam2_max = max(lam2); lam2_med = median(lam2);
end
