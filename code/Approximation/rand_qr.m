function [Qh,R,P] = rand_qr(A,tol,b)
%RAND_QR Martinsson's randomized blocked-QR
%% Uses a QB decomposition to reconstruct the QR decomposition
n = size(A,2);
Q = []; B = [];
for i=1:n
    Om = randn(n,b);
    Qi = qr(A*Om,0);
    Q = [Q,Qi];
    Bi = Qi.'*A;
    A = A - Qi*Bi;
    B = [B,Bi];
    if norm(A,'fro') < tol
        break;
    end
end
[Qtil,R,P] = qr(B,0);
Qh = Q*Qtil;

end

