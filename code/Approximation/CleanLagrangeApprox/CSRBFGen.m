function [el2,elinf,a_time,e_time,c_poly,cond_num,c_rbf, sparsity] = CSRBFGen(x,y,ell,xe,alph,rbf,ep,tree,ye_true)
% function [el2,elinf,a_time,e_time,c_poly,c_rbf, sparsity] = CSRBFGen(x,y,ell,xe,alph,rbf,ep,tree,ye_true)
%PLS General Interpolation with CSRBFs and polys

dim = size(x,2);
recurrence = @(N) jacobi_recurrence(N,alph,alph); 
a = total_degree_indices(dim,ell);    

xmax = max(x, [], 1);
xmin = min(x, [], 1);
c1 = 0.5 .* (xmax + xmin);
c2 = 2 ./ (xmax - xmin);

xc = (x - repmat(c1, [length(x) 1])) .* repmat(c2, [length(x) 1]); 
V = mpoly_eval(xc,a,recurrence);    


opts2.UT = true;
tic
rd = DistanceMatrixCSRBFwt(x,x,ep,tree);
A = rbf(ep,rd);
L_a = chol(A,'lower'); %get lower triangular cholesky factor            
dA = decomposition(full(L_a),'triangular','lower');
B = dA\V;
[Q,R] = qr(B,0);
g = dA\y;
c_poly = linsolve(R,Q'*g,opts2);
%c_rbf = dA'\(g - dA\(V*c_poly));         
c_rbf = (L_a.')\(g- B*c_poly);
a_time = toc;
density   = nnz(A)/numel(A);
sparsity  = 1 - density;  % degree of sparsity of the A matrix
cond_num = condest(L_a);  

tic
rde = DistanceMatrixCSRBFwt(xe,x,ep,tree);
xce = (xe - repmat(c1, [length(xe) 1])) .* repmat(c2, [length(xe) 1]);
Ae = rbf(ep,rde);
Ve = mpoly_eval(xce,a,recurrence);
ye = Ae*c_rbf + Ve*c_poly;
e_time = toc;


el2 = norm(ye - ye_true)./norm(ye_true);
elinf = norm(ye - ye_true,inf)./norm(ye_true,inf);

end

