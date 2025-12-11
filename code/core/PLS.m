function [el2,elinf,a_time,e_time,c_poly] = PLS(x,y,ell,xe,alph,ye_true)
%PLS Polynomial Least Squares approximation

dim = size(x,2);
recurrence = @(N) jacobi_recurrence(N,alph,alph); 
a = total_degree_indices(dim,ell);    

xmax = max(x, [], 1);
xmin = min(x, [], 1);
c1 = 0.5 .* (xmax + xmin);
c2 = 2 ./ (xmax - xmin);

xc = (x - repmat(c1, [length(x) 1])) .* repmat(c2, [length(x) 1]); 
V = mpoly_eval(xc,a,recurrence);    

tic
[Q,R] = qr(V,0); %use QR for least squares
opts.UT = true;
c_poly = linsolve(R,Q'*y,opts);       
a_time = toc;


tic
xce = (xe - repmat(c1, [length(xe) 1])) .* repmat(c2, [length(xe) 1]);
Ve = mpoly_eval(xce,a,recurrence);
ye = Ve*c_poly;
e_time = toc;


el2 = norm(ye - ye_true)./norm(ye_true);
elinf = norm(ye - ye_true,inf)./norm(ye_true,inf);

end

