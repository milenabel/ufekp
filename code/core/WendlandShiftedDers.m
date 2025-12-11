function [rbf,drbfor,lrbf,l2rbf] = WendlandShiftedDers(d,k,s)
%WENDLANDSHIFTEDDERS Generates shifted form of Wendland functions W(d,k,s)
%and derivative over r, and laplacian in s dims
rbf = Wendland(d,k);
syms e r;
rbf = subs(rbf,r,e*r);
drbfor = simplify(diff(rbf,r)./r);
lrbf = simplify(diff(rbf,r,2) + (s-1).*drbfor);
l2rbf = simplify(diff(lrbf,r,2) + (s-1).*(diff(lrbf,r)./r));
syms u;
rbf = simplify(subs(subs(rbf,e*r-1,u),e*r,1-u));
rbf = (subs(rbf,u,r));
drbfor = simplify(subs(subs(drbfor,e*r-1,u),e*r,1-u));
drbfor = (subs(drbfor,u,r));
lrbf = simplify(subs(subs(lrbf,e*r-1,u),e*r,1-u));
lrbf = (subs(lrbf,u,r));
l2rbf = simplify(subs(subs(l2rbf,e*r-1,u),e*r,1-u));
l2rbf = (subs(l2rbf,u,r));
end

