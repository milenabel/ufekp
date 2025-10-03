function [rbf,drbfor,lrbf] = WendlandShiftedDers(d,k)
%WENDLANDSHIFTEDDERS Generates shifted form of Wendland functions W(d,k)
%and derivative over r, and laplacian in d dims
rbf = Wendland(d,k);
syms e r;
rbf = subs(rbf,r,e*r);
drbfor = simplify(diff(rbf,r)./r);
lrbf = simplify(diff(rbf,r,2) + (d-1).*drbfor);
syms u;
rbf = simplify(subs(subs(rbf,e*r-1,u),e*r,1-u));
rbf = vectorize(subs(rbf,u,r));
drbfor = simplify(subs(subs(drbfor,e*r-1,u),e*r,1-u));
drbfor = vectorize(subs(drbfor,u,r));
lrbf = simplify(subs(subs(lrbf,e*r-1,u),e*r,1-u));
lrbf = vectorize(subs(lrbf,u,r));

end

