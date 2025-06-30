%% For a given Wendland RBF, returns its shifted form for derivatives
%%and Laplacian.

syms e r;
r1 = (1-e*r).^8; %power function part
r2 = (32*(e*r).^3 + 25*(e*r).^2 + 8*(e*r) + 1); %polynomial part
rbf = r1.*r2;

%% First derivative wrt r
drbfor = simplify(diff(r1).*r2 + r1.*diff(r2))./r;

%% Laplacian
dim = 3;
lrbf = simplify(diff(rbf,r,2) + (dim-1).*drbfor);

%% Shifted forms
mod_rbf = r1.*subs(r2,e*r,r)
syms u;
shifted_drbfor = simplify(subs(subs(drbfor,e*r-1,u),e*r,1-u));
shifted_lrbf = simplify(subs(subs(lrbf,e*r-1,u),e*r,1-u),10)

%% Now replace u with r and vectorize
shifted_drbfor = vectorize(subs(shifted_drbfor,u,r))
shifted_lrbf = vectorize(subs(shifted_lrbf,u,r))