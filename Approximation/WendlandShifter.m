%% For a given Wendland RBF, returns its shifted form for derivatives
%%and Laplacian.

syms e r;

%% W(3,2) = Wendland(3,1)
r1 = (1-e*r).^4;
r2 = (4*e*r + 1);

% %% W(3,4) = Wendland(3,2)
% r1 = (1-e*r).^6;
% r2 = (35*(e*r).^2 + 18*(e*r) + 3)/3;

% %% W(3,6)= Wendland(3,3)
% r1 = (1-e*r).^8; %power function part
% r2 = (32*(e*r).^3 + 25*(e*r).^2 + 8*(e*r) + 1); %polynomial part

rbf = r1.*r2;
%% First derivative wrt r
drbfor = simplify(diff(r1).*r2 + r1.*diff(r2))./r;

%% Laplacian
dim = 2;
lrbf = simplify(diff(rbf,r,2) + (dim-1).*drbfor);

%% Shifted forms
mod_rbf = r1.*subs(r2,e*r,r);
syms u;
shifted_rbf = simplify(subs(subs(mod_rbf,e*r-1,u),e*r,1-u));
shifted_drbfor = simplify(subs(subs(drbfor,e*r-1,u),e*r,1-u));
shifted_lrbf = simplify(subs(subs(lrbf,e*r-1,u),e*r,1-u),10);

%% Now replace u with r and vectorize
shifted_rbf = vectorize(subs(shifted_rbf,u,r))
shifted_drbfor = vectorize(subs(shifted_drbfor,u,r))
shifted_lrbf = vectorize(subs(shifted_lrbf,u,r))