clear 
close all

dim = 3;
N = 40;

x = cos(pi*rand([ N dim ]));

loi_info = loi_decomposition(x, 1e-3);

indices = zeros([0 dim]);
for q = 1:length(loi_info.ranks)
  indices = [indices; loi_levels(dim, q-1)];
end

% Make sure permutation info is correct
u = (x - repmat(loi_info.b, [N 1])) .* repmat(loi_info.a, [N 1]);
u = u(loi_info.row_permutation,:);
norm( u - loi_info.u )

% Make sure Vandermonde-like matrix is fine
V = chebyshev_eval( loi_info.u, indices)*loi_info.H.';
%V = V(loi_info.row_permutation,:);

V2 = loi_info.L * loi_info.U;

V3 = loi_vandermonde(loi_info, x);
V3 = V3(loi_info.row_permutation,:);

norm( V2 - V )
norm( V3 - V )
