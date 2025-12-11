function [el2,elinf,a_time,e_time,c_poly,cond_num,c_rbf,sparsity] = CSRBFGenWarp(x,y,ell,xe,alph,rbf,ep,ye_true)
 %--- 1) compute warp parameters ---
  [~, warp] = warpsites(x, y);

  %--- 2) warp centers and rebuild tree ---
  Z    = apply_warp(x, warp);           % N×d warped centers
  treeZ= KDTreeSearcher(Z);

  %--- 3) build kernel matrix on warped centers ---
  DM0  = DistanceMatrixCSRBFwt(Z, Z, ep, treeZ);
  A    = rbf(ep, DM0);

  %--- 4) build polynomial block on original sites ---
  N   = size(x,1); d = size(x,2);
  recurrence = @(n) jacobi_recurrence(n,alph,alph);
  a   = total_degree_indices(d, ell);
  xmin= min(x, [],1); xmax= max(x,[],1);
  c1  = 0.5*(xmax + xmin); c2 = 2./(xmax - xmin);
  xc  = (x - c1) .* c2;
  V   = mpoly_eval(xc, a, recurrence);

  %--- 5) solve interpolation system ---
  cond_num = condest(A);
  tic;
  L      = chol(A,'lower');
  dA     = decomposition(full(L),'triangular','lower');
  B      = dA \ V;
  g      = dA \ y;
  opts.UT= true;
  [Q,R]  = qr(B,0);
  c_poly = linsolve(R, Q' * g, opts);
  c_rbf  = (L.') \ (g - B * c_poly);
  a_time = toc;

  % sparsity
  sparsity = 1 - nnz(A)/numel(A);

  %--- 6) evaluate at warped evaluation points ---
  Xe = apply_warp(xe, warp);
  DM0_e = DistanceMatrixCSRBFwt(Xe, Z, ep, treeZ);
  Ae     = rbf(ep, DM0_e);
  xce    = (xe - c1) .* c2;
  Ve_e   = mpoly_eval(xce, a, recurrence);
  tic;
  ye     = Ae * c_rbf + Ve_e * c_poly;
  e_time = toc;

  % errors
  el2   = norm(ye - ye_true) / norm(ye_true);
  elinf = norm(ye - ye_true, Inf) / norm(ye_true, Inf);
end
