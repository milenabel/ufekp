function [TX,warp] = warpsites(X,F)
%WARP_SITES  Warp site locations focusing on high-gradient regions
%   TX = warp_sites(X) returns warped N×d points TX from original X (N×d).
%

  N = size(X,1);

  dim = size(X,2);
  alph = 0;
  recurrence = @(N) jacobi_recurrence(N,alph,alph); 
  a = total_degree_indices(dim,3);   
  P = mpoly_eval(X,a,recurrence); 
  d = P\F;

  if dim==1
     gpx =  mpoly_eval(X,a,recurrence,1);
     gradNorm = abs(gpx*d);
  elseif dim==2
     gpx =  mpoly_eval(X,a,recurrence,[1,0]);
     gpy =  mpoly_eval(X,a,recurrence,[0,1]);
     gradNorm = sqrt((gpx*d).^2 + (gpy*d).^2);
  elseif dim==3
     gpx =  mpoly_eval(X,a,recurrence,[1,0,0]);
     gpy =  mpoly_eval(X,a,recurrence,[0,1,0]);
     gpz =  mpoly_eval(X,a,recurrence,[0,0,1]);
     gradNorm = sqrt((gpx*d).^2 + (gpy*d).^2 + (gpz*d).^2);
  end

  %---- 2) Estimate feature center ----
  wgt = gradNorm / (sum(gradNorm) + eps);
  x0  = (wgt' * X)';   % D×1 weighted centroid

  %---- 3) Compute distance and warp factors ----
  diff = X - x0';            % N×D
  r    = sqrt(sum(diff.^2, 2));
  Rf   = quantile(r, 0.8);  % focus radius (e.g. 25th percentile)
  wmin = 0.8;                % minimum shrink at center
  w    = ones(N,1);

  % C^inf bump: psi(r) = exp(-Rf^2/(Rf^2 - r^2)), for r<Rf
  mask = r < Rf;
  psi0 = exp(-1);
  rr   = r(mask);
  psi  = exp(-Rf^2 ./ (Rf^2 - rr.^2));
  Q    = psi / psi0;         % normalized bump in (0,1]

  % warp factor: shrink inside focus, identity outside
  w(mask) = 1 - (1 - wmin) * Q;

  %---- 4) Apply local warp ----
  % new position = x0 + (X-x0)*w; points with w=1 => unchanged
  TX = x0' + diff .* w;

  %---- 5) Return warp parameters ----
 warp.x0   = x0;
  warp.Rf   = Rf;
  warp.wmin = wmin;
  warp.w    = w;     % <-- per-site shrink weights in (0,1]