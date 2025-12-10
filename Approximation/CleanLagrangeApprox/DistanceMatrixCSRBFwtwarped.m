function DM = DistanceMatrixCSRBFwtwarped(dsites, ctrs, ep, Tree, warp)
%DISTANCEMATRIXCSRBFWTWARPED  Wraps DistanceMatrixCSRBFwt + metric‐warp
%
%   DM = DistanceMatrixCSRBFwtwarped(dsites, ctrs, ep, Tree, warp)
%   calls the original DistanceMatrixCSRBFwt to get the sparse “u–distance”  
%       u0_ij = 1 – ep * ||dsites(i,:) – ctrs(j,:)||  
%   then rescales each entry by your C^∞ bump weights:
%       u_ij = 1 – (1 – u0_ij) * sqrt(w_ds(i) * w_cts(j)).
%
%   Inputs:
%     dsites  – N×d data (or eval) sites  
%     ctrs    – M×d RBF centers  
%     ep      – shape parameter  
%     Tree    – KDTreeSearcher built on whichever set DistanceMatrixCSRBFwt expects  
%     warp    – struct with fields x0, Rf, wmin (and warp.w if you like)  
%
%   Output:
%     DM      – N×M sparse, ready for A = rbf(ep,DM)

  % 1) get the original sparse matrix of u–values
  DM0 = DistanceMatrixCSRBFwt(dsites, ctrs, ep, Tree);

  % 2) extract triplets
  [i, j, u0] = find(DM0);

  % 3) compute per‐site bump weights
  w_ds = local_bump_weights(dsites, warp);   % N×1
  w_ct = local_bump_weights(ctrs,   warp);   % M×1

  % 4) build outer‐product sqrt(w_i w_j)
  Wij    = sqrt( w_ds(i) .* w_ct(j) );

  % 5) warp the u‐distance
  u_warp = 1 - (1 - u0) .* Wij;

  % 6) reassemble sparse matrix
  DM = sparse(i, j, u_warp, size(DM0,1), size(DM0,2));
end

function w = local_bump_weights(P, warp)
%LOCAL_BUMP_WEIGHTS  C^∞ bump weights for each point in P (N×d)
  dif  = P - warp.x0.';             % N×d
  r2   = sum(dif.^2, 2);
  mask = (r2 < warp.Rf^2);
  w    = ones(size(r2));
  if any(mask)
    rr   = sqrt(r2(mask));
    psi0 = exp(-1);
    psi  = exp(-warp.Rf^2 ./ (warp.Rf^2 - rr.^2));
    Q    = psi ./ psi0;             % in (0,1]
    w(mask) = 1 - (1 - warp.wmin) * Q;
  end
end
