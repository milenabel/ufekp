function YW = apply_warp(Y, warp)
%APPLY_WARP  Apply local C^\infty warp to new M×D points
%   YW = apply_warp(Y, warp) warps Y according to parameters in 'warp'.

  [M, D] = size(Y);
  x0   = warp.x0;
  Rf   = warp.Rf;
  wmin = warp.wmin;

  diff = Y - x0';
  r    = sqrt(sum(diff.^2, 2));
  w    = ones(M,1);

  mask = r < Rf;
  psi0 = exp(-1);
  rr   = r(mask);
  psi  = exp(-Rf^2 ./ (Rf^2 - rr.^2));
  Q    = psi / psi0;
  w(mask) = 1 - (1 - wmin) * Q;

    % apply only local warp, keep other points unchanged
  YW = Y;
  YW(mask, :) = x0' + diff(mask, :) .* w(mask);
end