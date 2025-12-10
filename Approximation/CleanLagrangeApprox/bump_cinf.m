function Z = bump_cinf(x, y, xc, yc, R,p)
 r2 = (x-xc).^2 + (y-yc).^2;
  Z  = zeros(size(r2));
  mask = (r2 < R^2);
  rr   = sqrt(r2(mask));
  psi  = exp(-p * R^2 ./ (R^2 - rr.^2));
  Z(mask) = psi / exp(-p);
end
