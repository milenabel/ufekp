function X = hemispherePts(N, p)
% HEMISPHERE_NODES  quasi‐uniform nodes on unit upper hemisphere
%   X = hemisphere_nodes(N)          % uniform area (p=1)
%   X = hemisphere_nodes(N, p<1)     % clusters toward boundary (z=0)
%
% Output X is N×3 with rows [x y z].

if nargin<2, p = 1; end
phi_const = (sqrt(5)-1)/2;       % golden ratio offset
X = zeros(N,3);

for k = 0:N-1
  t   = k/(N-1);
  z   = t^p;                      % clustering in z
  phi = 2*pi*phi_const*k;
  r   = sqrt(max(0,1 - z^2));
  X(k+1,:) = [r*cos(phi), r*sin(phi), z];
end
end
