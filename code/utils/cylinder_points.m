function X = cylinder_points(N, cluster)
% CYLINDER_NODES  quasi‐uniform nodes on finite cylinder x^2+y^2=1, z∈[0,1]
%   X = cylinder_nodes(N)              % default: uniform in z
%   X = cylinder_nodes(N, true)        % mild clustering at z=0 and z=1
%
%   X is approximately N×3 Cartesian points on the lateral surface.

% cylinder parameters
R = 1;
H = 1;

% approximate spacing h so that h^2 ≈ surface area/N
A  = 2*pi*R*H;
h  = sqrt(A/N);

% numbers in azimuth and axial
Nu = max(4, round(2*pi*R / h));
Nz = max(2, round(H      / h));

% azimuthal samples
u = linspace(0, 2*pi, Nu+1);
u(end) = [];                   % drop duplicate at 2π

% axial samples with optional boundary clustering
s = linspace(0,1,Nz);
if nargin>1 && cluster
  % Chebyshev‐type spacing clusters at both ends
  z = (1 - cos(pi*s))/2;
else
  z = s;
end

% build grid and map to Cartesian
[UU, ZZ] = meshgrid(u, z);
X = [ cos(UU(:)), sin(UU(:)), ZZ(:) ];
end
