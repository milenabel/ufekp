function X = half_torus_points(N, q)
% HALF_TORUS_NODES  quasi‐uniform nodes on the “cut” torus (u∈[0,π],v∈[0,2π))
%   X = half_torus_nodes(N)       % default q=2 clusters moderately toward u=0,π
%   X = half_torus_nodes(N, q)    % q>1 increases clustering at the cut‐edges
%
%   N ≈ total desired nodes; output X is #nodes×3 in Cartesian coords.

if nargin<2, q = 2; end

R = 1;   % major radius
r = 1/3; % minor radius

% choose grid dimensions so Nu*Nv ≈ N
Nu = round(sqrt(N));
Nv = ceil(N/Nu);

% parameter lines
s = linspace(0,1,Nu);             % uniform in [0,1]
u = pi* ( sin(pi*s/2).^q );       % sin^q–map clusters at u=0 and u=π
v = 2*pi*(0:(Nv-1))/Nv;            % uniform in [0,2π)

[U,V] = meshgrid(u,v);

X = zeros(Nu*Nv,3);
X(:,1) = (R + r*cos(V(:))) .* cos(U(:));
X(:,2) = (R + r*cos(V(:))) .* sin(U(:));
X(:,3) = r *          sin(V(:));
end
