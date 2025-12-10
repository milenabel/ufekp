function [x,nr,t,p] = computeHexagonTorusNodes(N)

% Torus parameters.
R = 1;  
r = 1/3;

p = linspace(0,2*pi,N+1); p = p(1:N);
h = p(2)-p(1);
p2 = p(1:end)+h/2;
t = linspace(0,2*pi,3*N+1); t = t(1:3*N);
h = t(2)-t(1);
t2 = t(1:end)+h/2;
[t,p] = meshgrid(t,p); t = t(:); p = p(:);
[t2,p2] = meshgrid(t2,p2); t2 = t2(:); p2 = p2(:);
t = [t;t2]; p = [p;p2];
x = [(R + r*cos(p)).*cos(t) (R + r*cos(p)).*sin(t) r*sin(p)];

nr(:,1) = r.*cos(p).*cos(t).*(R + r.*cos(p));
nr(:,2) = r.*cos(p).*sin(t).*(R + r.*cos(p));
nr(:,3) = r.*sin(p).*(R + r.*cos(p)).*cos(t).^2 + r.*sin(p).*(R + r.*cos(p)).*sin(t).^2;
nr=nr./repmat(sqrt(sum(nr.^2,2)),[1 3]);
end
