%% Solve Poisson's equation on an arbitrary 2D domain using CSRBFs

%% Type of boundary condition?
bctype=1; %1-Dirichlet, 2-Neumann, 3- Robin

%% weighted fekete nodes or not?
weighted_nodes = 1;

%% Setup anonymous functions for PDE solution and rhs for interior
%and boundary. Assume 2D.
syms x y;
syms nrx nry;
%u = sin(pi*x).*cos(pi*y); 
u = sin(pi*x).*cos(pi*y).*cos(pi*x);
f = (diff(u,x,2) + diff(u,y,2)); %lap u = f

if bctype==1
    g = u;
elseif bctype==2
    g = nrx.*diff(u,x) + nry.*diff(u,y); %n.grad(u) = g;
elseif bctype==3
    g = nrx.*diff(u,x) + nry.*diff(u,y) + u; %n.grad(u) = g;
end
u = matlabFunction(u);
f = matlabFunction(f);
g = matlabFunction(g);
clear x y nrx nry;

%% Load up the node set
st = load('DiskPoissonNodesLarge.mat');
%st = load('DomainPoissonNodes2.mat');
xi = st.fullintnodes{1};
xb = st.bdrynodes{1};
nr = st.normals{1};
x_candidates = [xi;xb];

%% Get approximate weighted fekete points
if weighted_nodes==1
    dim = 2;
    alph = 0; 
    recurrence = @(N) jacobi_recurrence(N,alph,alph);
    ell = 15;
    a = total_degree_indices(dim,ell);
    N = size(a,1);
    V = mpoly_eval(x_candidates, a, recurrence);
    weights = 1./sqrt(sum(V.^2, 2));
    V = repmat(weights, [1 N]).*V;
    [~,~,e] = qr(V', 0);
    ei = e((e(1:N)<=length(xi)));
    eb = e((e(1:N)>length(xi)));
    V = V(e(1:N),:);
    nr_ind = eb - length(xi);
    xi = x_candidates(ei,:);
    xb = x_candidates(eb,:);
    nr = nr(nr_ind,:);
    x = [xi;xb];
    weights = [weights(ei);weights(eb)]; %Not sure what to do with these!
else
    x = x_candidates;
end
%% Setup RBF approximation with compactly supported RBFs. We'll use Wendland's C^6 RBF.
%%Here we use the shifted form.
Ni = length(xi); Nb = length(xb);
tree = KDTreeSearcher(x);
rbf = @(e,r) r.^8.*(66*spones(r) - 154*r + 121*r.^2 - 32*r.^3);
drbfor = @(e,r) -22.*e.^2.*r.^7.*(39*e.*spones(r) + 16*e.^2.*r.^2 + 24);
lrbf = @(e,r) 44*e.^2.*r.^6.*(84*spones(r)-264*r+267*r.^2-88*r.^3);
ep = 1;
support = 1/ep;
rd = DistanceMatrixCSRBF(x,x,ep);
Lap = lrbf(ep,rd);

if bctype==1 %Dirichlet
    D = rbf(ep,rd(Ni+1:end,:));
elseif bctype==2 %Neumann
    
    %% Get gradients    
    Drbf = drbfor(ep,rd(Ni+1:end,:));
    if size(x,2)>=1
       [xj,xk] = ndgrid(xb(:,1),x(:,1));
       nr1 = repmat(nr(:,1),[1,size(xj,2)]);
       D = nr1.*(xj-xk).*Drbf;
    end
    if size(x,2)>=2
       [yj,yk] = ndgrid(xb(:,2),x(:,2));
       nr2 = repmat(nr(:,2),[1,size(xj,2)]);
       D = D + nr2.*(yj-yk).*Drbf;
    end
    if size(x,3)==3
       [zj,zk] = ndgrid(xb(:,3),x(:,3));
       nr3 = repmat(nr(:,3),[1,size(xj,2)]);
       D = D + nr(:,3).*(zj-zk).*Drbf; 
    end
    
elseif bctype==3 %Robin
    Drbf = drbfor(ep,rd(Ni+1:end,:));
    A = rbf(ep,rd(Ni+1:end,:));
    
    %% Get gradients  
    if size(x,2)>=1
       [xj,xk] = ndgrid(xb(:,1),x(:,1));
       nr1 = repmat(nr(:,1),[1,size(xj,2)]);
       D = nr1.*(xj-xk).*Drbf;
    end
    if size(x,2)>=2
       [yj,yk] = ndgrid(xb(:,2),x(:,2));
       nr2 = repmat(nr(:,2),[1,size(xj,2)]);
       D = D + nr2.*(yj-yk).*Drbf;
    end
    if size(x,3)==3
       [zj,zk] = ndgrid(xb(:,3),x(:,3));
       nr3 = repmat(nr(:,3),[1,size(xj,2)]);
       D = D + nr(:,3).*(zj-zk).*Drbf; 
    end  
    D = D + A;
end

%% Adjust Lap for BCs
Lii = Lap(1:Ni,1:Ni);
Lib = Lap(1:Ni,Ni+1:end);
lhs = [[Lii,Lib];D];


%% Get rhs for linear system
rhsi = f(xi(:,1),xi(:,2));

if bctype==1
    rhsb = g(xb(:,1),xb(:,2));
elseif bctype>=2
    rhsb = g(nr(:,1),nr(:,2),xb(:,1),xb(:,2));
end
%% Solve
c = lhs\[rhsi;rhsb];
sol = rbf(ep,rd)*c; sol = sol(1:Ni);
if bctype==2
    sol = sol - mean(sol);
end
true_sol = u(xi(:,1),xi(:,2));

%% plot stuff
plot(sol,'.')
hold on
plot(true_sol,'o');