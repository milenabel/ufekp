%% Solve Poisson's equation on an arbitrary 3D domain using Jacobi polynomials.
%%For this to work, we use weighted Fekete points. We naively form and solve
%%Kansa-style unsymmetric collocation systems. Stepping stone to fast algorithms. 

%% Spatial dimension
dim = 2;

%% Type of boundary condition?
bctype=1; %1-Dirichlet, 2-Neumann, 3- Robin

%% Setup anonymous functions for PDE solution and rhs for interior
%and boundary. 
if dim==2
    syms x y;
    syms nrx nry;
    u = 1 + sin(pi*x).*cos(pi*y);     
    f = (diff(u,x,2) + diff(u,y,2)); %lap u = f
    if bctype==1
        g = u;
    elseif bctype==2
        g = nrx.*diff(u,x) + nry.*diff(u,y); %n.grad(u) = g;
    elseif bctype==3
        g = nrx.*diff(u,x) + nry.*diff(u,y) + u; %n.grad(u) = g;
    end    
elseif dim==3
    syms x y z;
    syms nrx nry nrz;
    u = 1 + sin(pi*x).*cos(pi*y).*sin(pi*z);     
    f = (diff(u,x,2) + diff(u,y,2) + diff(u,z,2)); %lap u = f    
    if bctype==1
        g = u;
    elseif bctype==2
        g = nrx.*diff(u,x) + nry.*diff(u,y) + nrz.*diff(u,z); %n.grad(u) = g;
    elseif bctype==3
        g = nrx.*diff(u,x) + nry.*diff(u,y) + nrz.*diff(u,z) + u; %n.grad(u) = g;
    end      
end
u = matlabFunction(u);
f = matlabFunction(f);
g = matlabFunction(g);
if dim==2
    clear x y nrx nry;
elseif dim==3
    clear x y z nrx nry nrz;
end

%% Load up the node set
if dim==2
    st = load('DiskPoissonNodesLarge.mat');
    %st = load('DomainPoissonNodes2.mat');
elseif dim==3
    st = load('SpherePoissonNodes.mat');
end
if dim==2
    xi = st.fullintnodes{6};
    xb = st.bdrynodes{6};    
    nr = st.normals{6};    
elseif dim==3
    xi = st.fullintnodes{5};
    xb = st.bdrynodes{5};    
    nr = st.normals{5};
end
x_candidates = [xi;xb];

%% Get approximate weighted fekete points
alph = 0; 
recurrence = @(N) jacobi_recurrence(N,alph,alph);
ell = 60;
a = total_degree_indices(dim,ell);
%a = hyperbolic_cross_indices(dim,ell);
N = size(a,1);   
V = mpoly_eval(x_candidates, a, recurrence);
weights = 1./sqrt(sum(V.^2, 2));
V = repmat(weights, [1 N]).*V;
[~,~,e] = qr(V', 0);
ei = e((e(1:N)<=length(xi)));
eb = e((e(1:N)>length(xi)));
nr_ind = eb - length(xi);
xi = x_candidates(ei,:);
xb = x_candidates(eb,:);
nr = nr(nr_ind,:);
x = [xi;xb];
weights = [weights(ei);weights(eb)]; %Not sure what to do with these!


%% Setup differentiation matrices
V = mpoly_eval(x,a,recurrence);
Ni = length(xi); Nb = length(xb);

if dim==2
    Lv = mpoly_eval(x,a,recurrence,[2,0]) + mpoly_eval(x,a,recurrence,[0,2]);
elseif dim==3
    Lv = mpoly_eval(x,a,recurrence,[2,0,0])...
         + mpoly_eval(x,a,recurrence,[0,2,0])...
         + mpoly_eval(x,a,recurrence,[0,0,2]);    
end

if bctype==1 %Dirichlet    
    Dv = V(Ni+1:end,:); %get poly
elseif bctype==2 %Neumann
    
    %% Get polynomial gradients
    if dim==1
        Dvx = mpoly_eval(xb,a,recurrence,1);
        nr1 = repmat(nr(:,1),[1,size(Dvx,2)]);
        Dv = nr1.*Dvx;
    elseif dim==2
        Dvx =  mpoly_eval(xb,a,recurrence,[1,0]);
        Dvy =  mpoly_eval(xb,a,recurrence,[0,1]);
        nr1 = repmat(nr(:,1),[1,size(Dvx,2)]);
        nr2 = repmat(nr(:,2),[1,size(Dvy,2)]);
        Dv = nr1.*Dvx + nr2.*Dvy;
    elseif dim==3
        Dvx =  mpoly_eval(xb,a,recurrence,[1,0,0]);
        Dvy =  mpoly_eval(xb,a,recurrence,[0,1,0]);
        Dvz =  mpoly_eval(xb,a,recurrence,[0,0,1]);
        nr1 = repmat(nr(:,1),[1,size(Dvx,2)]);
        nr2 = repmat(nr(:,2),[1,size(Dvy,2)]);
        nr3 = repmat(nr(:,3),[1,size(Dvz,2)]);
        Dv = nr1.*Dvx + nr2.*Dvy + nr3.*Dvz;
    end
    
    
elseif bctype==3 %Robin   
    %% Get polynomial gradients
    if dim==1
        Dvx = mpoly_eval(xb,a,recurrence,1);
        nr1 = repmat(nr(:,1),[1,size(Dvx,2)]);
        Dv = nr1.*Dvx;
    elseif dim==2
        Dvx =  mpoly_eval(xb,a,recurrence,[1,0]);
        Dvy =  mpoly_eval(xb,a,recurrence,[0,1]);
        nr1 = repmat(nr(:,1),[1,size(Dvx,2)]);
        nr2 = repmat(nr(:,2),[1,size(Dvy,2)]);
        Dv = nr1.*Dvx + nr2.*Dvy;
    elseif dim==3
        Dvx =  mpoly_eval(xb,a,recurrence,[1,0,0]);
        Dvy =  mpoly_eval(xb,a,recurrence,[0,1,0]);
        Dvz =  mpoly_eval(xb,a,recurrence,[0,0,1]);
        nr1 = repmat(nr(:,1),[1,size(Dvx,2)]);
        nr2 = repmat(nr(:,2),[1,size(Dvy,2)]);
        nr3 = repmat(nr(:,3),[1,size(Dvz,2)]);
        Dv = nr1.*Dvx + nr2.*Dvy + nr3.*Dvz;
    end    
    
    
    %% form robin matrix for polys
    Dv = Dv + V(Ni+1:end,:);
end

%% Adjust Lap for BCs
lhs = [Lv(1:Ni,:);Dv];
z = spnull(lhs);
lhs = [[lhs,z];[z.',zeros(size(z,2))]];

%% Get rhs for linear system
if dim==2
    rhsi = f(xi(:,1),xi(:,2));
    if bctype==1
        rhsb = g(xb(:,1),xb(:,2));
    elseif bctype>=2
        rhsb = g(nr(:,1),nr(:,2),xb(:,1),xb(:,2));
    end
elseif dim==3
    rhsi = f(xi(:,1),xi(:,2),xi(:,3));
    if bctype==1
        rhsb = g(xb(:,1),xb(:,2),xb(:,3));
    elseif bctype>=2
        rhsb = g(nr(:,1),nr(:,2),nr(:,3),xb(:,1),xb(:,2),xb(:,3));
    end    
end
%% Solve
c = lhs\[rhsi;rhsb;zeros(size(z,2),1)];
c = c(1:end-size(z,2));
sol = V*c; sol_inner = sol(1:Ni);

if dim==2
    true_sol = u(x(:,1),x(:,2));
elseif dim==3
    true_sol = u(x(:,1),x(:,2),x(:,3));
end

if bctype==2
    sol = sol-mean(sol);
    true_sol = true_sol - mean(true_sol);
end

%% plot stuff
plot(sol,'.')
hold on
plot(true_sol,'o');
l2 = norm(sol-true_sol,2)./norm(true_sol,2)
linf = norm(sol-true_sol,inf)./norm(true_sol,inf)