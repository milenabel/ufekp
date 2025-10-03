%% Solve Poisson's equation on an arbitrary 2D domain using a coefficient
%%space spectral method and weighted fekete samples selected from some
%%precomputed node set. Uses Akil's Jacobi polynomial transformation tricks
%%for fast evaluation and differentiation.

%% Setup anonymous functions for PDE solution and rhs for interior
%and boundary. Assume 2D.
syms x y;
syms nrx nry;
%u = sin(pi*x).*cos(pi*y); 
u = x.^3 + y.^3;
f = -diff(u,x,2) - diff(u,y,2); %lap u = f
g = nrx.*diff(u,x) + nry.*diff(u,y); %n.grad(u) = g;
u = matlabFunction(u);
f = matlabFunction(f);
g = matlabFunction(g);

%% Setup Jacobi approximations
dim = 2;
alph = 0; 
recurrence = @(N) jacobi_recurrence(N,alph,alph);
ell = 12;
a = total_degree_indices(dim,ell);
N = size(a,1);

%% Load up the node set
st = load('DiskPoissonNodesLarge.mat');
xi = st.fullintnodes{1};
xb = st.bdrynodes{1};
nr = st.normals{1};
x_candidates = [xi;xb];

%% Get the Vandermonde-like matrix for the Jacobi polynomials
V = mpoly_eval(x_candidates, a, recurrence);
weights = 1./sqrt(sum(V.^2, 2));
V = repmat(weights, [1 N]).*V;
[~,~,e] = qr(V', 0);

%% Identify interpolation grid and update weights and V
%weights = weights(e(1:N));
ei = e((e(1:N)<=length(xi)));
eb = e((e(1:N)>length(xi)));
V = V(e(1:N),:);

%% Separate out points into interior and boundary, get the normals.
%%Reorder everything so it's always interior followed by boundary. This
%%will make it easier to do boundary bordering.
nr_ind = eb - length(xi);
xi = x_candidates(ei,:);
xb = x_candidates(eb,:);
nr = nr(nr_ind,:);
x = [xi;xb];
weights = [weights(ei);weights(eb)]; %Not sure what to do with these!


%% We'll use gmres to solve the linear system. For now, no preconditioner.
%One-time setup for differentiation
[a_structure, dimension_jumps] = mindex_derivative_analysis(a);
Cmat = jacobi_sparse_diffmat_inv(size(a_structure, 1), alph);
[lc,uc,pc,qc] = lu(Cmat);

%Get rhs for linear system
rhs = [f(xi(:,1),xi(:,2));g(nr(:,1),nr(:,2),xb(:,1),xb(:,2))];
%rhs = [f(xi(:,1),xi(:,2));g(xb(:,1),xb(:,2))];

%Run gmres using special wrapper
c = rungmres(weights,x,rhs,nr,alph,a,a_structure,dimension_jumps,recurrence,lc,uc,pc,qc,length(xi));

