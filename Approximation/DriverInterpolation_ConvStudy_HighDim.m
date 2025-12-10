%% Interpolation in 1D/2D/3D domains using CSRBFs combined with
%%Jacobi polynomials. Done efficiently using RBFs with very narrow
%%supports.


%% Spatial dimension
dim = 2;
dim = min([dim,100]);


%% Define RBF
%% Wendland C2 in 100d, pd in all lower dimensions
rbf = @(e,r) spones(r) - 49608.*(r - 1).^3 - 948753.*(r - 1).^4 - 12650040.*(r - 1).^5 - (4437073160935833123.*(r - 1).^6)./34359738368 - 1062603360.*(r - 1).^7 - 7283260530.*(r - 1).^8 - (11420181159456275685.*(r - 1).^9)./268435456 - 215376418530.*(r - 1).^10 - 957228526800.*(r - 1).^11 - (1978178352693903837.*(r - 1).^12)./524288 - 13298113226160.*(r - 1).^13 - 42189847318710.*(r - 1).^14 - (3970190120441610717.*(r - 1).^15)./32768 - 316423854890325.*(r - 1).^16 - 754453740287520.*(r - 1).^17 - (1687293876056355999.*(r - 1).^18)./1024 - 3305698625207160.*(r - 1).^19 - 6106359960452115.*(r - 1).^20 - (5328296450704281123.*(r - 1).^21)./512 - 16390755683318835.*(r - 1).^22 - 23890459629516480.*(r - 1).^23 - (8258859498591634095.*(r - 1).^24)./256 - 40396595373546048.*(r - 1).^25 - 46935146868302700.*(r - 1).^26 - (6479414408556917001.*(r - 1).^27)./128 - 50689958617766916.*(r - 1).^28 - 47129361269137056.*(r - 1).^29 - (10413344585180756655.*(r - 1).^30)./256 - 32577899494795200.*(r - 1).^31 - 24195877437280185.*(r - 1).^32 - (8525274321126850083.*(r - 1).^33)./512 - 10605783089206305.*(r - 1).^34 - 6244097403169080.*(r - 1).^35 - (6947680666114408149.*(r - 1).^36)./2048 - 1697520915646920.*(r - 1).^37 - 780512175396135.*(r - 1).^38 - (10776230326912945419.*(r - 1).^39)./32768 - 126569541956130.*(r - 1).^40 - 44327044087200.*(r - 1).^41 - (7373210223677276637.*(r - 1).^42)./524288 - 4020359812560.*(r - 1).^43 - 1029020666310.*(r - 1).^44 - (1962843636781547361.*(r - 1).^45)./8388608 - 46820960550.*(r - 1).^46 - 8146625760.*(r - 1).^47 - (10427121928199208483.*(r - 1).^48)./8589934592 - 151800480.*(r - 1).^49 - 15496299.*(r - 1).^50 - 1240200.*(r - 1).^51 - 72981.*(r - 1).^52 - 2808.*(r - 1).^53 - 53.*(r - 1).^54 - 1431.*(r - 1).^2;


%% This partially determines the number of polynomial terms
fac = 0.25*sqrt(dim);
ind_type = 1; %1 - total degree, 2 - hyperbolic cross

%% Setup anonymous function for function interpolation true solution
%We'll use Trefethen's function gallery from his Spectral methods book
%f = @(x) prod(abs(x.^3),2);
%f = @(x) exp( - sum(x.^(-2),2));
%f = @(x) prod(x.^(10),2);
f = @(x) 1./(1 + sum(x.^2,2));

%% Load up the Halton node set in dim dimensions
h = [0.7,0.6,0.5,0.4,0.3]';
%h = [0.5,0.25,0.125,0.06,0.03,0.01]';
start_nodes = 1;
end_nodes = size(h,1);

%% Get a fixed set of eval points
he = 0.2;
Ne = floor(1./he.^(dim));
pxe = haltonset(dim); %Quasi-random node set
pxe = scramble(pxe,'RR2');
xe = net(pxe,Ne);
xmax = max(xe, [], 1);
xmin = min(xe, [], 1);
c1 = 0.5 .* (xmax + xmin);
c2 = 2 ./ (xmax - xmin);
xe = (xe - repmat(c1, [size(xe,1) 1])) .* repmat(c2, [size(xe,1) 1]);    


%% Linear algebra solve
solve_type=2; %1 - QR, 2 - SVD

%% Start convergence study
for k=start_nodes:end_nodes
    x = net(pxe,floor(1./h(k).^dim));
    
    %% Pick family of orthogonal polynomials
    %%a = b = 0- Legendre
    %%a = b = -0.5 - Chebyshev first kind
    %%a = b = 0.5 - Chebyshev second kind
    %%a = b - All ultraspherical
    alph = 0; 
    recurrence = @(N) jacobi_recurrence(N,alph,alph);    
    ell = floor(fac*nthroot(size(x,1),dim));    
    ell = max([ell,1])
    
    if ind_type==1
        a = total_degree_indices(dim,ell);      
    else
        a = hyperbolic_cross_indices(dim,ell);
    end
    N = size(a,1);      
    
    %% Rescale points to lie on [-1,1]^d
    xmax = max(x, [], 1);
    xmin = min(x, [], 1);
    c1 = 0.5 .* (xmax + xmin);
    c2 = 2 ./ (xmax - xmin);

    x = (x - repmat(c1, [size(x,1) 1])) .* repmat(c2, [size(x,1) 1]);    
    
    %% We'll now select the support. By setting this to be smaller
    %%than the separation distance, we can guarantee an identity matrix
    %%for interpolation, and a zero matrix for differentiation at the
    %%interpolation nodes
    tree = KDTreeSearcher(x);
    [~,dist] = knnsearch(tree,x,'k',2);
    dist = dist(:,2);
    dist = min(dist); %separation radius
    ep = 1.002/dist; %support = 1/ep
    clear tree;

    %% Setup RBF approximation with compactly supported RBFs. We'll use Wendland's C^6 RBF.
    %%Since the support is chosen to be smaller than the sep radius, we
    %%have identity matrices for rd and A
    rd = speye(size(x,1),size(x,1)); 
    A = speye(size(rd)); 
    V = mpoly_eval(x,a,recurrence);        
    
    %% Get rhs for linear system
    rhs = f(x);
    
    %% Let's take the Schur complement approach. Then, we have two
    %%solves formally: 
    %%1. V^T A^{-1} V c_poly = V^T A^{-1} f
    %%2. A c_rbf = f - V c_poly
    %%If A is identity matrix, eq 1 above becomes
    %%V^T*V c_poly = V^T f, which is the set of normal equations
    %%for the least squares solution of the rectangular system
    %%V c_poly = f. eq 2 then becomes c_rbf = f - V c_poly.
    
    if solve_type==1 %QR
        tic
        [Q,R] = qr(V,0); %use QR for least squares
        opts.UT = true;
        c_poly = linsolve(R,Q'*rhs,opts);        
        c_rbf = rhs - V*c_poly;        
        toc
    else %SVD
        tic
        [Uv,Sv,Vv] = svd(V,0);
        c_poly = Vv*(Sv\(Uv.'*rhs));
        c_rbf = rhs - V*c_poly;        
        toc
    end
    lhs = A;
    
    
    %% Evaluate interpolant and derivative at x
    solx = rhs; %actually this, since c_rbf = rhs - V*c_poly

    %% Build RBF and poly evaluation matrices
    rde = DistanceMatrixCSRBF(xe,x,ep);
    Ae = rbf(ep,rde); Ve = mpoly_eval(xe,a,recurrence);
     
    %% Do RBF and poly evaluations
    sol = Ae*c_rbf + Ve*c_poly;
 
    %% Compute true solutions
    true_solx = f(x);
    true_sol = f(xe);

    %% Store a bunch of stuff
    sN(k) = nthroot(length(x),dim);
    eps(k) = ep;
    deg(k) = ell;
    conds(k) = cond(V);
    
    el2(k) = norm(sol-true_sol,2)./norm(true_sol,2)
    elinf(k) = norm(sol-true_sol,inf)./norm(true_sol,inf)
    
    %% Diagnostic errors
    poly_el2(k) =norm(Ve*c_poly-true_sol,2)./norm(true_sol,2);
    poly_elinf(k) =norm(Ve*c_poly-true_sol,inf)./norm(true_sol,inf);
    x_el2(k) = norm(solx-true_solx,2)./norm(true_solx,2);
    x_elinf(k) = norm(solx-true_solx,inf)./norm(true_solx,inf);
    x_poly_el2(k) =norm(V*c_poly-true_solx,2)./norm(true_solx,2);
    x_poly_elinf(k) =norm(V*c_poly-true_solx,inf)./norm(true_solx,inf);
    

end
fname = strcat('Interp_',num2str(dim),'d.mat');
save(fname,'el2','elinf','poly_el2','poly_elinf','sN','dim','deg','eps');

figure(1)
semilogy(deg,el2,'LineWidth',2);
hold on
semilogy(deg,el2,'blacko');
semilogy(deg,x_el2,'black','LineWidth',2);
semilogy(deg,x_el2,'blackx');
title('Convergence at Eval pts');
xlabel('Polynomial degree');
ylabel('Relative l2 error');
axis tight;

% figure(2)
% semilogy(sN,el2,'LineWidth',2);
% hold on
% semilogy(sN,el2,'blacko');
% semilogy(sN,x_el2,'black','LineWidth',2);
% semilogy(sN,x_el2,'blackx');
% title('Convergence at Eval pts');
% xlabel('d-th root(N)');
% ylabel('Relative l2 error');
% axis tight;

% 
% figure(3)
% semilogy(deg,poly_el2,'LineWidth',2);
% hold on
% semilogy(deg,poly_el2,'blacko');
% semilogy(deg,x_poly_el2,'black','LineWidth',2);
% hold on
% semilogy(deg,x_poly_el2,'blackx');
% title('Convergence for Polynomial part');
% xlabel('Polynomial degree');
% ylabel('Relative l2 error');
% axis tight;
% 
% figure(4)
% semilogy(sN,poly_el2,'LineWidth',2);
% hold on
% semilogy(sN,poly_el2,'blacko');
% semilogy(sN,x_poly_el2,'black','LineWidth',2);
% hold on
% semilogy(sN,x_poly_el2,'blackx');
% title('Convergence for Polynomial part');
% xlabel('d-th root(N)');
% ylabel('Relative l2 error');
% axis tight;
% 


