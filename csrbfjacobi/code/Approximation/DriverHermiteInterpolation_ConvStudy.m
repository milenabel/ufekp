%% Hermite Interpolation in 1D/2D/3D domains using CSRBFs combined with
%%Jacobi polynomials. Done efficiently using RBFs with very narrow
%%supports. Interpolate f in interior, lap f on boundary.


%% Spatial dimension
dim = 2;

%% Define RBF
% %% Wendland C4 in 3d
% rbf = @(e,r) (r.^6.*(35.*r.^2 - 88.*r + 56*spones(r)))./3;
% drbfor = @(e,r) -(56.*e.^2.*r.^5.*(5.*r - 6*spones(r)))./3;
% if dim==2
%     lrbf = @(e,r) (112.*e.^2.*r.^4.*(20.*r.^2 - 36.*r + 15*spones(r)))./3;   
%     l2rbf = @(e,r) 2240.*e.^4.*r.^2.*(12.*r.^2 - 17.*r + 3*spones(r));
% elseif dim==3    
%     lrbf = @(e,r) 56.*e.^2.*r.^4.*(15.*r.^2 - 26.*r + 10*spones(r));
%     l2rbf = @(e,r) 1680.*e.^4.*r.^2.*(21.*r.^2 - 30.*r + 4*spones(r));
% end

%% Wendland C6 in 3d
rbf = @(e,r) r.^8.*(66*spones(r) - 154*r + 121*r.^2 - 32*r.^3);
drbfor = @(e,r) 22.*e.^2.*r.^7.*(16.*r.^2 - 39.*r + 24*spones(r));
if dim==2
    lrbf = @(e,r) 44*e.^2.*r.^6.*(84*spones(r)-264*r+267*r.^2-88*r.^3);
    l2rbf = @(e,r) -1056.*e.^4.*r.^4.*(483.*r - 679.*r.^2 + 297.*r.^3 - 105*spones(r));
elseif dim==3
    lrbf = @(e,r) -66.*e.^2.*r.^6.*((r - spones(r)).^2 - 6.*r + 64.*(r - spones(r)).^3 + 7*spones(r));
    l2rbf = @(e,r) -7920.*e.^4.*r.^4.*(70.*r - 105.*r.^2 + 48.*r.^3 - 14*spones(r));
end

% %% Wendland C8 in 3d
% rbf = @(e,r) (6552290047271679.*r.^10.*(4134.*r.^2 - 3536.*r - 2166.*r.^3 + 429.*r.^4 + 1144*spones(r)))./32761450236358396;
% drbfor = @(e,r) -(85179770614531827.*e.^2.*r.^9.*(1056.*r - 852.*r.^2 + 231.*r.^3 - 440*spones(r)))./16380725118179198;
% if dim==2
%     lrbf = @(e,r) (85179770614531827.*e.^2.*r.^8.*(11022.*r.^2 - 7700.*r - 6924.*r.^3 + 1617.*r.^4 + 1980*spones(r)))./8190362559089599;
% elseif dim==3
%     lrbf = @(e,r) (1277696559217977405.*e.^2.*r.^8.*(1540.*r.^2 - 1056.*r - 980.*r.^3 + 231.*r.^4 + 264*spones(r)))./16380725118179198;
% end


%% This partially determines the number of polynomial terms
if dim==1
    fac = 0.08;
elseif dim==2
    fac = 0.2;
elseif dim==3
    fac = 2.5;
end

%% Setup anonymous function for function interpolation true solution
%%We'll use Trefethen's function gallery from his Spectral methods book
%%with 2d and 3d analogues
if dim==1
    syms x;       
    %f = abs(x.^3);
    %f = exp(-x.^(-2));    
    f = 1./(1 + x.^2);
    %f = x.^(10);
    lf = diff(f,x,2);
elseif dim==2
    syms x y;        
    %f = abs(x.^3).*abs(y.^3);
    %f = exp(-x.^(-2)).*exp(-y.^(-2));    
    %f = 1./(1 + x.^2 + y.^2);
    f = x.^(5).*y.^(5);
    lf = diff(f,x,2) + diff(f,y,2);
elseif dim==3
    syms x y z;      
    %f = abs(x.^3).*abs(y.^3).*abs(z.^3);
    %f = exp(-x.^(-2)).*exp(-y.^(-2)).*exp(-z.^(-2));    
    f = 1./(1 + x.^2 + y.^2 + z.^2);
    %f = x.^(10).*y.^(10).*z.^(10);    
    lf = diff(f,x,2) + diff(f,y,2) + diff(f,z,2);
end
f = matlabFunction(f);
if dim==1
    lf = matlabFunction(lf);
    clear x;
elseif dim==2
    lf = matlabFunction(lf);
    clear x y;
elseif dim==3
    lf = matlabFunction(lf);
    clear x y z;    
end

%% Load up the node set
if dim==1
    N = 2.^(4:12); N = N';
    for k=1:length(N)
        xi = linspace(-1,1,N(k)-2).';
        xi = xi(2:end-1,:);
        xb = [-1;1];
        st.fullintnodes{k,1} = xi;
        st.bdrynodes{k,1} = xb;
    end
    clear N xb xi k;
    xe = linspace(-1,1,2^14).';    
elseif dim==2
    %% Get evaluation nodes
    st = load('DiskPoissonNodesLarge.mat');
    xe =  [st.fullintnodes{7}; st.bdrynodes{7}];
   
    st = load('DiskPoissonNodes.mat');
    %st = load('DiskPoissonNodesClustered.mat');
    %st = load('StarPoissonNodesClustered.mat');
    %st = load('DomainPoissonNodes2.mat');
elseif dim==3
    st = load('SpherePoissonNodesLarge.mat');
    xe = [st.fullintnodes{1}; st.bdrynodes{1}];
    
    st = load('SpherePoissonNodes.mat');
    
    %st = load('RBCPoissonNodesClustered.mat');
    %xe = [st.fullintnodes{7}; st.bdrynodes{7}];
    %st = load('BumpySpherePoissonNodes.mat');
    %st = load('BumpySpherePoissonNodesClustered.mat');
end
start_nodes = 1;
end_nodes = size(st.fullintnodes,1);


%% Start convergence study
for k=start_nodes:end_nodes
    xi = st.fullintnodes{k};
    xb = st.bdrynodes{k};
    x = [xi;xb];    

    %% Lower half is Lagrange, upper Hermite
    
    %xl = x(x(:,2)<=0,:);
    %xh = x(x(:,2)>0,:);    
    xl = xi; xh = xb;
    x = [xl;xh];
    Nl = length(xl); Nh = length(xh);
    
    %% Pick family of orthogonal polynomials
    %%a = b = 0- Legendre
    %%a = b = -0.5 - Chebyshev first kind
    %%a = b = 0.5 - Chebyshev second kind
    %%a = b - All ultraspherical
    alph = 0; 
    recurrence = @(N) jacobi_recurrence(N,alph,alph);    
    ell = floor(nthroot(fac*length(x),dim));    
    ell = max([ell,1])
    a = total_degree_indices(dim,ell);      
    N = size(a,1);      
    

   
    %if k==1
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
       dflag = 1;
    %end
    
    
    
    %% Setup RBF approximation with compactly supported RBFs. We'll use Wendland's C^6 RBF.
    %%Since the support is chosen to be smaller than the sep radius, we
    %%have identity matrices for rd and A
    if dflag==1
        A11 = rbf(ep,1)*speye(Nl,Nl);
        A12 = sparse(Nl,Nh);
        A21 = A12.';
        A22 = l2rbf(ep,1)*speye(Nh,Nh);
    else
        rd = DistanceMatrixCSRBF(x,x,ep);
        A11 = rbf(ep,rd(1:Nl,1:Nl));
        A12 = lrbf(ep,rd(1:Nl,Nl+1:end));
        A21 = A12.';
        A22 = l2rbf(ep,rd(Nl+1:end,Nl+1:end));
    end
    A = [[A11,A12];[A21,A22]];
    Vl = mpoly_eval(xl,a,recurrence); %polynomial in lagrange  
    Vh = mpoly_eval(xh,a,recurrence);
    Vf = [Vl;Vh]; 
    %Laplacian of poly on hermite
    if dim==1
        Lvh = mpoly_eval(xh,a,recurrence,2);
    elseif dim==2
        Lvh = mpoly_eval(xh,a,recurrence,[2,0]) + mpoly_eval(xh,a,recurrence,[0,2]);
    elseif dim==3
        Lvh = mpoly_eval(xh,a,recurrence,[2,0,0]) +...
            mpoly_eval(xh,a,recurrence,[0,2,0]) +...
            mpoly_eval(xh,a,recurrence,[0,0,2]);
    end
    V = [Vl;Lvh]; %concatenate        
    
    %% Get rhs for linear system
    if dim==1
        rhs = [f(xl(:,1));lf(xh(:,1))];
    elseif dim==2
        rhs = [f(xl(:,1),xl(:,2));lf(xh(:,1),xh(:,2))];
    elseif dim==3
        rhs = [f(xl(:,1),xl(:,2),xl(:,3));lf(xh(:,1),xh(:,2),xh(:,3))];
    end
    
    %% Let's take the Schur complement approach. Then, we have two
    %%solves formally: 
    %%1. V^T A^{-1} V c_poly = V^T A^{-1} f
    %%2. A c_rbf = f - V c_poly
    [la,ua,pa,qa] = lu(A);        
    tic 
    lhs_schur = full(V.'*(qa*(ua\(la\(pa*V)))));
    rhs_schur = full(V.'*(qa*(ua\(la\(pa*rhs)))));     
%     [Us,Ss,Vs] =svd(lhs_schur,0);
%     Ss = diag(Ss); cutoff = Ss>1e-13;
%     Us = Us(:,cutoff); Ss = Ss(cutoff); Vs = Vs(:,cutoff);
%     c_poly = Vs*((Us.'*rhs_schur)./Ss);
    c_poly = lhs_schur\rhs_schur;
    c_rbf = qa*(ua\(la\(pa*(rhs - V*c_poly))));
    toc     
   
    %% Evaluate interpolant and derivative at x
    rd = DistanceMatrixCSRBF(x,x,ep);
    Ae = [rbf(ep,rd(:,1:Nl)),lrbf(ep,rd(:,Nl+1:end))]; 
    solx = Ae*c_rbf + Vf*c_poly;

    %% Build RBF and poly evaluation matrices
    rd = DistanceMatrixCSRBF(xe,x,ep);
    Ae = [rbf(ep,rd(:,1:Nl)),lrbf(ep,rd(:,Nl+1:end))]; 
    Ve = mpoly_eval(xe,a,recurrence);    

    %% Do RBF and poly evaluations
    sol = Ae*c_rbf + Ve*c_poly;


    %% Compute true solutions
    if dim==1
        true_solx = f(x(:,1));
        true_sol = f(xe(:,1));
    elseif dim==2
        true_solx = f(x(:,1),x(:,2));
        true_sol = f(xe(:,1),xe(:,2));
    elseif dim==3
        true_solx = f(x(:,1),x(:,2),x(:,3));
        true_sol = f(xe(:,1),xe(:,2),xe(:,3));
    end


    %% Store a bunch of stuff
    sN(k) = nthroot(length(xi),dim);
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
    x_poly_el2(k) =norm(Vf*c_poly-true_solx,2)./norm(true_solx,2);
    x_poly_elinf(k) =norm(Vf*c_poly-true_solx,inf)./norm(true_solx,inf);
    

end
fname = strcat('HermInterp_',num2str(dim),'d.mat');
save(fname,'el2','elinf','poly_el2','poly_elinf','sN','dim','deg','eps');

figure(1)
semilogy(deg,el2,'LineWidth',2);
hold on
semilogy(deg,el2,'blacko');
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


