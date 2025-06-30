%% Interpolation in 1D/2D/3D domains using CSRBFs combined with
%%Jacobi polynomials. Done efficiently using RBFs with very narrow
%%supports. Uses Halton points on [-1,1]^d


%% Spatial dimension
s_dim = 2;

%% Define RBF
%% Wendland C2 in 3d, pd in all lower dimensions
rbf = @(e,r) -r.^4.*(4.*r - 5*spones(r));
drbfor = @(e,r) 20.*e.^2.*r.^3;
if s_dim==1
    lrbf = @(e,r) -20.*e.^2.*r.^2.*(2.*r - 3*spones(r));
elseif s_dim==2
    lrbf = @(e,r) -20.*e.^2.*r.^2.*(r - 3*spones(r));
elseif s_dim==3
    lrbf = @(e,r) 60.*e.^2.*r.^2;
end

% %% Wendland C4 in 3d
% rbf = @(e,r) (r.^6.*(35.*r.^2 - 88.*r + 56*spones(r)))./3;
% drbfor = @(e,r) -(56.*e.^2.*r.^5.*(5.*r - 6*spones(r)))./3;
% if s_dim==2
%     lrbf = @(e,r) (112.*e.^2.*r.^4.*(20.*r.^2 - 36.*r + 15*spones(r)))./3;   
% elseif s_dim==3    
%     lrbf = @(e,r) 56.*e.^2.*r.^4.*(15.*r.^2 - 26.*r + 10*spones(r));
% end

% %% Wendland C6 in 3d
% rbf = @(e,r) r.^8.*(66*spones(r) - 154*r + 121*r.^2 - 32*r.^3);
% drbfor = @(e,r) 22.*e.^2.*r.^7.*(16.*r.^2 - 39.*r + 24*spones(r));
% % if s_dim==2
% %     lrbf = @(e,r) 44*e.^2.*r.^6.*(84*spones(r)-264*r+267*r.^2-88*r.^3);
% % elseif s_dim==3
% %     lrbf = @(e,r) -66.*e.^2.*r.^6.*((r - spones(r)).^2 - 6.*r + 64.*(r - spones(r)).^3 + 7*spones(r));
% % end

% %% Wendland C8 in 3d
% rbf = @(e,r) (6552290047271679.*r.^10.*(4134.*r.^2 - 3536.*r - 2166.*r.^3 + 429.*r.^4 + 1144*spones(r)))./32761450236358396;
% drbfor = @(e,r) -(85179770614531827.*e.^2.*r.^9.*(1056.*r - 852.*r.^2 + 231.*r.^3 - 440*spones(r)))./16380725118179198;
% if s_dim==2
%     lrbf = @(e,r) (85179770614531827.*e.^2.*r.^8.*(11022.*r.^2 - 7700.*r - 6924.*r.^3 + 1617.*r.^4 + 1980*spones(r)))./8190362559089599;
% elseif s_dim==3
%     lrbf = @(e,r) (1277696559217977405.*e.^2.*r.^8.*(1540.*r.^2 - 1056.*r - 980.*r.^3 + 231.*r.^4 + 264*spones(r)))./16380725118179198;
% end


%% This partially determines the number of polynomial terms
if s_dim==1
    fac = 0.08;
elseif s_dim==2
    fac = 0.3;
elseif s_dim==3
    fac = 3;
end

%% Setup anonymous function for function interpolation true solution
%%We'll use Trefethen's function gallery from his Spectral methods book
%%with 2d and 3d analogues
if s_dim==1
    syms x;       
    %f = abs(x.^3);
    %f = exp(-x.^(-2));    
    f = 1./(1 + x.^2);
    %f = x.^(10);
    dfx = diff(f,x);
elseif s_dim==2
    syms x y;    
    f = abs(x.^3).*abs(y.^3);
    %f = exp(-x.^(-2)).*exp(-y.^(-2));    
    %f = 1./(1 + 25*(x.^2 + y.^2));
    %f = x.^(10).*y.^(10);
    dfx = diff(f,x); dfy = diff(f,y);
elseif s_dim==3
    syms x y z;      
    %f = abs(x.^3).*abs(y.^3).*abs(z.^3);
    %f = exp(-x.^(-2)).*exp(-y.^(-2)).*exp(-z.^(-2));    
    f = 1./(1 + 25*(x.^2 + y.^2 + z.^2));
    %f = x.^(10).*y.^(10).*z.^(10);    
    dfx = diff(f,x); dfy = diff(f,y); dfz = diff(f,z);
end
f = matlabFunction(f);
if s_dim==1
    dfx = matlabFunction(dfx);
    clear x;
elseif s_dim==2
    dfx = matlabFunction(dfx);
    dfy = matlabFunction(dfy);
    clear x y;
elseif s_dim==3
    dfx = matlabFunction(dfx);
    dfy = matlabFunction(dfy);
    dfz = matlabFunction(dfz);
    clear x y z;
end

%% Load up the node set
pxe = haltonset(s_dim,'Skip',1e3,'Leap',1e2); %Quasi-random node set
pxe = scramble(pxe,'RR2');
N = (4:4:70).^s_dim; N = N';
for k=1:length(N)
    p = net(pxe,N(k));
    p = p*2 - 1;
    st.nodes{k,1} = p;    
end
xe = net(pxe,1000);
xe = xe*2 - 1;
clear p N pxe;
start_nodes = 1;
end_nodes = size(st.nodes,1);



%% Start convergence study
for k=start_nodes:end_nodes
    x = st.nodes{k};
    
    %% Pick family of orthogonal polynomials
    %%a = b = 0- Legendre
    %%a = b = -0.5 - Chebyshev first kind
    %%a = b = 0.5 - Chebyshev second kind
    %%a = b - All ultraspherical
    alph = 0; 
    recurrence = @(N) jacobi_recurrence(N,alph,alph);    
    ell = floor(nthroot(fac*length(x),s_dim));    
    ell = max([ell,1])
    a = total_degree_indices(s_dim,ell);      
    N = size(a,1);      
    
    
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
    V = mpoly_eval(x,a,recurrence);    
    %Ni = length(xi); Nb = length(xb);    
    
    %% Get rhs for linear system
    if s_dim==1
        rhs = f(x(:,1));
    elseif s_dim==2
        rhs = f(x(:,1),x(:,2));
    elseif s_dim==3
        rhs = f(x(:,1),x(:,2),x(:,3));
    end
    
    %% Let's take the Schur complement approach. Then, we have two
    %%solves formally: 
    %%1. V^T A^{-1} V c_poly = V^T A^{-1} f
    %%2. A c_rbf = f - V c_poly
    %%If A is identity matrix, eq 1 above becomes
    %%V^T*V c_poly = V^T f, which is the set of normal equations
    %%for the least squares solution of the rectangular system
    %%V c_poly = f. We can use QR decompositions to solve this.
    %%eq 2 then becomes c_rbf = f - V c_poly.
    tic
    [Q,R] = qr(V,0); %use QR for least squares
    opts.UT = true;
    c_poly = linsolve(R,Q'*rhs,opts);        
    c_rbf = rhs - V*c_poly;        
    toc
    
    %% Setup fast Jacobi diff
    [a_structure, dimension_jumps] = mindex_derivative_analysis(a);
    Cmat = jacobi_sparse_diffmat_inv(size(a_structure, 1), alph);
    [lc,uc,pc,qc] = lu(Cmat);

    
    %% Evaluate interpolant and derivative at x
%     %solx = rhs; %actually this, since c_rbf = rhs - V*c_poly
%     if s_dim==1
%         cdx = mjacobi_symm_faster_differentiation(c_poly, a, alph, [1], a_structure, dimension_jumps, lc,uc,pc,qc);
%     elseif s_dim==2
%         cdx = mjacobi_symm_faster_differentiation(c_poly, a, alph, [1 0], a_structure, dimension_jumps, lc,uc,pc,qc);
%         cdy = mjacobi_symm_faster_differentiation(c_poly, a, alph, [0 1], a_structure, dimension_jumps, lc,uc,pc,qc);        
%     elseif s_dim==3
%         cdx = mjacobi_symm_faster_differentiation(c_poly, a, alph, [1 0 0], a_structure, dimension_jumps, lc,uc,pc,qc);
%         cdy = mjacobi_symm_faster_differentiation(c_poly, a, alph, [0 1 0], a_structure, dimension_jumps, lc,uc,pc,qc);
%         cdz = mjacobi_symm_faster_differentiation(c_poly, a, alph, [0 0 1], a_structure, dimension_jumps, lc,uc,pc,qc);
%     end
    
    %% Build RBF and poly evaluation matrices
    rde = DistanceMatrixCSRBF(xe,x,ep);
    Ae = rbf(ep,rde); Ve = mpoly_eval(xe,a,recurrence);
%     Drbf = drbfor(ep,rde); %derivative of rbf divided by r
%     if s_dim==1
%         [xj,xk] = ndgrid(xe(:,1),x(:,1));
%         Dxe = (xj-xk).*Drbf; %x der
%     elseif s_dim==2
%         [xj,xk] = ndgrid(xe(:,1),x(:,1));
%         [yj,yk] = ndgrid(xe(:,2),x(:,2));
%         Dxe = (xj-xk).*Drbf; %x der
%         Dye = (yj-yk).*Drbf; %y der
%     elseif s_dim==3
%         [xj,xk] = ndgrid(xe(:,1),x(:,1));
%         [yj,yk] = ndgrid(xe(:,2),x(:,2));
%         [zj,zk] = ndgrid(xe(:,3),x(:,3));
%         Dxe = (xj-xk).*Drbf; %x der
%         Dye = (yj-yk).*Drbf; %y der
%         Dze = (zj-zk).*Drbf; %z der
%     end    
    
    %% Do RBF and poly evaluations
    sol = Ae*c_rbf + Ve*c_poly;
%     dxe = Dxe*c_rbf + Ve*cdx;
%     if s_dim>=2
%         dye = Dye*c_rbf + Ve*cdy;
%     end
%     if s_dim>=3
%         dze = Dze*c_rbf + Ve*cdz;
%     end

    %% Compute true solutions
    if s_dim==1
        true_solx = f(x(:,1));
%        tdx = dfx(x(:,1));
        true_sol = f(xe(:,1));
%         tdxe = dfx(xe(:,1));
    elseif s_dim==2
        true_solx = f(x(:,1),x(:,2));
        true_sol = f(xe(:,1),xe(:,2));
%         tdx = dfx(x(:,1),x(:,2)); tdxe = dfx(xe(:,1),xe(:,2));
%         tdy = dfy(x(:,1),x(:,2)); tdye = dfy(xe(:,1),xe(:,2));
    elseif s_dim==3
        true_solx = f(x(:,1),x(:,2),x(:,3));
        true_sol = f(xe(:,1),xe(:,2),xe(:,3));
%         tdx = dfx(x(:,1),x(:,2),x(:,3)); tdxe = dfx(xe(:,1),xe(:,2),xe(:,3));
%         tdy = dfy(x(:,1),x(:,2),x(:,3)); tdye = dfy(xe(:,1),xe(:,2),xe(:,3));
%         tdz = dfz(x(:,1),x(:,2),x(:,3)); tdze = dfz(xe(:,1),xe(:,2),xe(:,3));
    end


    %% Store a bunch of stuff
    sN(k) = nthroot(length(x),s_dim);
    eps(k) = ep;
    deg(k) = ell;
    conds(k) = cond(V);
    
    el2(k) = norm(sol-true_sol,2)./norm(true_sol,2)
    elinf(k) = norm(sol-true_sol,inf)./norm(true_sol,inf)
    
    %% Diagnostic errors
    poly_el2(k) =norm(Ve*c_poly-true_sol,2)./norm(true_sol,2);
    poly_elinf(k) =norm(Ve*c_poly-true_sol,inf)./norm(true_sol,inf);
    

end
fname = strcat('Interp_',num2str(s_dim),'d.mat');
save(fname,'el2','elinf','poly_el2','poly_elinf','sN','s_dim','deg','eps');

figure(1)
semilogy(deg,el2,'LineWidth',2);
hold on
semilogy(deg,el2,'blacko');
%semilogy(deg,x_el2,'black','LineWidth',2);
%semilogy(deg,x_el2,'blackx');
title('Convergence at Eval pts');
xlabel('Polynomial degree');
ylabel('Relative l2 error');
axis tight;

% figure(2)
% semilogy(deg,dx_el2,'LineWidth',2);
% hold on
% semilogy(deg,dx_el2,'blacko');
% title('x-der convergence at Eval pts');
% xlabel('Polynomial degree');
% ylabel('Relative l2 error');
% axis tight;
% 
% 
% if s_dim>=2
%     figure(3)
%     semilogy(deg,dy_el2,'LineWidth',2);
%     hold on
%     semilogy(deg,dy_el2,'blacko');
%     title('y-der convergence at Eval pts');
%     xlabel('Polynomial degree');
%     ylabel('Relative l2 error');
%     axis tight;
% end
% 
% if s_dim>=3
%     figure(4)
%     semilogy(deg,dz_el2,'LineWidth',2);
%     hold on
%     semilogy(deg,dz_el2,'blacko');
%     title('z-der convergence at Eval pts');
%     xlabel('Polynomial degree');
%     ylabel('Relative l2 error');
%     axis tight;
% end

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


