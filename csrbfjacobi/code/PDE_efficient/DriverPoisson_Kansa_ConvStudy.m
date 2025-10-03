%% Solve Poisson's equation on an arbitrary 2D/3D domain using CSRBFs combined with
%%Jacobi polynomials. We form and solve Kansa-style unsymmetric
%%collocation systems. 


%% Spatial dimension
dim = 3;

%% Define RBF

% % Wendland C2 in 3d, pd in all lower dimensions
% rbf = @(e,r) -r.^4.*(4.*r - 5*spones(r));
% drbfor = @(e,r) 20.*e.^2.*r.^3;
% if dim==1
%     lrbf = @(e,r) -20.*e.^2.*r.^2.*(2.*r - 3*spones(r));
% elseif dim==2
%     lrbf = @(e,r) -20.*e.^2.*r.^2.*(r - 3*spones(r));
% elseif dim==3
%     lrbf = @(e,r) 60.*e.^2.*r.^2;
% end

% %% Wendland C4 in 3d
% rbf = @(e,r) (r.^6.*(35.*r.^2 - 88.*r + 56*spones(r)))./3;
% drbfor = @(e,r) -(56.*e.^2.*r.^5.*(5.*r - 6*spones(r)))./3;
% if dim==2
%     lrbf = @(e,r) (112.*e.^2.*r.^4.*(20.*r.^2 - 36.*r + 15*spones(r)))./3;   
% elseif dim==3    
%     lrbf = @(e,r) 56.*e.^2.*r.^4.*(15.*r.^2 - 26.*r + 10*spones(r));
% end

% %% Wendland C6 in 3d
% rbf = @(e,r) r.^8.*(66*spones(r) - 154*r + 121*r.^2 - 32*r.^3);
% drbfor = @(e,r) 22.*e.^2.*r.^7.*(16.*r.^2 - 39.*r + 24*spones(r));
% if dim==2
%     lrbf = @(e,r) 44*e.^2.*r.^6.*(84*spones(r)-264*r+267*r.^2-88*r.^3);
% elseif dim==3
%     lrbf = @(e,r) -66.*e.^2.*r.^6.*((r - spones(r)).^2 - 6.*r + 64.*(r - spones(r)).^3 + 7*spones(r));
% end

%% Wendland C8 in 3d
rbf = @(e,r) (6552290047271679.*r.^10.*(4134.*r.^2 - 3536.*r - 2166.*r.^3 + 429.*r.^4 + 1144*spones(r)))./32761450236358396;
drbfor = @(e,r) -(85179770614531827.*e.^2.*r.^9.*(1056.*r - 852.*r.^2 + 231.*r.^3 - 440*spones(r)))./16380725118179198;
if dim==2
    lrbf = @(e,r) (85179770614531827.*e.^2.*r.^8.*(11022.*r.^2 - 7700.*r - 6924.*r.^3 + 1617.*r.^4 + 1980*spones(r)))./8190362559089599;
elseif dim==3
    lrbf = @(e,r) (1277696559217977405.*e.^2.*r.^8.*(1540.*r.^2 - 1056.*r - 980.*r.^3 + 231.*r.^4 + 264*spones(r)))./16380725118179198;
end

%% Type of boundary condition?
bctype=1; %1-Dirichlet, 2-Neumann, 3- Robin

%% This partially determines the number of polynomial terms
if dim==1
    fac = 0.1;
elseif dim==2
    fac = 0.1;
elseif dim==3
    fac = 2;
end

%% Setup anonymous functions for PDE solution and rhs for interior
%and boundary. 
if dim==1
    syms x;
    syms nrx;
    u = 1 + sin(pi*x);
    f = -diff(u,x,2);
    if bctype==1
        g = u;
    elseif bctype==2
        g = nrx.*diff(u,x);
    elseif bctype==3
        g = nrx.*diff(u,x)+ u; %n.grad(u) = g;
    end 
elseif dim==2
    syms x y;
    syms nrx nry;
    u = 1 + sin(pi*x).*cos(pi*y);     
    %u = 1 + x.^5 + y.^2 + x.^2.*y.^2; %verify polynomial reproduction
    f = -(diff(u,x,2) + diff(u,y,2)); %lap u = f
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
    %u = x.^3 + y.^4 + z.^5;
    f = -(diff(u,x,2) + diff(u,y,2) + diff(u,z,2)); %lap u = f    
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
if dim==1
    clear x nrx;
elseif dim==2
    clear x y nrx nry;
elseif dim==3
    clear x y z nrx nry nrz;
end

%% Load up the node set
if dim==1
    N = 2.^linspace(4,15,11);
    for k=1:11
        xi = linspace(-1,1,N(k)-2).';
        xi = xi(2:end-1,:);
        xb = [-1;1];
        nr = [-1;1];
        st.fullintnodes{k,1} = xi;
        st.bdrynodes{k,1} = xb;
        st.normals{k,1} = nr;
    end
    clear N xb xi k;    
elseif dim==2
    %st = load('DiskPoissonNodes.mat');
    st = load('DiskPoissonNodesClustered.mat');
    %st = load('StarPoissonNodesClustered.mat');
    %st = load('DomainPoissonNodes2.mat');
elseif dim==3
    st = load('SpherePoissonNodes.mat');
    %st = load('RBCPoissonNodesClustered.mat');
    %st = load('BumpySpherePoissonNodes.mat');
    %st = load('BumpySpherePoissonNodesClustered.mat');
end
start_nodes = 1;
end_nodes = min([6,size(st.fullintnodes,1)]);

%% Start convergence study
%ep = 2;
for k=start_nodes:end_nodes    
    xi = st.fullintnodes{k};
    xb = st.bdrynodes{k};
    nr = st.normals{k};
    x = [xi;xb];
    
    %% Pick family of orthogonal polynomials
    %%a = b = 0- Legendre
    %%a = b = -0.5 - Chebyshev first kind
    %%a = b = 0.5 - Chebyshev second kind
    %%a = b - Ultraspherical
    alph = 0; 
    recurrence = @(N) jacobi_recurrence(N,alph,alph);    
    ell = floor(nthroot(fac*length(xi),dim)); 
    ell = max([ell,1])
    a = total_degree_indices(dim,ell);      
    N = size(a,1);      
    
    %% We'll now select the support. By setting this to be smaller
    %%than the separation distance, we can guarantee an identity matrix
    %%for interpolation, and a zero matrix for differentiation at the
    %%interpolation nodes
    tree = KDTreeSearcher(x);

    
    if bctype~=2 %%Dirichlet and Robin need different tech from Neumann
        [~,dist] = knnsearch(tree,x,'k',2);
        dist = dist(:,2);
        dist = min(dist); %separation radius        
        ep = 1.2/dist; %support = 1/ep
        clear tree;
        
        %% Setup RBF approximation with compactly supported RBFs. We'll use Wendland's C^6 RBF.
        %%Here we use the shifted form.
        V = mpoly_eval(x,a,recurrence);
        rd = speye(length(x),length(x)); %diagonal matrix due to choice of support
        lconst = -lrbf(ep,1); %get the number on the diagonal
        Lap = lconst*speye(length(x),length(x)); %Laplacian of RBF, also diagonal
        Ni = length(xi); Nb = length(xb);

        if dim==1
            Lv = -mpoly_eval(x,a,recurrence,2);
        elseif dim==2
            Lv = -mpoly_eval(x,a,recurrence,[2,0]) - mpoly_eval(x,a,recurrence,[0,2]);
        elseif dim==3
            Lv = -mpoly_eval(x,a,recurrence,[2,0,0])...
                 - mpoly_eval(x,a,recurrence,[0,2,0])...
                 - mpoly_eval(x,a,recurrence,[0,0,2]);    
        end

        if bctype==1 %Dirichlet
            D = [sparse(Nb,Ni),speye(Nb,Nb)];
            Dv = V(Ni+1:end,:); %get poly
        elseif bctype==3 %Robin
            dconst = drbfor(ep,1);
            Drbf = [sparse(Nb,Ni),dconst*speye(Nb,Nb)];      
            A = [sparse(Nb,Ni),speye(Nb,Nb)];     

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

            %% Form Robin matrix for RBFs
            D = D + A;

            %% form robin matrix for polys
            Dv = Dv + V(Ni+1:end,:);
        end

        %% Get rhs for linear system
        if dim==1
            rhsi = f(xi(:,1));
            if bctype==1
                rhsb = g(xb(:,1));
            elseif bctype>=2
                rhsb = g(nr(:,1),xb(:,1));
            end
        elseif dim==2
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

        %% Efficient Solve!
        rhs = [rhsi;rhsb];        
        lhs_A = [Lap(1:Ni,:);D];
        lhs_B = [Lv(1:Ni,:);Dv];

        %% This is a block-symmetric system [[A B];[B' O]] [c_rbf;c_poly] = [rhs;zeros]
        %% We now have the following system for the polynomial
        %%coefficients:
        %% B^T A^{-1} B c_poly = B^T A^{-1} rhs
        %%Since A is diagonal, this is actually a weighted set of normal
        %%equations. We will first get the square root of the diagonal
        %%of A^{-1} and fold it into B. B = W^{1/2}B.
        Ainv_d = 1./diag(lhs_A);
        Ainv = spdiags(Ainv_d,0,Ni+Nb,Ni+Nb);        
        Wh = spdiags(sqrt(Ainv_d),0,Ni+Nb,Ni+Nb);        
        B = Wh*lhs_B;

        %%Now, qr decompose B.
        [Q,R] = qr(B,0);
        opts.UT = true;
        c_poly = linsolve(R,Q'*(Wh*rhs),opts);        
        c_rbf = Ainv*(rhs - lhs_B*c_poly);        
        c = [c_rbf; c_poly];
        
    else % different solve for Neumann case
        tree = KDTreeSearcher(xi);
        [~,dist] = knnsearch(tree,xb,'k',1);
        %dist = dist(:,2);
        dist = min(dist); %separation radius
        %ep = 1.02/dist; %support = 1/ep
        ep = 0.9/dist;
        clear tree;
        
        %% Setup RBF approximation with compactly supported RBFs. We'll use Wendland's C^6 RBF.
        %%Here we use the shifted form.
        V = mpoly_eval(x,a,recurrence);
        %rd = speye(length(x),length(x)); %diagonal matrix due to choice of support
        rd = DistanceMatrixCSRBF(x,x,ep);        
        %Lap = lconst*speye(length(x),length(x)); %Laplacian of RBF, also diagonal
        Lap = lrbf(ep,rd);
        Ni = length(xi); Nb = length(xb);
        
        %% Get polynomial laplacian
        if dim==1
            Lv = mpoly_eval(x,a,recurrence,2);
        elseif dim==2
            Lv = mpoly_eval(x,a,recurrence,[2,0]) + mpoly_eval(x,a,recurrence,[0,2]);
        elseif dim==3
            Lv = mpoly_eval(x,a,recurrence,[2,0,0])...
                 + mpoly_eval(x,a,recurrence,[0,2,0])...
                 + mpoly_eval(x,a,recurrence,[0,0,2]);    
        end
        
        %% Get RBF gradients            
        %rdb =  DistanceMatrixCSRBF(xb,x,ep); %not diagonal matrix 
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

        %% Get polynomial gradients
        if dim==1
            Dvx =  mpoly_eval(xb,a,recurrence,1);
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
        
        if dim==1
            rhsi = f(xi(:,1));
            rhsb = g(nr(:,1),xb(:,1));
        elseif dim==2
            rhsi = f(xi(:,1),xi(:,2));
            rhsb = g(nr(:,1),nr(:,2),xb(:,1),xb(:,2));
        elseif dim==3
            rhsi = f(xi(:,1),xi(:,2),xi(:,3));
            rhsb = g(nr(:,1),nr(:,2),nr(:,3),xb(:,1),xb(:,2),xb(:,3));
        end        
        
        rhs = [rhsi;rhsb];
        Vt = [Lv(1:Ni,:).', Dv.'];        
        lhsa = [[Lap(1:Ni,:),Lv(1:Ni,:)];[D,Dv];[Vt,sparse(N,N)]];
        z = spnull(lhsa);
        if size(z,2)>1
            here = 1;
        end
        lhs = [[lhsa,z];[z.',sparse(size(z,2),size(z,2))]];
        rhs = [rhs;zeros(N,1);zeros(size(z,2),1)];
        c = lhs\rhs;
        c_rbf = c(1:Ni+Nb,:);
        c_poly = c(Ni+Nb+1:Ni+Nb+N,:);
    end

    %% Evaluate interpolant.
    %rd = [rd(1:Ni,:);rdb];
    sol = rbf(ep,rd)*c_rbf + V*c_poly;
    sol_inner = sol(1:Ni);

    if dim==1
        true_sol = u(x(:,1));
    elseif dim==2
        true_sol = u(x(:,1),x(:,2));
    elseif dim==3
        true_sol = u(x(:,1),x(:,2),x(:,3));
    end

    if bctype==2
        sol = sol-mean(sol);
        true_sol = true_sol - mean(true_sol);
    end

    sN(k) = nthroot(length(xi),dim);
    eps(k) = ep;
    deg(k) = ell;
    
    el2(k) = norm(sol-true_sol,2)./norm(true_sol,2)
    elinf(k) = norm(sol-true_sol,inf)./norm(true_sol,inf)

end
%save('Poisson_PolyReprod.mat','el2','elinf','sN','dim','deg','eps');
figure(1)
semilogy(deg,el2,'LineWidth',2);
hold on
semilogy(deg,el2,'blacko');
title('Convergence for Poisson equation');
xlabel('Polynomial degree');
ylabel('Relative l2 error');
axis tight;

% figure(2)
% semilogy(sN,el2,'LineWidth',2);
% hold on
% semilogy(sN,el2,'blacko');
% title('Convergence for Poisson equation');
% xlabel('d-th root(N)');
% ylabel('Relative l2 error');
% axis tight;
% 
% figure(3)
% semilogy(sN,nonzeros,'LineWidth',2);
% hold on
% semilogy(sN,nonzeros,'blacko');
% title('Percent non-zeros in system matrix');
% xlabel('d-th root(N)');
% ylabel('Percent non-zeros');
% axis tight;
% 
% figure(4)
% spy(lhs);
% title('System matrix');
% 
% figure(5)
% semilogy(sN,conds,'LineWidth',2);
% hold on
% semilogy(sN,conds,'blacko');
% title('System matrix');


