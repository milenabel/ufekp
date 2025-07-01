%% Solve Poisson's equation on an arbitrary 2D domain using CSRBFs combined with
%%Jacobi polynomials. We naively form and solve Kansa-style unsymmetric
%%collocation systems. Stepping stone to fast algorithms. 

%% Solve Poisson's equation on an arbitrary 2D/3D domain using CSRBFs combined with
%%Jacobi polynomials. We form and solve Kansa-style unsymmetric
%%collocation systems. 


%% Spatial dimension
dim = 2;

%% Define RBF

% %% Wendland C4
% rbf = @(e,r) (r.^6.*(35.*r.^2 - 88.*r + 56*spones(r)))./3;
% drbfor = @(e,r) -(56.*e.^2.*r.^5.*(5.*r - 6*spones(r)))./3;
% if dim==2
%     lrbf = @(e,r) (112.*e.^2.*r.^4.*(20.*r.^2 - 36.*r + 15*spones(r)))./3;   
% elseif dim==3    
%     lrbf = @(e,r) 56.*e.^2.*r.^4.*(15.*r.^2 - 26.*r + 10*spones(r));
% end

%% Wendland C6
rbf = @(e,r) r.^8.*(66*spones(r) - 154*r + 121*r.^2 - 32*r.^3);
drbfor = @(e,r) 22.*e.^2.*r.^7.*(16.*r.^2 - 39.*r + 24*spones(r));
if dim==2
    lrbf = @(e,r) 44*e.^2.*r.^6.*(84*spones(r)-264*r+267*r.^2-88*r.^3);
elseif dim==3
    lrbf = @(e,r) -66.*e.^2.*r.^6.*((r - spones(r)).^2 - 6.*r + 64.*(r - spones(r)).^3 + 7*spones(r));
end

% %% Wendland C8
% rbf = @(e,r) (6552290047271679.*r.^10.*(4134.*r.^2 - 3536.*r - 2166.*r.^3 + 429.*r.^4 + 1144*spones(r)))./32761450236358396;
% drbfor = @(e,r) -(85179770614531827.*e.^2.*r.^9.*(1056.*r - 852.*r.^2 + 231.*r.^3 - 440*spones(r)))./16380725118179198;
% if dim==2
%     lrbf = @(e,r) (85179770614531827.*e.^2.*r.^8.*(11022.*r.^2 - 7700.*r - 6924.*r.^3 + 1617.*r.^4 + 1980*spones(r)))./8190362559089599;
% elseif dim==3
%     lrbf = @(e,r) (1277696559217977405.*e.^2.*r.^8.*(1540.*r.^2 - 1056.*r - 980.*r.^3 + 231.*r.^4 + 264*spones(r)))./16380725118179198;
% end

%% Type of boundary condition?
bctype=2; %1-Dirichlet, 2-Neumann, 3- Robin

%% This partially determines the number of polynomial terms
if dim<3
    fac = 0.1;
else
    fac = 4;
end

%% Type of linear algebra technique?
schur_method = 0; %1 - use schur complement, 0-use full system

%% Setup anonymous functions for PDE solution and rhs for interior
%and boundary. 
if dim==2
    syms x y;
    syms nrx nry;
    u = x.^4 + x.^3.*sin(pi*x).*cos(pi*y);     
    %u = 1 + x.^5 + y.^2 + x.^2.*y.^2; %verify polynomial reproduction
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
    %st = load('DiskPoissonNodes.mat');
    st = load('DiskPoissonNodesClustered.mat');
    %st = load('StarPoissonNodesClustered.mat');
    %st = load('DomainPoissonNodes2.mat');
elseif dim==3
    %st = load('SpherePoissonNodes.mat');
    %st = load('RBCPoissonNodesClustered.mat');
    %st = load('BumpySpherePoissonNodes.mat');
    st = load('BumpySpherePoissonNodesClustered.mat');
end
start_nodes = 1;
end_nodes = min([6,size(st.fullintnodes,1)]);

%% Start convergence study
ep = 1;
for k=start_nodes:end_nodes    
    xi = st.fullintnodes{k};
    xb = st.bdrynodes{k};
    nr = st.normals{k};
    
    tree = KDTreeSearcher(xi);
    [~,dist] = knnsearch(tree,xb,'k',1);    
    dist_scalar = mean(dist);

    %Generate ghost points    
    xg = xb + repmat(dist,[1,dim]).*nr;
    %xg = xb + 0.5*dist_scalar.*nr;
    clear tree;

    x = [xi;xb;xg];
    alph = 0; %0-legendre, -0.5- Chebyshev 1, 0.5- Chebyshev 2
    recurrence = @(N) jacobi_recurrence(N,alph,alph);    
    
    ell = floor(nthroot(fac*length(xi),dim));    
    a = total_degree_indices(dim,ell);      
    N = size(a,1);      
    if k>1
        ep = ep*2*(nthroot(length(xi)/length(st.fullintnodes{k-1}),dim));
    end


    %% Setup RBF approximation with compactly supported RBFs. We'll use Wendland's C^6 RBF.
    %%Here we use the shifted form.
    V = mpoly_eval(x,a,recurrence);
     
    rd = DistanceMatrixCSRBF(x,x,ep);
    Lap = lrbf(ep,rd); %Laplacian of RBF
    Ni = length(xi); Nb = length(xb);

    if dim==2
        Lv = mpoly_eval(x,a,recurrence,[2,0]) + mpoly_eval(x,a,recurrence,[0,2]);
    elseif dim==3
        Lv = mpoly_eval(x,a,recurrence,[2,0,0])...
             + mpoly_eval(x,a,recurrence,[0,2,0])...
             + mpoly_eval(x,a,recurrence,[0,0,2]);    
    end

    if bctype==1 %Dirichlet
        D = rbf(ep,rd(Ni+1:Ni+Nb,:)); %Get RBF
        Dv = V(Ni+1:Ni+Nb,:); %get poly
    elseif bctype==2 %Neumann

        %% Get RBF gradients    
        Drbf = drbfor(ep,rd(Ni+1:Ni+Nb,:));
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


    elseif bctype==3 %Robin
        Drbf = drbfor(ep,rd(Ni+1:Ni+Nb,:));
        A = rbf(ep,rd(Ni+1:Ni+Nb,:));

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
        Dv = Dv + V(Ni+1:Ni+Nb,:);
    end

    %% Adjust Lap for BCs
    Lii_rbf = Lap(1:Ni,1:Ni);
    Libg_rbf = Lap(1:Ni,Ni+1:end);
    Lbi_rbf = Lap(Ni+1:Ni+Nb,1:Ni);
    Lbbg_rbf = Lap(Ni+1:Ni+Nb,Ni+1:end);

    %% Get nullspace (if any)
    Vt = [Lv(1:Ni+Nb,:).', Dv.']; 
    lhsa = [[Lii_rbf,Libg_rbf,Lv(1:Ni,:)];[Lbi_rbf,Lbbg_rbf,Lv(Ni+1:Ni+Nb,:)];[D,Dv];[Vt,sparse(N,N)]];
    
    %% Form lhs based on type of solve. Also account for possible nullspace
    if bctype==2
        z = spnull(lhsa);
        lhs = [[lhsa,z];[z.',sparse(size(z,2),size(z,2))]];
    else
        z = [];
        lhs = lhsa;
    end

    %% Get rhs for linear system
    if dim==2
        rhsi = [f(xi(:,1),xi(:,2));f(xb(:,1),xb(:,2))];
        if bctype==1
            rhsb = g(xb(:,1),xb(:,2));
        elseif bctype>=2
            rhsb = g(nr(:,1),nr(:,2),xb(:,1),xb(:,2));
        end
    elseif dim==3
        rhsi = [f(xi(:,1),xi(:,2),xi(:,3));f(xb(:,1),xb(:,2),xb(:,3))];
        if bctype==1
            rhsb = g(xb(:,1),xb(:,2),xb(:,3));
        elseif bctype>=2
            rhsb = g(nr(:,1),nr(:,2),nr(:,3),xb(:,1),xb(:,2),xb(:,3));
        end    
    end

    c = lhs\[rhsi;rhsb;zeros(N,1);zeros(size(z,2),1)];
    c = c(1:end-size(z,2));
    sol = [rbf(ep,rd),V]*c; 
    sol = sol(1:Ni+Nb,:);
    sol_inner = sol(1:Ni);


    xib = [xi;xb];
    if dim==2
        true_sol = u(xib(:,1),xib(:,2));
    elseif dim==3
        true_sol = u(xib(:,1),xib(:,2),xib(:,3));
    end

    if bctype==2
        sol = sol-mean(sol);
        true_sol = true_sol - mean(true_sol);
    end

    sN(k) = nthroot(length(xi),dim);
    eps(k) = ep;
    deg(k) = ell;
    
    if schur_method==1
        nonzeros(k) = nnz(lhs_A)*100/numel(lhs_A);
        conds(k) = condest(lhs_A);        
    elseif schur_method==0
        nonzeros(k) = nnz(lhs)*100/numel(lhs);
        conds(k) = condest(lhs);
    end
    el2(k) = norm(sol-true_sol,2)./norm(true_sol,2)
    elinf(k) = norm(sol-true_sol,inf)./norm(true_sol,inf)
end
if schur_method==0
    save('Poisson_PolyReprod_Ghost.mat','el2','elinf','sN','dim','nonzeros','deg','eps');
else
    save('Poisson_PolyReprod_schur.mat','el2','elinf','sN','dim','nonzeros','deg','eps');
end
figure(1)
semilogy(deg,el2,'LineWidth',2);
hold on
semilogy(deg,el2,'blacko');
title('Convergence for Poisson equation');
xlabel('Polynomial degree');
ylabel('Relative l2 error');
axis tight;

figure(2)
semilogy(sN,el2,'LineWidth',2);
hold on
semilogy(sN,el2,'blacko');
title('Convergence for Poisson equation');
xlabel('d-th root(N)');
ylabel('Relative l2 error');
axis tight;

figure(3)
semilogy(sN,nonzeros,'LineWidth',2);
hold on
semilogy(sN,nonzeros,'blacko');
title('Percent non-zeros in system matrix');
xlabel('d-th root(N)');
ylabel('Percent non-zeros');
axis tight;

figure(4)
spy(lhs);
title('System matrix');

figure(5)
semilogy(sN,conds,'LineWidth',2);
hold on
semilogy(sN,conds,'blacko');
title('System matrix');

