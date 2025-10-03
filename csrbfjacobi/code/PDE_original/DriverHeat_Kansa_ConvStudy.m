%% Solve Heat equation on an arbitrary domain using Kansa's
%%method, but collocating with CSRBFs + jacobi polynomials.

%% Spatial dimension
dim = 2;

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
bctype=1; %1-Dirichlet, 2-Neumann, 3- Robin

%% This partially determines the number of polynomial terms
if dim<3
    fac = 0.3;
else
    fac = 1;
end

%% PDE parameters
nu = 1;

%% Setup anonymous functions for PDE solution and rhs for interior
%and boundary. 
if dim==2
    syms x y t;
    syms nrx nry;
    u = 1 + sin(pi*x).*cos(pi*y).*exp(-pi*t);     
    %u = 1 + x.^5 + y.^2 + x.^2.*y.^2; %verify polynomial reproduction
    f = diff(u,t) - nu*(diff(u,x,2) + diff(u,y,2)); 
    if bctype==1
        g = u;
    elseif bctype==2
        g = nu*(nrx.*diff(u,x) + nry.*diff(u,y)); %n.grad(u) = g;
    elseif bctype==3
        g = nu*(nrx.*diff(u,x) + nry.*diff(u,y)) + u; %n.grad(u) = g;
    end    
elseif dim==3
    syms x y z t;
    syms nrx nry nrz;
    u = 1 + sin(pi*x).*cos(pi*y).*sin(pi*z).*exp(-pi*t);     
    f = diff(u,t) - nu*(diff(u,x,2) + diff(u,y,2) + diff(u,z,2)); %lap u = f    
    if bctype==1
        g = u;
    elseif bctype==2
        g = nu*(nrx.*diff(u,x) + nry.*diff(u,y) + nrz.*diff(u,z)); %n.grad(u) = g;
    elseif bctype==3
        g = nu*(nrx.*diff(u,x) + nry.*diff(u,y) + nrz.*diff(u,z)) + u; %n.grad(u) = g;
    end      
end
u = matlabFunction(u);
f = matlabFunction(f);
g = matlabFunction(g);
if dim==2
    clear x y t nrx nry;
elseif dim==3
    clear x y z t nrx nry nrz;
end

%% Load up the node set
if dim==2
    %st = load('DiskPoissonNodes.mat');
    st = load('DiskPoissonNodesClustered.mat');
    %st = load('StarPoissonNodesClustered.mat');
    %st = load('DomainPoissonNodes2.mat');
elseif dim==3
    st = load('SpherePoissonNodes.mat');
end
start_nodes = 1;
end_nodes = min([6,size(st.fullintnodes,1)]);

%% Final time
T = 0.2; 


%% Start convergence study
ep = 2;
for k=start_nodes:end_nodes    
    xi = st.fullintnodes{k};
    xb = st.bdrynodes{k};
    nr = st.normals{k};
    x = [xi;xb];
    alph = 0; %0-legendre, -0.5- Chebyshev 1, 0.5- Chebyshev 
    recurrence = @(N) jacobi_recurrence(N,alph,alph);    
    
    ell = floor(nthroot(fac*length(xi),dim));
    a = total_degree_indices(dim,ell);  
    %a = hyperbolic_cross_indices(dim,ell);
    N = size(a,1);      
    if k>1
        ep = ep*2*nthroot(length(xi)/length(st.fullintnodes{k-1}),dim);  
    end

    %% Set dt based on h
    error_estimate = (1./nthroot(N,dim)).^(ell+1);
    dt = sqrt(error_estimate); %BDF2
    if dt<1e-5
        dt = 1e-5;
    end
    
    %% Now fix dt so we have integer number of steps
    nsteps = floor(T/dt);
    dt = T/nsteps;

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
        D = rbf(ep,rd(Ni+1:end,:)); %Get RBF
        Dv = V(Ni+1:end,:); %get poly
    elseif bctype==2 %Neumann

        %% Get RBF gradients    
        Drbf = drbfor(ep,rd(Ni+1:end,:));
        if size(x,2)>=1
           [xj,xk] = ndgrid(xb(:,1),x(:,1));
           nr1 = repmat(nr(:,1),[1,size(xj,2)]);
           D = nu*nr1.*(xj-xk).*Drbf;
        end
        if size(x,2)>=2
           [yj,yk] = ndgrid(xb(:,2),x(:,2));
           nr2 = repmat(nr(:,2),[1,size(xj,2)]);
           D = D + nu*nr2.*(yj-yk).*Drbf;
        end
        if size(x,3)==3
           [zj,zk] = ndgrid(xb(:,3),x(:,3));
           nr3 = repmat(nr(:,3),[1,size(xj,2)]);
           D = D + nu*nr(:,3).*(zj-zk).*Drbf; 
        end

        %% Get polynomial gradients
        if dim==1
            Dvx =  mpoly_eval(xb,a,recurrence,1);
            nr1 = repmat(nr(:,1),[1,size(Dvx,2)]);
            Dv = nu*nr1.*Dvx;
        elseif dim==2
            Dvx =  mpoly_eval(xb,a,recurrence,[1,0]);
            Dvy =  mpoly_eval(xb,a,recurrence,[0,1]);
            nr1 = repmat(nr(:,1),[1,size(Dvx,2)]);
            nr2 = repmat(nr(:,2),[1,size(Dvy,2)]);
            Dv = nu*nr1.*Dvx + nu*nr2.*Dvy;
        elseif dim==3
            Dvx =  mpoly_eval(xb,a,recurrence,[1,0,0]);
            Dvy =  mpoly_eval(xb,a,recurrence,[0,1,0]);
            Dvz =  mpoly_eval(xb,a,recurrence,[0,0,1]);
            nr1 = repmat(nr(:,1),[1,size(Dvx,2)]);
            nr2 = repmat(nr(:,2),[1,size(Dvy,2)]);
            nr3 = repmat(nr(:,3),[1,size(Dvz,2)]);
            Dv = nu*(nr1.*Dvx + nr2.*Dvy + nr3.*Dvz);
        end
    elseif bctype==3 %Robin
        Drbf = drbfor(ep,rd(Ni+1:end,:));
        A = rbf(ep,rd(Ni+1:end,:));

        %% Get gradients  
        if size(x,2)>=1
           [xj,xk] = ndgrid(xb(:,1),x(:,1));
           nr1 = repmat(nr(:,1),[1,size(xj,2)]);
           D = nu*nr1.*(xj-xk).*Drbf;
        end
        if size(x,2)>=2
           [yj,yk] = ndgrid(xb(:,2),x(:,2));
           nr2 = repmat(nr(:,2),[1,size(xj,2)]);
           D = D + nu*nr2.*(yj-yk).*Drbf;
        end
        if size(x,3)==3
           [zj,zk] = ndgrid(xb(:,3),x(:,3));
           nr3 = repmat(nr(:,3),[1,size(xj,2)]);
           D = D + nu*nr(:,3).*(zj-zk).*Drbf; 
        end  

        %% Get polynomial gradients
        if dim==1
            Dvx = mpoly_eval(xb,a,recurrence,1);
            nr1 = repmat(nr(:,1),[1,size(Dvx,2)]);
            Dv = nu*nr1.*Dvx;
        elseif dim==2
            Dvx =  mpoly_eval(xb,a,recurrence,[1,0]);
            Dvy =  mpoly_eval(xb,a,recurrence,[0,1]);
            nr1 = repmat(nr(:,1),[1,size(Dvx,2)]);
            nr2 = repmat(nr(:,2),[1,size(Dvy,2)]);
            Dv = nu*(nr1.*Dvx + nr2.*Dvy);
        elseif dim==3
            Dvx =  mpoly_eval(xb,a,recurrence,[1,0,0]);
            Dvy =  mpoly_eval(xb,a,recurrence,[0,1,0]);
            Dvz =  mpoly_eval(xb,a,recurrence,[0,0,1]);
            nr1 = repmat(nr(:,1),[1,size(Dvx,2)]);
            nr2 = repmat(nr(:,2),[1,size(Dvy,2)]);
            nr3 = repmat(nr(:,3),[1,size(Dvz,2)]);
            Dv = nu*(nr1.*Dvx + nr2.*Dvy + nr3.*Dvz);
        end    

        %% Form Robin matrix for RBFs
        D = D + A;

        %% form robin matrix for polys
        Dv = Dv + V(Ni+1:end,:);
    end
    
    %% Get the RBF matrix
    A = rbf(ep,rd); %RBF
    
    %% Get polynomial conditions
    if bctype==1 || bctype==3 %If it's Dirichlet or Robin, it should be fine
        Vt = [Lv(1:Ni,:).', Dv.'];        
    else %If it's Neumann, we'll enforce a Robin-style reproduction
        Vt = [Lv(1:Ni,:).', (Dv+V(Ni+1:end,:)).'];        
    end
    
    %Vt = [V(1:Ni,:).',Dv.'];
    %rt = rank(Vt);
    
    %% Get matrix blocks
    Lif_rbf = Lap(1:Ni,:);    
    Aif_rbf = A(1:Ni,:);    
    Lif_poly = Lv(1:Ni,:);
    Vif_poly = V(1:Ni,:);
   
    %% Get BDF2 matrix blocks
    lhs_11 = Aif_rbf - (2/3)*nu*dt*Lif_rbf;
    lhs_12 = Vif_poly - (2/3)*nu*dt*Lif_poly;        
    lhsBDF2 = [[lhs_11,lhs_12];[D,Dv];[Vt,sparse(N,N)]];
    z = spnull(lhsBDF2);
    lhsBDF2 = [[lhsBDF2,z];[z.',zeros(size(z,2),size(z,2))]];
    [L,U,P,Q] = lu(lhsBDF2);
    
    %% Get initial condition and first two steps
    if dim==2
        u0f = u(0,x(:,1),x(:,2));
        u1f = u(dt,x(:,1),x(:,2));
    elseif dim==3
        u0f = u(0,x(:,1),x(:,2),x(:,3));
        u1f = u(dt,x(:,1),x(:,2),x(:,3));
    end
    Aeval = [rbf(ep,rd),V];
    

    %% Now do BDF2
    for nt=2:nsteps
        %%Get BCs
        if dim==2        
            if bctype==1
                rhsb = g(nt*dt,xb(:,1),xb(:,2));
            elseif bctype>=2
                rhsb = g(nr(:,1),nr(:,2),nt*dt,xb(:,1),xb(:,2));
            end
        elseif dim==3        
            if bctype==1
                rhsb = g(nt*dt,xb(:,1),xb(:,2),xb(:,3));
            elseif bctype>=2
                rhsb = g(nr(:,1),nr(:,2),nr(:,3),nt*dt,xb(:,1),xb(:,2),xb(:,3));
            end    
        end    

        %%Assemble rhs
        rhs = [((4/3)*u1f(1:Ni) - (1/3)*u0f(1:Ni,:)+ (2/3)*dt*f(nt*dt,xi(:,1),xi(:,2)));rhsb;zeros(N,1);zeros(size(z,2),1)];
        c = Q*(U\(L\(P*rhs))); ch = c(1:end-size(z,2),:);
        sol = Aeval*ch;
        if norm(sol,inf)> 4
            error('Blew up');            
        end
        
        u0f = u1f; u1f = sol;
    end

    if dim==2
        true_sol = u(T,x(:,1),x(:,2));
    elseif dim==3
        true_sol = u(T,x(:,1),x(:,2),x(:,3));
    end

%     if bctype==2
%         sol = sol-mean(sol);
%         true_sol = true_sol - mean(true_sol);
%     end

    sN(k) = nthroot(length(xi),dim);
    eps(k) = ep;
    deg(k) = ell;
    nonzeros(k) = nnz(lhsBDF2)*100/size(lhsBDF2,2).^2;
    el2(k) = norm(sol-true_sol,2)./norm(true_sol,2)
    elinf(k) = norm(sol-true_sol,inf)./norm(true_sol,inf)
end
save('Heat_PolyDerReprod_schur.mat','el2','elinf','sN','dim','nonzeros','deg','eps');
figure(1)
semilogy(deg,el2,'LineWidth',2);
hold on
semilogy(deg,el2,'blacko');
title('Convergence for Heat equation');
xlabel('Polynomial degree');
ylabel('Relative l2 error');
axis tight;

figure(2)
semilogy(sN,el2,'LineWidth',2);
hold on
semilogy(sN,el2,'blacko');
title('Convergence for Heat equation');
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


