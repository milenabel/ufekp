%% Interpolation in 1D/2D/3D domains using CSRBFs combined with
%%Jacobi polynomials. 


%% Spatial dimension
dim = 2;

%% Load up the node set
if dim==1
    N = 2.^(2:8); N = N';
    for k=1:length(N)
        X = chebspace2(-1,1,N(k));
        xi = X(2:end-1,:);
        xb = [X(1,:); X(end,:)];
        st.fullintnodes{k,1} = xi;
        st.bdrynodes{k,1} = xb;
    end
    clear N xb xi k;
    xe = linspace(-1,1,2^14).';    
elseif dim==2
    %% Get evaluation nodes
    st = load('DiskPoissonNodesLarge.mat');
    xe =  [st.fullintnodes{7}; st.bdrynodes{7}];

    %st = load('DiskPoissonNodes.mat');
    st = load('DiskPoissonNodesClustered.mat');
    %st = load('StarPoissonNodesClustered.mat');
    %st = load('DomainPoissonNodes2.mat');
%     pxe = haltonset(dim,'Skip',1e3,'Leap',1e2); %Quasi-random node set
%     pxe = scramble(pxe,'RR2');
%     N = (4:5:60).^2; N = N';
%     for k=1:length(N)
%         p = net(pxe,N(k));
%         p = p*2 - 1;
%         st.fullintnodes{k,1} = p;    
%     end
%     xe = net(pxe,10000);
%     xe = xe*2 - 1;
    clear p N pxe;
elseif dim==3
     st = load('SpherePoissonNodesLarge.mat');
     xe = [st.fullintnodes{2}; st.bdrynodes{2}];
% 
     st = load('SpherePoissonNodesClustered.mat');

%    st = load('RBCPoissonNodesClustered.mat');
%    xe = [st.fullintnodes{7}; st.bdrynodes{7}];
%    st = load('BumpySpherePoissonNodes.mat');
%    st = load('BumpySpherePoissonNodesClustered.mat');
%     pxe = haltonset(dim,'Skip',1e3,'Leap',1e2); %Quasi-random node set
%     pxe = scramble(pxe,'RR2');
%     N = (4:5:20).^3; N = N';
%     for k=1:length(N)
%         p = net(pxe,N(k));
%         p = p*2 - 1;
%         st.fullintnodes{k,1} = p;    
%     end
%     xe = net(pxe,10000);
%     xe = xe*2 - 1;
%     clear p N pxe;
end
start_nodes = 1;
end_nodes = size(st.fullintnodes,1);

%% This partially determines the number of polynomial terms
%slightly smaller fac helps with increasing accuracy for both rough and
%smooth functions. There seems to be a notion of the "best" polynomial
%degree for a given number of nodes
if dim==1
    fac = 1;
elseif dim==2
    fac = 0.8; 
elseif dim==3
    fac = 1.0;
end

%% Setup anonymous function for function interpolation true solution
%%We'll use Trefethen's function gallery from his Spectral methods book
%%with 2d and 3d analogues
if dim==1
    syms x;       
    %f = abs(x);
    %f = exp(-x.^(-2));    
    f = 1./(1 + 16*x.^2);
    %f = x.^(10);
    dfx = diff(f,x);
elseif dim==2
    syms x y;    
    %f = abs(x).^3.*abs(y).^3;
    %f = exp(-x.^(-2)).*exp(-y.^(-2));    
    %f = 1./(1 + 25*(x.^2 + y.^2));
    f = exp(-10*((x-.3).^(-2)+y.^(-2)));    
    %f = exp(-10*((x-.3).^2+y.^2));
    %f = x.^(8).*y.^(8);
    dfx = diff(f,x); dfy = diff(f,y);
elseif dim==3
    syms x y z;      
    %f = abs(x).^3.*abs(y).^3.*abs(z).^3;
    %f = exp(-10*((x-.3).^(-2)+y.^(-2) + z.^(-2)));  
    f = exp(-10*((x-.3).^2+y.^2 + z.^2));  
    %f = 1./(1 + 16*(x.^2 + y.^2 + z.^2));
    %f = x.^(4).*y.^(2).*z.^(2);    
    dfx = diff(f,x); dfy = diff(f,y); dfz = diff(f,z);
end
f = matlabFunction(f);
if dim==1
    dfx = matlabFunction(dfx);
    clear x;
elseif dim==2
    dfx = matlabFunction(dfx);
    dfy = matlabFunction(dfy);
    clear x y;
elseif dim==3
    dfx = matlabFunction(dfx);
    dfy = matlabFunction(dfy);
    dfz = matlabFunction(dfz);
    clear x y z;
end
%% All possible tests:
%% 1. Different smoothness for CSRBF
%% 2. Different shape parameter strats
%% 3. Different interp techniques

%% Get the standard polynomial least squares stuff out of the way
for k=start_nodes:end_nodes
    xi = st.fullintnodes{k};
    xb = st.bdrynodes{k};
    x = [xi;xb];
    alph = 0; %legendre
    ell = floor(fac*nthroot(length(x),dim));            
    ell = max([ell,1]); 
    if dim==1
        y = f(x(:,1));
        ye_true = f(xe(:,1));
    elseif dim==2
        y = f(x(:,1),x(:,2));
        ye_true = f(xe(:,1),xe(:,2));
    elseif dim==3
        y = f(x(:,1),x(:,2),x(:,3));
        ye_true = f(xe(:,1),xe(:,2),xe(:,3));
    end    
    [el2_poly(k,1),elinf_poly(k,1),a_time_poly(k,1),e_time_poly(k,1),c_poly{k,1}] = PLS(x,y,ell,xe,alph,ye_true);   
    sN(k,1) = nthroot(length(x),dim);
end


%% Next, get the diagonal approximation stuff out of the way for different smoothnesses
for smoothness=1:3
    if smoothness==1
        %% Wendland C2 in 3d, pd in all lower dimensions
        rbf = @(e,r) -r.^4.*(4.*r - 5*spones(r));
        drbfor = @(e,r) 20.*e.^2.*r.^3;
        if dim==1
            lrbf = @(e,r) -20.*e.^2.*r.^2.*(2.*r - 3*spones(r));
        elseif dim==2
            lrbf = @(e,r) -20.*e.^2.*r.^2.*(r - 3*spones(r));
        elseif dim==3
            lrbf = @(e,r) 60.*e.^2.*r.^2;
        end
    elseif smoothness==2
        %% Wendland C4 in 3d
        rbf = @(e,r) (r.^6.*(35.*r.^2 - 88.*r + 56*spones(r)))./3;
        drbfor = @(e,r) -(56.*e.^2.*r.^5.*(5.*r - 6*spones(r)))./3;
        if dim==2
            lrbf = @(e,r) (112.*e.^2.*r.^4.*(20.*r.^2 - 36.*r + 15*spones(r)))./3;   
        elseif dim==3    
            lrbf = @(e,r) 56.*e.^2.*r.^4.*(15.*r.^2 - 26.*r + 10*spones(r));
        end
    elseif smoothness==3
        %% Wendland C6 in 3d
        rbf = @(e,r) r.^8.*(66*spones(r) - 154*r + 121*r.^2 - 32*r.^3);
        drbfor = @(e,r) 22.*e.^2.*r.^7.*(16.*r.^2 - 39.*r + 24*spones(r));
        if dim==2
            lrbf = @(e,r) 44*e.^2.*r.^6.*(84*spones(r)-264*r+267*r.^2-88*r.^3);
        elseif dim==3
            lrbf = @(e,r) -66.*e.^2.*r.^6.*((r - spones(r)).^2 - 6.*r + 64.*(r - spones(r)).^3 + 7*spones(r));
        end
    else
        %% Wendland C8 in 3d
        rbf = @(e,r) (6552290047271679.*r.^10.*(4134.*r.^2 - 3536.*r - 2166.*r.^3 + 429.*r.^4 + 1144*spones(r)))./32761450236358396;
        drbfor = @(e,r) -(85179770614531827.*e.^2.*r.^9.*(1056.*r - 852.*r.^2 + 231.*r.^3 - 440*spones(r)))./16380725118179198;
        if dim==2
            lrbf = @(e,r) (85179770614531827.*e.^2.*r.^8.*(11022.*r.^2 - 7700.*r - 6924.*r.^3 + 1617.*r.^4 + 1980*spones(r)))./8190362559089599;
        elseif dim==3
            lrbf = @(e,r) (1277696559217977405.*e.^2.*r.^8.*(1540.*r.^2 - 1056.*r - 980.*r.^3 + 231.*r.^4 + 264*spones(r)))./16380725118179198;
        end
    end

    for k=start_nodes:end_nodes
        xi = st.fullintnodes{k};
        xb = st.bdrynodes{k};
        x = [xi;xb];
        alph = 0; %legendre
        ell = floor(fac*nthroot(length(x),dim));            
        ell = max([ell,1]); 
        if dim==1
            y = f(x(:,1));
            ye_true = f(xe(:,1));
        elseif dim==2
            y = f(x(:,1),x(:,2));
            ye_true = f(xe(:,1),xe(:,2));
        elseif dim==3
            y = f(x(:,1),x(:,2),x(:,3));
            ye_true = f(xe(:,1),xe(:,2),xe(:,3));
        end    
        tree = KDTreeSearcher(x);
        [el2_diag(k,smoothness),elinf_diag(k,smoothness),a_time_diag(k,smoothness),e_time_diag(k,smoothness),c_poly_diag{k,smoothness}] = CSRBFDiag(x,y,ell,xe,alph,rbf,tree,ye_true);    
    end
end

%% Now do a couple of different fixed shape parameter strategies: a large one and a small one.
%% Again, different smoothnesses
for smoothness=1:3
    if smoothness==1
        %% Wendland C2 in 3d, pd in all lower dimensions
        rbf = @(e,r) -r.^4.*(4.*r - 5*spones(r));
        drbfor = @(e,r) 20.*e.^2.*r.^3;
        if dim==1
            lrbf = @(e,r) -20.*e.^2.*r.^2.*(2.*r - 3*spones(r));
        elseif dim==2
            lrbf = @(e,r) -20.*e.^2.*r.^2.*(r - 3*spones(r));
        elseif dim==3
            lrbf = @(e,r) 60.*e.^2.*r.^2;
        end
    elseif smoothness==2
        %% Wendland C4 in 3d
        rbf = @(e,r) (r.^6.*(35.*r.^2 - 88.*r + 56*spones(r)))./3;
        drbfor = @(e,r) -(56.*e.^2.*r.^5.*(5.*r - 6*spones(r)))./3;
        if dim==2
            lrbf = @(e,r) (112.*e.^2.*r.^4.*(20.*r.^2 - 36.*r + 15*spones(r)))./3;   
        elseif dim==3    
            lrbf = @(e,r) 56.*e.^2.*r.^4.*(15.*r.^2 - 26.*r + 10*spones(r));
        end
    elseif smoothness==3
        %% Wendland C6 in 3d
        rbf = @(e,r) r.^8.*(66*spones(r) - 154*r + 121*r.^2 - 32*r.^3);
        drbfor = @(e,r) 22.*e.^2.*r.^7.*(16.*r.^2 - 39.*r + 24*spones(r));
        if dim==2
            lrbf = @(e,r) 44*e.^2.*r.^6.*(84*spones(r)-264*r+267*r.^2-88*r.^3);
        elseif dim==3
            lrbf = @(e,r) -66.*e.^2.*r.^6.*((r - spones(r)).^2 - 6.*r + 64.*(r - spones(r)).^3 + 7*spones(r));
        end
    else
        %% Wendland C8 in 3d
        rbf = @(e,r) (6552290047271679.*r.^10.*(4134.*r.^2 - 3536.*r - 2166.*r.^3 + 429.*r.^4 + 1144*spones(r)))./32761450236358396;
        drbfor = @(e,r) -(85179770614531827.*e.^2.*r.^9.*(1056.*r - 852.*r.^2 + 231.*r.^3 - 440*spones(r)))./16380725118179198;
        if dim==2
            lrbf = @(e,r) (85179770614531827.*e.^2.*r.^8.*(11022.*r.^2 - 7700.*r - 6924.*r.^3 + 1617.*r.^4 + 1980*spones(r)))./8190362559089599;
        elseif dim==3
            lrbf = @(e,r) (1277696559217977405.*e.^2.*r.^8.*(1540.*r.^2 - 1056.*r - 980.*r.^3 + 231.*r.^4 + 264*spones(r)))./16380725118179198;
        end
    end

    for k=start_nodes:end_nodes
        xi = st.fullintnodes{k};
        xb = st.bdrynodes{k};
        x = [xi;xb];
        alph = 0; %legendre
        ell = floor(fac*nthroot(length(x),dim));            
        ell = max([ell,1]); 
        if dim==1
            y = f(x(:,1));
            ye_true = f(xe(:,1));
        elseif dim==2
            y = f(x(:,1),x(:,2));
            ye_true = f(xe(:,1),xe(:,2));
        elseif dim==3
            y = f(x(:,1),x(:,2),x(:,3));
            ye_true = f(xe(:,1),xe(:,2),xe(:,3));
        end    
        tree = KDTreeSearcher(x);
        [~,dist] = knnsearch(tree,x,'k',2);
        dist = dist(:,2);
        dist = min(dist); %separation radius       
        if k==1 %fix under refinement
            ep1 = 0.9/dist; %support = 1/ep
        end
        [el2_fs1(k,smoothness),elinf_fs1(k,smoothness),a_time_fs1(k,smoothness),e_time_fs1(k,smoothness),c_poly_fs1{k,smoothness}] = CSRBFGen(x,y,ell,xe,alph,rbf,ep1,tree,ye_true);    
        
        if k==1 %fix under refinement
            ep2 = 0.1/dist; %support = 1/ep
        end
        [el2_fs2(k,smoothness),elinf_fs2(k,smoothness),a_time_fs2(k,smoothness),e_time_fs2(k,smoothness),c_poly_fs2{k,smoothness}] = CSRBFGen(x,y,ell,xe,alph,rbf,ep2,tree,ye_true);            
        
        
        if k==1 %fix under refinement
            ep3 = 0.05/dist; %support = 1/ep
        end
        [el2_fs3(k,smoothness),elinf_fs3(k,smoothness),a_time_fs3(k,smoothness),e_time_fs3(k,smoothness),c_poly_fs3{k,smoothness}] = CSRBFGen(x,y,ell,xe,alph,rbf,ep3,tree,ye_true);                  
    end
end

%% Now get the variable shape parameter (aka fixed spatial support) strategies out of the way.
%% Again, different smoothnesses
for smoothness=1:3
    if smoothness==1
        %% Wendland C2 in 3d, pd in all lower dimensions
        rbf = @(e,r) -r.^4.*(4.*r - 5*spones(r));
        drbfor = @(e,r) 20.*e.^2.*r.^3;
        if dim==1
            lrbf = @(e,r) -20.*e.^2.*r.^2.*(2.*r - 3*spones(r));
        elseif dim==2
            lrbf = @(e,r) -20.*e.^2.*r.^2.*(r - 3*spones(r));
        elseif dim==3
            lrbf = @(e,r) 60.*e.^2.*r.^2;
        end
    elseif smoothness==2
        %% Wendland C4 in 3d
        rbf = @(e,r) (r.^6.*(35.*r.^2 - 88.*r + 56*spones(r)))./3;
        drbfor = @(e,r) -(56.*e.^2.*r.^5.*(5.*r - 6*spones(r)))./3;
        if dim==2
            lrbf = @(e,r) (112.*e.^2.*r.^4.*(20.*r.^2 - 36.*r + 15*spones(r)))./3;   
        elseif dim==3    
            lrbf = @(e,r) 56.*e.^2.*r.^4.*(15.*r.^2 - 26.*r + 10*spones(r));
        end
    elseif smoothness==3
        %% Wendland C6 in 3d
        rbf = @(e,r) r.^8.*(66*spones(r) - 154*r + 121*r.^2 - 32*r.^3);
        drbfor = @(e,r) 22.*e.^2.*r.^7.*(16.*r.^2 - 39.*r + 24*spones(r));
        if dim==2
            lrbf = @(e,r) 44*e.^2.*r.^6.*(84*spones(r)-264*r+267*r.^2-88*r.^3);
        elseif dim==3
            lrbf = @(e,r) -66.*e.^2.*r.^6.*((r - spones(r)).^2 - 6.*r + 64.*(r - spones(r)).^3 + 7*spones(r));
        end
    else
        %% Wendland C8 in 3d
        rbf = @(e,r) (6552290047271679.*r.^10.*(4134.*r.^2 - 3536.*r - 2166.*r.^3 + 429.*r.^4 + 1144*spones(r)))./32761450236358396;
        drbfor = @(e,r) -(85179770614531827.*e.^2.*r.^9.*(1056.*r - 852.*r.^2 + 231.*r.^3 - 440*spones(r)))./16380725118179198;
        if dim==2
            lrbf = @(e,r) (85179770614531827.*e.^2.*r.^8.*(11022.*r.^2 - 7700.*r - 6924.*r.^3 + 1617.*r.^4 + 1980*spones(r)))./8190362559089599;
        elseif dim==3
            lrbf = @(e,r) (1277696559217977405.*e.^2.*r.^8.*(1540.*r.^2 - 1056.*r - 980.*r.^3 + 231.*r.^4 + 264*spones(r)))./16380725118179198;
        end
    end

    for k=start_nodes:end_nodes
        xi = st.fullintnodes{k};
        xb = st.bdrynodes{k};
        x = [xi;xb];
        alph = 0; %legendre
        ell = floor(fac*nthroot(length(x),dim));            
        ell = max([ell,1]); 
        if dim==1
            y = f(x(:,1));
            ye_true = f(xe(:,1));
        elseif dim==2
            y = f(x(:,1),x(:,2));
            ye_true = f(xe(:,1),xe(:,2));
        elseif dim==3
            y = f(x(:,1),x(:,2),x(:,3));
            ye_true = f(xe(:,1),xe(:,2),xe(:,3));
        end    
        tree = KDTreeSearcher(x);
        [~,dist] = knnsearch(tree,x,'k',2);
        dist = dist(:,2);
        dist = min(dist); %separation radius               
        ep1 = 0.9/dist; %support = 1/ep        
        [el2_vs1(k,smoothness),elinf_vs1(k,smoothness),a_time_vs1(k,smoothness),e_time_vs1(k,smoothness),c_poly_vs1{k,smoothness}] = CSRBFGen(x,y,ell,xe,alph,rbf,ep1,tree,ye_true);    
        
        ep2 = 0.1/dist; %support = 1/ep        
        [el2_vs2(k,smoothness),elinf_vs2(k,smoothness),a_time_vs2(k,smoothness),e_time_vs2(k,smoothness),c_poly_vs2{k,smoothness}] = CSRBFGen(x,y,ell,xe,alph,rbf,ep2,tree,ye_true);         
        
        ep2 = 0.05/dist; %support = 1/ep        
        [el2_vs3(k,smoothness),elinf_vs3(k,smoothness),a_time_vs3(k,smoothness),e_time_vs3(k,smoothness),c_poly_vs3{k,smoothness}] = CSRBFGen(x,y,ell,xe,alph,rbf,ep3,tree,ye_true);           
    end
end

%% Plot everything