%% Interpolation in 1D/2D/3D domains using CSRBFs combined with
%%Jacobi polynomials. 


%% Spatial dimension
dim = 1;

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

% Sparsity storage initialization 
sparsity_fs1 = zeros(end_nodes, 3);
sparsity_fs2 = zeros(end_nodes, 3);
sparsity_fs3 = zeros(end_nodes, 3);
sparsity_vs1 = zeros(end_nodes, 3);
sparsity_vs2 = zeros(end_nodes, 3);
sparsity_vs3 = zeros(end_nodes, 3);

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
    f = abs(x);                function_name = 'abs_1d';
    %f_sym = exp(-x.^(-2));     function_name = 'gauss_1d';
    %f_sym = 1./(1 + 16*x.^2);  function_name = 'rational_1d';
    %f_sym = x.^(10);           function_name = 'poly10_1d';
    dfx = diff(f,x);
elseif dim==2
    syms x y;    
    f = abs(x).^3 .* abs(y).^3;              function_name = 'abs3x3_2d';
    %f = exp(-x.^(-2)).*exp(-y.^(-2));        function_name = 'gauss_2d';
    %f = 1./(1 + 25*(x.^2 + y.^2));           function_name = 'rational_2d';
    %%f = exp(-10*((x-.3).^(-2)+y.^(-2)));    function_name = 'gauss10_2d';
    %f = exp(-10*((x-.3).^2+y.^2));           function_name = 'gauss10inv_2d';
    %f = x.^8 .* y.^8;                        function_name = 'poly8x8_2d';
    dfx = diff(f,x); dfy = diff(f,y);
elseif dim==3
    syms x y z;      
    f = abs(x).^3.*abs(y).^3.*abs(z).^3;              function_name = 'abs333_3d';
    %f = exp(-10*((x-.3).^(-2)+y.^(-2) + z.^(-2)));    function_name = 'gaussinv_3d'
    %f = exp(-10*((x-.3).^2+y.^2 + z.^2));             function_name = 'gauss_3d';
    %f = 1./(1 + 16*(x.^2 + y.^2 + z.^2));             function_name = 'rational_3d';
    %f = x.^(4).*y.^(2).*z.^(2);                       function_name = 'poly422_3d';
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
        % dist = 0.5*min(dist); %separation radius
        dist = 0.25*median(dist);
        % max_d = median(dist);
        % frac = [0.25, 0.50, 0.75];   % 25%, 50%, 75%
        % r = frac*max_d;
        % eps = 1 ./ r;
        % ep1 = eps(1); %support = 1/ep
        % ep3 = eps(3); %support = 1/ep
        % ep2 = eps(2); %support = 1/ep
        % 
        % [el2_fs1(k,smoothness),elinf_fs1(k,smoothness),a_time_fs1(k,smoothness),e_time_fs1(k,smoothness),c_poly_fs1{k,smoothness}, cond_fs1(k,smoothness), ~, sparsity_fs1(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep1,tree,ye_true);       
        % [el2_fs2(k,smoothness),elinf_fs2(k,smoothness),a_time_fs2(k,smoothness),e_time_fs2(k,smoothness),c_poly_fs2{k,smoothness}, cond_fs2(k,smoothness), ~, sparsity_fs2(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep2,tree,ye_true); 
        % [el2_fs3(k,smoothness),elinf_fs3(k,smoothness),a_time_fs3(k,smoothness),e_time_fs3(k,smoothness),c_poly_fs3{k,smoothness}, cond_fs3(k,smoothness), ~, sparsity_fs3(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep3,tree,ye_true);

        if k==1 %fix under refinement
            ep1 = 2.5/dist; %support = 1/ep
        end
        [el2_fs1(k,smoothness),elinf_fs1(k,smoothness),a_time_fs1(k,smoothness),e_time_fs1(k,smoothness),c_poly_fs1{k,smoothness}, cond_fs1(k,smoothness), ~, sparsity_fs1(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep1,tree,ye_true);    

        if k==1 %fix under refinement
            ep2 = 1.75/dist; %support = 1/ep
        end
        [el2_fs2(k,smoothness),elinf_fs2(k,smoothness),a_time_fs2(k,smoothness),e_time_fs2(k,smoothness),c_poly_fs2{k,smoothness}, cond_fs2(k,smoothness), ~, sparsity_fs2(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep2,tree,ye_true); 


        if k==1 %fix under refinement
            ep3 = 1.25/dist; %support = 1/ep
        end
        [el2_fs3(k,smoothness),elinf_fs3(k,smoothness),a_time_fs3(k,smoothness),e_time_fs3(k,smoothness),c_poly_fs3{k,smoothness}, cond_fs3(k,smoothness), ~, sparsity_fs3(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep3,tree,ye_true);
    end
end

%% Now get the variable shape parameter (aka fixed level of sparsity) strategies out of the way.
%% Again, different smoothnesses
%s_targets = [0.1, 0.95, 0.99];    % three levels of sparsity to try
 s_targets = [0.40, 0.55, 0.70];    % 1 - s_target(i) is target sparsity value
tolS = 1e-5;                        % your tight tolerance
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

        Nnodes = size(x,1);

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
        % tree = KDTreeSearcher(x);
        % [~,dist] = knnsearch(tree,x,'k',2);
        % dist = dist(:,2);
        % dist = 0.5*min(dist); %separation radius (smallest neighbor spacing (q))
        % diam = norm(max(x)-min(x));  % largest inter‐point span (w)
        % % diam = dist;
        % 
        % tolS = 1e-5;   % tolerance on sparsity
        % ep1 = findEpForSparsity(x, tree, rbf, s_targets(1), tolS, diam);
        % [el2_vs1(k,smoothness),elinf_vs1(k,smoothness),a_time_vs1(k,smoothness),e_time_vs1(k,smoothness),c_poly_vs1{k,smoothness}, ~, sparsity_vs1(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep1,tree,ye_true);
        % 
        % ep2 = findEpForSparsity(x, tree, rbf, s_targets(2), tolS, diam);
        % [el2_vs2(k,smoothness),elinf_vs2(k,smoothness),a_time_vs2(k,smoothness),e_time_vs2(k,smoothness),c_poly_vs2{k,smoothness}, ~, sparsity_vs2(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep2,tree,ye_true);
        % 
        % ep3 = findEpForSparsity(x, tree, rbf, s_targets(3), tolS, diam);       
        % [el2_vs3(k,smoothness),elinf_vs3(k,smoothness),a_time_vs3(k,smoothness),e_time_vs3(k,smoothness),c_poly_vs3{k,smoothness}, ~, sparsity_vs3(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep3,tree,ye_true); 

        % assume: x, y, ell, xe, alph, rbf, ye_true, s_targets = [0.25,0.50,0.75]
        tree   = KDTreeSearcher(x);
        
        % precompute “natural” bracket in ε = 1/radius space
        [~,D2] = knnsearch(tree,x,'k',2); % N×2, second col = nearest neighbor dist
        r = min(D2(:,2)); % separation radius
        diam = max(pdist(x)); % diameter
        ep_min = 1/diam;
        ep_max = 1/r;
        
        eps_vs = zeros(1,3);
        for j = 1:3
            s_t = s_targets(j);

            ep_lo = ep_min;      % reset for each j
            ep_hi = ep_max;
            
            % objective: sparsity(ε)−target = 0
            fun = @(ep)(sum(cellfun(@numel, rangesearch(tree,x,1/ep)))/(Nnodes^2)) - (1-s_t);
            
            % bracket locally
            f_lo = fun(ep_lo);
            f_hi = fun(ep_hi);
            tries = 0;
            while f_lo * f_hi > 0 && tries < 5 % tried 10, 50, 100, 200;  none worked without error
                ep_lo = ep_lo/2;
                ep_hi = ep_hi*2;
                f_lo = fun(ep_lo);
                f_hi = fun(ep_hi);
                tries = tries + 1;
            end
            
            if f_lo * f_hi <= 0
                % fzero
                eps_vs(j) = fzero(fun, [ep_lo, ep_hi], optimset('TolX',1e-4));
            else
                warning('Couldn’t bracket $\epsilon$ for sparsity=%.2f; falling back to nearest endpoint', s_t, "latex");
                % pick whichever endpoint gives smaller |fun|
                if abs(f_lo) < abs(f_hi)
                    eps_vs(j) = ep_lo;
                else
                    eps_vs(j) = ep_hi;
                end
            end
            % clamp to a PD epsilon
            eps_vs(j) = clampToPD(eps_vs(j), x, tree, rbf);
        end

        ep1 = eps_vs(1);
        ep2 = eps_vs(2);
        ep3 = eps_vs(3);
        
        [el2_vs1(k,smoothness), elinf_vs1(k,smoothness), a_time_vs1(k,smoothness), e_time_vs1(k,smoothness), c_poly_vs1{k,smoothness}, cond_vs1(k,smoothness), ~, sparsity_vs1(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep1,tree,ye_true);
        
        [el2_vs2(k,smoothness), elinf_vs2(k,smoothness), a_time_vs2(k,smoothness),  e_time_vs2(k,smoothness), c_poly_vs2{k,smoothness}, cond_vs2(k,smoothness), ~, sparsity_vs2(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep2,tree,ye_true);
        
        [el2_vs3(k,smoothness), elinf_vs3(k,smoothness), a_time_vs3(k,smoothness), e_time_vs3(k,smoothness), c_poly_vs3{k,smoothness}, cond_vs3(k,smoothness), ~, sparsity_vs3(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep3,tree,ye_true);
    end
end

function ep = clampToPD(ep_guess, x,tree,rbf)  % does not work for 5, 10, 100, 1000
  max_tries = 1000; 
  for t=1:max_tries
    rd = DistanceMatrixCSRBFwt(x,x,ep_guess,tree);
    A  = rbf(ep_guess, rd);
    % try the Cholesky 
    [L,p] = chol(A,'lower');
    if p==0
      ep = ep_guess; 
      return
    end
    % otherwise reduce guess:
    ep_guess = ep_guess * 0.8;  
  end
  error('Couldn’t find a PD ε after %d tries', max_tries);
end


%% Bisection for sparsity
function ep = findEpForSparsity(x, tree, rbf, target, tol, sepR)
% findEpForSparsity  binary‐searches for ep=1/r so that nnz(A)/numel(A) ≈ target
%
% Inputs:
%   x      : N×d array of node coordinates
%   tree   : KDTreeSearcher built on x
%   rbf    : function handle @(ep,distances)->RBF values
%   target : desired fraction of nonzero A‐entries (0<target<1)
%   tol    : allowable error in achieved sparsity (e.g. 1e-5)
%   sepR   : maximum search radius (the separation radius)
%
% Output:
%   ep     : shape parameter (1/r_mid) that yields sparsity ≈ target

    if nargin<5, tol = 1e-5; end
    if nargin<6
        error('You must pass sepR (the separation radius) as the 6th argument.');
    end
    
    N = size(x,1);
    
    % clamp search range to [0, sepR]
    r_low  = 0;
    r_high = sepR;
    
    for iter = 1:200
    % while true
        r_mid = 0.5*(r_low + r_high);
        % count nnz entries: how many neighbors within r_mid (including self)
        neigh = rangesearch(tree, x, r_mid);
        total_nz = sum(cellfun(@numel, neigh));  
        spars = 1 - (total_nz )/ (N^2);
        
        if abs(spars - target) < tol
        break;
        elseif spars > target
            r_high = r_mid;
        else
            r_low  = r_mid;
        end
        ep = 1 / r_mid;
    end
end

%% Save results 
% Timestamp for uniqueness
timestamp = datestr(datetime('now'), 'yyyyMMdd_HHmmss');

% Construct folder and filename
results_dir = fullfile('CleanLagrangeApprox/results/', sprintf('%s', function_name));
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end
results_filename = fullfile(results_dir, sprintf('results_%s.mat', function_name));


% Save everything
save(results_filename, ...
    'el2_poly', 'elinf_poly', 'a_time_poly', 'e_time_poly', 'c_poly', ...
    'el2_diag', 'elinf_diag', 'a_time_diag', 'e_time_diag', 'c_poly_diag', ...
    'el2_fs1', 'elinf_fs1', 'a_time_fs1', 'e_time_fs1', 'c_poly_fs1', 'sparsity_fs1', ...
    'el2_fs2', 'elinf_fs2', 'a_time_fs2', 'e_time_fs2', 'c_poly_fs2', 'sparsity_fs2', ...
    'el2_fs3', 'elinf_fs3', 'a_time_fs3', 'e_time_fs3', 'c_poly_fs3', 'sparsity_fs3', ...
    'el2_vs1', 'elinf_vs1', 'a_time_vs1', 'e_time_vs1', 'c_poly_vs1', 'sparsity_vs1', ...
    'el2_vs2', 'elinf_vs2', 'a_time_vs2', 'e_time_vs2', 'c_poly_vs2', 'sparsity_vs2', ...
    'el2_vs3', 'elinf_vs3', 'a_time_vs3', 'e_time_vs3', 'c_poly_vs3', 'sparsity_vs3', ...
    'sN', 'dim', 'function_name', 'timestamp');

%% Plot everything
% add flotting using export_fig
% first sort
[ sNs, I ] = sort(sN);

% reorder each error vector/matrix accordingly
el2_poly_s  = el2_poly  (I);
el2_diag_s  = el2_diag  (I,1);
el2_fs1_s   = el2_fs1   (I,1);
el2_fs2_s   = el2_fs2   (I,1);
el2_fs3_s   = el2_fs3   (I,1);
el2_vs1_s   = el2_vs1   (I,1);
el2_vs2_s   = el2_vs2   (I,1);
el2_vs3_s   = el2_vs3   (I,1);

e_fs1 = e_time_fs1(I,1);
e_fs2 = e_time_fs2(I,1);
e_fs3 = e_time_fs3(I,1);
e_vs1 = e_time_vs1(I,1);
e_vs2 = e_time_vs2(I,1);
e_vs3 = e_time_vs3(I,1);

a_fs1 = a_time_fs1(I,1);
a_fs2 = a_time_fs2(I,1);
a_fs3 = a_time_fs3(I,1);
a_vs1 = a_time_vs1(I,1);
a_vs2 = a_time_vs2(I,1);
a_vs3 = a_time_vs3(I,1);


%% ERROR vs N^{1/d}
h = figure; 
hold on; 
grid on;

mark = {'-o','-s','-^','--o','--s','--^','-.x','-.+','-.*'};

semilogy(sNs, el2_poly_s,  mark{1}, 'LineWidth',1.2);
semilogy(sNs, el2_diag_s,  mark{2}, 'LineWidth',1.2);
semilogy(sNs, el2_fs1_s,   mark{3}, 'LineWidth',1.2);
semilogy(sNs, el2_fs2_s,   mark{4}, 'LineWidth',1.2);
semilogy(sNs, el2_fs3_s,   mark{5}, 'LineWidth',1.2);
semilogy(sNs, el2_vs1_s,   mark{6}, 'LineWidth',1.2);
semilogy(sNs, el2_vs2_s,   mark{7}, 'LineWidth',1.2);
semilogy(sNs, el2_vs3_s,   mark{8}, 'LineWidth',1.2);

xlabel(sprintf('N^{1/%d}', dim), 'FontSize', 14);
ylabel('Relative $\ell_2$ error',     'FontSize', 14,'interpreter','latex');
legend({'PLS','Diag','FS1','FS2','FS3','VS1','VS2','VS3'}, ...
       'Location','best','FontSize',10);
title('$\ell_2$ error vs. $N^{1/d}$', 'FontSize',16,'Interpreter','latex');
set(gca,'FontSize',12);

export_fig(h, fullfile(results_dir,'error_vs_N.png'), '-png','-r300');

%% BUILD TIMES vs N^{1/d}
h = figure; hold on; grid on;
% plot just a_time for each strategy:
semilogy(sN, a_time_fs1(:,1), mark{3}, 'LineWidth',1.2);
semilogy(sN, a_time_fs2(:,1), mark{4}, 'LineWidth',1.2);
semilogy(sN, a_time_fs3(:,1), mark{5}, 'LineWidth',1.2);
semilogy(sN, a_time_vs1(:,1), mark{6}, 'LineWidth',1.2);
semilogy(sN, a_time_vs2(:,1), mark{7}, 'LineWidth',1.2);
semilogy(sN, a_time_vs3(:,1), mark{8}, 'LineWidth',1.2);
xlabel(sprintf('N^{1/%d}',dim),'FontSize',14)
ylabel('Assembly time (s)','FontSize',14)
legend({'FS1','FS2','FS3','VS1','VS2','VS3'},'Location','northwest')
title('Assembly time vs. N^{1/d}','FontSize',16)
set(gca,'FontSize',12)
export_fig(h, fullfile(results_dir,'assembly_time_vs_N.png'), '-png','-r300');


%% EVAL TIMES vs N^{1/d}
h = figure; hold on; grid on;
semilogy(sN, e_time_fs1(:,1), mark{3}, 'LineWidth',1.2);
semilogy(sN, e_time_fs2(:,1), mark{4}, 'LineWidth',1.2);
semilogy(sN, e_time_fs3(:,1), mark{5}, 'LineWidth',1.2);
semilogy(sN, e_time_vs1(:,1), mark{6}, 'LineWidth',1.2);
semilogy(sN, e_time_vs2(:,1), mark{7}, 'LineWidth',1.2);
semilogy(sN, e_time_vs3(:,1), mark{8}, 'LineWidth',1.2);
xlabel(sprintf('N^{1/%d}',dim),'FontSize',14)
ylabel('Evaluation time (s)','FontSize',14)
legend({'FS1','FS2','FS3','VS1','VS2','VS3'},'Location','northwest')
title('Evaluation time vs. N^{1/d}','FontSize',16)
set(gca,'FontSize',12)
export_fig(h, fullfile(results_dir,'eval_time_vs_N.png'), '-png','-r300');

%% ERROR vs SPARSITY (one panel per smoothness, same axes, shared legend)
strategies = {'fs1','fs2','fs3','vs1','vs2','vs3'};
labels     = {'FS@0.025','FS@0.05','FS@0.10','VS@0.025','VS@0.05','VS@0.10'};
markers    = {'-o','-s','-^','-x','-+','-*'};
colors     = lines(6);

allS = []; allE = [];
for sm = 1:3
  for j = 1:6
    S = eval(sprintf('sparsity_%s(:,%d)', strategies{j}, sm));
    E = eval(sprintf('el2_%s(:,%d)',      strategies{j}, sm));
    allS = [allS; S(:)];
    allE = [allE; E(:)];
  end
end
xlims = [0, max(allS)*1.05];
ylims = [0, max(allE)*1.05];

figure('Position',[100 100 600 800]);
t = tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact');

nStrat = numel(strategies);
hStrat = gobjects(nStrat+1,1);  % +1 for the PLS line

for sm = 1:3
  ax(sm) = nexttile;
  hold on, grid on

  for j = 1:nStrat
    S = eval(sprintf('sparsity_%s(:,%d)',strategies{j},sm));
    E = eval(sprintf('el2_%s(:,%d)',     strategies{j},sm));
    n = min(numel(S),numel(E));
    hh = plot(S(1:n), E(1:n), markers{j}, ...
              'LineWidth',1.2, 'Color',colors(j,:));
    % store handles only on the first smoothness
    if sm==1
      hStrat(j) = hh;
    end
  end

  % draw PLS baseline
  plh = yline(mean(el2_poly), '--k','PLS','LineWidth',1);
  if sm==1
    hStrat(end) = plh;
  end

  % unify axes across tiles
  xlim(xlims);  ylim(ylims);

  if sm==3
    xlabel('Achieved sparsity','FontSize',12)
  end
  ylabel('Rel. L^2 error','FontSize',12)
  title(sprintf('Wendland smoothness = C%d',sm*2),'FontSize',14)
end

% now make one big legend *on the first axes*:
lg = legend(ax(1), hStrat, [labels,'PLS'], ...
            'Orientation','horizontal', ...
            'Location','northoutside', ...
            'NumColumns',4, ...
            'FontSize',10);

% give a bit more breathing‐room:
t.Padding     = 'loose';
t.TileSpacing = 'loose';

export_fig(gcf, fullfile(results_dir,'error_vs_sparsity.png'), '-png','-r300');