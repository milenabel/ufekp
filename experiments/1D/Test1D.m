%% Interpolation in 1D/2D/3D domains using CSRBFs combined with
%%Jacobi polynomials. 


%% Spatial dimension
dim = 1;

%% Load up the node set
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

start_nodes = 1;
end_nodes = size(st.fullintnodes,1);

% Sparsity storage initialization 
sparsity_fs1 = zeros(end_nodes, 3);
sparsity_fs2 = zeros(end_nodes, 3);
sparsity_fs3 = zeros(end_nodes, 3);
sparsity_fc1 = zeros(end_nodes, 3);
sparsity_fc2 = zeros(end_nodes, 3);
sparsity_fc3 = zeros(end_nodes, 3);

%% This partially determines the number of polynomial terms
%slightly smaller fac helps with increasing accuracy for both rough and
%smooth functions. There seems to be a notion of the "best" polynomial
%degree for a given number of nodes
fac = 1;

%% Setup anonymous function for function interpolation true solution
syms x;       
f = abs(x);                function_name = 'abs_1d';
%f = 1./(1 + 25*x.^2);     function_name = 'rk_1d';
dfx = diff(f,x);

f = matlabFunction(f);
dfx = matlabFunction(dfx);
clear x;

%% All possible tests:
%% 1. Different smoothness for CSRBF
%% 2. Different shape parameter strats
%% 3. Different interp techniques

%% Get the standard polynomial least squares results
for k=start_nodes:end_nodes
    xi = st.fullintnodes{k};
    xb = st.bdrynodes{k};
    x = [xi;xb];
    alph = 0; %legendre
    ell = floor(fac*nthroot(length(x),dim));            
    ell = max([ell,1]); 
    y = f(x(:,1));
    ye_true = f(xe(:,1));   
    [el2_poly(k,1),elinf_poly(k,1),a_time_poly(k,1),e_time_poly(k,1),c_poly{k,1}] = PLS(x,y,ell,xe,alph,ye_true);   
    sN(k,1) = nthroot(length(x),dim);
end


%% Next, get the diagonal approximation results for different smoothnesses
for smoothness=1:3
    if smoothness==1
        %% Wendland C2 in 3d, pd in all lower dimensions
        rbf = @(e,r) -r.^4.*(4.*r - 5*spones(r));
        drbfor = @(e,r) 20.*e.^2.*r.^3;
        lrbf = @(e,r) -20.*e.^2.*r.^2.*(2.*r - 3*spones(r));
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
        y = f(x(:,1));
        ye_true = f(xe(:,1));  
        tree = KDTreeSearcher(x);
        [el2_diag(k,smoothness),elinf_diag(k,smoothness),a_time_diag(k,smoothness),e_time_diag(k,smoothness),c_poly_diag{k,smoothness}] = CSRBFDiag(x,y,ell,xe,alph,rbf,tree,ye_true);    
    end
end

%% Now find shape parameters that induce a target condition number on the finest node set
%% Use those on the coarser node sets, and it looks like the supports are increasing
%% Again, different smoothnesses
for smoothness=1:3
    if smoothness==1
        %% Wendland C2 in 3d, pd in all lower dimensions
        rbf = @(e,r) -r.^4.*(4.*r - 5*spones(r));
        drbfor = @(e,r) 20.*e.^2.*r.^3;
        lrbf = @(e,r) -20.*e.^2.*r.^2.*(2.*r - 3*spones(r));
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

    %% Three picked condition‐number targets
    K_targets = [1e12, 1e8, 1e4];

    %% Find the epsilons for each K_target on the finest node set
    k_finest = end_nodes;

    xi = st.fullintnodes{k_finest};
    xb = st.bdrynodes{k_finest};
    x_finest = [xi;xb];
    tree_finest = KDTreeSearcher(x_finest);
    [~,dist] = knnsearch(tree_finest,x_finest,'k',2);
    sep_dist = 0.5*mean(dist(:,2));

    p = 2*smoothness + 3 + 1; 
    eps_guess = K_targets.^(-1./p)/sep_dist; 
    eps_fs = zeros(size(eps_guess));  
    options.TolX = 1e-2; %loose tolerance for speed
    for kit=1:3        
        eps0 = eps_guess(kit);
        k_target = K_targets(kit);
        ep_func = @(ep) log10( cond(full(rbf(ep,DistanceMatrixCSRBFwt(x_finest,x_finest,ep,tree_finest))))) - log10(k_target);
        eps_fs(kit) = fzero(ep_func,[eps0*0.01,eps0*10],options); % a bracketed search
    end
    all_eps_fs(smoothness,:) = eps_fs; 

    for k=start_nodes:end_nodes
        xi = st.fullintnodes{k};
        xb = st.bdrynodes{k};
        x = [xi;xb];
        alph = 0; %legendre
        ell = floor(fac*nthroot(length(x),dim));            
        ell = max([ell,1]); 
        y = f(x(:,1));
        ye_true = f(xe(:,1));

        tree = KDTreeSearcher(x);

        % For K_target = 1e12 (j=1)
        ep1 = eps_fs(1);
        [el2_fs1(k,smoothness), elinf_fs1(k,smoothness), a_time_fs1(k,smoothness), e_time_fs1(k,smoothness), c_poly_fs1{k,smoothness}, cond_fs1(k,smoothness), ~, sparsity_fs1(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep1,tree,ye_true);

        % For K_target = 1e8 (j=2)
        ep2 = eps_fs(2);
        [el2_fs2(k,smoothness),elinf_fs2(k,smoothness),a_time_fs2(k,smoothness),e_time_fs2(k,smoothness),c_poly_fs2{k,smoothness}, cond_fs2(k,smoothness), ~, sparsity_fs2(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep2,tree,ye_true); 

        % For K_target = 1e4 (j=3)
        ep3 = eps_fs(3);
        [el2_fs3(k,smoothness),elinf_fs3(k,smoothness),a_time_fs3(k,smoothness),e_time_fs3(k,smoothness),c_poly_fs3{k,smoothness}, cond_fs3(k,smoothness), ~, sparsity_fs3(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep3,tree,ye_true);
    end
end


%% Now, we do fixed condition number strategies
%% Again, different smoothnesses
K_targets = [1e12, 1e8, 1e4]; 

for smoothness=1:3
    if smoothness==1
        %% Wendland C2 in 3d, pd in all lower dimensions
        rbf = @(e,r) -r.^4.*(4.*r - 5*spones(r));
        drbfor = @(e,r) 20.*e.^2.*r.^3;
        lrbf = @(e,r) -20.*e.^2.*r.^2.*(2.*r - 3*spones(r));
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

    options.TolX = 1e-4;

    for k=start_nodes:end_nodes
        xi = st.fullintnodes{k};
        xb = st.bdrynodes{k};
        x = [xi;xb];

        Nnodes = size(x,1);

        alph = 0; %legendre
        ell = floor(fac*nthroot(length(x),dim));            
        ell = max([ell,1]); 
        y = f(x(:,1));
        ye_true = f(xe(:,1));
        tree  = KDTreeSearcher(x);    
        [~,dist] = knnsearch(tree,x,'k',2);
        sep_dist = 0.5*mean(dist(:,2));
        p = 2*smoothness + 3 + 1; 
        eps_guess = K_targets.^(-1./p)/sep_dist; 
        eps_fc = zeros(size(eps_guess));  
        options.TolX = 1e-6;
        for kit=1:3        
            eps0 = eps_guess(kit);
            k_target = K_targets(kit);
            ep_func = @(ep) log10( cond(full(rbf(ep,DistanceMatrixCSRBFwt(x,x,ep,tree))))) - log10(k_target);
            eps_fc(kit) = fzero(ep_func,[eps0*0.001,eps0*100],options); % a bracketed search
        end

        all_eps_fc(smoothness,:) = eps_fc; 
        ep1 = eps_fc(1);
        ep2 = eps_fc(2);
        ep3 = eps_fc(3);

        [el2_fc1(k,smoothness), elinf_fc1(k,smoothness), a_time_fc1(k,smoothness), e_time_fc1(k,smoothness), c_poly_fc1{k,smoothness}, cond_fc1(k,smoothness), ~, sparsity_fc1(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep1,tree,ye_true);
        
        [el2_fc2(k,smoothness), elinf_fc2(k,smoothness), a_time_fc2(k,smoothness),  e_time_fc2(k,smoothness), c_poly_fc2{k,smoothness}, cond_fc2(k,smoothness), ~, sparsity_fc2(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep2,tree,ye_true);
        
        [el2_fc3(k,smoothness), elinf_fc3(k,smoothness), a_time_fc3(k,smoothness), e_time_fc3(k,smoothness), c_poly_fc3{k,smoothness}, cond_fc3(k,smoothness), ~, sparsity_fc3(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep3,tree,ye_true);
    end
end

%% Save results 
% Timestamp for uniqueness
timestamp = datestr(datetime('now'), 'yyyyMMdd_HHmmss');

% Construct folder and filename
results_dir = fullfile('results/', sprintf('%s', function_name),'/high');
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
    'el2_fc1', 'elinf_fc1', 'a_time_fc1', 'e_time_fc1', 'c_poly_fc1', 'sparsity_fc1', ...
    'el2_fc2', 'elinf_fc2', 'a_time_fc2', 'e_time_fc2', 'c_poly_fc2', 'sparsity_fc2', ...
    'el2_fc3', 'elinf_fc3', 'a_time_fc3', 'e_time_fc3', 'c_poly_fc3', 'sparsity_fc3', ...
    'sN', 'dim', 'function_name', 'timestamp');