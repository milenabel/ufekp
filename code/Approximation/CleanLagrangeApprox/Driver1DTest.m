%% Test unified interp and PLS on KTE points
dim = 1;
xe = linspace(-1,1,2^14).';
N = 2.^(2:8); N = N';
alpha = 0.8;
for k=1:length(N)
    s = chebspace2(-1,1,N(k));
    X = asin(alpha*s)./asin(alpha);
    xi = X(2:end-1,:);
    xb = [X(1,:); X(end,:)];
    st.fullintnodes{k,1} = xi;
    st.bdrynodes{k,1} = xb;
end
clear xb xi k;
f = @(x) 1./sqrt(1 + 25*x.^2);
start_nodes = 1;
end_nodes = size(st.fullintnodes,1);
fac = 1;

%% Get the standard polynomial least squares stuff out of the way
for k=start_nodes:end_nodes
    xi = st.fullintnodes{k};
    xb = st.bdrynodes{k};
    x = [xi;xb];
    alph = 0; %legendre
    ell = floor(fac*nthroot(length(x),dim));
    ell = max([ell,1]); 
    ell_poly(k,1) = ell;
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
for smoothness=1:1
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
        ell_diag(k,smoothness) = ell;

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

%% Again, different smoothnesses
for smoothness=1:1
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

    %% Three picked condition‐number targets
    %K_targets = [1e12, 1e8, 1e4];
    K_targets = [1e4];

    %% Find the epsilons for each K_target on the finest node set
    k_finest = end_nodes;

    xi = st.fullintnodes{k_finest};
    xb = st.bdrynodes{k_finest};
    x_finest = [xi;xb];
    tree_finest = KDTreeSearcher(x_finest);
    [~,dist] = knnsearch(tree_finest,x_finest,'k',2);
    sep_dist = 0.5*min(dist(:,2));
    if dim==1
        sep_dist = 0.5*mean(dist(:,2));
    end

    %% Milena: this works
    % guess an eigenvalue. 3 here is the dim of the Wendland kernel (always 3 in this code)
    % the fourier transform of the Wendland kernel
    % decays at this rate (p), so the eigenvalues must satisfy this on a given
    % point set. This should give us a good guess and possibly a bracket
    p = 2*smoothness + 3 + 1; 
    eps_guess = K_targets.^(-1./p)/sep_dist; 
    eps_fs = zeros(size(eps_guess));  
    options.TolX = 1e-2; %loose tolerance for speed
    for kit=1:length(K_targets)        
        eps0 = eps_guess(kit);
        k_target = K_targets(kit);
        ep_func = @(ep) log10( cond(full(rbf(ep,DistanceMatrixCSRBFwt(x_finest,x_finest,ep,tree_finest))))) - log10(k_target);
        eps_fs(kit) = fzero(ep_func,[eps0*0.01,eps0*10],options); % a bracketed search
        % eps_fs(kit) = fzero(ep_func,[eps0,45*eps0],options); % a bracketed search
    end
    all_eps_fs(smoothness,:) = eps_fs; 
    supports_fs(smoothness,:) = 1./all_eps_fs(smoothness,:);

    for k=start_nodes:end_nodes
        xi = st.fullintnodes{k};
        xb = st.bdrynodes{k};
        x = [xi;xb];
        alph = 0; %legendre
        ell = floor(fac*nthroot(length(x),dim));
        ell = max([ell,1]); 
        ell_fs1(k,smoothness) = ell;  % For K=1e12
        ell_fs2(k,smoothness) = ell;  % For K=1e8
        ell_fs3(k,smoothness) = ell;  % For K=1e4
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

        % For K_target = 1e4 (j=1)
        ep1 = eps_fs(1);
        [el2_fs1(k,smoothness), elinf_fs1(k,smoothness), a_time_fs1(k,smoothness), e_time_fs1(k,smoothness), c_poly_fs1{k,smoothness}, cond_fs1(k,smoothness), ~, sparsity_fs1(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep1,tree,ye_true);
    end
end
