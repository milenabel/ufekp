%% Interpolation in 1D/2D/3D domains using CSRBFs combined with
%%Jacobi polynomials. 
clearvars; close all; clc;

project_root = fileparts(mfilename('fullpath'));
results_root  = fullfile(project_root, 'results');

%% Define all test cases
test_cases = {
    % Dimension 1
    % struct('dim', 1, 'functions', {
    %     struct('f', 'abs(x)', 'name', 'abs_1d');
    %     struct('f', 'exp(-x.^(-2))', 'name', 'gauss_1d');
    %     struct('f', '1./(1 + 16*x.^2)', 'name', 'rational_1d');
    %     struct('f', 'x.^(10)', 'name', 'poly10_1d');
    % }),

    % Dimension 2
    struct('dim', 2, 'functions', {
        % struct('f', 'abs(x).^3 .* abs(y).^3', 'name', 'abs3x3_2d');
        struct('f', 'exp(-x.^(-2)).*exp(-y.^(-2))', 'name', 'gauss_2d');
        struct('f', '1./(1 + 25*(x.^2 + y.^2))', 'name', 'rational_2d');
        struct('f', 'exp(-10*((x-.3).^(-2)+y.^(-2)))', 'name', 'gauss10_2d');
        struct('f', 'exp(-10*((x-.3).^2+y.^2))', 'name', 'gauss10inv_2d');
        struct('f', 'x.^8 .* y.^8', 'name', 'poly8x8_2d');
    }),

    % % Dimension 3
    % struct('dim', 3, 'functions', {
    %     struct('f', 'abs(x).^3.*abs(y).^3.*abs(z).^3', 'name', 'abs333_3d');
    %     struct('f', 'exp(-10*((x-.3).^(-2)+y.^(-2) + z.^(-2)))', 'name', 'gaussinv_3d');
    %     struct('f', 'exp(-10*((x-.3).^2+y.^2 + z.^2))', 'name', 'gauss_3d');
    %     struct('f', '1./(1 + 16*(x.^2 + y.^2 + z.^2))', 'name', 'rational_3d');
    %     struct('f', 'x.^(4).*y.^(2).*z.^(2)', 'name', 'poly422_3d');
    % })
};

%% Main loop over dimensions and functions
for dim_case = 1:length(test_cases)
    current_dim = test_cases{dim_case}.dim;
    dim_functions = test_cases{dim_case}.functions;

    fprintf('\n=== Processing Dimension %d (%d functions) ===\n', current_dim, numel(dim_functions));

    for func_idx = 1:length(dim_functions)
        current_func = dim_functions(func_idx);
        fprintf('-- Testing function: %s\n', current_func.name);

        %% Run the interpolation tests for this dimension and function
        run_interpolation_tests(current_dim, current_func.f, current_func.name, results_root);
    end
end

function run_interpolation_tests(dim, func_str, function_name, results_root)

results_dir = fullfile(results_root, function_name);
    if ~exist(results_dir,'dir')
        mkdir(results_dir)
    end
    fprintf('saving into: %s\n', results_dir);

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
        st = load('DiskPoissonNodesClustered.mat');
        clear p N pxe;
    elseif dim==3
        st = load('SpherePoissonNodesLarge.mat');
        xe = [st.fullintnodes{2}; st.bdrynodes{2}];
        st = load('SpherePoissonNodesClustered.mat');
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

    if dim == 1
        syms x;
        f = eval(func_str);
        dfx = diff(f,x);
    elseif dim == 2
        syms x y;
        f = eval(func_str);
        dfx = diff(f,x); dfy = diff(f,y);
    elseif dim == 3
        syms x y z;
        f = eval(func_str);
        dfx = diff(f,x); dfy = diff(f,y); dfz = diff(f,z);
    end

    %% Create symbolic function
    if dim==1
        syms x
        f = str2sym(func_str);
        dfx = diff(f,x);
        f = matlabFunction(f, 'Vars', x);
        dfx = matlabFunction(dfx, 'Vars', x);
        clear x;
    elseif dim==2
        syms x y
        f = str2sym(func_str);
        dfx = diff(f,x); dfy = diff(f,y);
        f = matlabFunction(f, 'Vars', [x y]);
        dfx = matlabFunction(dfx, 'Vars', [x y]);
        dfy = matlabFunction(dfy, 'Vars', [x y]);
        clear x y;
    elseif dim==3
        syms x y z
        f = str2sym(func_str);
        dfx = diff(f,x); dfy = diff(f,y); dfz = diff(f,z);
        f = matlabFunction(f, 'Vars', [x y z]);
        dfx = matlabFunction(dfx, 'Vars', [x y z]);
        dfy = matlabFunction(dfy, 'Vars', [x y z]);
        dfz = matlabFunction(dfz, 'Vars', [x y z]);
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

    %% Now find shape parameters that induce a target condition number on the finest node set
    %% Use those on the coarser node sets, and it looks like the supports are increasing
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

        %% Three picked condition‐number targets
        K_targets = [1e12, 1e8, 1e4];

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

    %% Now do fixed condition number strategies
    %% Again, different smoothnesses
    K_targets = [1e12, 1e8, 1e4]; 

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

        options.TolX = 1e-4;

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
            tree  = KDTreeSearcher(x);    
            [~,dist] = knnsearch(tree,x,'k',2);
            % if dim==1
            %     sep_dist = 0.5*mean(dist(:,2)); % Same fix as in fs section
            % else
            %     sep_dist = 0.5*min(dist(:,2));
            % end
            sep_dist = 0.5*min(dist(:,2));


            %% Milena: CHECK if this works
            % guess an eigenvalue. 3 here is the dim of the Wendland kernel (always 3 in this code)
            % the fourier transform of the Wendland kernel
            % decays at this rate (p), so the eigenvalues must satisfy this on a given
            % point set. This should give us a good guess and possibly a bracket
            p = 2*smoothness + 3 + 1; 
            eps_guess = K_targets.^(-1./p)/sep_dist; 
            eps_vs = zeros(size(eps_guess));  
            options.TolX = 1e-2; %loose tolerance for speed
            for kit=1:3        
                eps0 = eps_guess(kit);
                k_target = K_targets(kit);
                ep_func = @(ep) log10( cond(full(rbf(ep,DistanceMatrixCSRBFwt(x,x,ep,tree))))) - log10(k_target);
                eps_vs(kit) = fzero(ep_func,[eps0*0.01,eps0*10],options); % a bracketed search
            end

            all_eps_vs(smoothness,:) = eps_vs; 
            ep1 = eps_vs(1);
            ep2 = eps_vs(2);
            ep3 = eps_vs(3);

            [el2_vs1(k,smoothness), elinf_vs1(k,smoothness), a_time_vs1(k,smoothness), e_time_vs1(k,smoothness), c_poly_vs1{k,smoothness}, cond_vs1(k,smoothness), ~, sparsity_vs1(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep1,tree,ye_true);

            [el2_vs2(k,smoothness), elinf_vs2(k,smoothness), a_time_vs2(k,smoothness),  e_time_vs2(k,smoothness), c_poly_vs2{k,smoothness}, cond_vs2(k,smoothness), ~, sparsity_vs2(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep2,tree,ye_true);

            [el2_vs3(k,smoothness), elinf_vs3(k,smoothness), a_time_vs3(k,smoothness), e_time_vs3(k,smoothness), c_poly_vs3{k,smoothness}, cond_vs3(k,smoothness), ~, sparsity_vs3(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep3,tree,ye_true);
        end
    end


    %% Save results 
    % Timestamp for uniqueness
    timestamp = datestr(datetime('now'), 'yyyyMMdd_HHmmss');

    % % Construct folder and filename
    % results_dir = fullfile('CleanLagrangeApprox/results/', sprintf('%s', function_name));
    % if ~exist(results_dir, 'dir')
    %     mkdir(results_dir);
    % end
    % results_filename = fullfile(results_dir, sprintf('results_%s.mat', function_name));

    save_file    = fullfile(results_dir, sprintf('results_%s_%s.mat', function_name, timestamp));
    fprintf('saving final MAT to:\n %s\n\n', save_file);
    save(save_file, '-v7.3');    


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

end