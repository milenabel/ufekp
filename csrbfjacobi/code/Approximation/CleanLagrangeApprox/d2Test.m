%% Test our hybrid CSRBF + polynomial method on a flat disk domain.
%% Pick of a spatial dimension.
d = 2;

%% Get the evaluation nodes.
start_nodes = 1;

% Evaluation nodes (fixed for all tests)
st = load('DiskPoissonNodesLarge.mat');
xe =  [st.fullintnodes{7}; st.bdrynodes{7}];

% Numbers of interpolation node sets for convergence
st = load('DiskPoissonNodesClustered.mat');
end_nodes = size(st.fullintnodes,1);

%% Polynomial‐degree scaling factor.
fac = 0.8;

%% Target function
%f = exp( ((x(:,1) + x(:,2) + x(:,3)).^(2))/0.1 );
f = @(x) (x(:,1).^2 + x(:,2).^2).^(3/2);    

%% Define a Wendland RBF.
rbf = @(e,r) -r.^4.*(4.*r - 5*spones(r)); %C2
%rbf = @(e,r) (r.^6.*(35.*r.^2 - 88.*r + 56*spones(r)))./3; %C4
%rbf = @(e,r) r.^8.*(66*spones(r) - 154*r + 121*r.^2 - 32*r.^3); %C6


%% Convergence study (2D disk, PLS vs. CSRBF + poly FS)
for k=start_nodes:end_nodes
    %% Build interpolation node set X = {interior; boundary}.
    xi = st.fullintnodes{k};
    xb = st.bdrynodes{k};
    x = [xi;xb];

    alph = 0; %legendre
    
    %% Get the polynomial degree.
    sN(k,1) = nthroot(length(x),dim);
    ell = floor(fac*sN);     
    ell = max([ell,1]); 
    ell_poly(k,1) = ell;
    ell_diag(k,1) = ell;
    ell_fs(k,1) = ell;

    %% Evaluate target at X and xe.
    y = f(x(:,1),x(:,2));
    ye_true = f(xe(:,1),xe(:,2));
    
    %% Build KD‐tree for CSRBF routines.
    tree = KDTreeSearcher(x);
    
    %% Pure polynomial least‐squares (PLS).
    [el2_poly(k,1),elinf_poly(k,1),a_time_poly(k,1),e_time_poly(k,1),c_poly{k,1}] = PLS(x,y,ell,xe,alph,ye_true);   
    
    %% “Diagonal” CSRBF+poly hybrid (CSRBFDiag)
    [el2_diag(k,smoothness),elinf_diag(k,smoothness),a_time_diag(k,smoothness),e_time_diag(k,smoothness),c_poly_diag{k,smoothness}] = CSRBFDiag(x,y,ell,xe,alph,rbf,tree,ye_true);    
    
    %% Fixed‐shape CSRBF+poly (CSRBFGen) with epsilon = 10
    ep = 10;
    [el2_fs(k,smoothness), elinf_fs(k,smoothness), a_time_fs(k,smoothness), e_time_fs(k,smoothness), c_poly_fs{k,smoothness}, cond_fs(k,smoothness), ~, sparsity_fs(k,smoothness)] = CSRBFGen(x,y,ell,xe,alph,rbf,ep1,tree,ye_true);
end