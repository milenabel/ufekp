function [ c ] = rungmres(weights,X,b,nr,alph,a,a_structure,dimension_jumps,recurrence,lc,uc,pc,qc,Ni)

    %% Call gmres    
    maxit = 600;
    restart = 5;
    tol = 1e-8;
    V = mpoly_eval(X, a, recurrence); %evaluate polynomial basis
    V = repmat(weights, [1 size(a,1)]).*V;
    c = gmres(@afun,b,restart,tol,maxit,[],[]);
    
    %% Helper function that applies system matrix to a vector using
    %fast algorithms for differentiation.
    function y = afun(x)
        %% Get coefficients for derivatives
        d2x = mjacobi_symm_faster_differentiation(x, a, alph, [2 0], a_structure, dimension_jumps, lc,uc,pc,qc);
        d2y = mjacobi_symm_faster_differentiation(x, a, alph, [0 2], a_structure, dimension_jumps, lc,uc,pc,qc);        
        dx = mjacobi_symm_faster_differentiation(x, a, alph, [1 0], a_structure, dimension_jumps, lc,uc,pc,qc);
        dy = mjacobi_symm_faster_differentiation(x, a, alph, [0 1], a_structure, dimension_jumps, lc,uc,pc,qc);
        
        %% Multiply interp mat by diff coeffients        
        lap = V*(d2x + d2y);
        grad_x = V*dx;
        grad_y = V*dy;
        neu = nr(:,1).*grad_x(Ni+1:end) +  nr(:,2).*grad_y(Ni+1:end);
        
        %% Boundary bordering
        y = weights.*[lap(1:Ni);neu];
    end

end

