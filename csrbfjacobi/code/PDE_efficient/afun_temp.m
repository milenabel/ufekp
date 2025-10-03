function[y] = afun_temp(x_geom,x_coeff,nr,alph,a,a_structure,dimension_jumps,recurrence,lc,uc,pc,qc,Ni)

V = mpoly_eval(x_geom, a, recurrence); %evaluate polynomial basis

%% Get coefficients for derivatives
d2x = mjacobi_symm_faster_differentiation(x_coeff, a, alph, [2 0], a_structure, dimension_jumps, lc,uc,pc,qc);
d2y = mjacobi_symm_faster_differentiation(x_coeff, a, alph, [0 2], a_structure, dimension_jumps, lc,uc,pc,qc);        
dx = mjacobi_symm_faster_differentiation(x_coeff, a, alph, [1 0], a_structure, dimension_jumps, lc,uc,pc,qc);
dy = mjacobi_symm_faster_differentiation(x_coeff, a, alph, [0 1], a_structure, dimension_jumps, lc,uc,pc,qc);

%% Multiply interp mat by diff coeffients        
lap = V*(d2x + d2y);
grad_x = V*dx;
grad_y = V*dy;
neu = nr(:,1).*grad_x(Ni+1:end) +  nr(:,2).*grad_y(Ni+1:end);

%% Boundary bordering
y = [lap(1:Ni);neu];
