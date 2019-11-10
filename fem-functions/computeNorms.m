%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute displacement and energy norm
% this script must be used for problems having
% analytical solution.
% For infinite plate with a circular hole: use line 71
% For Timoshenko beam use line 72
% It does not work for extended IGA!!!
% Vinh Phu Nguyen
% Delft University of Technology, The Netherlands

Cinv           = inv(C);
energy_norm    = 0;
disp_norm      = 0;
disp_normDenom = 0;

[W,Q] = quadrature(4, 'GAUSS', 2); 

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
                      
    sctr   = element(e,:);          %  element scatter vector
    sctrB  = [sctr sctr+noCtrPts]; %  vector that scatters a B matrix
    nn     = length(sctr);
    pts    = controlPts(sctr,:);
    
    B      = zeros(3,2*nn);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE, pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        % compute basis and derivatives
        
        [N dRdxi dRdeta] = NURBS2DBasisDers([Xi; Eta],p,q,uKnot,vKnot,weights');
                            
        jacob      = pts' * [dRdxi' dRdeta'];
        strPoint   = N * pts;
        J1         = det(jacob);
        
        % Jacobian inverse and spatial derivatives
        
        invJacob   = inv(jacob);
        dRdx       = [dRdxi' dRdeta'] * invJacob;
        
        % B matrix
        
        B(1,1:nn)       = dRdx(:,1)';
        B(2,nn+1:2*nn)  = dRdx(:,2)';
        B(3,1:nn)       = dRdx(:,2)';
        B(3,nn+1:2*nn)  = dRdx(:,1)';
        
        % Numerical strain and stress and displacements
        
        strain          = B*U(sctrB);                        
        num_stress      = C*strain;
        num_disp        = N*[Ux(sctr) Uy(sctr)];
              
        % Exact stresses and displacements
            
        [exact_stress exact_disp] = exact_solution_hole(strPoint,a,E,nu);
        
%         [exact_stress exact_disp]= exact_solution_timoshenkobeam(...
%                                       strPoint,E0,nu0,I,L,D,P);
                                 
        errorInNorm = num_stress - exact_stress';  
        errorInDisp = num_disp   - exact_disp;
        
        fac            = wt*J1*J2;
        energy_norm    = energy_norm    + (Cinv*errorInNorm)'*errorInNorm*fac;        
        disp_normDenom = disp_normDenom + exact_disp*exact_disp'  *fac;        
        disp_norm      = disp_norm      + errorInDisp*errorInDisp'*fac;        
    end
end

disp_norm   = sqrt(disp_norm)
%disp_norm   = sqrt(disp_norm/disp_normDenom)
energy_norm = sqrt(0.5*energy_norm)
