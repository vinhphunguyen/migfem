
disp_norm1 = 0 ;
disp_norm2 = 0 ;

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
                          
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE,pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        % compute basis and derivatives
        
        [N dRdxi dRdeta] = NURBS2DBasisDers([Xi; Eta],p,q,uKnot,vKnot,weights');
                    
        pts        = controlPts(sctr,:);
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
        
        % Numerical strain and stress
        
        strain          = B*U(sctrB);                        
        num_stress      = C*strain;
        %num_stress(3)   = 0.5*num_stress(3); shear stress
        
        % Exact stress
        
        x = strPoint(1,1);
        y = strPoint(1,2);
        
        exact_stress(1) = (1/I)*P*(L-x)*y;
        exact_stress(2) = 0;
        exact_stress(3) = -0.5*(P/I)*(0.25*D^2 - y^2);
                      
        errorInNorm = num_stress - exact_stress';  
        %(Cinv*errorInNorm)'*errorInNorm*wt*J1*J2
        energy_norm = energy_norm + ...
            (Cinv*errorInNorm)'*errorInNorm*wt*J1*J2;        
    end
end


energy_norm = sqrt(0.5*energy_norm);
