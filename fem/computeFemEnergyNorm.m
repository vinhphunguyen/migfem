energy_norm = 0;

for e=1:numelem                          % start of element loop
    sctr=element(e,:);           % element scatter vector
    nn=length(sctr);
    for q=1:size(W,1)
        pt=Q(q,:);
        wt=W(q);
        [N,dNdxi]=lagrange_basis(elemType,pt); % element shape functions
        J0=node(sctr,:)'*dNdxi;                % element Jacobian matrix
        invJ0=inv(J0);
        dNdx=dNdxi*invJ0;
        
        B(1,1:nn)       = dNdx(:,1)';
        B(2,nn+1:2*nn)  = dNdx(:,2)';
        B(3,1:nn)       = dNdx(:,2)';
        B(3,nn+1:2*nn)  = dNdx(:,1)';
        
        strain          = B*U(sctrB);
        num_stress      = C*strain;
        %num_stress(3)   = 0.5*num_stress(3);
        % Exact stress
        
        x = N'*node(sctr,1);
        y = N'*node(sctr,2);
        
        exact_stress(1) = (1/I0)*P*(L-x)*y;
        exact_stress(2) = 0;
        exact_stress(3) = -0.5*(P/I0)*(c^2 - y^2);
        
        errorInNorm = num_stress - exact_stress';
        energy_norm = energy_norm + ...
            ((inv(C)*errorInNorm)'*errorInNorm)*wt*det(J0);
        
    end  % of quadrature loop
end    % of element loop

energy_norm = sqrt(0.5*energy_norm)