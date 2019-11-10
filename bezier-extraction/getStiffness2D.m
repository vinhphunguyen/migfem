function [K,fint] = getStiffness2D (u)

global element C weights controlPts W Q shapes derivs noElems De noBasis noDofs

K    = zeros(noDofs,noDofs);  % global stiffness matrix
fint = zeros(noDofs,1);

% Loop over elements

for e=1:noElems
    sctr   = element(e,:);         %  element scatter vector
    nn     = length(sctr);
    
    sctrB(1,1:2:2*nn) = 2*sctr-1;
    sctrB(1,2:2:2*nn) = 2*sctr;
    
    elemDisp = u(sctrB);
        
    B      = zeros(3,2*nn);
    
    Ce     = C(:,:,e);             % element Bezier extraction operator
    we     = diag(weights(sctr));  % element weights
    pts    = controlPts(sctr,:);   % element nodes
    Wb     = Ce'*weights(sctr);    % element Bezier weights
    
    % loop over Gauss points
    
    for gp=1:size(W,1)        
        wt      = W(gp);
        
        %% Bernstein basis and derivatives at GP gp
        Be      = shapes(gp,:)';
        dBedxi  = reshape(derivs(gp,:,:),noBasis,2);
        
        %% Bezier weight functions (denomenator of NURBS)
        wb        = dot(Be,Wb);            % Be(I)*Wb(I)
        dwbdxi(1) = dot(dBedxi(:,1),Wb);   % Be(I)_{,xi} * Wb(I)
        dwbdxi(2) = dot(dBedxi(:,2),Wb);   % Be(I)_{,et} * Wb(I)
        %% Shape function and derivatives        
        dRdxi(:,1) = we*Ce*(dBedxi(:,1)/wb-dwbdxi(1)*Be/(wb*wb));
        dRdxi(:,2) = we*Ce*(dBedxi(:,2)/wb-dwbdxi(2)*Be/(wb*wb));
        
        %% Jacobian matrix
        dxdxi = pts'*dRdxi;
        
        dxidx = inv(dxdxi);
        dRdx  = dRdxi*dxidx;
        detJ  = det(dxdxi);
        
        % B matrix
        B(1,1:2:2*nn)  = dRdx(:,1)';
        B(2,2:2:2*nn)  = dRdx(:,2)';
        B(3,1:2:2*nn)  = dRdx(:,2)';
        B(3,2:2:2*nn)  = dRdx(:,1)';
        
        sigma = De * B * elemDisp;
        
        % compute elementary stiffness matrix and
        % assemble it to the global matrix
        
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * De * B * detJ * wt;
        fint(sctrB)    = fint(sctrB)    + B' * sigma * detJ * wt;
    end   
end
