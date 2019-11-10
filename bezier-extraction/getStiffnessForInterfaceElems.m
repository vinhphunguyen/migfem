function [K,fint] = getStiffnessForInterfaceElems (K,fint,u)

global noElemsU Cxi Cet iElements weights controlPts p noDofs
global W1 Q1 damage damage0

N       = zeros(2,2*(p+1));

% elemMat = zeros(2*(p+1),2*(p+1));
% elemVec = zeros(2*(p+1),1);


for e=1:noElemsU
    
    % node numbers at the lower and upper faces
    
    sctrL  = iElements(e,1:p+1);
    sctrU  = iElements(e,p+2:end);
    
    % global positions
    
    dofsL(1:2:2*(p+1))  = sctrL*2-1;
    dofsL(2:2:2*(p+1))  = sctrL*2;
    dofsU(1:2:2*(p+1))  = sctrU*2-1;
    dofsU(2:2:2*(p+1))  = sctrU*2;
    
    dispU  = u(dofsU);
    dispL  = u(dofsL);
    
    dJump  = dispU - dispL;
    
    Ce     = Cxi(:,:,e);            % element Bezier extraction operator
    we     = diag(weights(sctrL));  % element weights
    Wb     = Ce'*weights(sctrL);    % element Bezier weights
    pts    = controlPts(sctrL,:);   % element nodes
    
    elemMat = 0.;
    elemVec = 0.;
    
    % loop over Gauss points
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        [Be,dBe]= getShapeGradBernstein(p,xi);
        wb      = dot(Be,Wb);            % Be(I)*Wb(I)
        dwbdxi  = dot(dBe,Wb);           % Be(I)_{,xi} * Wb(I)
        R       = we*Ce*Be/wb;
        dRdxi   = we*Ce*(dBe/wb-dwbdxi*Be/(wb*wb));
        
        jacob   = dRdxi' * pts;
        J       = norm (jacob);
        
        N(1,1:2:2*(p+1)) = R;
        N(2,2:2:2*(p+1)) = R;
        
        % jump at GP
        
        jump(1) = dot(R,dJump(1:2:2*p+2));
        jump(2) = dot(R,dJump(2:2:2*p+2));
        
        % compute traction at GP
       
        [t,T] = getTraction (jump,gp);
        
        % compute stiffness matrix and internal force
        
        elemMat = elemMat + N'*T*N*wt*J;
        elemVec = elemVec + N'*t*wt*J;
    end
    
    % assemble them to the global matrix and vector
    
    fint(dofsL) = fint(dofsL) + elemVec;
    fint(dofsU) = fint(dofsU) - elemVec;
    
    K(dofsL,dofsL) = K(dofsL,dofsL) + elemMat;
    K(dofsU,dofsU) = K(dofsU,dofsU) + elemMat;
    
    elemMat = -elemMat;
    
    K(dofsL,dofsU) = K(dofsL,dofsU) + elemMat;
    K(dofsU,dofsL) = K(dofsU,dofsL) + elemMat;
end



