% -------------------------------------------------------------------------
% Compute the Stress Intensity Factors
% Using the Interaction integral method

% Steps :
% 1- detection of the elements on which we integrate
% 2- loop over these elements
% 3- loop over Gauss points
% 4- computation of stress, strain... in local coordinates !!!   ATTENTION
% 5- computation of the auxilliary fields: AuxStress and AuxEps and AuxGradDisp
% 6- computation of I1 and I2

% Determine J domain and weight function

[Jdomain1,qnode1,radius1] = jIntegrationDomain(tip_elem(1),xTip(1,:),node,elementV);
[Jdomain2,qnode2,radius2] = jIntegrationDomain(tip_elem(2),xTip(2,:),node,elementV);

for it=1:2
    
    if it==1
        Jdomain = Jdomain1;
        qnode   = qnode1;
        radius  = radius1;
        xtip    = xTip(1,:);
        alpha   = alpha1;
    else
        Jdomain = Jdomain2;
        qnode   = qnode2;
        radius  = radius2;
        xtip    = xTip(2,:);
        alpha   = alpha2;
    end
    
    QT  =[cos(alpha) sin(alpha); 
         -sin(alpha) cos(alpha)];
 
    JIntegralOneTip
    
    if ( strcmp(stressState,'PLANE_STRAIN') )
        Eb = E0/(1-nu0^2);
    else
        Eb = E0;
    end
    
    Knum(it,:) = 0.5*Eb*I;
end

Knum



