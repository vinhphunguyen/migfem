function C = elasticityMatrix(E0,nu0, stressState)
%
% Elasticity matrix for isotropic elastic materials.
%
% VP Nguyen
% Cardiff University, UK

if     ( strcmp(stressState,'PLANE_STRESS')  )      % Plane Stress case
    C=E0/(1-nu0^2)*[  1      nu0         0;
                      nu0     1          0;
                      0       0  (1-nu0)/2  ];
elseif ( strcmp(stressState,'PLANE_STRAIN')  )      % Plane Strain case
    C=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0   nu0        0;
                             nu0    1-nu0       0;
                             0      0  1/2-nu0 ];
else                                                 % 3D
    C=zeros(6,6);
    C(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
                                      nu0 1-nu0 nu0;
                                      nu0 nu0 1-nu0];
    C(4:6,4:6)=E0/2/(1+nu0)*eye(3);
end