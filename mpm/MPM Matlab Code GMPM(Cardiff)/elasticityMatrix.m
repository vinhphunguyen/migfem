function C = elasticityMatrix(E0,nu0, stressState)

if ( strcmp(stressState,'PLANE_STRESS')  )      % Plane Strain case
  C=E0/(1-nu0^2)*[  1      nu0          0;
                  nu0        1          0;
                    0        0  (1-nu0)/2  ];
else                                            % Plane Strain case
  C=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0      nu0        0;
                            nu0    1-nu0        0;
                              0        0  1/2-nu0 ];
end