function Ne = computeNMatrix(N)
%
% Compute the shape function matrix for a rotation free thin shell
% element.
%
% VP Nguyen
% Cardiff University, UK, Cardiff

nn=length(N);

Ne=zeros(3,3*nn);

Ne(1,1:3:3*nn) = N;
Ne(2,2:3:3*nn) = N;
Ne(3,3:3:3*nn) = N;
