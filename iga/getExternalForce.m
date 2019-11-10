function [fext] = getExternalForce(bndMesh,bndPoints,D,func,t)

%
% Compute external force vector applied on an 2D edge.
%
% VP Nguyen
% Cardiff University

global element controlPts p q uKnot vKnot weights index elRangeU elRangeV 
global Q W Q1 W1 rho C noDofs noCtrPts elConnV I

fext = zeros(noDofs,1);        % external force vector

noElemsV = size(bndMesh,1);

for e=1:noElemsV
    xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
    conn  = elConnV(e,:);
    sctry = 2*bndMesh(e,:);
    pts   = bndPoints(conn,:);
    
    % loop over Gauss points
    
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
        J2      = 0.5 * ( xiE(2) - xiE(1) );
        
        [N dNdxi] = NURBS1DBasisDers(Xi,q,vKnot,weights);
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        x        = N    *pts; % global coord of GP
        jacob1   = dNdxi*pts;
        J1       = norm (jacob1);
        
        gt       = func(t);         
        f        = 1000 * gt;
        ty       = -f/(2*I)*((D*D)/4-x(1,2)^2);
              

        fext(sctry) = fext(sctry) + N' * ty * J1 * J2 * wt;
    end
end
