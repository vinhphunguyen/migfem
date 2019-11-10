%---------------------------------------------
% Compute interaction integral

I1 = 0;
I2 = 0;
I  = [zeros(2,1)];
Ideb=[];

globGP=[];

% ---------------------------
% Starting LOOP over ELEMENTS
%----------------------------

for iel = 1 : size(Jdomain,2)    
    e      = Jdomain(iel) ; % current element
    sctr   = element(e,:);
    sctrQ4 = elementV(e,:);
    nn     = length(sctr);
    
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    % Choose Gauss quadrature rule
    
    if (ismember(e,split_elem))     % split element
        order = 12; % 13 GPs for each sub-triangle
        %phi   = ls(sctr,1);
        %[W,Q] = discontQ4quad(order,phi);
        [W,Q] = quadrature(order,'GAUSS',2);
    else
        order = 8 ; 
        [W,Q] = quadrature(order,'GAUSS',2);
    end
    
    % -----------------------------
    % start loop over Gauss points
    % -----------------------------
    
    for igp = 1:size(W,1)        
        pt = Q(igp,:);
        wt = W(igp);
        
        % Q4 element for weighting q
        [N,dNdxi] = lagrange_basis('Q4',pt);
        J0    = node(sctrQ4,:)'*dNdxi;
        invJ0 = inv(J0);
        dNdx  = dNdxi*invJ0;
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE,pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
    
        % compute derivative of basis functions w.r.t parameter coord

        [N dRdxi dRdeta] = NURBS2DBasisDers([Xi; Eta],p,q,uKnot,vKnot,weights');

        pts        = controlPts(sctr,:);
        jacob      = pts'*[dRdxi' dRdeta'];
        J1         = det(jacob);
        invJacob   = inv(jacob);
        dRdx       = [dRdxi' dRdeta'] * invJacob;
        Gpt        = N*pts;     % GP in global coord
        
        globGP     = [globGP; Gpt]; % for plotting GPs only
        
        % +++++++++++++++++++++++++ 
        % Gradient of displacement
        % +++++++++++++++++++++++++ 
        
        % need to compute u,x u,y v,x v,y, stored in matrix H
        
        [B,J1] = BMatrixXIGA1(Xi,Eta,e,enrich_node,...
                             xCr,xtip,alpha,N,dRdxi,dRdeta);
        leB = size(B,2);
        
        % nodal displacement of current element
        % taken from the total nodal parameters U
        
        elemDisp = element_disp(e,pos,enrich_node,U);
        
        % compute derivatives of u w.r.t x and y
        
        H(1,1) = B(1,1:2:leB)*elemDisp(1:2:leB);    % u,x
        H(1,2) = B(2,2:2:leB)*elemDisp(1:2:leB);    % u,y
        H(2,1) = B(1,1:2:leB)*elemDisp(2:2:leB);    % v,x
        H(2,2) = B(2,2:2:leB)*elemDisp(2:2:leB);    % v,y
        
        % +++++++++++++++++++
        % Gradient of weight
        % +++++++++++++++++++ 
        
        weight  = qnode(iel,:);
        %gradq   = weight*dRdx;
        gradq   = weight*dNdx;

        % ++++++++++++++
        % Stress at GPs
        % ++++++++++++++ 
        
        epsilon = B*elemDisp ;
        sigma   = C*epsilon;
        
        % +++++++++++++++++++++++++++++++++++
        % Transformation to local coordinate
        % +++++++++++++++++++++++++++++++++++ 
        
        voit2ind    = [1 3;3 2];
        gradqloc    = QT*gradq';
        graddisploc = QT*H*QT';
        
        stressloc   = QT*sigma(voit2ind)*QT';
        
        epsilon(3)  = 0.5 * epsilon(3);
        strainloc   = QT*epsilon(voit2ind)*QT';

        % ++++++++++++++++++
        %  Auxiliary fields
        % ++++++++++++++++++ 
        
        xp    = QT*(Gpt-xtip)';           % local coordinates
        r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
        theta = atan2(xp(2),xp(1));

        K1 = 1.0 ;
        K2 = K1  ;

        mu = E0/(2.+ nu0 + nu0);
        
        if ( strcmp(stressState,'PLANE_STRAIN')  )
            kappa = 3-4*nu0;   
        else
            kappa = (3-nu0)/(1+nu0);   
        end
        
        SQR  = sqrt(r);
        CT   = cos(theta);
        ST   = sin(theta);
        CT2  = cos(theta/2);
        ST2  = sin(theta/2);
        C3T2 = cos(3*theta/2);
        S3T2 = sin(3*theta/2);

        drdx = CT;
        drdy = ST;
        dtdx = -ST/r;
        dtdy = CT/r;

        FACStress1 = sqrt(1/(2*pi));
        FACStress2 = FACStress1;

        FACDisp1 = sqrt(1/(2*pi))/(2*mu);
        FACDisp2 = FACDisp1;

        AuxStress   = zeros(2,2);
        AuxGradDisp = zeros(2,2);
        AuxEps      = zeros(2,2);

        for mode = 1:2            
            if     (mode == 1)           
               AuxiliaryFieldModeI          
            elseif (mode == 2)                
               AuxiliaryFieldModeII 
            end
            
            % +++++++++++++++
            %   J integral
            % +++++++++++++++
            
            I1= (stressloc(1,1) * AuxGradDisp(1,1) + stressloc(2,1) * AuxGradDisp(2,1) ) * gradqloc(1) + ...
                (stressloc(1,2) * AuxGradDisp(1,1) + stressloc(2,2) * AuxGradDisp(2,1) ) * gradqloc(2);

            I2= (AuxStress(1,1) * graddisploc(1,1) + AuxStress(2,1) * graddisploc(2,1) ) * gradqloc(1) + ...
                (AuxStress(2,1) * graddisploc(1,1) + AuxStress(2,2) * graddisploc(2,1) ) * gradqloc(2);

            StrainEnergy = 0;
            
            for i=1:2 %size(AuxEpsm1,1)
                for j=1:2  %size(AuxEpsm1,2)
                   % StrainEnergy = StrainEnergy +  stressloc(i,j)*AuxEps(i,j);
                   StrainEnergy = StrainEnergy +  AuxStress(i,j) * strainloc(i,j);
                end
            end          
            
            % Interaction integral I
            I(mode,1) = I(mode,1) + ...
                 (I1 + I2 - StrainEnergy*gradqloc(1))*J1*J2*wt;
        end   %loop on mode

    end       % of quadrature loop
end           % end of element loop

% plot
figure
hold on
% plot the circle
theta = -pi:0.1:pi;
xo = xtip(1) + radius*cos(theta) ;
yo = xtip(2) + radius*sin(theta) ;
plot(xo,yo,'k-');
plot_mesh(node,elementV,'Q4','b-',1.2)
plot_mesh(node,elementV(Jdomain,:),'Q4','r-',1.2)
cr = plot(xCr(:,1),xCr(:,2),'k-');
set(cr,'LineWidth',2);
cr = plot(globGP(:,1),globGP(:,2),'b*');
% -------------------------------------