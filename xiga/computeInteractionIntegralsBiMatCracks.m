% -------------------------------------------------------------------------
% Compute the Stress Intensity Factors for bi-material cracks
% Using the Interaction integral method
% Some codes are copied from MXFEM written by Mathiew Pais

% Steps :
% 1- detection of the elements on which we integrate
% 2- loop over these elements
% 3- loop over Gauss points
% 4- computation of stress, strain... in local coordinates !!!   ATTENTION
% 5- computation of the auxilliary fields: AuxStress and AuxEps and AuxGradDisp
% 6- computation of I1 and I2

% Determine J domain and weight function
[Jdomain,qnode,radius] = jIntegrationDomain(tip_elem,xTip,node,elementV);

%---------------------------------------------
% Compute interaction integral


I  = zeros(2,1);

globGP=[]; % GPs in global coords. 

% ---------------------------
% Starting LOOP over ELEMENTS
%----------------------------

for iel = 1 : size(Jdomain,2)
    e      = Jdomain(iel) ; % current element
    sctr   = element(e,:);
    sctrQ4 = elementV(e,:);
    nn     = length(sctr);
    
    pts    = controlPts(sctr,:);
    levelS = CHI(1,sctr);
    
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    % Choose Gauss quadrature rule
    
    if     (ismember(e,splitElems))     % split element
        %[W,Q] = discontQ4quad(7,levelSets(1,sctrQ4,1));  
        [W,Q] = quadrature(16,'GAUSS',2);
    elseif (any(intersect(tip_nodes,sctr) ) ~= 0) || ...
           (any(intersect(iTip_nodes,sctr)) ~= 0) || ...
           (any(intersect(iTIp_nodes,sctr)) ~= 0) || ...
           (any(intersect(itip_nodes,sctr)) ~= 0) % having tip enriched nodes
        [W,Q] = quadrature(16,'GAUSS',2);
%     elseif (ismember(e,splitElems))     % split element by mat. interface
%         [W,Q] = discontQ4quad(7,chi(1,sctrQ4));
    else
        [W,Q] = quadrature(12,'GAUSS',2);
    end
    
    % nodal displacement of current element
    % taken from the total nodal parameters U
    
    elemDisp = element_disp(e,pos,enrich_node,U);
        
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
        Xi      = parent2ParametricSpace(xiE, pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping   (xiE,etaE);
        
        % compute derivative of basis functions w.r.t parameter coord
        
        [R dRdxi dRdeta] = NURBS2DBasisDers([Xi;Eta],p,q,uKnot,vKnot,weights');
        
        jacob      = pts'*[dRdxi' dRdeta'];
        J1         = det(jacob);
        invJacob   = inv(jacob);
        dRdx       = [dRdxi' dRdeta'] * invJacob;
        Gpt        = R*pts;     % GP in global coord
        
        globGP     = [globGP; Gpt]; % for plotting GPs only
        
        % +++++++++++++++++++++++++
        % Gradient of displacement
        % +++++++++++++++++++++++++
        
        % need to compute u,x u,y v,x v,y, stored in matrix H
        
        [B,J1] = BMatrixXIGA(Xi,Eta,e,enrich_node,R,dRdxi,dRdeta);
        %[B,J1] = BMatrixXIGAGen(Xi,Eta,e,enrich_node,R,dRdxi,dRdeta,N,dNdxi);
        
        leB    = size(B,2);
               
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
        
        levelset = dot(R,levelS);
        %levelset = Gpt(2) - y0;
        
        if levelset >= 0
            C = Cm;
            Eo = E1; vo = nu1; Go = G1; k = k1;
        else
            C = Ci;
            Eo = E2; vo = nu2; Go = G2; k = k2;
        end
        
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
        
        xp    = QT*(Gpt-xTip)';           % local coordinates
        r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
        theta = atan2(xp(2),xp(1));
        
        if ( theta > pi || theta < -pi)
            disp (['something wrong with angle ',num2str(thet)]);
        end
            
        % Constants (appendix B of Sukumar's IJNME paper), checked
        
        if ((theta >= 0) && (theta < pi)) ||...
                ((theta < -pi) && (theta >= -2*pi))
            d = exp(-(pi-theta)*ep);                                 % Constant, delta, upper half plane
        elseif ((theta >= pi) && (theta < 2*pi)) ||...
                ((theta < 0) && (theta >= -pi))
            d = exp((pi+theta)*ep);                                  % Constant, delta, lower half plane
        end
        
        lgr       = log(r);
        epLogr    = ep*lgr;
        temp      = 0.25+ep^2;
        cosEpLogr = cos(epLogr);
        sinEpLogr = sin(epLogr);
        
        BTA  = (0.5*cosEpLogr + ep*sinEpLogr)/temp;      % Constant, beta
        BTAP = (0.5*sinEpLogr - ep*cosEpLogr)/temp;      % Constant, beta prime
        GMA  = k*d-1/d;                                           % Constant, gamma
        GMAP = k*d+1/d;                                           % Constant, gamma prime
        PHI  = epLogr+theta/2;                                    % Constant, phi
        
        st   = sin(theta);
        ct   = cos(theta);
        st2  = sin(theta/2);
        ct2  = cos(theta/2);
        
        A = 1/(4*Go*cosh(ep*pi));                     % Constant, A
        B = sqrt(r/2/pi);                             % Constant, B
        C = BTAP*GMA *ct2 - BTA *GMAP*st2;            % Constant, C
        D = BTA *GMA *ct2 + BTAP*GMAP*st2;            % Constant, D
        E = BTAP*GMAP*ct2 - BTA *GMA *st2;            % Constant, E
        F = BTA *GMAP*ct2 + BTAP*GMA *st2;            % Constant, F
        
        dCdr =  ep*D/r;                               % Derivative of constant C with respect to r
        dCdt = -F/2 + ep*E;                           % Derivative of constant C with respect to theta
        dDdr = -ep*C/r;                               % Derivative of constant D with respect to r
        dDdt =  E/2 + ep*F;                           % Derivative of constant D with respect to theta
       
        
        T1 = 2*d*st*sin(PHI);                         % Constant, T1
        T2 = 2*d*st*cos(PHI);                         % Constant, T2
        T3 = 2*d*ct*sin(PHI);                         % Constant, T3
        T4 = 2*d*ct*cos(PHI);                         % Constant, T4
        
        dT1dr =  ep*T2/r;                             % Derivative of constant T1 with respect to r
        dT1dt =  ep*T1+T2/2+T3;                       % Derivative of constant T1 with respect to theta
        dT2dr = -ep*T1/r;                             % Derivative of constant T2 with respect to r
        dT2dt =  ep*T2-T1/2+T4;                       % Derivative of constant T2 with respect to theta
        
        drdx =  ct;                                   % Derivative of r with respect to x
        drdy =  st;                                         % Derivative of r with respect to y
        dtdx = -st/r;                                       % Derivative of theta with respect to x
        dtdy =  ct/r;
                      
        AuxStress   = zeros(2,2);
        AuxGradDisp = zeros(2,2);
        AuxStrain   = zeros(2,2);
        
        inv4piB     = 1/(4*pi*B);
        
        for mode = 1:2
            if mode == 1                                            % K1 = 1.0 and K2 = 0.0
                f1 =  D+T1;                                         % Function, f1
                f2 = -C-T2;                                         % Function, f2
                
                df1dr = dDdr+dT1dr;                                 % Derivative of function f1 with respect to r
                df1dt = dDdt+dT1dt;                                 % Derivative of function f1 with respect to theta
                df2dr = -dCdr-dT2dr;                                % Derivative of function f2 with respect to r
                df2dt = -dCdt-dT2dt;                                % Derivative of function f2 with respect to theta                           
            elseif mode == 2                                        % K1 = 0.0 and K2 = 1.0
                f1 = -C+T2;                                         % Function, f1
                f2 = -D+T1;                                         % Function, f2
                
                df1dr = -dCdr+dT2dr;                                % Derivative of f1 with respect to r
                df1dt = -dCdt+dT2dt;                                % Derivative of f1 with respect to theta
                df2dr = -dDdr+dT1dr;                                % Derivative of f2 with respect to t
                df2dt = -dDdt+dT1dt;                                % Derivative of f2 with respect to theta                         
            end
            
            df1dx = df1dr*drdx+df1dt*dtdx;                      % Derivative of f1 with respect to x
            df1dy = df1dr*drdy+df1dt*dtdy;                      % Derivative of f1 with respect to y
            df2dx = df2dr*drdx+df2dt*dtdx;                      % Derivative of f2 with respect to x
            df2dy = df2dr*drdy+df2dt*dtdy;                      % Deriva
                
            AuxGradDisp(1,1) = A*(B*df1dx+drdx*f1*inv4piB);    % Auxiliary displacement gradient of u1 with respect to x
            AuxGradDisp(1,2) = A*(B*df1dy+drdy*f1*inv4piB);    % Auxiliary displacement gradient of u1 with respect to y
            AuxGradDisp(2,1) = A*(B*df2dx+drdx*f2*inv4piB);    % Auxiliary displacement gradient of u2 with respect to x
            AuxGradDisp(2,2) = A*(B*df2dy+drdy*f2*inv4piB);    % Auxiliary displacement grad
                
            AuxStrain(1,1) = AuxGradDisp(1,1);
            AuxStrain(1,2) = 0.5*(AuxGradDisp(1,2)+AuxGradDisp(2,1));
            AuxStrain(2,1) = AuxStrain(1,2);
            AuxStrain(2,2) = AuxGradDisp(2,2);
            
            if strcmp(stressState,'PLANE_STRESS')
                AuxStress(1,1) = Eo/(1-vo^2)*(AuxStrain(1,1)+vo*AuxStrain(2,2));
                AuxStress(1,2) = Eo/(1+vo)*AuxStrain(1,2);
                AuxStress(2,1) = AuxStress(1,2);
                AuxStress(2,2) = Eo/(1-vo^2)*(vo*AuxStrain(1,1)+AuxStrain(2,2));
            else                                 % Plane strain
                AuxStress(1,1) = Eo/(1+vo)/(1-2*vo)*((1-vo)*AuxStrain(1,1)+vo*AuxStrain(2,2));
                AuxStress(1,2) = Eo/(1+vo)*AuxStrain(1,2);
                AuxStress(2,1) = AuxStress(1,2);
                AuxStress(2,2) = Eo/(1+vo)/(1-2*vo)*(vo*AuxStrain(1,1)+(1-vo)*AuxStrain(2,2));
            end
            
            I1= (stressloc(1,1) * AuxGradDisp(1,1) + stressloc(2,1) * AuxGradDisp(2,1) ) * gradqloc(1) + ...
                (stressloc(1,2) * AuxGradDisp(1,1) + stressloc(2,2) * AuxGradDisp(2,1) ) * gradqloc(2);
            
            I2= (AuxStress(1,1) * graddisploc(1,1) + AuxStress(2,1) * graddisploc(2,1) ) * gradqloc(1) + ...
                (AuxStress(2,1) * graddisploc(1,1) + AuxStress(2,2) * graddisploc(2,1) ) * gradqloc(2);
                        
            StrainEnergy = stressloc(1,1)*AuxStrain(1,1) + ...
                           stressloc(2,2)*AuxStrain(2,2) + ...
                           2*stressloc(1,2)*AuxStrain(1,2);
            
            % Interaction integral I
            
            I(mode,1) = I(mode,1) + (I1 + I2 - StrainEnergy*gradqloc(1))*J1*J2*wt;
            
        end   %loop on mode
        
    end       % of quadrature loop
end           % end of element loop

% Find the effective modulus

if strcmp(stressState,'PLANE_STRESS')
    Emeff = E1;
    Efeff = E2;
else                                      % Plane strain
    Emeff = E1/(1-nu1^2);
    Efeff = E2/(1-nu2^2);
end

Eeff = Emeff*Efeff/(Emeff+Efeff);

% Solve for the mixed-mode stress intensity factors
Kcalc   = I*Eeff*cosh(ep*pi)*cosh(ep*pi);            % Solve for Mode I and Mode II SIF
KI      = Kcalc(1);                                    % Mode I  SIF
KII     = Kcalc(2);


figure
hold on
%plot the circle
theta = -pi:0.1:pi;
xo = xTip(1) + radius*cos(theta) ;
yo = xTip(2) + radius*sin(theta) ;
plot(xo,yo,'k-');
plot_mesh(node,elementV,'Q4','b-')
plot_mesh(node,elementV(Jdomain,:),'Q4','r-')
cr = plot(xCr(:,1),xCr(:,2),'k-');
set(cr,'LineWidth',2);
cr = plot(globGP(:,1),globGP(:,2),'b*');


