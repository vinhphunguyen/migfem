%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extended Isogeometric analysis for two dimensional linear elastic
% fracture mechanics problems.
%
% Two edge cracks plate under tension
%
% Vinh Phu Nguyen,
% Johns Hopkins University
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util/
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/
addpath ../nurbs-geopdes/inst/
addpath ../nurbs-util/

clc
clear all

global p q controlPts weights element xCrack xTips crack_node jDomainFac

E0          = 1e3;  % Young modulus
nu0         = 0.3;  % Poisson ratio
stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
sigma0      = 1;

jDomainFac = 3;

vtuFile     = 'twoEdgeCracks';

% COMPUTE ELASTICITY MATRIX

if ( strcmp(stressState,'PLANE_STRESS')  )      % Plane Strain case
    C=E0/(1-nu0^2)*[  1      nu0          0;
        nu0        1          0;
        0        0  (1-nu0)/2  ];
else                                            % Plane Strain case
    C=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0      nu0        0;
        nu0    1-nu0        0;
        0        0  1/2-nu0 ];
end

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%twoEdgeCracksC2Data
%twoEdgeCracksC3Data
twoEdgeCracksKRefData

noGPs         = p+1; % # of Gauss points along one direction
noGP1         = p+1;

% find boundary nodes for boundary conditions

fixedYNodes  =  find(controlPts(:,2)==0); %bottom edge
fixedXNodes  =  1;
forcedNodes  =  find(controlPts(:,2)==D); %top edge

% build connectivity ...

generateIGA2DMesh

% build boundary mesh for force vector computation

bndPoints      = controlPts(forcedNodes,:);
rightEdgeMesh  = zeros(noElemsU,p+1);

for i=1:noElemsU
    rightEdgeMesh(i,:) = forcedNodes(i:i+p);
end

% Each split node is enriched by ONE function, H(x)
% Each tip node is enriched by FOUR functions, B_i(x)
% then, the total dofs is :
% total dof = numnode*nsd + numsplitnode*1*nsd + numstipnode*4*nsd
% here, two dimension, nsd = 2

noDofs = noDofs + size(split_nodes,1)*1*2 + ...
    size(tip_nodes,1)*4*2;

K = sparse(noDofs,noDofs);  % global stiffness matrix
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

% We use fictitious nodes/control points to handle these
% additional dofs. At a H(x) enriched node, we add one fantom node and
% at tip enriched node, four fantom nodes are added. These fictitious nodes
% are numbered from the total number of true nodes, ie, from numnode+1 ...

pos    = zeros(noCtrPts,1);
nsnode = 0 ;
ntnode = 0 ;

for i = 1 : noCtrPts
    if (enrich_node(i) == 1)
        pos(i) = (noCtrPts + nsnode*1 + ntnode*4) + 1 ;
        nsnode = nsnode + 1 ;
    elseif (enrich_node(i) == 2)
        pos(i) = (noCtrPts + nsnode*1 + ntnode*4) + 1 ;
        ntnode = ntnode + 1 ;
    end
end

jacob = zeros(2,2);

% essential boundary conditions

uFixed     = zeros(size(fixedXNodes))';
vFixed     = zeros(size(fixedYNodes))';

udofs      = fixedXNodes;
vdofs      = 2*fixedYNodes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])

% Loop over elements (knot spans)

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);          %  element scatter vector
    nn     = length(sctr);
    
    % -----------------------------------------------
    % Choose Gauss quadrature rules for elements
    
    if     (ismember(e,split_elem))     % split element
        [W,Q] = quadrature(12,'GAUSS',2);
    elseif (ismember(e,tip_elem))       % tip element
        [W,Q] = quadrature(12,'GAUSS',2);
    elseif ( any(intersect(tip_nodes,sctr)) ~= 0)% having tip enriched nodes
        [W,Q] = quadrature(10,'GAUSS',2);
    else
        [W,Q] = quadrature(noGPs,'GAUSS',2);
    end
    
    % Determine the position in the global matrix K
    
    sctrB = assembly(e,enrich_node,pos);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE,pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        % compute derivative of basis functions w.r.t parameter coord
        
        [N dRdxi dRdeta] = NURBS2DBasisDers([Xi; Eta],p,q,uKnot,vKnot,weights');
        
        % B matrix
        
        [B,J1] = BMatrixXIGA(Xi,Eta,e,enrich_node,N,dRdxi,dRdeta);
        
        % Stiffness matrix
        
        K(sctrB,sctrB) = K(sctrB,sctrB) + B'*C*B*W(gp)*J1*J2;
    end
end

% Computing external force

[W1,Q1] = quadrature(noGP1, 'GAUSS', 1 );

% Loop over elements along top edge

for e=1:noElemsU
    xiE   = elRangeU(e,:); % [xi_i,xi_i+1]
    conn  = elConnU(e,:);
    pts   = bndPoints(conn,:);
    
    sctry = 2*rightEdgeMesh(e,:);
    
    % loop over Gauss points
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
        J2      = 0.5 * ( xiE(2) - xiE(1) );

        [N dNdxi] = NURBS1DBasisDers(Xi,p,uKnot,weights);
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        jacob1   = dNdxi * pts;
        J1       = norm (jacob1);
        
        f(sctry) = f(sctry) + N' * sigma0 * J1 * J2 * wt;
    end
end

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])
bcwt=mean(diag(K)); % a measure of the average  size of an element in K
% used to keep the  conditioning of the K matrix

applyBC


% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ux = U(1:2:noDofs);
Uy = U(2:2:noDofs);

% compute stresses and displacements at mesh vertices
disp([num2str(toc),'  COMPUTING STRESSES '])

plotStressXIGA

% plot stress

stressComp=2;
figure
clf
plot_field(node,elementV,'Q4',stress(:,:,stressComp));
hold on
colorbar
title('Stress in x direction')
axis off
%plot_mesh(node,elementV,'Q4','g.-');

figure
clf
plot_field(node,elementV,'Q4',disp(:,:,2));
hold on
colorbar
title('Displacement in y direction')
axis off

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)

clear disp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([num2str(toc),'  SIFs COMPUTATION'])


if ( strcmp(stressState,'PLANE_STRAIN') )
    Eb = E0/(1-nu0^2);
else
    Eb = E0;
end

for iCr=1:noCracks
    
    xTip    = xTips(iCr,:);
    xCr     = reshape(xCrack(iCr,:,:),2,2);
    seg     = xCr(2,:) - xCr(1,:);   % tip segment
    alpha   = atan2(seg(2),seg(1));  % inclination angle
    QT      = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
    
    [Jdomain,qnode,radius] = jIntegrationDomain(tip_elem(iCr),xTip,node,elementV);
    
    % plot
    figure
    hold on
    % plot the circle
    theta = -pi:0.1:pi;
    xo = xTip(1) + radius*cos(theta) ;
    yo = xTip(2) + radius*sin(theta) ;
    plot(xo,yo,'k-');
    plot_mesh(node,elementV,'Q4','b-',1)
    plot_mesh(node,elementV(Jdomain,:),'Q4','r-',1)
    cr = plot(xCr(:,1),xCr(:,2),'k-');
    set(cr,'LineWidth',2);
    % -------------------------------------
    
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
        
        pts    = controlPts(sctr,:);
        
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
            
            [B,J1] = BMatrixXIGA(Xi,Eta,e,enrich_node,N,dRdxi,dRdeta);
            leB    = size(B,2);
            
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
            
            xp    = QT*(Gpt-xTip)';           % local coordinates
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
    
    % Compute SIFs from I integral
    
    Knum = 0.5*Eb*I
end % end of crack loop


% Compute the exact SIFs

L      = 1; % width of the plate
F      = 1.12 + 0.43*(a/L) - 4.79*(a/L)^2 + 15.46*(a/L)^3;
Kexact = F * sigma0*sqrt(pi*a)



