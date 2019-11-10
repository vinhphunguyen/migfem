%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for Kirchoff-Love shell problems.
%
% Rotation-free thin shells. Fully clamped or simply supported 
% rectangular plates. Used for verifying the implementation.
%
% Vinh Phu Nguyen,
% Cardiff University
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util/
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../meshing/
addpath ../nurbs-util/
addpath ../nurbs-geopdes/inst/

clc
clear all

global p q

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

squareShellData

% constitutive matrix

memStiff = E*t/(1-nu^2);
benStiff = E*t^3/12/(1-nu^2);

%C  = [1 nu 0;nu 1 0; 0 0 0.5*(1-nu)];

% find boundary nodes for boundary conditions

EPS = 1e-8;
bottomNodes  =  find(abs(controlPts(:,2))  <EPS);
topNodes     =  find(abs(controlPts(:,2)-b)<EPS);
leftNodes    =  find(abs(controlPts(:,1))  <EPS);
rightNodes   =  find(abs(controlPts(:,1)-a)<EPS);

fixedNodes  =  unique([bottomNodes;rightNodes;topNodes;leftNodes]);

% if clamped plate, then fixing two rows of control points.
if clamped
    nextToBotNodes = noPtsX+2:2*noPtsX-1;
    nextToRgtNodes = 2*noPtsX-1:noPtsX:noPtsX*(noPtsY-1)-1;
    nextToTopNodes = noPtsX*(noPtsY-2)+2:noPtsX*(noPtsY-1)-1;
    nextToLefNodes = noPtsX+2:noPtsX:noPtsX*(noPtsY-2)+2;
    
    nextNodes      = unique([nextToBotNodes';nextToRgtNodes';...
                             nextToTopNodes';nextToLefNodes']);
    
    fixedNodes     = [fixedNodes; nextNodes(:)];
end

% build connectivity ...

generateIGA2DMesh

noCtrPts       = noPtsX   * noPtsY;
noDofs         = noCtrPts * 3;   % three displacement dofs per node

% initialization

K = zeros(noDofs,noDofs);  % global stiffness matrix
f = zeros(noDofs,1);        % external force vector

% essential boundary conditions

udofs      = 3*fixedNodes-2;
vdofs      = 3*fixedNodes-1;
wdofs      = 3*fixedNodes-0;


uFixed = zeros(size(fixedNodes))';
vFixed = zeros(size(fixedNodes))';
wFixed = zeros(size(fixedNodes))';

% udofs=[];
% vdofs=[];
% uFixed = [];
% vFixed = [];
%% fast assembly using the triple sparse matrix

nElNod = size(element,2);
nElDof = nElNod*3;
nElmLK = nElDof^2;
nSprGK = nElmLK*noElems;

jSprLK = 1:nElmLK;
vIdGrd = ones(1,nElDof);

% values, row indices, columns indices of the global K matrix
vSprGK = zeros(nSprGK,1);
jSprRw = zeros(nSprGK,1);
jSprCl = zeros(nSprGK,1);

Ke0    = zeros(nElDof,nElDof); % element Ke

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Gauss quadrature rule
noGPs = p+1;
noGpEle = noGPs^2;
[W,Q]=quadrature(  noGPs, 'GAUSS', 2 ); % noGPs x noGPs point quadrature

% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM']);

% Loop over elements (knot spans)

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);          %  element connectivity
    nn     = length(sctr);
    pts    = controlPts(sctr,:);
    
    nn3      = 3*nn;
    sctrB = zeros(1,nn3);
    
    sctrB(1:3:nn3) = 3*sctr-2;
    sctrB(2:3:nn3) = 3*sctr-1;
    sctrB(3:3:nn3) = 3*sctr;
    
    sctrf          = 3*sctr; % scatter for distributed force 

    Ke = Ke0;
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE, pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        % shape functions, first and second derivatives w.r.t natural coords
        
        [R dRdxi dRdeta dR2dxi dR2det dR2dxe] = ...
            NURBS2DBasis2ndDers([Xi; Eta],p,q,uKnot,vKnot,weights');
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        jacob  = [dRdxi; dRdeta]          * pts; % 2x3 matrix
        jacob2 = [dR2dxi; dR2det; dR2dxe] * pts; % 3x2 matrix
                                          
        % a1, a2 and a3 vectors (surface basis vectors)
        % and its derivatives
  
        a1    = jacob(1,:);
        a2    = jacob(2,:);
        a3    = cross(a1,a2); 
        norma = norm(a3);
        inorma = 1/norma;
        a3    = a3*inorma; J1    = norma;
        
        a11   = jacob2(1,:);
        a22   = jacob2(2,:);
        a12   = jacob2(3,:);
        
        % dot products of ai and ei
        
        a1e1  = a1(1); a1e2  = a1(2); a1e3  = a1(3);
        a2e1  = a2(1); a2e2  = a2(2); a2e3  = a2(3);
        
        % R_I,2*a1 + R_I,1*a2 for all shape functions
        
        noBasis = length(R);
        dRIa    = zeros(3,noBasis);
        for i=1:noBasis
          dRIa(:,i) = dRdeta(i)*a1 + dRdxi(i)*a2;
        end
        
        % compute the constitutive matrix C
        a_11 = dot(a1,a1);
        a_12 = dot(a1,a2);
        a_21 = dot(a2,a1); 
        a_22 = dot(a2,a2);
                
        au11 =   a_22 / ( a_11 * a_22 - a_21 * a_12 );
        au12 = - a_12 / ( a_11 * a_22 - a_21 * a_12 );
        au22 =   a_11 / ( a_11 * a_22 - a_21 * a_12 );
        
        C = [au11^2 nu*au11*au22+(1-nu)*au12^2 au11*au12;
             nu*au11*au22+(1-nu)*au12^2 au22^2 au22*au12;
             au11*au12 au22*au12 0.5*((1-nu)*au11*au22+(1+nu)*au12^2)];
        
        % membrane and bending B matrices
        
        a11xa2 = cross(a11,a2);
        a12xa2 = cross(a12,a2);
        a1xa11 = cross(a1,a11);
        a1xa22 = cross(a1,a22);
        a1xa12 = cross(a1,a12);
        a2xa3  = cross(a2,a3);
        a3xa1  = cross(a3,a1);
        a22xa2 = cross(a22,a2);

        a3a22  = dot(a3,a22);
        a3a11  = dot(a3,a11);
        a3a12  = dot(a3,a12); 
        
        Bmem = zeros(3,noBasis*3);
        Bben = zeros(3,noBasis*3);
        for i = 1:noBasis
            dRIdx = dRdxi (i);
            dRIdy = dRdeta(i);
            
            id    = (i-1)*3+1:3*i;
            
            Bmem(:,id)=[dRIdx*a1e1 dRIdx*a1e2 dRIdx*a1e3;
                        dRIdy*a2e1 dRIdy*a2e2 dRIdy*a2e3;
                        dRIa(1,i)  dRIa(2,i)  dRIa(3,i)];
            
            BI1 = -dR2dxi(i)*a3 + inorma*(dRIdx*a11xa2 + dRIdy*a1xa11 + ...
                          a3a11*(dRIdx*a2xa3 + dRIdy*a3xa1)); 
                      
            BI2 = -dR2det(i)*a3 + inorma*(dRIdx*a22xa2 + dRIdy*a1xa22 + ...
                          a3a22*(dRIdx*a2xa3 + dRIdy*a3xa1)); 
                      
            BI3 = -dR2dxe(i)*a3 + inorma*(dRIdx*a12xa2 + dRIdy*a1xa12 + ...
                          a3a12*(dRIdx*a2xa3 + dRIdy*a3xa1)); 
                      
            Bben(:,id)=[BI1;BI2;2*BI3];                       
        end
        
        % compute elementary stiffness matrix 
        
        Ke = Ke + memStiff * Bmem' * C * Bmem * J1 * J2 * wt + ... 
                  benStiff * Bben' * C * Bben * J1 * J2 * wt ;
                      
        f(sctrf)      = f(sctrf)      + q0 * R' * J1 * J2 * wt;
    end
    
    %----------------------------------------------------------------------
    % Global Assembly (Sparse)
    %----------------------------------------------------------------------
    
    vElDof = sctrB';
    mRwGrd = vElDof(:,vIdGrd);
    mClGrd = mRwGrd';
    
    jSprRw(jSprLK) = mRwGrd(:);
    jSprCl(jSprLK) = mClGrd(:);
    vSprGK(jSprLK) = Ke(:);
    
    jSprLK         = jSprLK + nElmLK; % move to the next position 
    
end

% Here comes the total stiffness matrix in one shot!!!

K = sparse(jSprRw,jSprCl,vSprGK,noDofs,noDofs);

%%
disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS']);

[K,f]=applyDirichletBCs3D(K,f,udofs,vdofs,wdofs,uFixed,vFixed,wFixed);

% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM']);
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([num2str(toc),'  POST-PROCESSING']);

%% Visualization using a Q4 visualization mesh

buildVMeshShell;
%plot_mesh(node,elementV,'Q4','g.-',1);

disp   = zeros(noElems,size(elementV,2),3);

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);         %  element scatter vector
    pts    = controlPts(sctr,:);
       
    sctrUx = 3*sctr-2;
    sctrUy = 3*sctr-1;
    sctrUz = 3*sctr;
    
    elemDisp = [U(sctrUx) U(sctrUy) U(sctrUz)];
    
    uspan = FindSpan(noPtsX-1,p,xiE(1), uKnot);
    vspan = FindSpan(noPtsY-1,q,etaE(1),vKnot);
    
    % loop over Gauss points
    
    gp = 1;
    for iv=1:2
        if (iv==2)
            xiE = sort(xiE,'descend');
        end
        for iu=1:2
            Xi  = xiE(iu);
            Eta = etaE(iv);
            [N dRdxi dRdeta] = NURBS2DBasisDersSpecial([Xi; Eta],...
                p,q,uKnot,vKnot,weights',[uspan;vspan]);
            
            disp(e,gp,:)    = N * elemDisp;
            
            gp = gp +1;
        end
    end
end

X = zeros(4,noElemsV);
Y = zeros(4,noElemsV);
Z = zeros(4,noElemsV);
C0 = zeros(4,noElemsV);

component = 3; % 1 and 2 are other options
for i = 1:size(elementV,1)
    sctr   = elementV(i,:);
    X(:,i) = node(sctr,1);
    Y(:,i) = node(sctr,2);
    Z(:,i) = node(sctr,3);
    C0(:,i) = disp(i,:,component);
end

figure
fill3(X,Y,C0,C0);
colorbar
title('deflection')
axis on



opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)


wbar = min(C0(:))*benStiff*1000/q0/a^4
wext = 1.26532  % CCCC
%wext = 4.06235; % SSSS








