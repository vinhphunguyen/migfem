% This file implements the Nitsche method to join two mechanical models:
% a 3D continuum model and a Mindlin plate model using the embedded mesh
% approach. The entire geometry is discretised by a plate model and
% overlaped by a 3D model on the clamped region.
%
% Discretisation: NURBS elements.
%
% Problem: Timoshenko beam in bending.
% 
% Vinh Phu Nguyen
% Cardiff University, UK
% 11 July 2013

addpath ../../fem_util/
addpath ../../gmshFiles/
addpath ../../post-processing/
addpath ../../fem-functions/
addpath ../../analytical-solutions/
addpath ../../meshing/
addpath ../../nurbs-util/
addpath ../../nurbs-geopdes/inst/
addpath ../../delamination/
addpath ../../C_files/

clear all
colordef white
state = 0;
tic;

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% MATERIAL PROPERTIES-----------------------------------------------------
E0  = 1000;  % Young modulus
nu0 = 0.3;   % Poisson ratio
k   = 5/6;   % shear correction factor 

% BEAM PROPERTIES
L     = 320;     % length of the beam
width = 80;
t     = 20;     % thicknes

% TIP LOAD
P = -10; % the peak magnitude of the traction at the right edge
q  = 0;    % uniformly distributed loads

% Constitutive matrices: continuum and plate
C=zeros(6,6);
C(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
                                  nu0 1-nu0 nu0;
                                  nu0 nu0 1-nu0];
C(4:6,4:6)=E0/2/(1+nu0)*eye(3);

D   = E0*t^3/(12*(1-nu0^2));
Cb  = D*[1 nu0 0;nu0 1 0; 0 0 0.5*(1-nu0)];
Cs  = E0*t*k/2/(1+nu0)*[1 0;0 1];

% penalty parameter

alpha = 5000;

%% GENERATE FINITE ELEMENT MESH
%
plotMesh  = 1;
disp([num2str(toc),'   GENERATING MESH'])

% MESH PROPERTIES

%% domain 1 ---------------------------------------------------------------
L1    = 160;

controlPts          = zeros(4,2,2,2);

controlPts(1:3,1,1,1) = [0;0;0];
controlPts(1:3,2,1,1) = [L1;0;0];
controlPts(1:3,1,2,1) = [0;width;0];
controlPts(1:3,2,2,1) = [L1;width;0];

controlPts(1:3,1,1,2) = [0;0;t];
controlPts(1:3,2,1,2) = [L1;0;t];
controlPts(1:3,1,2,2) = [0;width;t];
controlPts(1:3,2,2,2) = [L1;width;t];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];
wKnot = [0 0 1 1];

solid1 = nrbmak(controlPts,{uKnot vKnot wKnot});

% evaluate order

solid1 = nrbdegelev(solid1,[2 2 2]);

uKnot     = cell2mat(solid1.knots(1));
vKnot     = cell2mat(solid1.knots(2));

refineCountX = 5;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    newKnotsY = [];
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {newKnotsX [] []};
    solid1     = nrbkntins(solid1,newKnots);
    uKnot      = cell2mat(solid1.knots(1));
    vKnot      = cell2mat(solid1.knots(2));
end

solid1    = nrbkntins(solid1,{[] [0.25 0.5 0.75] [0.1 0.3 0.6 0.9]});
mesh1     = buildIGA3DMesh (solid1);

node1     = mesh1.controlPts;
element1  = mesh1.locElems;


numx1 = mesh1.noElemsU;
numy1 = mesh1.noElemsV;
numz1 = mesh1.noElemsW;

%% domain 2 (plate)--------------------------------------------------------

controlPts          = zeros(4,2,2);

controlPts(1:3,1,1) = [0;0;t/2];
controlPts(1:3,2,1) = [L;0;t/2];
controlPts(1:3,1,2) = [0;width;t/2];
controlPts(1:3,2,2) = [L;width;t/2];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

solid2 = nrbmak(controlPts,{uKnot vKnot});

% evaluate order

solid2 = nrbdegelev(solid2,[2 2]);

uKnot     = cell2mat(solid2.knots(1));
vKnot     = cell2mat(solid2.knots(2));
solid2    = nrbkntins(solid2,{[] [0.5] });
refineCountX = 5;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {newKnotsX []};
    solid2     = nrbkntins(solid2,newKnots);
    uKnot      = cell2mat(solid2.knots(1));
    vKnot      = cell2mat(solid2.knots(2));
end

mesh2     = buildIGA2DMeshForSurface (solid2);

numx2 = mesh2.noElemsU;
numy2 = mesh2.noElemsV;

node2     = mesh2.controlPts;
element2  = mesh2.locElems;
mesh2.globElems=mesh2.globElems+size(node1,1);

vMesh=buildVisualizationMesh2D(solid2);

%% plot mesh---------------------------------------------------------------

figure
hold on
nrbkntplot(solid1);
nrbkntplot(solid2);
plot3(mesh1.controlPts(:,1),mesh1.controlPts(:,2),mesh1.controlPts(:,3),'o','MarkerEdgeColor','k',...
    'MarkerFaceColor','r','MarkerSize',4.5);
plot3(mesh2.controlPts(:,1),mesh2.controlPts(:,2),mesh2.controlPts(:,3),'s','MarkerEdgeColor','k',...
    'MarkerFaceColor','b','MarkerSize',4.5);
axis off
%plot3(node1(fixedNode,1),node1(fixedNode,2),node1(fixedNode,3),'r*');
%    plot3(node2(rightNode,1),node2(rightNode,2),node2(rightNode,3),'rs');
%axis([0 L -c c])

%% determine cut elements, void elememnts and standard elements------------

linInterfaceGeo   = [L1 0; L1 width];

% Compute level sets

nNode = size(vMesh.node,1);
chi   = zeros(1,nNode);


x0   = linInterfaceGeo(1,1); y0 = linInterfaceGeo(1,2);
x1   = linInterfaceGeo(2,1); y1 = linInterfaceGeo(2,2);

for i = 1 : nNode
    x                  = vMesh.node(i,1);
    y                  = vMesh.node(i,2);    
    phi                = x - x0;
    chi(1,i) = phi;
end

elem2State = zeros(size(element2,1),1);

% loop over elements

for iel = 1 : size(element2,1)
    sctr    = element2(iel,:);
    sctrV   = vMesh.element(iel,:);
        
    phi  = chi(1,sctrV);
        
    if     ( max(phi)*min(phi) < 0 )
        elem2State(iel) = 1;             % cut elements
    elseif ( max(phi) < 0 )  
        elem2State(iel) = 0;             % void elements   
    else
        elem2State(iel) = 2;             % normal elements 
    end
end

voidElems = find(elem2State==0);
cutElems  = find(elem2State==1);

% find inactive plate dofs

inactiveNodes1 = [];
inactiveNodes2 = [];

for i=1:size(element2,1)
    sctr = element2(i,:);
    if     elem2State(i) == 0
        inactiveNodes1 = [inactiveNodes1;sctr];
    elseif elem2State(i) == 1
        inactiveNodes2 = [inactiveNodes2;sctr];
    end
end

inactiveNodes1 = unique(inactiveNodes1);
inactiveNodes2 = unique(inactiveNodes2);

inactivePlateNodes = setdiff(inactiveNodes1,inactiveNodes2);
inactivePlateNodes = inactivePlateNodes + size(node1,1);

      
%% boundary edges----------------------------------------------------------

bndNodes  = find(abs(node1(:,1)-L1)<1e-14);  

[bndMesh1,botIndex]     = surfaceMesh (mesh1.vKnot,mesh1.wKnot,bndNodes,mesh1.q,mesh1.r,...
    mesh1.noPtsY,mesh1.noPtsZ,mesh1.rangeV,mesh1.rangeW,mesh1.elConnV,mesh1.elConnW);


map = [];

for iz=1:numz1
    for iy=1:numy1
        map = [map;numx1*iy + numx1*numy1*(iz-1)];
    end
end

map2 = [];

for iz=1:numz1
    for iy=1:numy1
        if iy <= numy1/numy2
            map2 = [map2;1];
        else
            map2 = [map2;2];
        end
    end
end

bndMesh2 = find(elem2State==1);

%% GPs for coupling terms -------------------------------------------------

ngp = (mesh1.q+1) * (mesh1.q+1);

[W1,Q1]=quadrature( mesh1.q+1, 'GAUSS', 2 ); % two point quadrature

GP1 = [];
GP2 = [];

for e=1:size(bndMesh1,1)    
    be = map (e);
    pe = bndMesh2(map2(e));
    sctrS    = bndMesh1(e,:);    
    sctrV    = element1(be,:);    
    pts1     = node1(sctrV,:);    
    pts2     = node2(element2(pe,:),[1 2]);    
    ptsS     = node1(sctrS,:);
    
    idv    = botIndex(e,1);
    idw    = botIndex(e,2);
    
    etaE   = mesh1.rangeV(idv,:); % [eta_j,eta_j+1]
    zetaE  = mesh1.rangeW(idw,:); % [zeta_k,zeta_k+1]
    
    xiE1     = mesh1.rangeU(end,:);                  % [xi_i,xi_i+1]
    etE1     = mesh1.rangeV(mesh1.index(be,2),:);   % [et_j,eta_j+1]
    zetE1    = mesh1.rangeW(mesh1.index(be,3),:);   % [zet_j,zeta_j+1]
           
    xiE2     = mesh2.elRangeU(mesh2.index(pe,1),:);   % [xi_i,xi_i+1]
    etE2     = mesh2.elRangeV(mesh2.index(pe,2),:);   % [et_j,eta_j+1]
        
    for q=1:size(W1)
        pt=Q1(q,:);
        wt=W1(q);
 
        Xi      = parent2ParametricSpace(etaE, pt(1));
        Eta     = parent2ParametricSpace(zetaE,pt(2));
        J2      = jacobianPaPaMapping(etaE,zetaE);
        
        % compute derivatives of shape functions
        [R dRdxi dRdeta] = NURBS2DBasisDers([Xi; Eta],mesh1.q,mesh1.r,mesh1.vKnot,mesh1.wKnot,mesh1.weights');
        J0       = [dRdxi; dRdeta]*ptsS;
        a1       = J0(1,:);
        a2       = J0(2,:);
        a3       = cross(a1,a2);
        norma    = norm(a3); a3 = a3/norma;
        
        x   = R*ptsS;
        xP  = x([1 2]);
        X1  = global2LocalMapNURBS3D(x,pts1,xiE1,etE1,zetE1,mesh1.p,mesh1.q,mesh1.r,...
            mesh1.uKnot,mesh1.vKnot,mesh1.wKnot,mesh1.weights);                              
        X2  =  global2LocalMapNURBS2D(xP,pts2,xiE2,etE2,mesh2.p,mesh2.q,...
            mesh2.uKnot,mesh2.vKnot,mesh2.weights);    
        GP1 = [GP1;X1 wt*norma*J2 a3];
        GP2 = [GP2;X2 x(3)-t/2];
    end
end

%% GET NODES ON DISPLACEMENT BOUNDARY--------------------------------------
%      Here we get the nodes on the essential boundaries

fixedNode = find(node1(:,1)==0);  
rightNode = find(node2(:,1)==L); 


bndPoints   = node2(rightNode,:);

bndMesh = zeros(mesh2.noElemsV,mesh2.q+1);

for i=1:mesh2.noElemsV     
    bndMesh(i,:) = rightNode(i:i+mesh2.q);
end

fixedNode = [fixedNode;inactivePlateNodes];

uFixed    = zeros(1,length(fixedNode));  % a vector of u_x for the nodes
vFixed    = zeros(1,length(fixedNode));     
wFixed    = zeros(1,length(fixedNode)); 

udofs     = 3*fixedNode-2;
vdofs     = 3*fixedNode-1;
wdofs     = 3*fixedNode-0;

numnodes = size(node1,1) + size(node2,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM DATA STRUCTURES

disp([num2str(toc),' INITIALIZING DATA STRUCTURES'])

f=zeros(3*numnodes,1);          % external load vector
K=zeros(3*numnodes,3*numnodes);  % stiffness matrix

% ******************************************************************************
% ***                          P R O C E S S I N G                           ***
% ******************************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% COMPUTE STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'   COMPUTING STIFFNESS MATRIX'])

[W,Q]=quadrature( mesh1.p+1, 'GAUSS', 3 ); % 2x2x2 Gaussian quadrature

%% for domain 1 (continuum)

for e=1:size(element1,1)                          % start of element loop
    sctr            = element1(e,:);              % element scatter vector
    nn              = length(sctr);
    sctrB(1:3:3*nn) = 3*sctr-2;
    sctrB(2:3:3*nn) = 3*sctr-1;
    sctrB(3:3:3*nn) = 3*sctr-0;
    pts             = node1(sctr,:);
    
    idu    = mesh1.index(e,1);
    idv    = mesh1.index(e,2);
    idw    = mesh1.index(e,3);
    
    xiE    = mesh1.rangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = mesh1.rangeV(idv,:); % [eta_j,eta_j+1]
    zetaE  = mesh1.rangeW(idw,:); % [zeta_k,zeta_k+1]
    
    for q=1:size(W,1)
        pt=Q(q,:);
        wt=W(q);
        
        Xi      = parent2ParametricSpace(xiE,  pt(1));
        Eta     = parent2ParametricSpace(etaE, pt(2));
        Zeta    = parent2ParametricSpace(zetaE,pt(3));
        J2      = jacobianPaPaMapping3d (xiE,etaE,zetaE);
        % compute derivative of basis functions w.r.t parameter coord
        [N dRdxi dRdeta dRdzeta] = NURBS3DBasisDers([Xi;Eta;Zeta],...
            mesh1.p,mesh1.q,mesh1.r,mesh1.uKnot,mesh1.vKnot,mesh1.wKnot,mesh1.weights');
        jacob      = pts'*[dRdxi' dRdeta' dRdzeta'];
        J1         = det(jacob);
        invJacob   = inv(jacob);
        dRdx       = [dRdxi' dRdeta' dRdzeta'] * invJacob;        
        B          = getBmatrix3D(nn,dRdx);        
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * C * B * J1 * J2 * wt;
    end  % of quadrature loop
end

%% for domain 2 (plate)----------------------------------------------------
clear sctrB;

[W,Q]=quadrature(  mesh2.p+1, 'GAUSS', 2 ); 

for e=1:mesh2.noElems
    % skip void elements
    if elem2State(e) == 0
        continue
    end
    idu    = mesh2.index(e,1);
    idv    = mesh2.index(e,2);
    xiE    = mesh2.elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = mesh2.elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element2(e,:);          %  element scatter vector
    sctrg  = sctr + size(node1,1);
    nn     = length(sctr);
    pts    = node2(sctr,[1 2]);
            
    nn3    = nn*3;
    sctrw  = 3*sctr-2;
    sctrB  = zeros(1,nn3);
    
    sctrB(1:3:nn3) = 3*sctrg-2; % deflection
    sctrB(2:3:nn3) = 3*sctrg-1; % rotation 1
    sctrB(3:3:nn3) = 3*sctrg;   % rotation 2
    
    Bb     = zeros(3,3*nn);
    Bs     = zeros(2,3*nn);
    
    % cut element, special integration rule
    if elem2State(e) == 1
        [W,Q]=quadrature(  (mesh2.p+1)*3, 'GAUSS', 2 ); % 2 x2point quadrature
        x=zeros(size(W,1),1);
        for gp=1:size(W,1)
            pt      = Q(gp,:);    
            Xi      = parent2ParametricSpace(xiE,  pt(1));
            Eta     = parent2ParametricSpace(etaE, pt(2));
            [R dRdxi dRdeta ] = ...
            NURBS2DBasisDers([Xi; Eta],mesh2.p,mesh2.q,mesh2.uKnot,mesh2.vKnot,mesh2.weights');
            x = R*pts(:,1);
        end
        id = find(x<L1);
        W(id)=[];
        Q(id)=[];
    end
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE,pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        % shape functions, first derivatives w.r.t natural coords
        
        [R dRdxi dRdeta ] = ...
            NURBS2DBasisDers([Xi; Eta],mesh2.p,mesh2.q,mesh2.uKnot,mesh2.vKnot,mesh2.weights');
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        jacob  = [dRdxi; dRdeta]  * pts; % 2x2 matrix               
        J1     = det(jacob);
                   
        % Jacobian inverse and spatial 1st derivatives
        
        invJacob   = inv(jacob);
        dRdx       = invJacob*[dRdxi;dRdeta];        
        
        % bending and shear B matrices
        
        Bb(1,2:3:3*nn) = dRdx(1,:);
        Bb(2,3:3:3*nn) = dRdx(2,:);
        Bb(3,2:3:3*nn) = dRdx(2,:);
        Bb(3,3:3:3*nn) = dRdx(1,:);
        
        Bs(1,2:3:3*nn) = -R;
        Bs(1,1:3:3*nn) = dRdx(1,:);
        Bs(2,1:3:3*nn) = dRdx(2,:);
        Bs(2,3:3:3*nn) = -R;
        
        % compute elementary stiffness matrix and               
        K(sctrB,sctrB) = K(sctrB,sctrB) + ...
                         (Bb' * Cb * Bb + Bs' * Cs * Bs) * J1 * J2 * wt;                             
    end
end

%% interface integrals

Cplate  = zeros(5,5);

D   = E0/(1-nu0^2);
Cb  = D*[1 nu0 0;nu0 1 0; 0 0 0.5*(1-nu0)];
Cs  = E0*k/2/(1+nu0)*[1 0;0 1];

Cplate(1:3,1:3) = Cb;
Cplate(4:5,4:5) = Cs;

Csolid = C([1 2 4 5 6],:);

for i=1:size(bndMesh1,1)                     % start of element loop
    e1     = map (i);  
    e2     = bndMesh2(map2(i));  
    sctr1  = element1(e1,:);      
    sctr2  = element2(e2,:); 
    sctr2n = sctr2 + size(node1,1);
    
    nn1    = length(sctr1);
    nn2    = length(sctr2);
    
    sctrB1(1:3:3*nn1) = 3*sctr1-2;
    sctrB1(2:3:3*nn1) = 3*sctr1-1;
    sctrB1(3:3:3*nn1) = 3*sctr1-0;
    
    sctrB2(1:3:3*nn2) = 3*sctr2n-2;
    sctrB2(2:3:3*nn2) = 3*sctr2n-1;
    sctrB2(3:3:3*nn2) = 3*sctr2n-0;
                
    pts1 = node1(sctr1,:);
    pts2 = node2(sctr2,[1 2]);
    
    Kp11 = zeros(nn1*3,nn1*3);
    Kp12 = zeros(nn1*3,nn2*3);
    Kp22 = zeros(nn2*3,nn2*3);
    
    Kd11 = zeros(nn1*3,nn1*3);
    Kd12 = zeros(nn1*3,nn2*3);
    Kd21 = zeros(nn2*3,nn1*3);
    Kd22 = zeros(nn2*3,nn2*3);
    
    for q=ngp*(i-1)+1:ngp*(i-1)+ngp
        pt1    = GP1(q,1:3);
        pt2    = GP2(q,1:2);
        wt1    = GP1(q,4);
        y      = GP2(q,3);
        normal = GP1(q,5:end);
        n = [normal(1) 0 normal(2) 0 normal(3);...
             0 normal(2) normal(1) normal(3) 0;...
             0 0 0 normal(2) normal(1)];
        n=-n;
        
        [N1 dN1dxi dN1deta dN1dzeta] = NURBS3DBasisDers(pt1,...
            mesh1.p,mesh1.q,mesh1.r,mesh1.uKnot,mesh1.vKnot,mesh1.wKnot,mesh1.weights');
        [N2 dN2dxi dN2deta] = NURBS2DBasisDers(pt2,...
            mesh2.p,mesh2.q,mesh2.uKnot,mesh2.vKnot,mesh2.weights');
        
        J1      = pts1'*[dN1dxi' dN1deta' dN1dzeta'];
        J2      = pts2'*[dN2dxi' dN2deta'];
        
        dN1dx   = [dN1dxi' dN1deta' dN1dzeta']*inv(J1);
        dN2dx   = [dN2dxi' dN2deta']*inv(J2);
        
        B1      = getBmatrix3D(nn1,dN1dx);
        
        B2(1,2:3:3*nn2)  = -y*dN2dx(:,1)';
        B2(2,3:3:3*nn2)  = -y*dN2dx(:,2)';
        B2(3,2:3:3*nn2)  = -y*dN2dx(:,2)';
        B2(3,3:3:3*nn2)  = -y*dN2dx(:,1)';
        B2(4,1:3:3*nn2)  =    dN2dx(:,2)';
        B2(4,3:3:3*nn2)  =   -N2;
        B2(5,1:3:3*nn2)  =    dN2dx(:,1)';
        B2(5,2:3:3*nn2)  =   -N2;
        
        % N matrices
        Nm1(1,1:3:3*nn1)  = N1';
        Nm1(2,2:3:3*nn1)  = N1';
        Nm1(3,3:3:3*nn1)  = N1';
        
        Nm2(1,2:3:3*nn2)  = -y*N2';
        Nm2(2,3:3:3*nn2)  = -y*N2';
        Nm2(3,1:3:3*nn2)  =    N2';
        
        % COMPUTE ELEMENT STIFFNESS AT QUADRATURE POINT
        
        Kp11 = Kp11 + alpha*(Nm1'*Nm1)*wt1;
        Kp12 = Kp12 + alpha*(Nm1'*Nm2)*wt1;
        Kp22 = Kp22 + alpha*(Nm2'*Nm2)*wt1;
        
        Kd11 = Kd11 + 0.5 * Nm1'* n * Csolid * B1 *wt1;
        Kd12 = Kd12 + 0.5 * Nm1'* n * Cplate * B2 *wt1;
        Kd21 = Kd21 + 0.5 * Nm2'* n * Csolid * B1 *wt1;
        Kd22 = Kd22 + 0.5 * Nm2'* n * Cplate * B2 *wt1;
               
    end  % of quadrature loop
   
    K(sctrB1,sctrB1)  = K(sctrB1,sctrB1)  + Kd11 + Kd11' + Kp11;    
    K(sctrB1,sctrB2)  = K(sctrB1,sctrB2)  + Kd12 - Kd21' - Kp12;
    K(sctrB2,sctrB1)  = K(sctrB2,sctrB1)  - Kd21 + Kd12' - Kp12';
    K(sctrB2,sctrB2)  = K(sctrB2,sctrB2)  - Kd22 - Kd22' + Kp22;
end    

%% External force

[W1,Q1]=quadrature( mesh2.q+1, 'GAUSS', 1 ); 

for e=1:mesh2.noElemsV
    xiE   = mesh2.elRangeV(e,:); % [xi_i,xi_i+1]
    conn  = mesh2.elConnV(e,:);    
    sctr  = bndMesh(e,:);           % element scatter vector
    sctrg = sctr + size(node1,1);
    sctry = 3*sctrg-2;
    pts   = bndPoints(conn,1:2);
    
    % loop over Gauss points
    
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
        J2      = 0.5 * ( xiE(2) - xiE(1) );
        
        [N dNdxi] = NURBS1DBasisDers(Xi,mesh2.q,mesh2.vKnot,mesh2.weights);        
        jacob1   = dNdxi*pts;
        J1       = norm (jacob1);                        
        f(sctry) = f(sctry) + N' * P * J1 * J2 * wt;
    end
end  

%%%%%%%%%%%%%%%%%%% END OF STIFFNESS MATRIX COMPUTATION %%%%%%%%%%%%%%%%%%%


% APPLY ESSENTIAL BOUNDARY CONDITIONS
disp([num2str(toc),'   APPLYING BOUNDARY CONDITIONS'])

[K,f]=applyDirichletBCs3D(K,f,udofs,vdofs,wdofs,uFixed,vFixed,wFixed);

% SOLVE SYSTEM
disp([num2str(toc),'   SOLVING SYSTEM'])
U=K\f;

%******************************************************************************
%*** POST - PROCESSING *** 
%***************************************************

%%

disp([num2str(toc),'   POST-PROCESSING'])

numnode1 = size(node1,1);
numnode2 = size(node2,1);

index1   = 1:numnode1;
index2   = numnode1+1:numnodes;

Ux1      = U(3*index1-2);
Uy1      = U(3*index1-1);
Uz1      = U(3*index1-0);

Ux2      = U(3*index2-2);
Uy2      = U(3*index2-1);
Uz2      = U(3*index2-0);

Uz2(:) = 0;
Uy2(:) = 0;

%% write to VTK
%
% stresses=[sigmaXX sigmaYY sigmaZZ sigmaXY sigmaYZ sigmaZX sigmaVM]

vtuFile1    = 'nitsche-iga-plate3D-continuumEmbed';
vtuFile2    = 'nitsche-iga-plate3D-plateEmbed';
vtuFileName = 'nitsche-iga-plate3DEmbed';

displacement = [U(1:3:end) U(2:3:end) U(3:3:end)];
damage       = zeros(length(U),1);

matMap1 = ones(mesh1.noElems,1);
matMap2 = ones(mesh2.noElems,1);
vMesh1=buildVisualizationMesh3D(solid1);
vMesh2=buildVisualizationMesh2D(solid2);
meshes{1}  = mesh1;
meshes{2}  = mesh2;
vmeshes{1} = vMesh1;
vmeshes{2} = vMesh2;
meshData.mesh    = meshes;
meshData.vmesh   = vmeshes;
meshData.matMap{1}=matMap1;
meshData.matMap{2}=matMap2;
material.stiffMat=C;
materials{1}      = material;
meshData.voids    = voidElems;

%figure; hold on;
ok      = plotStress3DForPatch   (meshData,1,vtuFile1,displacement,damage,materials);
ok      = plotStressShellForPatch(meshData,2,vtuFile2,displacement,damage,materials);

pvdFile = fopen(strcat(vtuFileName,'.pvd'), 'wt');

fprintf(pvdFile,'<VTKFile byte_order="LittleEndian" type="Collection" version="0.1">\n');
fprintf(pvdFile,'<Collection>\n');

vtuFile = sprintf('%s%s',vtuFile1,'.vtu');
fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''%d'' timestep=''%d''/>\n',...
    vtuFile,0,1);

vtuFile = sprintf('%s%s',vtuFile2,'.vtu');
fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''%d'' timestep=''%d''/>\n',...
    vtuFile,0,1);


fprintf(pvdFile,'</Collection>\n');
fprintf(pvdFile,'</VTKFile>\n');

fclose(pvdFile);



