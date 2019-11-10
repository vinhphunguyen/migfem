%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for Kirchoff-Love shell problems.
%
% ScordelisLo shell problem.
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

scordelisLoData

% constitutive matrix

memStiff = E*t/(1-nu^2);
benStiff = E*t^3/12/(1-nu^2);

%% Dirichlet BCs (symmetry conditions)
%   z
%   | 
% (D)------- (C)
%   |      |
%   |    L | 
%   |  R   | 
% (A)------- (B) --->x 
%
% rigid diaphram: AB: u_x, u_y = 0
% Symmetry conditions on CD and AD
% BC is free  

% boundary nodes
EPS = 1e-8;
nodesOnAB  = find(abs(controlPts(:,3))<EPS)'; 
nodesOnCD  = find(abs(controlPts(:,3)-L)<EPS)';
nodesOnAD  = find(abs(controlPts(:,1))<EPS)';

% nodes (control points) right next to boundary nodes
nextToADNodes = noPtsX-1:noPtsX:noPtsX*(noPtsY)-1;
nextToCDNodes = noPtsX*(noPtsY-2)+1:noPtsX*(noPtsY-1);

xConsNodes = unique([nodesOnAB nodesOnAD]);
yConsNodes = unique([nodesOnAB]);                
zConsNodes = unique([nodesOnCD]);

uFixed     = zeros(size(xConsNodes));
vFixed     = zeros(size(yConsNodes));
wFixed     = zeros(size(zConsNodes));

udofs      = 3*xConsNodes-2; % global indecies  of the fixed x disps
vdofs      = 3*yConsNodes-1; % global indecies  of the fixed y disps
wdofs      = 3*zConsNodes;   % global indecies  of the fixed z disps

% build connectivity ...

generateIGA2DMesh

noCtrPts       = noPtsX   * noPtsY;
noDofs         = noCtrPts * 3;   % three displacement dofs per node


%%
figure 
hold on
nrbkntplot (solid)
nrbctrlplot(solid)
view(3)

plot3(controlPts(xConsNodes,1),controlPts(xConsNodes,2),controlPts(xConsNodes,3),'bs');
plot3(controlPts(yConsNodes,1),controlPts(yConsNodes,2),controlPts(yConsNodes,3),'cs');
plot3(controlPts(zConsNodes,1),controlPts(zConsNodes,2),controlPts(zConsNodes,3),'ys');

%% initialization

K = zeros(noDofs,noDofs);  % global stiffness matrix
f = zeros(noDofs,1);        % external force vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Gauss quadrature rule
noGPs = 5;
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
    
    nn3   = nn*3;
    sctrB = zeros(1,nn3);
    
    sctrB(1:3:nn3) = 3*sctr-2;
    sctrB(2:3:nn3) = 3*sctr-1;
    sctrB(3:3:nn3) = 3*sctr;
    
    sctrf          = 3*sctr-1; % scatter for distributed force (y dir)

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
        
        jacob  = [dRdxi; dRdeta]          * pts; % 2x2 matrix
        jacob2 = [dR2dxi; dR2det; dR2dxe] * pts; % 3x2 matrix
                       
        dxdxi = jacob(1,1); dydxi = jacob(1,2);
        dxdet = jacob(2,1); dydet = jacob(2,2);
        
        j33   = [dxdxi^2     dydxi^2     2*dxdxi*dydxi;
            dxdet^2     dydet^2     2*dxdet*dydet;
            dxdxi*dxdet dydxi*dydet dxdxi*dydet+dxdet*dydxi];
        
        % Jacobian inverse and spatial 1st and 2nd derivatives
        
%         invJacob   = inv(jacob);
%         dRdx       = invJacob*[dRdxi;dRdeta];
%         dR2dx      = inv(j33)*([dR2dxi; dR2det; dR2dxe]-jacob2*dRdx);
        
        % a1, a2 and a3 vectors (surface basis vectors)
        % and its derivatives
        
        a1    = jacob(1,:);
        a2    = jacob(2,:);
        a3    = cross(a1,a2); 
        norma = norm(a3);
        a3    = a3/norma; J1    = norma;
        
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
        a_11 = dot(a1,a1); a_12 = dot(a1,a2);
        a_21 = dot(a2,a1); a_22 = dot(a2,a2);
        
        aa1 = [a_11 a_21; a_12 a_22 ] \ [1;0];
        aa2 = [a_11 a_21; a_12 a_22 ] \ [0;1];
        
        au11 = aa1(1);
        au12 = aa1(2);
        au22 = aa2(2);
        
        C = [au11^2 nu*au11*au22+(1-nu)*au12^2 au11*au12;
             nu*au11*au22+(1-nu)*au12^2 au22^2 au22*au12;
             au11*au12 au22*au12 0.5*((1-nu)*au11*au22+(1+nu)*au12^2)];
         
        % membrane and bending B matrices
        
        Bmem = zeros(3,noBasis*3);
        Bben = zeros(3,noBasis*3);
        for i = 1:noBasis
            dRIdx = dRdxi (i);
            dRIdy = dRdeta(i);
            
            id    = (i-1)*3+1:3*i;
            
            Bmem(:,id)=[dRIdx*a1e1 dRIdx*a1e2 dRIdx*a1e3;
                        dRIdy*a2e1 dRIdy*a2e2 dRIdy*a2e3;
                        dRIa(1,i)  dRIa(2,i)  dRIa(3,i)];
            
            BI1 = -dR2dxi(i)*a3 + 1/norma*(dRIdx*cross(a11,a2) + dRIdy*cross(a1,a11) + ...
                          dot(a3,a11)*(dRIdx*cross(a2,a3) + dRIdy*cross(a3,a3))); 
                      
            BI2 = -dR2det(i)*a3 + 1/norma*(dRIdx*cross(a22,a2) + dRIdy*cross(a1,a22) + ...
                          dot(a3,a22)*(dRIdx*cross(a2,a3) + dRIdy*cross(a3,a3))); 
                      
            BI3 = -dR2dxe(i)*a3 + 1/norma*(dRIdx*cross(a12,a2) + dRIdy*cross(a1,a12) + ...
                          dot(a3,a12)*(dRIdx*cross(a2,a3) + dRIdy*cross(a3,a3))); 
                      
            Bben(:,id)=[BI1;BI2;2*BI3];                       
        end
        
        % compute elementary stiffness matrix and
        % assemble it to the global matrix
        
        K(sctrB,sctrB) = K(sctrB,sctrB) + ...
                          memStiff * Bmem' * C * Bmem * J1 * J2 * wt + ... 
                          benStiff * Bben' * C * Bben * J1 * J2 * wt ;
         
         % volume force g            
        f(sctrf)      = f(sctrf)      + g * R' * J1 * J2 * wt;
    end
end

%% external force

%% Enforcing symmetry BCs

w     = 1e7;
penaltyStiffness = w*[1 -1;-1 1];

for i=1:length(nodesOnCD)
    sctr  = [nodesOnCD(i) nextToCDNodes(i)];
    sctrx = 3*sctr-2;
    sctry = 3*sctr-1;
    sctrz = 3*sctr-0;
    
    K(sctrx,sctrx) = K(sctrx,sctrx) + penaltyStiffness;
    K(sctry,sctry) = K(sctry,sctry) + penaltyStiffness;
    K(sctrz,sctrz) = K(sctrz,sctrz) + penaltyStiffness;
end

for i=1:length(nodesOnAD)
    sctr  = [nodesOnAD(i) nextToADNodes(i)];
    sctrx = 3*sctr-2;
    sctry = 3*sctr-1;
    sctrz = 3*sctr-0;
    
    K(sctrx,sctrx) = K(sctrx,sctrx) + penaltyStiffness;
    K(sctry,sctry) = K(sctry,sctry) + penaltyStiffness;
    K(sctrz,sctrz) = K(sctrz,sctrz) + penaltyStiffness;
end

%%
disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS']);
bcwt=mean(diag(K)); % a measure of the average  size of an element in K
% used to keep the  conditioning of the K matrix

% modify the force vector
f        = f-K(:,udofs)*uFixed';  
f        = f-K(:,vdofs)*vFixed';  
f        = f-K(:,wdofs)*wFixed';  
f(udofs) = bcwt*uFixed;
f(vdofs) = bcwt*vFixed;
f(wdofs) = bcwt*wFixed;

% modify the stiffness matrix
K(udofs,:)    = 0;  K(:,udofs)    = 0;
K(vdofs,:)    = 0;  K(:,vdofs)    = 0;
K(wdofs,:)    = 0;  K(:,wdofs)    = 0;

K(udofs,udofs)=bcwt*speye(length(udofs));  
K(vdofs,vdofs)=bcwt*speye(length(vdofs)); 
K(wdofs,wdofs)=bcwt*speye(length(wdofs)); 

% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM']);
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vertical displacement at C
% The exact value is: 0.3024

U(3*nodesOnCD(1)-1)

disp([num2str(toc),'  POST-PROCESSING']);

%% Visualization using a Q4 visualization mesh

buildVMeshShell;
plot_mesh(node,elementV,'Q4','g.-',1);
view(3)

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
Ux = zeros(4,noElemsV);
Uy = zeros(4,noElemsV);
Uz = zeros(4,noElemsV);

for i = 1:size(elementV,1)
    sctr   = elementV(i,:);
    X(:,i) = node(sctr,1);
    Y(:,i) = node(sctr,2);
    Z(:,i) = node(sctr,3);
    Ux(:,i) = disp(i,:,1);
    Uy(:,i) = disp(i,:,2);
    Uz(:,i) = disp(i,:,3);
end

factor=1e6;
figure
fill3(X+factor*Ux,Y+factor*Uy,Z+factor*Uz,Uy);
colorbar
title('deflection')
axis on

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)


%%

dispX = zeros(size(node,1),1);
dispY = zeros(size(node,1),1);

for e=1:size(elementV,1)
    connect = elementV(e,:);
    for in=1:4
        nid = connect(in);
        
        dispX(nid) = disp(e,in,1);
        dispY(nid) = disp(e,in,2);
        dispZ(nid) = disp(e,in,3);
    end
end

dim      = size(node,2);
numNodes = size(node,1);
numCells = size(elementV,1);
dof      = 3;
x        = node;
connect  = elementV;

% Output files

outfileVTU  = '../results/scordelisLoShell.vtu';
results_vtu = fopen(outfileVTU, 'wt');

numVertexesPerCell = 4;
VTKCellCode        = 9;
dof_per_vertex     = 3;


%% Write headers
fprintf(results_vtu, '<VTKFile type="UnstructuredGrid"  version="0.1"   > \n');
fprintf(results_vtu, '<UnstructuredGrid> \n');
fprintf(results_vtu, '<Piece  NumberOfPoints="  %g" NumberOfCells=" %g"> \n', numNodes, numCells);

%% Write point data
fprintf(results_vtu, '<Points> \n');

fprintf(results_vtu, '<DataArray  type="Float64"  NumberOfComponents="3"  format="ascii" > \n');


for i=1:numNodes    
    fprintf(results_vtu, '%f ',  x(i,1:3));
    fprintf(results_vtu, '\n');
end

fprintf(results_vtu, '</DataArray> \n');
fprintf(results_vtu, '</Points> \n');

%% Print cells
fprintf(results_vtu, '<Cells> \n');

%% Print cell connectivity
fprintf(results_vtu, '<DataArray  type="Int32"  Name="connectivity"  format="ascii"> \n');

for i=1:numCells
    fprintf(results_vtu, '%g ',  connect(i,1:numVertexesPerCell)-1 );
    fprintf(results_vtu, '\n');
end

fprintf(results_vtu, '</DataArray> \n');

% Print cell offsets
fprintf(results_vtu, '<DataArray  type="Int32"  Name="offsets"  format="ascii"> \n');

offset = 0;
for i=1:numCells
    offset = offset + numVertexesPerCell;
    fprintf(results_vtu, '%g ', offset);
    fprintf(results_vtu, '\n');
end

fprintf(results_vtu, '</DataArray> \n');

% Print cell types
fprintf(results_vtu, '<DataArray  type="UInt8"  Name="types"  format="ascii"> \n');

for i=1:numCells
    fprintf(results_vtu, '%g ', VTKCellCode);
    fprintf(results_vtu, '\n');
end

fprintf(results_vtu, '</DataArray> \n');
fprintf(results_vtu, '</Cells> \n');

%% Print result data

fprintf(results_vtu, '<PointData  Vectors="sigma"> \n');

% print displacement field

fprintf(results_vtu, '<DataArray  type="Float64"  Name="U" NumberOfComponents="3" format="ascii"> \n');
for i=1:numNodes
    fprintf(results_vtu, '%f   ', dispX(i) );
    fprintf(results_vtu, '%f   ', dispY(i) );
    fprintf(results_vtu, '%f   ', dispZ(i) );    
    fprintf(results_vtu, '\n');
end
fprintf(results_vtu, '</DataArray> \n');
fprintf(results_vtu, '</PointData> \n');

% end of VTK file

fprintf(results_vtu, '</Piece> \n');
fprintf(results_vtu, '</UnstructuredGrid> \n');
fprintf(results_vtu, '</VTKFile> \n');

fclose(results_vtu);







