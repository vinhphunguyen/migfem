% This file implements the non-symmetric Nitsche method to join two mechanical models:
% a 3D continuum model and a Mindlin plate model.
%
% Discretisation: standard Lagrange finite elements.
%
% Problem: Timoshenko beam in bending.
% 
% Vinh Phu Nguyen
% Cardiff University, UK
% 4 June 2013

addpath ../fem_util/
addpath ../gmshFiles/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/
addpath ../meshing/
addpath ../nurbs-util/
addpath ../nurbs-geopdes/inst/
addpath ../delamination/

clear all
colordef white
state = 0;
tic;

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


% MATERIAL PROPERTIES
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

noGPs = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE FINITE ELEMENT MESH
%
plotMesh  = 1;
disp([num2str(toc),'   GENERATING MESH'])

% MESH PROPERTIES

elemType  = 'B8'; % the element type for solid;
elemTypeB = 'Q4'; % the element type for plate;

% domain 1
L1    = 100;

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

solid = nrbmak(controlPts,{uKnot vKnot wKnot});

% evaluate order

%solid = nrbdegelev(solid,[0 0]);


uKnot     = cell2mat(solid.knots(1));
vKnot     = cell2mat(solid.knots(2));

refineCountX = 5;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    newKnotsY = [];
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {newKnotsX [] []};
    solid     = nrbkntins(solid,newKnots);
    uKnot      = cell2mat(solid.knots(1));
    vKnot      = cell2mat(solid.knots(2));
end

solid     = nrbkntins(solid,{[] [0.25 0.5 0.75] [0.1 0.3 0.6 0.9]});
mesh = buildIGA3DMesh (solid);

node1     = mesh.controlPts;
element1  = mesh.locElems;
element1(:,[3 4 7 8]) = element1(:,[4 3 8 7]);

numx1 = mesh.noElemsU;
numy1 = mesh.noElemsV;
numz1 = mesh.noElemsW;

% domain 2 (plate)
numx2 = 20;
numy2 = 2;

nnx=numx2+1;
nny=numy2+1;
node2t=square_node_array([L1 0],[L 0],[L width],[L1 width],nnx,nny);
inc_u=1;
inc_v=nnx;
node_pattern=[ 1 2 nnx+2 nnx+1 ];
element2 =make_elem(node_pattern,numx2,numy2,inc_u,inc_v);

node2 = zeros(size(node2t,1),3);
node2(:,[1 2]) = node2t;
node2(:,3) = t/2;
      
% boundary edges

bndNodes  = find(abs(node1(:,1)-L1)<1e-14);  

noElems  = numy1 * numz1;
bndMesh1 = zeros(noElems,4);
elConnU  = zeros(numy1,2);
elConnV  = zeros(numz1,2);

for i=1:numy1
    elConnU(i,:) = i:i+1;
end

for i=1:numz1
    elConnV(i,:) = i:i+1;
end

chan = reshape(bndNodes,numy1+1,numz1+1)';

e = 1;
for v=1:numz1
    vConn = elConnV(v,:);
    for u=1:numy1
        c = 1;
        uConn = elConnU(u,:);
        for i=1:length(vConn)
            for j=1:length(uConn)
              bndMesh1(e,c) = chan(vConn(i),uConn(j));
              c = c + 1;
            end
        end
        e = e + 1;
    end        
end

% note that bndMesh1 nodes are numbered as IGA elements
% for standard Lagrange elements need to be renumbered

bndMesh1(:,[3 4]) = bndMesh1(:,[4 3]);

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

bndMesh2 = [1 numx2+1];

[W1,Q1]=quadrature( 2, 'GAUSS', 2 ); % two point quadrature

GP1 = [];
GP2 = [];

for e=1:noElems    
    be = map (e);
    pe = bndMesh2(map2(e));
    sctrS    = bndMesh1(e,:);    
    sctrV    = element1(be,:);    
    pts1     = node1(sctrV,:);    
    pts2     = node2(element2(pe,:),[1 2]);    
    ptsS     = node1(sctrS,:);
    for q=1:size(W1)
        pt=Q1(q,:);
        wt=W1(q);
        [N,dNdxi]= lagrange_basis(elemTypeB,pt);  % element shape functions                
        J0       = dNdxi'*ptsS;
        a1       = J0(1,:);
        a2       = J0(2,:);
        a3       = cross(a1,a2); 
        norma    = norm(a3); a3 = a3/norma;
        
        x   = N'*ptsS;
        xP  = x([1 2]);
        X1  = global2LocalMap3D(x,pts1);        
        X2  = global2LocalMap  (xP,pts2,elemTypeB);     
        GP1 = [GP1;X1 wt*norma a3];
        GP2 = [GP2;X2 x(3)-t/2];
    end
end

% GET NODES ON DISPLACEMENT BOUNDARY
%      Here we get the nodes on the essential boundaries

fixedNode = find(node1(:,1)==0);  
rightNode = find(node2(:,1)==L); 

bndEdge  = zeros(numy2,2);

for i=1:numy2
    bndEdge(i,:) = rightNode(i:i+1);
end

uFixed    = zeros(1,length(fixedNode));  % a vector of u_x for the nodes
vFixed    = zeros(1,length(fixedNode));     
wFixed    = zeros(1,length(fixedNode)); 

udofs     = 3*fixedNode-2;
vdofs     = 3*fixedNode-1;
wdofs     = 3*fixedNode-0;

numnodes = size(node1,1) + size(node2,1);

% PLOT MESH
if ( plotMesh )  
    clf
    hold on
    plot_mesh(node1,element1,elemType,'g.-',2.4);
    plot_mesh(node2,element2,elemTypeB,'r.-',1.7);
    plot3(node1(fixedNode,1),node1(fixedNode,2),node1(fixedNode,3),'r*');
%    plot3(node2(rightNode,1),node2(rightNode,2),node2(rightNode,3),'rs');    
    axis on
    %axis([0 L -c c])
end


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

[W,Q]=quadrature( 2, 'GAUSS', 3 ); % 2x2x2 Gaussian quadrature

%% for domain 1 (continuum)

for e=1:size(element1,1)                          % start of element loop
    sctr = element1(e,:);           % element scatter vector
    nn   = length(sctr);
    sctrB(1:3:3*nn) = 3*sctr-2;
    sctrB(2:3:3*nn) = 3*sctr-1;
    sctrB(3:3:3*nn) = 3*sctr-0;
    pts = node1(sctr,:);
    for q=1:size(W,1)
        pt=Q(q,:);
        wt=W(q);
        [N,dNdxi]=lagrange_basis(elemType,pt); % element shape functions
        J0    = pts'*dNdxi;                % element Jacobian matrix
        invJ0 = inv(J0);
        dNdx  = dNdxi*invJ0;        
        B     = getBmatrix3D(nn,dNdx);                
        K(sctrB,sctrB)=K(sctrB,sctrB)+B'*C*B*W(q)*det(J0);
    end  % of quadrature loop
end    

%% for domain 2 (plate)
clear sctrB;

[W,Q]=quadrature(  noGPs, 'GAUSS', 2 ); % 2 x2point quadrature

for e=1:size(element2,1)  
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
    
    Bb     = zeros(3,nn3);
    Bs     = zeros(2,nn3);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
               
        [R,dRdxi]= lagrange_basis(elemTypeB,pt); % element shape functions
        J0       = pts'*dRdxi;                   % element Jacobian matrix
        invJ0    = inv(J0);
        dRdx     = dRdxi*invJ0;
                      
        % bending and shear B matrices
        
        Bb(1,2:3:3*nn) = dRdx(:,1);
        Bb(2,3:3:3*nn) = dRdx(:,2);
        Bb(3,2:3:3*nn) = dRdx(:,2);
        Bb(3,3:3:3*nn) = dRdx(:,1);
        
        Bs(1,2:3:3*nn) = -R;
        Bs(1,1:3:3*nn) = dRdx(:,1);
        Bs(2,1:3:3*nn) = dRdx(:,2);
        Bs(2,3:3:3*nn) = -R;
        
        % compute elementary stiffness matrix 
                
        K(sctrB,sctrB) = K(sctrB,sctrB)+(Bb' * Cb * Bb + Bs' * Cs * Bs) * det(J0) * wt;
        %f(sctrw)       = f(sctrw) + q0 * R' * det(J0) * wt;
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

for i=1:noElems                     % start of element loop
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
    
    for q=4*(i-1)+1:4*(i-1)+4
        pt1    = GP1(q,1:3);
        pt2    = GP2(q,1:2);
        wt1    = GP1(q,4);
        y      = GP2(q,3);
        normal = GP1(q,5:end);
        n = [normal(1) 0 normal(2) 0 normal(3);...
             0 normal(2) normal(1) normal(3) 0;...
             0 0 0 normal(2) normal(1)];
        n=-n;
        
        [N1,dN1dxi] = lagrange_basis(elemType, pt1); 
        [N2,dN2dxi] = lagrange_basis(elemTypeB,pt2); 
        
        J1 = pts1'*dN1dxi;                
        J2 = pts2'*dN2dxi;                
                       
        dN1dx   = dN1dxi*inv(J1);        
        dN2dx   = dN2dxi*inv(J2);
        
        % B matrices
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
           
        Kd11 = Kd11 + 0.5 * Nm1'* n * Csolid * B1 *wt1;
        Kd12 = Kd12 + 0.5 * Nm1'* n * Cplate * B2 *wt1;
        Kd21 = Kd21 + 0.5 * Nm2'* n * Csolid * B1 *wt1;
        Kd22 = Kd22 + 0.5 * Nm2'* n * Cplate * B2 *wt1;
               
    end  % of quadrature loop
   
    K(sctrB1,sctrB1)  = K(sctrB1,sctrB1)  + Kd11 - Kd11';    
    K(sctrB1,sctrB2)  = K(sctrB1,sctrB2)  + Kd12 + Kd21';
    K(sctrB2,sctrB1)  = K(sctrB2,sctrB1)  - Kd21 - Kd12';
    K(sctrB2,sctrB2)  = K(sctrB2,sctrB2)  - Kd22 + Kd22';
end    

%% External force

[W,Q]=quadrature( 2, 'GAUSS', 1 ); 

for e=1:size(bndEdge,1)                          % start of element loop
    sctr  = bndEdge(e,:);           % element scatter vector
    sctrg = sctr + size(node1,1);
    nn    = length(sctr);
    sctrf = 3*sctrg-2;
    pts   = node2(sctr,:);
    for q=1:size(W,1)
        pt=Q(q,:);
        wt=W(q);
        [N,dNdxi]=lagrange_basis('L2',pt); % element shape functions
        J0       = dNdxi'*pts;                % element Jacobian matrix               
        f(sctrf)=f(sctrf)+N*P*W(q)*norm(J0);
    end  % of quadrature loop
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

% Here we plot the stresses and displacements of the solution.
disp([num2str(toc),'   POST-PROCESSING'])

dispNorm=L/max(sqrt(Ux1.^2+Uy1.^2));
scaleFact=0.1*dispNorm;


colordef black
figure
clf
hold on
%plot_field(node1+scaleFact*[Ux1 Uy1],element1,elemType,Ux1);
%plot_field(node2+scaleFact*[Ux2 Uy2],element2,elemType,Ux2);
plot_mesh(node1+scaleFact*[Ux1 Uy1 Uz1],element1,elemType,'w.-',1);
plot_mesh(node2+scaleFact*[Uz2 Uy2 Ux2],element2,elemTypeB,'w.-',1);
%colorbar
%axis off
%title('DEFORMED DISPLACEMENT IN Y-DIRECTION')

%% write to VTK

stress=zeros(size(element1,1),size(element1,2),6);


stressPoints=[-1 -1 -1;1 -1 -1;1 1 -1;-1 1 -1; ...
              -1 -1  1;1 -1  1;1 1 1; -1 1 1];

for e=1:size(element1,1)
    sctr=element1(e,:);
    nn=length(sctr);
    sctrB(1:3:3*nn) = 3*sctr-2;
    sctrB(2:3:3*nn) = 3*sctr-1;
    sctrB(3:3:3*nn) = 3*sctr-0;    
    for q=1:nn
        pt=stressPoints(q,:);
        [N,dNdxi]=lagrange_basis(elemType,pt);
        J0=node1(sctr,:)'*dNdxi;
        invJ0=inv(J0);
        dNdx=dNdxi*invJ0;
        B      = getBmatrix3D(nn,dNdx);  
        strain=B*U(sctrB);
        stress(e,q,:)=C*strain;
    end
end   % of element loop

vtuFile1   = '../results/nitscheNS-plate-3D-continuum';
vtuFile2    = '../results/nitscheNS-plate-3D-plate';
vtuFileName = '../results/nitscheNS-plate-3D';

numNode  = size(node1,1);
numNode2 = size(node2,1);

% normal stresses
sigmaXX = zeros(numNode,2);
sigmaYY = zeros(numNode,2);
sigmaZZ = zeros(numNode,2);

% shear stresses
sigmaXY = zeros(numNode,2);
sigmaYZ = zeros(numNode,2);
sigmaZX = zeros(numNode,2);
% von Mises stress
sigmaVM = zeros(numNode,1);

for e=1:size(element1,1)
    connect = element1(e,:);
    for in=1:8
        nid = connect(in);
        sigmaXX(nid,:) = sigmaXX(nid,:) + [stress(e,in,1) 1];
        sigmaYY(nid,:) = sigmaYY(nid,:) + [stress(e,in,2) 1];
        sigmaZZ(nid,:) = sigmaZZ(nid,:) + [stress(e,in,3) 1];
        sigmaXY(nid,:) = sigmaXY(nid,:) + [stress(e,in,4) 1];
        sigmaYZ(nid,:) = sigmaYZ(nid,:) + [stress(e,in,5) 1];
        sigmaZX(nid,:) = sigmaZX(nid,:) + [stress(e,in,6) 1];
    end
end

% Average nodal stress values (learned from Mathiew Pais XFEM code)
sigmaXX(:,1) = sigmaXX(:,1)./sigmaXX(:,2); sigmaXX(:,2) = [];
sigmaYY(:,1) = sigmaYY(:,1)./sigmaYY(:,2); sigmaYY(:,2) = [];
sigmaZZ(:,1) = sigmaZZ(:,1)./sigmaZZ(:,2); sigmaZZ(:,2) = [];
sigmaXY(:,1) = sigmaXY(:,1)./sigmaXY(:,2); sigmaXY(:,2) = [];
sigmaYZ(:,1) = sigmaYZ(:,1)./sigmaYZ(:,2); sigmaYZ(:,2) = [];
sigmaZX(:,1) = sigmaZX(:,1)./sigmaZX(:,2); sigmaZX(:,2) = [];

VTKPostProcess3d(node1,element1,'B8',vtuFile1,...
    [sigmaXX sigmaYY sigmaZZ sigmaXY sigmaYZ sigmaZX sigmaVM],[Ux1 Uy1 Uz1]);


VTKPostProcess(node2,element2,3,'Quad4',vtuFile2,...
             [zeros(numNode2,1) zeros(numNode2,1) zeros(numNode2,1)],[Uz2 Uy2 Ux2]);

pvdFile = fopen(strcat('../results/',vtuFileName,'.pvd'), 'wt');

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



