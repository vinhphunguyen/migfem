%
%
% Vinh Phu Nguyen
% Cardiff University, UK
% 5 July 2013

addpath ../fem_util/
addpath ../gmshFiles/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/
addpath ../meshing/
addpath ../nurbs-util/
addpath ../nurbs-geopdes/inst/
addpath ../delamination/
addpath ../C_files/

clear all
colordef white
state = 0;
tic;

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


E0  = 1000;  % Young modulus
nu0 = 0.3;   % Poisson ratio


% BEAM PROPERTIES
L     = 320;     % length of the beam
width = 80;
t     = 20;     % thicknes

% TIP LOAD
P = -10/20; % the thickness = 20, to make it equivalent to 3D-Plate model

% Constitutive matrices:
C=zeros(6,6);
C(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
    nu0 1-nu0 nu0;
    nu0 nu0 1-nu0];
C(4:6,4:6)=E0/2/(1+nu0)*eye(3);

noGPs = 2;

% penalty parameter

alpha = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE FINITE ELEMENT MESH
%
plotMesh  = 1;
disp([num2str(toc),'   GENERATING MESH'])

% domain 1 ----------------------------------------------------------------

controlPts          = zeros(4,2,2,2);

controlPts(1:3,1,1,1) = [0;0;0];
controlPts(1:3,2,1,1) = [L;0;0];
controlPts(1:3,1,2,1) = [0;width;0];
controlPts(1:3,2,2,1) = [L;width;0];

controlPts(1:3,1,1,2) = [0;0;t];
controlPts(1:3,2,1,2) = [L;0;t];
controlPts(1:3,1,2,2) = [0;width;t];
controlPts(1:3,2,2,2) = [L;width;t];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];
wKnot = [0 0 1 1];

solid1 = nrbmak(controlPts,{uKnot vKnot wKnot});

% evaluate order

solid1 = nrbdegelev(solid1,[2 2 2]);


uKnot     = cell2mat(solid1.knots(1));
vKnot     = cell2mat(solid1.knots(2));
wKnot     = cell2mat(solid1.knots(3));

refineCountX = 6;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    uKnotVectorW = unique(wKnot);
    newKnotsY = [];
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {newKnotsX [] []};
    solid1     = nrbkntins(solid1,newKnots);
    uKnot     = cell2mat(solid1.knots(1));
    vKnot     = cell2mat(solid1.knots(2));
    wKnot     = cell2mat(solid1.knots(3));
end


% refineCountX = 2;
% for i=1:refineCountX
%     uKnotVectorU = unique(uKnot);
%     uKnotVectorV = unique(vKnot);
%     uKnotVectorW = unique(wKnot);
%     newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
%     newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
%     newKnotsZ = uKnotVectorW(1:end-1) + 0.5*diff(uKnotVectorW);
%     newKnots  = {[] newKnotsY newKnotsZ};
%     solid1     = nrbkntins(solid1,newKnots);
%     uKnot     = cell2mat(solid1.knots(1));
%     vKnot     = cell2mat(solid1.knots(2));
%     wKnot     = cell2mat(solid1.knots(3));
% end

solid1    = nrbkntins(solid1,{[] [0.25 0.5 0.75] [0.1 0.3 0.6 0.9]});

solid1     = nrbkntins(solid1,newKnots);

mesh1      = buildIGA3DMesh (solid1);

node1     = mesh1.controlPts;
element1  = mesh1.locElems;

numx1 = mesh1.noElemsU;
numy1 = mesh1.noElemsV;
numz1 = mesh1.noElemsW;


%% ------------------------------------------------------------------------
% boundary nodes

fixedNode = find(node1(:,1)==0);
rightNode = find(node1(:,1)==L);


uFixed    = zeros(1,length(fixedNode));  % a vector of u_x for the nodes
vFixed    = zeros(1,length(fixedNode));
wFixed    = zeros(1,length(fixedNode));

udofs     = 3*fixedNode-2;
vdofs     = 3*fixedNode-1;
wdofs     = 3*fixedNode-0;

numnodes = size(node1,1);


[bndMesh,botIndex]     = surfaceMesh (mesh1.vKnot,mesh1.wKnot,rightNode,mesh1.q,mesh1.r,...
    mesh1.noPtsY,mesh1.noPtsZ,mesh1.rangeV,mesh1.rangeW,mesh1.elConnV,mesh1.elConnW);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM DATA STRUCTURES

f=zeros(3*numnodes,1);          % external load vector
K=sparse(3*numnodes,3*numnodes);  % stiffness matrix

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

%%%%%%%%%%%%%%%%%%% END OF STIFFNESS MATRIX COMPUTATION %%%%%%%%%%%%%%%%%%%

[W1,Q1]=quadrature( mesh1.q+1, 'GAUSS', 2 ); % 2x2x2 Gaussian quadrature

for e=1:size(bndMesh,1);
    idu    = botIndex(e,1);
    idv    = botIndex(e,2);
    xiE    = mesh1.rangeV(idu,:); % [xi_i,xi_i+1]
    etaE   = mesh1.rangeW(idv,:); % [eta_j,eta_j+1]
    
    sctr   = bndMesh(e,:);          %  element scatter vector
    nn     = length(sctr);
    
    sctrf  = 3*sctr;
    pts    = node1(sctr,:);
    
    for gp=1:size(W1,1)
        pt      = Q1(gp,:);
        wt      = W1(gp);
        
        Xi      = parent2ParametricSpace(xiE, pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE, etaE);
        [R dRdxi dRdeta] = NURBS2DBasisDers([Xi; Eta],mesh1.q,mesh1.r,mesh1.vKnot,mesh1.wKnot,mesh1.weights);
                                        
        point   = R*pts;
        % basis vectors of the loaded surface
        jacob   = [dRdxi; dRdeta] * pts; % 2x3 matrix
        a1      = jacob(1,:);
        a2      = jacob(2,:);
        a3      = cross(a1,a2); % normal to the surface
        J1      = norm(a3);
        
        f(sctrf) = f(sctrf) + R'*P*J1*J2*wt;
    end
end

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

displacement = [U(1:3:end) U(2:3:end) U(3:3:end)];
damage       = zeros(length(U),1);

matMap1 = ones(mesh1.noElems,1);

vMesh1=buildVisualizationMesh3D(solid1);

meshes{1}  = mesh1;
vmeshes{1} = vMesh1;

meshData.mesh    = meshes;
meshData.vmesh   = vmeshes;
% meshData.noElems = meshes{1}.noElems;
% meshData.noPts   = meshes{1}.noPts;
meshData.matMap{1}=matMap1;
material.stiffMat=C;
materials{1}      = material;

vtuFileName = '../results/iga-3DPlate';
for ip=1:1
    vtuFile = strcat(vtuFileName,'-mesh',num2str(ip));
    %figure; hold on;
    ok      = plotStress3DForPatch(meshData,ip,vtuFile,displacement,damage,materials);
end

pvdFile = fopen(strcat('../results/',vtuFileName,'.pvd'), 'wt');

fprintf(pvdFile,'<VTKFile byte_order="LittleEndian" type="Collection" version="0.1">\n');
fprintf(pvdFile,'<Collection>\n');


for im=1:length(meshes)
    vtuFile = sprintf('%s%s%d%s',vtuFileName,'-mesh',im,'.vtu');
    fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''%d'' timestep=''%d''/>\n',...
        vtuFile,im,1);
end


fprintf(pvdFile,'</Collection>\n');
fprintf(pvdFile,'</VTKFile>\n');

fclose(pvdFile);


