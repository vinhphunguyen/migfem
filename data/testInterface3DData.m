% data file for the 3D double cantilever beam problem

addpath ../C_files/
addpath ../fem_util/;
addpath ../nurbs-geopdes/inst/
addpath ../nurbs-util/
addpath ../meshing/
addpath ../fem-functions/
addpath ../post-processing/
addpath ../xiga/
addpath ../nurbs-geopdes/inst/
addpath ../delamination/
addpath ../xml_toolbox/

L = 1;
W = 1;
b = 2; % thickness


% initial mesh- linear x linear x linear

controlPts          = zeros(4,2,2,2);

controlPts(1:3,1,1,1) = [0;0;0];
controlPts(1:3,2,1,1) = [L;0;0];
controlPts(1:3,1,2,1) = [0;W;0];
controlPts(1:3,2,2,1) = [L;W;0];

controlPts(1:3,1,1,2) = [0;0;b];
controlPts(1:3,2,1,2) = [L;0;b];
controlPts(1:3,1,2,2) = [0;W;b];
controlPts(1:3,2,2,2) = [L;W;b];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];
wKnot = [0 0 1 1];

%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot wKnot});

% evaluate order 

solid = nrbdegelev(solid,[1 1 1]); 

% h-refinement in Y direction to make sure it is C^{-1} along the 
% delamination path. The multiplicity must be p+1.
  
solid     = nrbkntins(solid,{[] [] [0.5 0.5 0.5]});

uKnot     = cell2mat(solid.knots(1));
vKnot     = cell2mat(solid.knots(2));
wKnot     = cell2mat(solid.knots(3));

refineCountX = 1;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    wKnotVectorV = unique(wKnot);
   
    newKnotsY = [];
    newKnotsZ = [];
    
    % new knots along two directions (uniform)
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    
    newKnots  = {newKnotsX newKnotsY newKnotsZ};
    
    % h-refinement
    
    solid     = nrbkntins(solid,newKnots);
    
    uKnot     = cell2mat(solid.knots(1));
    vKnot     = cell2mat(solid.knots(2));
    wKnot     = cell2mat(solid.knots(3));
end

refineCountY = 1;
for i=1:refineCountY
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    uKnotVectorW = unique(wKnot);
   
    newKnotsZ = [];
    newKnotsX = [];
    
    % new knots along two directions (uniform)
    
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    
    newKnots  = {newKnotsX newKnotsY newKnotsZ};
    
    % h-refinement
    
    solid     = nrbkntins(solid,newKnots);
    
    uKnot     = cell2mat(solid.knots(1));
    vKnot     = cell2mat(solid.knots(2));
    wKnot     = cell2mat(solid.knots(3));
end

nrbplot(solid,[20 20 20]);

%% convert to IGA format

convert3DNurbs
noCtrPts       = noPtsX * noPtsY * noPtsZ;
noDofs         = noCtrPts * 3;
    
%% generate element connectivity ...

generateIGA3DMesh

mesh.p          = p;
mesh.q          = q;
mesh.r          = r;
mesh.uKnot      = uKnot;
mesh.vKnot      = vKnot;
mesh.wKnot      = wKnot;
mesh.noPtsX     = noPtsX;
mesh.noPtsY     = noPtsY;
mesh.noPtsZ     = noPtsZ;
mesh.weights    = weights;
mesh.controlPts = controlPts(:,1:3);
mesh.elements   = element;
mesh.rangeU     = elRangeU;
mesh.rangeV     = elRangeV;
mesh.rangeW     = elRangeW;
mesh.index      = index;
mesh.noElemsU   = noElemsU;
mesh.noElemsV   = noElemsV;
mesh.noElemsW   = noElemsW;
mesh.elConnU    = elConnU;
mesh.elConnV    = elConnV;
mesh.elConnW    = elConnW;
mesh.elemCount  = noElems;



% plot the mesh

vMesh=buildVMesh3DInterface(solid);

% Plot mesh and enriched nodes to check
figure
hold on
plot_mesh(vMesh.node,vMesh.element,'B8','b-',1.1);
n5 = plot3(controlPts(:,1),controlPts(:,2),controlPts(:,3),'go');
set(n5,'MarkerSize',6, 'MarkerFaceColor','g','MarkerEdgeColor','k','LineWidth',0.9);
axis off

% vMesh.node([9:12],3) = vMesh.node([9:12],3) + 0.03;
% 
% figure
% hold on
% plot_mesh(vMesh.node,vMesh.element,'B8','b-',1.1);
% n5 = plot3(controlPts(:,1),controlPts(:,2),controlPts(:,3),'go');
% set(n5,'MarkerSize',6, 'MarkerFaceColor','g','MarkerEdgeColor','k','LineWidth',0.9);
% axis off

%% Bezier extraction operators

[C,Cxi,Cet,Cze] = bezierExtraction3D(uKnot,vKnot,wKnot,p,q,r);
[Cxe]           = bezierExtraction2D(uKnot,vKnot,p,q);
    
eps = 1e-10;
fixedNodes    = find(abs(controlPts(:,3)-0) < eps);
forcedNodesU  = find(abs(controlPts(:,3)-b) < eps);
                                                 
% build mesh of interface elements

delaminationNodes  =  find(abs(controlPts(:,3) -b/2  ) <1e-10);
mm                 = 0.5*length(delaminationNodes);
lowerNodes         = delaminationNodes(1:mm);
upperNodes         = delaminationNodes(mm+1:end);

iElements   = zeros(noElemsU*noElemsV,2*(p+1)*(q+1));
iElementS   = buildIGA2DMesh (uKnot,vKnot,noPtsX,noPtsY,p,q);

for e=1:noElemsU*noElemsV
    iElements(e,1:(p+1)*(q+1))     = lowerNodes(iElementS(e,:));
    iElements(e,(p+1)*(q+1)+1:end) = upperNodes(iElementS(e,:));
end

%% write to jive mesh 


noElems  = size(C,3);

fileName = '~/code/jive/bezier/delamination/test/dcb3d.mesh';

file = fopen(fileName, 'wt');

fprintf(file, '<Nodes>\n');

for i=1:length(controlPts)
   fprintf(file, '  %1d %2.6f %2.6f', i, ...
            controlPts(i,1),controlPts(i,2),controlPts(i,3));
   fprintf(file, ';\n');
end

fprintf(file, '</Nodes>\n\n');

fprintf(file, '<Elements>\n');

% write solid elements

for i=1:noElems
   fprintf(file, '  %1d %1d', i-1, element(i,:) );
   fprintf(file, ';\n');
end

% write interface elements

for i=1:size(iElements,1)
   fprintf(file, '  %1d %1d', noElems+i-1, iElements(i,:) );
   fprintf(file, ';\n');
end

fprintf(file, '</Elements>\n\n');

% write Bezier extractors 
% first for solid elements
% then for interface elements

fprintf(file, '<ElementDatabase name="C">\n');

fprintf(file, ' <Column name = "irows" type = "int">\n');

for e=1:noElems
    Ce = C(:,:,e);
    [row,col] = find(Ce);
    fprintf(file, '  %1d ', e-1);
    for i=1:length(row)
        fprintf(file, '%1d ', row(i)-1);
    end
    fprintf(file, ';\n');
end

for e=1:size(iElements,1)
    Ce = Cxe(:,:,e);
    [row,col] = find(Ce);
    fprintf(file, '  %1d ', noElems+e-1);
    for i=1:length(row)
        fprintf(file, '%1d ', row(i)-1);
    end
    fprintf(file, ';\n');
end

fprintf(file, ' </Column>\n');

fprintf(file, ' <Column name = "jcols" type = "int">\n');

for e=1:noElems
    Ce = C(:,:,e);
    [row,col] = find(Ce);
    
    fprintf(file, '  %d ', e-1);
    for i=1:length(row)
        fprintf(file, '%1d ', col(i)-1);
    end
    fprintf(file, ';\n');
end

for e=1:size(iElements,1)
    Ce = Cxe(:,:,e);
    [row,col] = find(Ce);
    
    fprintf(file, '  %d ', noElems+e-1);
    for i=1:length(row)
        fprintf(file, '%1d ', col(i)-1);
    end
    fprintf(file, ';\n');
end

fprintf(file, ' </Column>\n');

fprintf(file, ' <Column name = "values" type = "float">\n');
for e=1:noElems
    Ce = C(:,:,e);
    [row,col,val] = find(Ce);
    
    fprintf(file, '  %d ', e-1);
    for i=1:length(row)
        fprintf(file, '%2.1f ', val(i));
    end
    fprintf(file, ';\n');
end
for e=1:size(iElements,1)
    Ce = Cxe(:,:,e);
    [row,col,val] = find(Ce);
    
    fprintf(file, '  %d ', noElems+e-1);
    for i=1:length(row)
        fprintf(file, '%2.1f ', val(i));
    end
    fprintf(file, ';\n');
end
fprintf(file, ' </Column>\n');

% write weights 

fprintf(file, ' <Column name = "weights" type = "float">\n');

for e=1:noElems
    w = weights(element(e,:));
    fprintf(file, '  %1d ',e-1);
    for j=1:length(w)
        fprintf(file, '%2.4f ', w(j));
    end
    fprintf(file, ';\n');
end

for e=1:size(iElements,1)
    sctr = iElements(e,:);
    w = weights(sctr(1:0.5*end));
    fprintf(file, '  %1d ',noElems+e-1);
    for j=1:length(w)
        fprintf(file, '%2.4f ', w(j));
    end
    fprintf(file, ';\n');
end

fprintf(file, ' </Column>\n');
fprintf(file, '</ElementDatabase>\n\n');

% write element groups

fprintf(file, '<ElementGroup name="solid">\n{');

for i=1:noElems
   fprintf(file, '  %1d', i-1);
  
end
fprintf(file, '}\n');
fprintf(file, '</ElementGroup>\n');

fprintf(file, '<ElementGroup name="interface">\n{');

for i=1+noElems:noElems+size(iElements,1)
   fprintf(file, '  %1d', i-1);
  
end
fprintf(file, '}\n');
fprintf(file, '</ElementGroup>\n');


% write node groups

fprintf(file, '<NodeGroup name="gr1">\n{');

for i=1:length(fixedNodes)
   fprintf(file, '  %1d', fixedNodes(i));
  
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="gr2">\n{');

for i=1:length(forcedNodesU)
   fprintf(file, '  %1d', forcedNodesU(i));
  
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');


fclose(file);

pause

%%
% read jem-jive result and do post-processing

E0           = 100;  % Young modulus
nu0          = 0.;  % Poisson’s ratio

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE COMPLIANCE MATRIX
D=zeros(6,6);
D(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
                                  nu0 1-nu0 nu0;
                                  nu0 nu0 1-nu0];
D(4:6,4:6)=E0/2/(1+nu0)*eye(3);

%vMesh=buildVMesh2DInterface(solid);

vtuFileName = '../results/testInterface3d';
resultFile = '~/code/jive/bezier/delamination/test/dcb3d.dat';
V = xml_parseany(fileread(resultFile));

noStep = size(V,2);

for it=1:noStep
    vtuFile = sprintf('../results/%s%d',vtuFileName,it);
    cc=V{it}.Section{1};
    disp = str2num(cc.CONTENT);
    U    = disp(:,2:end);
    ok=writeVTKForJive3D(mesh,vMesh,vtuFile,U,D);
end

%% write a PVD file that collects all VTU files
% An example is given
% <VTKFile byte_order='LittleEndian' type='Collection' version='0.1'>
% <Collection>
% <DataSet file='cantilever8-0.vtu' groups='' part='0' timestep='0'/>
% <DataSet file='cantilever8-1.vtu' groups='' part='0' timestep='1'/>
% </Collection>
% </VTKFile>

pvdFile = fopen(strcat('../results/',vtuFileName,'.pvd'), 'wt');

fprintf(pvdFile,'<VTKFile byte_order="LittleEndian" type="Collection" version="0.1">\n');
fprintf(pvdFile,'<Collection>\n');

for i = 1:noStep
    vtuFile = sprintf('%s%d%s',vtuFileName,i,'.vtu');
    fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''0'' timestep=''%d''/>\n',vtuFile,i);
end

fprintf(pvdFile,'</Collection>\n');
fprintf(pvdFile,'</VTKFile>\n');

fclose(pvdFile);
