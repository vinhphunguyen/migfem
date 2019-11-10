% data file for the double cantilever beam problem
% VP Nguyen, H Nguyen-Xuan High order Bsplines for delamination paper.

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
addpath ../C_files/

clear all

global noElems p q index elRangeU elRangeV noPtsX noPtsY uKnot vKnot De element
global controlPts weights

E           = 150e3;  % Young modulus
nu          = 0.25;  % Poisson’s ratio
stressState = 'PLANE_STRESS';
Gc = 0.352;
tensile=80;
cohLength = E*Gc/(tensile^2);

% Elasticity matrix

De = elasticityMatrix(E,nu,stressState);

eps = 1e-10;

L = 100;
W = 3;
a0 = 20;

p = 1;
q = 1;

% initial mesh- linear x linear

controlPts          = zeros(4,2,2);

controlPts(1:2,1,1) = [0;0];
controlPts(1:2,2,1) = [L;0];
controlPts(1:2,1,2) = [0;W];
controlPts(1:2,2,2) = [L;W];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot});

% evaluate order 

solid = nrbdegelev(solid,[1 1]); 

% h-refinement in Y direction to make sure it is C^{-1} along the 
% delamination path. The multiplicity must be q+1.
  
solid     = nrbkntins(solid,{[0.5 0.5] [0.5 0.5 0.5]});

% insert for intitial crack 

solid     = nrbkntins(solid,{[1-a0/L 1-a0/L] []});

% insert to make the mesh more square
solid     = nrbkntins(solid,{[0.25] []});

% figure
% hold on
% nrbctrlplot(solid);
% axis equal
% axis off
% 
% pause

uKnot     = cell2mat(solid.knots(1));
vKnot     = cell2mat(solid.knots(2));
refineCountX = 5;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
    id = find( uKnotVectorU == 1-a0/L);
    subKnot = uKnotVectorU(1:id);
    newKnotsY = [];
    newKnotsZ = [];
    
    % new knots along two directions (uniform)
    
    newKnotsX = subKnot(1:end-1) + 0.5*diff(subKnot);
    
    newKnots  = {newKnotsX []};
    
    % h-refinement
    
    solid     = nrbkntins(solid,newKnots);
    
    uKnot     = cell2mat(solid.knots(1));
    vKnot     = cell2mat(solid.knots(2));    
end

refineCountX = 5;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
    id = find( uKnotVectorU == 1-a0/L);
    subKnot   = uKnotVectorU(id:end);
    newKnotsX = subKnot(1:end-1) + 0.5*diff(subKnot);
    
    newKnots  = {newKnotsX newKnotsY};
    
    % h-refinement
    
    solid     = nrbkntins(solid,newKnots);
    
    uKnot      = cell2mat(solid.knots(1));
    vKnot      = cell2mat(solid.knots(2));
end

%% convert to IGA format

convert2DNurbs
noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;
    

%% generate element connectivity ...

generateIGA2DMesh

% plot the mesh

vMesh=buildVMesh2DInterface(solid);

% Plot mesh and enriched nodes to check
figure
hold on
plot_mesh(vMesh.node,vMesh.element,'Q4','b-',1.1);
n5 = plot(controlPts(:,1),controlPts(:,2),'go');
set(n5,'MarkerSize',4, 'MarkerFaceColor','g','MarkerEdgeColor','k','LineWidth',1.01);
%axis off

% vMesh.node(2*(noElemsU+1)+1:3*(noElemsU+1),2) = vMesh.node(2*(noElemsU+1)+1:3*(noElemsU+1),2) + 0.05;
% figure
% hold on
% plot_mesh(vMesh.node,vMesh.element,'Q4','b-',1.1);
% n5 = plot(vMesh.node(:,1),vMesh.node(:,2),'go');
% set(n5,'MarkerSize',3, 'MarkerFaceColor','g','MarkerEdgeColor','k','LineWidth',1.01);
% %axis off

%% Bezier extraction operators

[C,Cxi,Cet]  = bezierExtraction2D(uKnot,vKnot,p,q);
    
% build mesh of interface elements
[ielements] = buildGA1DMesh (uKnot,p);
delaminationNodes  =  find(abs(controlPts(:,2) -W/2  ) <1e-10);
mm                 = 0.5*length(delaminationNodes);
lowerNodes         = delaminationNodes(1:mm);
upperNodes         = delaminationNodes(mm+1:end);

iElements   = zeros(noElemsU,2*(p+1));

for i=1:noElemsU
    sctr = ielements.element(i,:);
    iElements(i,1:p+1)   = lowerNodes(sctr);
    iElements(i,p+2:end) = upperNodes(sctr);
end

%% write to jive mesh 

fixedNodes   =  [1 noPtsX];
forcedNodes  = intersect(find(abs(controlPts(:,1)-L/2)<eps),...
                         find(abs(controlPts(:,2)-W  )<eps));

forcedNodes  = [forcedNodes noCtrPts]; 
                     
% xi     = 0.5;
% span   = findspan(noPtsX-1,p,xi,uKnot);
% ielem  = span - 1;
% elem   = element(ielem,:);
% forcedNodes=elem(end-2:end);

noElems  = size(C,3);

fileName = '~/code/jive/bezier/delamination/mmb/mmb.mesh';

file = fopen(fileName, 'wt');

fprintf(file, '<Nodes>\n');

for i=1:length(controlPts)
   fprintf(file, '  %1d %2.6f %2.6f', i, controlPts(i,1),controlPts(i,2));
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

for i=1:noElemsU
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

for e=1:noElemsU
    Ce = Cxi(:,:,e);
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

for e=1:noElemsU
    Ce = Cxi(:,:,e);
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
for e=1:noElemsU
    Ce = Cxi(:,:,e);
    [row,col,val] = find(Ce);
    
    fprintf(file, '  %d ', noElems+e-1);
    for i=1:length(row)
        fprintf(file, '%2.1f ', val(i));
    end
    fprintf(file, ';\n');
end
fprintf(file, ' </Column>\n');

%write weights (not necessary for weights=1)

fprintf(file, ' <Column name = "weights" type = "float">\n');

for e=1:noElems
    W = weights(element(e,:));
    fprintf(file, '  %1d ',e-1);
    for j=1:length(W)
        fprintf(file, '%2.4f ', W(j));
    end
    fprintf(file, ';\n');
end

for e=1:noElemsU
    sctr = iElements(e,:);
    W = weights(sctr(1:0.5*end));
    fprintf(file, '  %1d ',noElems+e-1);
    for j=1:length(W)
        fprintf(file, '%2.4f ', W(j));
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

for i=1+noElems:noElems+noElemsU
   fprintf(file, '  %1d', i-1);
  
end
fprintf(file, '}\n');
fprintf(file, '</ElementGroup>\n');

id = find(unique(uKnot)==1-a0/L);

contacts  = id:noElemsU;
cohesives = setdiff(1:noElemsU, contacts);
contacts  = contacts - 1;
cohesives = cohesives - 1;

fprintf(file, '<ElementGroup name="contacts">\n{');

for i=1:length(contacts)
   fprintf(file, '  %1d', contacts(i) + noElems);
  
end
fprintf(file, '}\n');
fprintf(file, '</ElementGroup>\n');

fprintf(file, '<ElementGroup name="cohesives">\n{');

for i=1:length(cohesives)
   fprintf(file, '  %1d', cohesives(i) + noElems);
  
end
fprintf(file, '}\n');
fprintf(file, '</ElementGroup>\n');

% write node groups

fprintf(file, '<NodeGroup name="gr1">\n{');
fprintf(file, '  %1d', fixedNodes(1));
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="gr2">\n{');
fprintf(file, '  %1d', fixedNodes(2));
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="gr3">\n{');

for i=1:length(forcedNodes)
   fprintf(file, '  %1d', forcedNodes(i));
  
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fclose(file);

disp('writing mesh file for jive is done!!!')

pause

%%
% read jem-jive result and do post-processing

%vMesh=buildVMesh2DInterface(solid);

vtuFileName = '../results/mmb';
resultFile = '~/code/jive/bezier/delamination/mmb/mmb2x128-B5x2.dat';
V = xml_parseany(fileread(resultFile));

noStep = size(V,2)/2;

for it=1:noStep
    vtuFile = sprintf('../results/%s%d',vtuFileName,it);
    uData=V{2*it-1}.Section{1};   % displacement
    dData=V{2*it}.Section{1}; % damage
    disp = str2num(uData.CONTENT);
    dam  = str2num(dData.CONTENT);
    damage = zeros(noCtrPts,1);
    damage(dam(:,1),1) = dam(:,2);
    U    = disp(:,2:3);
    ok   = writeVTKForJive(solid,U,damage,vtuFile,vMesh);
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

