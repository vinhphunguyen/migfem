% data file for the 3D double cantilever beam problem


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
clc

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%%
l = 85/2;
L = 80;
t = 10.1;
w = 40;
h = 20;
no = 45; % number of plies
 
%%
uKnot = [0 0 0 1 2 3 4 5 5 5];
%uKnot = uKnot/max(uKnot);

controlPts          = zeros(4,7);

controlPts(1:2,1) = [0;0];
controlPts(1:2,2) = [10;0];

controlPts(1:2,3) = [0.5*L-8;h-3];
controlPts(1:2,4) = [0.5*L;h];
controlPts(1:2,5) = [0.5*L+8;h-3];

controlPts(1:2,6) = [L-10;0];
controlPts(1:2,7) = [L;0];

controlPts(4,:,:)   = 1;

curve = nrbmak(controlPts,uKnot);


%% offset curve

alpha = 0.1;
beta  = 0.7;

eps1   = 1e-2;
eps2   = 1e-2;

maxIter = 50;

%[oCurve,offsetPts]= offsetCurveNormal(curve,t,alpha,beta,eps1,maxIter);
[oCurve,offsetPts]= offsetCurveNormal(curve,t,alpha,beta,eps1,maxIter);

%oCurve.coefs(1,end) = 2*l+L;

figure
hold on
nrbctrlplot(curve);
nrbctrlplot(oCurve);
%nrbctrlplot(oCurve1);
%nrbctrlplot(solid1);
n5 = plot(offsetPts(:,1),offsetPts(:,2),'r--');
axis equal
axis off

% cross section from the two curves
surf  = surfaceFromTwoCurves(curve, oCurve);

% solid from surface by extrusion
solid = nrbextrude(surf, [0,0,w]);

% order elevation => quartic-quadratic-linear

solid = nrbdegelev(solid,[2 1 0]);

% h-refinement in Y direction to make sure it is C^{-1} along the
% delamination path. The multiplicity must be q+1.

knots=[];
for ip=1:no-1
    dd = ip/no;
    knots = [knots dd dd dd];
end

solid     = nrbkntins(solid,{[] knots []});

% h-refinement
ods      = solid.order -1;
refCount = [3 0 2];

pnew = ods(1);
qnew = ods(2);
knew = ods(3);

refCountX = refCount(1);
refCountY = refCount(2);
refCountZ = refCount(3);

orders = solid.order - 1;
p      = orders(1);
q      = orders(2);
k      = orders(3);

solid  = nrbdegelev(solid,[pnew-p qnew-q knew-k]); 

uKnot  = cell2mat(solid.knots(1));
vKnot  = cell2mat(solid.knots(2));
wKnot  = cell2mat(solid.knots(3));

for i=1:refCountX
    uKnotVectorU = unique(uKnot);          
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);   
    newKnots  = {newKnotsX [] []};    
    % h-refinement    
    solid     = nrbkntins(solid,newKnots);    
    uKnot     = cell2mat(solid.knots(1));
end

for i=1:refCountY
    vKnotVectorU = unique(vKnot);  
    %newKnotsY = [];    
    % new knots along two directions (uniform)    
    newKnotsY = vKnotVectorU(1:end-1) + 0.5*diff(vKnotVectorU);   
    newKnots  = {[] newKnotsY []};    
    % h-refinement    
    solid     = nrbkntins(solid,newKnots);    
    vKnot     = cell2mat(solid.knots(2));
end

for i=1:refCountZ
    wKnotVectorU = unique(wKnot);  
    %newKnotsY = [];    
    % new knots along two directions (uniform)    
    newKnotsZ = wKnotVectorU(1:end-1) + 0.5*diff(wKnotVectorU);   
    newKnots  = {[] [] newKnotsZ};    
    % h-refinement    
    solid     = nrbkntins(solid,newKnots);    
    wKnot     = cell2mat(solid.knots(3));
end


% FE mesh from NURBS solids
mesh = buildIGA3DMesh(solid);

% plot the mesh

vMesh=buildVMesh3DInterface(solid);

% Plot mesh and enriched nodes to check
figure
hold on
plot_mesh(vMesh.node,vMesh.element,'B8','b-',1.1);
n5 = plot3(controlPts(:,1),controlPts(:,2),controlPts(:,3),'go');
set(n5,'MarkerSize',6, 'MarkerFaceColor','g','MarkerEdgeColor','k','LineWidth',0.9);
axis off

%% Bezier extraction operators

[C,Cxi,Cet,Cze] = bezierExtraction3D(mesh.uKnot,mesh.vKnot,mesh.wKnot,mesh.p,mesh.q,mesh.r);
[Cxz]           = bezierExtraction2D(mesh.uKnot,mesh.wKnot,mesh.p,mesh.r);


uKnot = mesh.uKnot;
vKnot = mesh.vKnot;
wKnot = mesh.wKnot;
noPtsX = mesh.noPtsX;
noPtsY = mesh.noPtsY;
noPtsZ = mesh.noPtsZ;
p      = mesh.p;
q      = mesh.q;
r      = mesh.r;
noElems    = mesh.elemCount;
noElemsU   = mesh.noElemsU;
noElemsV   = mesh.noElemsV;
noElemsW   = mesh.noElemsW;
controlPts = mesh.controlPts;
element    = mesh.globElems;
mesh.elements = element;
weights    = mesh.weights;
noCtrPts   = size(mesh.controlPts,1);

% build mesh of interface elements
iElementS   = buildIGA2DMesh (uKnot,wKnot,noPtsX,noPtsZ,p,r);
iElements   = zeros(noElemsU*noElemsW*(no-1),2*(p+1)*(r+1));

e = 1;
eps = 1e-9;
% loop over plies
for ip=1:no-1
    y0 = t/no*ip;
    dd = find(abs(controlPts(:,2)-y0)<eps);
    start0 = dd(1);
    start1 = start0 + noPtsX;
    step   = find(dd==start1);
    
    id = ones(length(dd)/(step-1),1);
    id(1) = 1;
    for i=2:length(dd)/(step-1)
        id(i) = id(i-1) + step-1;
    end
    lowerNodes=[];
    upperNodes=[];
    for ii=1:noPtsZ
      start1  = dd(id(2*ii-1));
      start2  = dd(id(2*ii-0));
      lowerNodes         = [lowerNodes; start1:start1+noPtsX-1];
      upperNodes         = [upperNodes; start2:start2+noPtsX-1];
    end
    
    for i=1:noElemsU*noElemsW
        sctr = iElementS(i,:);        
        iElements(e,1:(p+1)*(r+1))     = lowerNodes(sctr);
        iElements(e,(p+1)*(r+1)+1:end) = upperNodes(sctr);
        e = e + 1;
    end
    
    %n5 = plot(controlPts(lowerNodes,1),controlPts(lowerNodes,2),'go');
    %set(n5,'MarkerSize',6, 'MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',0.9);
end

%% element sets

matMap = ones(noElems,1);
elemset0 = [];
count    = 1;
for iw=1:noElemsW
    for iv=1:noElemsV
        for iu=1:noElemsU
            if mod(iv,2) ~= 0
                elemset0 = [elemset0;count];
            end
            count = count + 1;
        end
    end
end

elemset90 = setdiff(1:noElems,elemset0);

matMap(elemset90,1) = 2;


%% node sets

eps = 1e-10;
fixedNodes    = find(abs(controlPts(:,1)) < eps);
forcedNodes   = find(abs(controlPts(:,1)-(2*l+L)) < eps);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write to jive mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noElems  = size(C,3);

fileName = '~/code/jive/bezier/delamination/curved3d/curved3d.mesh';

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

e= 1;
for il=1:no-1
    for ie=1:noElemsU*noElemsW
        Ce = Cxz(:,:,ie);
        [row,col] = find(Ce);
        fprintf(file, '  %1d ', noElems+e-1);
        for i=1:length(row)
            fprintf(file, '%1d ', row(i)-1);
        end
        fprintf(file, ';\n');
        e=e+1;
    end
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

e= 1;
for il=1:no-1
    for ie=1:noElemsU*noElemsW
        Ce = Cxz(:,:,ie);
        [row,col] = find(Ce);
        
        fprintf(file, '  %d ', noElems+e-1);
        for i=1:length(row)
            fprintf(file, '%1d ', col(i)-1);
        end
        fprintf(file, ';\n');
        e=e+1;
    end
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

e= 1;
for il=1:no-1
    for ie=1:noElemsU*noElemsW
        Ce = Cxz(:,:,ie);
        [row,col,val] = find(Ce);
        
        fprintf(file, '  %d ', noElems+e-1);
        for i=1:length(row)
            fprintf(file, '%2.4f ', val(i));
        end
        fprintf(file, ';\n');
        e=e+1;
    end
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

fprintf(file, '<ElementGroup name="solid0">\n{');

for i=1:length(elemset0)
    fprintf(file, '  %1d', elemset0(i)-1);
    
end
fprintf(file, '}\n');
fprintf(file, '</ElementGroup>\n');

fprintf(file, '<ElementGroup name="solid90">\n{');

for i=1:length(elemset90)
    fprintf(file, '  %1d', elemset90(i)-1);
    
end
fprintf(file, '}\n');
fprintf(file, '</ElementGroup>\n');

fprintf(file, '<ElementGroup name="interface">\n{');

for i=1+noElems:noElems+size(iElements,1)
    fprintf(file, '  %1d', i-1);
    
end
fprintf(file, '}\n');
fprintf(file, '</ElementGroup>\n');

% % differentiate contact elements and cohesive elements
% 
% id = find(unique(uKnot)==a0/L) - 1;
% 
% surf = zeros(noElemsV,noElemsU);
% count = 1;
% for j=1:noElemsV
%     for i=1:noElemsU
%         surf(j,i) = count;
%         count = count + 1;
%     end
% end
% 
% contacts  = surf(:,1:id);
% cohesives = surf(:,id+1:end);
% 
% contacts  = contacts - 1;
% cohesives = cohesives - 1;
% 
% contacts  = contacts(:);
% cohesives = cohesives(:);
% 
% fprintf(file, '<ElementGroup name="contacts">\n{');
% 
% for i=1:length(contacts)
%     fprintf(file, '  %1d', contacts(i) + noElems);
%     
% end
% fprintf(file, '}\n');
% fprintf(file, '</ElementGroup>\n');
% 
% fprintf(file, '<ElementGroup name="cohesives">\n{');
% 
% for i=1:length(cohesives)
%     fprintf(file, '  %1d', cohesives(i) + noElems);
%     
% end
% fprintf(file, '}\n');
% fprintf(file, '</ElementGroup>\n');

% write node groups

fprintf(file, '<NodeGroup name="fix">\n{');

for i=1:length(fixedNodes)
    fprintf(file, '  %1d', fixedNodes(i));
    
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="force">\n{');

for i=1:length(forcedNodes)
    fprintf(file, '  %1d', forcedNodes(i));
    
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fclose(file);

disp('writing jive mesh done!')

pause

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read jem-jive result and do post-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Materials

e1 = 115e3;
e2 = 8.5e3;
e3 = 8.5e3;
nu12 = 0.29;
nu23 = 0.29;
nu31 = 0.29;
g12  = 4.5e3;
g23  = e2/(2+2*nu23);
g31  = 4.5e3;

dirs = [0,90];

for i =1:length(dirs)
    material = createOrthotropicMaterial(dirs(i),3,e1,e2,e3,nu12,nu23,nu31,g12,g23,g31);
    materials{i} = material;
end

%vMesh=buildVMesh2DInterface(solid);

vtuFileName = '../results/curved3d';
resultFile = '~/code/jive/bezier/delamination/curved3d/curved3d.dat';
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
    U    = disp(:,2:end);
    ok   = writeVTKForJive3D(mesh,vMesh,vtuFile,U,damage,materials,matMap);
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

