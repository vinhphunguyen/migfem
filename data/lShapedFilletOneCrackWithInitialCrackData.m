% L-shaped fillet composite specimen.
% Only one crack between ply 5 and 6 and there is one initial crack.
% VP Nguyen
% Cardiff University
% 5 May, 2013

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
addpath ../examples/

clear all
clc

global noElems p q index elRangeU elRangeV noPtsX noPtsY uKnot vKnot De element
global controlPts weights

E           = 150e3;  % Young modulus
nu          = 0.25;  % Poisson’s ratio
stressState = 'PLANE_STRAIN';


% Elasticity matrix

De = elasticityMatrix(E,nu,stressState);

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

%%
H  = 6.4;
R  = 2.55;
R0 = 2.25;
no = 15;   % number of plies
t  = R0/no;

% initial mesh- quadratic x linear

uKnot = [0 0 0 1 1 2 2 3 3 3];
uKnot = uKnot/max(uKnot);
vKnot = [0 0 1 1];
controlPts      = zeros(4,7,2);

controlPts(1:2,1,1) = [H+R;0];
controlPts(1:2,2,1) = [(H+R)/2;0];
controlPts(1:2,3,1) = [R;0];
controlPts(1:2,4,1) = [0;0];
controlPts(1:2,5,1) = [0;R];
controlPts(1:2,6,1) = [0;(H+R)/2];
controlPts(1:2,7,1) = [0;H+R];

controlPts(1:2,1,2) = [H+R;-R0];
controlPts(1:2,2,2) = [(H+R)/2;-R0];
controlPts(1:2,3,2) = [R;-R0];
controlPts(1:2,4,2) = [-R0;-R0];
controlPts(1:2,5,2) = [-R0;R];
controlPts(1:2,6,2) = [-R0;(H+R)/2];
controlPts(1:2,7,2) = [-R0;H+R];

fac = 1/sqrt(2);
controlPts(4,:) = 1;
controlPts(4,4,1) = fac;
controlPts(4,4,2) = fac;
controlPts(1:2,4,1) = fac*controlPts(1:2,4,1);
controlPts(1:2,4,2) = fac*controlPts(1:2,4,2);

% build NURBS object
solid = nrbmak(controlPts,{uKnot,vKnot});

%% evaluate order => bi-quadratic

solid = nrbdegelev(solid,[0 1]);

% h-refinement in Y direction to make sure it is C^{-1} along the
% delamination path. The multiplicity must be q+1.

knots=[];
for ip=1:no-1
    dd = t*ip/R0;
    knots = [knots dd dd];
end

solid     = nrbkntins(solid,{[] knots});

solid     = nrbkntins(solid,{[] 5*t/R0});

%%
% find position of intial crack

convert2DNurbs

tol = 1e-3;
f   = 1;
p = 2;
mesh1d = buildGA1DMesh (uKnot,p);
e = 2;
xiE   = mesh1d.range(e,:); % [xi_i,xi_i+1]
conn  = mesh1d.element(e,:);


y0 = -t*5;
dd = find(abs(controlPts(:,2)-y0)<1e-10);
nodes= dd(1):dd(1)+noPtsX-1;
pts  = controlPts(nodes,:);
wts  = weights(nodes);

[sCurve,iKnot] = plot2DNURBSCurve (uKnot, pts, p, wts,200);

figure
hold on
nrbctrlplot(solid);
view([0 90])
% n5 = plot(pts(conn,1),pts(conn,2),'go');
% set(n5,'MarkerSize',6, 'MarkerFaceColor','g','MarkerEdgeColor','k','LineWidth',0.9);

pt1    = xiE(1);
x0    = 1.;
while abs(f) > tol    
    [N dNdxi] = NURBS1DBasisDers(pt1,p,uKnot,wts);
    f = x0 - N*pts(conn,1);
    pt1 = pt1 + f/(dNdxi*pts(conn,1));
end

pt2    = xiE(1);
x0    = 0.25;
f=1;
while abs(f) > tol    
    [N dNdxi] = NURBS1DBasisDers(pt2,p,uKnot,wts);
    f = x0 - N*pts(conn,1);
    pt2 = pt2 + f/(dNdxi*pts(conn,1));
end

%initCracks=[pt1;pt2];

initCracks=[0.3;0.7];
initCrackPts=[];
ss=[1 3];
for i=1:length(initCracks)
    e=ss(i);
    conn  = mesh1d.element(e,:);
    [N dNdxi] = NURBS1DBasisDers(initCracks(i),p,uKnot,wts);
    initCrackPts = [initCrackPts;N*pts(conn,:)];
end

figure
hold on
plot(sCurve(1,:),sCurve(2,:),'b-');
n5 = plot(initCrackPts(:,1),initCrackPts(:,2),'go');
set(n5,'MarkerSize',6, 'MarkerFaceColor','g','MarkerEdgeColor','k','LineWidth',0.9);
axis off
axis equal

% insert initial cracks knot 

solid =nrbkntins(solid,{[initCracks(1) initCracks(1) initCracks(2) initCracks(2)] []});

%%

uKnot     = cell2mat(solid.knots(1));
vKnot     = cell2mat(solid.knots(2));

% h-refinement along xi direction

refineCountX = 6;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
    newKnotsY = [];
    
    % new knots along two directions (uniform)
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    %newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    
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
plot_mesh(vMesh.node,vMesh.element,'Q4','b-',0.9);
n5 = plot(initCrackPts(:,1),initCrackPts(:,2),'go');
set(n5,'MarkerSize',6, 'MarkerFaceColor','g','MarkerEdgeColor','k','LineWidth',0.9);
axis off


%% Bezier extraction operators

[C,Cxi,Cet]  = bezierExtraction2D(uKnot,vKnot,p,q);

%% build mesh of interface elements

[ielements] = buildGA1DMesh (uKnot,p);
iElements   = zeros(noElemsU*(no-1),2*(p+1));
e = 1;
eps = 1e-9;
ip = 5;
y0 = -t*ip;
dd = find(abs(controlPts(:,2)-y0)<eps);

lowerNodes         = dd(1):dd(1)+noPtsX-1;
upperNodes         = dd(1+length(dd)/2):dd(1+length(dd)/2)+noPtsX-1;

for i=1:noElemsU
    sctr = ielements.element(i,:);
    iElements(e,1:p+1)   = lowerNodes(sctr);
    iElements(e,p+2:end) = upperNodes(sctr);
    e = e + 1;
end

%n5 = plot(controlPts(lowerNodes,1),controlPts(lowerNodes,2),'go');
%set(n5,'MarkerSize',6, 'MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',0.9);


%% check connectivity of interface elements

for ie=1:noElemsU
    sctr = iElements(ie,:);
    uNodes = controlPts(sctr(1:end/2));
    lNodes = controlPts(sctr(end/2+1:end));
    err    = uNodes - lNodes;
    if norm(err) > eps
        error ('WRONG')
    end
end

%% element sets

elemset0 = [];
count    = 1;
for iv=1:noElemsV
    for iu=1:noElemsU
        if mod(iv,2) ~= 0
            elemset0 = [elemset0;count];
        end
        count = count + 1;
    end
end

elemset90 = setdiff(1:noElems,elemset0);

%% node sets
eps=1e-10;
topNodes   = find(abs(controlPts(:,2)-H-R) < eps);
fixedNodes = find(abs(controlPts(:,2)+R0)  < eps);
fixedNodes(end) = [];

%% write to jive mesh

noElems  = size(C,3);

fileName = '~/code/jive/bezier/delamination/lshape/lshape-1crack-with1.mesh';

file = fopen(fileName, 'wt');

fprintf(file, '<Nodes>\n');

for i=1:length(controlPts)
    fprintf(file, '  %1d %2.4f %2.4f', i, controlPts(i,1),controlPts(i,2));
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

for e=1:noElemsU
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

%% write element groups

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

for i=1+noElems:noElems+noElemsU
    fprintf(file, '  %1d', i-1);
    
end
fprintf(file, '}\n');
fprintf(file, '</ElementGroup>\n');

% write element groups for interface elements
% cohesive elements and initial cracks

id1 = find( unique(uKnot)== initCracks(1) );
id2 = find( unique(uKnot)== initCracks(2) );

contacts  = id1:id2-1;
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

%% write node groups

fprintf(file, '<NodeGroup name="gr1">\n{');

for i=1:length(fixedNodes)
    fprintf(file, '  %1d', fixedNodes(i));
    
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="gr2">\n{');

for i=1:length(topNodes)
    fprintf(file, '  %1d', topNodes(i));
    
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fclose(file);

%%

disp('writing jive mesh done!!!')

pause

%%
% read jem-jive result and do post-processing

%vMesh=buildVMesh2DInterface(solid);

vtuFileName = '../results/lshape-1crack-with1';
resultFile = '~/code/jive/bezier/delamination/lshape/lshape-1crack-with1.dat';
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

