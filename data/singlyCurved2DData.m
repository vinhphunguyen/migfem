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

global noElems p q index elRangeU elRangeV noPtsX noPtsY uKnot vKnot De element
global controlPts weights

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

% controlPts(1:2,1) = [0;0];
% controlPts(1:2,2) = [0.5*l;0];
controlPts(1:2,1) = [0;0];
controlPts(1:2,2) = [10;0];

controlPts(1:2,3) = [0.5*L-8;h-3];
controlPts(1:2,4) = [0.5*L;h];
controlPts(1:2,5) = [0.5*L+8;h-3];

controlPts(1:2,6) = [L-10;0];
controlPts(1:2,7) = [L;0];
% controlPts(1:2,10) = [1.5*l+L;0];
% controlPts(1:2,11) = [2*l+L;0];

controlPts(4,:,:)   = 1;

curve = nrbmak(controlPts,uKnot);


%%

solid1 = nrbextrude(curve, [0,t]);

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

solid  = surfaceFromTwoCurves(curve, oCurve);

figure
hold on
nrbctrlplot(solid)
axis off


% evaluate order => bi-quadratic or cubic-quadratic

solid = nrbdegelev(solid,[2 1]);

% figure
% hold on
% nrbctrlplot(solid)
% axis off

% h-refinement in Y direction to make sure it is C^{-1} along the
% delamination path. The multiplicity must be q+1.

knots=[];
for ip=1:no-1
    dd = ip/no;
    knots = [knots dd dd dd];
end

solid     = nrbkntins(solid,{[] knots});

% initial crack

initCracks=[0.4;0.6];

% insert initial cracks knot 
%solid =nrbkntins(solid,{[initCracks(1) initCracks(1) initCracks(1)...
%                          initCracks(2) initCracks(2) initCracks(2)] []});

figure
hold on
nrbctrlplot(solid)
axis off

pause;

% nrbplot(solid,[40 2 2])
% axis off

% h-refinement

refCount = [3 0];

refCountX = refCount(1);
refCountY = refCount(2);


uKnot  = cell2mat(solid.knots(1));
vKnot  = cell2mat(solid.knots(2));

% uKnotVectorU = [0 0.1429];   
% for i=1:refCountX      
%     newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);   
%     newKnots  = {newKnotsX [] []};    
%     % h-refinement    
%     solid     = nrbkntins(solid,newKnots);    
%     uKnotVectorU     = unique([uKnotVectorU newKnotsX]);
% end
% 
% uKnotVectorU = [0.8571 1];   
% for i=1:refCountX      
%     newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);   
%     newKnots  = {newKnotsX [] []};    
%     % h-refinement    
%     solid     = nrbkntins(solid,newKnots);    
%     uKnotVectorU     = unique([uKnotVectorU newKnotsX]);
% end
% 
% uKnotVectorU = [ 0.1429 0.2857  0.7143 0.8571]; 
% for i=1:2        
%     newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);   
%     newKnots  = {newKnotsX [] []};    
%     % h-refinement    
%     solid     = nrbkntins(solid,newKnots);    
%     uKnotVectorU     = unique([uKnotVectorU newKnotsX]);
% end
% 
% uKnotVectorU = [0.2857    0.4286    0.5714    0.7143]; 
% for i=1:3        
%     newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);   
%     newKnots  = {newKnotsX [] []};    
%     % h-refinement    
%     solid     = nrbkntins(solid,newKnots);    
%     uKnotVectorU     = unique([uKnotVectorU newKnotsX]);
% end

for i=1:refCountX
    uKnotVectorU = unique(uKnot);          
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);   
    newKnots  = {newKnotsX []};    
    % h-refinement    
    solid     = nrbkntins(solid,newKnots);    
    uKnot     = cell2mat(solid.knots(1));
end

for i=1:refCountY
    vKnotVectorU = unique(vKnot);      
    % new knots along two directions (uniform)    
    newKnotsY = vKnotVectorU(1:end-1) + 0.5*diff(vKnotVectorU);   
    newKnots  = {[] newKnotsY};    
    % h-refinement    
    solid     = nrbkntins(solid,newKnots);    
    vKnot     = cell2mat(solid.knots(2));
end


% FE mesh from NURBS solids
convert2DNurbs
noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;
generateIGA2DMesh

% plot the mesh

vMesh=buildVMesh2DInterface(solid);

% Plot mesh and enriched nodes to check
figure
hold on
plot_mesh(vMesh.node,vMesh.element,'Q4','b-',1.1);
%$n5 = plot3(controlPts(:,1),controlPts(:,2),'go');
%set(n5,'MarkerSize',6, 'MarkerFaceColor','g','MarkerEdgeColor','k','LineWidth',0.9);
axis off

%% Bezier extraction operators


[C,Cxi,Cet]  = bezierExtraction2D(uKnot,vKnot,p,q);

% build mesh of interface elements
iElementS   = buildGA1DMesh (uKnot,p);
iElements   = zeros(noElemsU*(no-1),2*(p+1));

e = 1;
eps = 1e-9;
% loop over plies

for ip=1:no-1
    y0 = t/no*ip;
    dd = find(abs(controlPts(:,2)-y0)<eps);
    
    lowerNodes         = dd(1):dd(1)+noPtsX-1;
    upperNodes         = dd(1+length(dd)/2):dd(1+length(dd)/2)+noPtsX-1;
    
    for i=1:noElemsU
        sctr = iElementS.element(i,:);
        iElements(e,1:p+1)   = lowerNodes(sctr);
        iElements(e,p+2:end) = upperNodes(sctr);
        e = e + 1;
    end
    
    %n5 = plot(controlPts(lowerNodes,1),controlPts(lowerNodes,2),'go');
    %set(n5,'MarkerSize',6, 'MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',0.9);
end
  

%% element sets

matMap = ones(noElems,1);
elemset0 = [];
count    = 1;

for iv=1:noElemsV
    for iu=1:noElemsU
        if ismember(iv,[1:20])
            elemset0 = [elemset0;count];
        end
        count = count + 1;
    end
end


elemset90 = setdiff(1:noElems,elemset0);

matMap(elemset90,1) = 2;


%% node sets

eps = 1e-10;
fixedNodes    = find(abs(controlPts(:,1)) < eps);
forcedNodes   = find(abs(controlPts(:,1)-(L)) < eps);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write to jive mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noElems  = size(C,3);

fileName = '~/code/jive/bezier/delamination/curved3d/curved2d.mesh';

file = fopen(fileName, 'wt');

fprintf(file, '<Nodes>\n');

for i=1:length(controlPts)
    fprintf(file, '  %1d %2.6f %2.6f', i, ...
        controlPts(i,1),controlPts(i,2));
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
    for ie=1:noElemsU
        Ce = Cxi(:,:,ie);
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
    for ie=1:noElemsU
        Ce = Cxi(:,:,ie);
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
    for ie=1:noElemsU
        Ce = Cxi(:,:,ie);
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

% write element groups for interface elements
% cohesive elements and initial cracks
% Initial crack between ply 3 and 4

id1 = find( unique(uKnot)== initCracks(1) );
id2 = find( unique(uKnot)== initCracks(2) );

contacts  = 22*noElemsU+ [id1:id2-1];
cohesives = setdiff(1:(no-1)*noElemsU, contacts);
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

e1 = 110e3;
e2 = 10e3;
e3 = 10e3;
nu12 = 0.27;
nu23 = 0.30;
nu31 = 0.27;
g12  = 5e3;
g23  = e2/(2+2*nu23);
g31  = 5e3;

dirs = [0,0];

for i =1:length(dirs)
    material = createOrthotropicMaterial(dirs(i),2,e1,e2,e3,nu12,nu23,nu31,g12,g23,g31);
    materials{i} = material;
end

%vMesh=buildVMesh2DInterface(solid);

vtuFileName = '../results/curved2d';
resultFile = '~/code/jive/bezier/delamination/curved3d/curved2d.dat';
V = xml_parseany(fileread(resultFile));

noStep = size(V,2)/2;
step   = 1;
start  = 102;

for it=start:step:noStep
    vtuFile = sprintf('../results/%s%d',vtuFileName,it);
    uData=V{2*it-1}.Section{1};   % displacement
    dData=V{2*it}.Section{1}; % damage
    disp = str2num(uData.CONTENT);
    dam  = str2num(dData.CONTENT);
    damage = zeros(noCtrPts,1);
    damage(dam(:,1),1) = dam(:,2);
    U    = disp(:,2:end);
    ok   = writeVTKForJive(solid,vMesh,vtuFile,U,damage,materials,matMap);    
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

for i = 1:step:noStep
    vtuFile = sprintf('%s%d%s',vtuFileName,i,'.vtu');
    fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''0'' timestep=''%d''/>\n',vtuFile,i);
end

fprintf(pvdFile,'</Collection>\n');
fprintf(pvdFile,'</VTKFile>\n');

fclose(pvdFile);

