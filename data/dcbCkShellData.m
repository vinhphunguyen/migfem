% data file for the 3D double cantilever beam problem
% modelled by Kirchhoff Love shell elements.
% VP Nguyen, 7 May 2013

addpath ../fem_util/;
addpath ../nurbs-geopdes/inst/
addpath ../nurbs-util/
addpath ../meshing/
addpath ../fem-functions/
addpath ../post-processing/
addpath ../xiga/
addpath ../nurbs-geopdes/inst/
addpath ../delamination/
addpath ../structural-mechanics/
addpath ../xml_toolbox/
addpath ../C_files/

clear all

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

%%
L = 100;
W = 20;
b = 3; % thickness
a0= 30;  % initial crack length

%% shell surface

% initial mesh- linear x linear

controlPts          = zeros(4,2,2);

controlPts(1:3,1,1) = [0;0;0];
controlPts(1:3,2,1) = [L;0;0];
controlPts(1:3,1,2) = [0;W;0];
controlPts(1:3,2,2) = [L;W;0];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];


%% build NURBS object
solid = nrbmak(controlPts,{uKnot vKnot});

% evaluate order

solid = nrbdegelev(solid,[1 1]);


% initial crack

solid     = nrbkntins(solid,{1/L*[a0] []});

uKnot     = cell2mat(solid.knots(1));
vKnot     = cell2mat(solid.knots(2));

% h-refinement

% refine the initial crack
refineCountX = 3;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    id = find( uKnotVectorU == a0/L);
    subKnot = uKnotVectorU(1:id);
    
    newKnotsY = [];
    
    % new knots along two directions (uniform)
    
    newKnotsX = subKnot(1:end-1) + 0.5*diff(subKnot);
    newKnots  = {newKnotsX newKnotsY };
    % h-refinement
    solid     = nrbkntins(solid,newKnots);
    uKnot     = cell2mat(solid.knots(1));
    vKnot     = cell2mat(solid.knots(2));
end

refineCountX = 7;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    id = find( uKnotVectorU == a0/L);
    subKnot = uKnotVectorU(id:end);
    newKnotsY = [];
    
    % new knots along two directions (uniform)
    
    newKnotsX = subKnot(1:end-1) + 0.5*diff(subKnot);
    newKnots  = {newKnotsX newKnotsY};
    
    % h-refinement
    solid     = nrbkntins(solid,newKnots);
    uKnot     = cell2mat(solid.knots(1));
    vKnot     = cell2mat(solid.knots(2));
end


refineCountY = 2;
for i=1:refineCountY
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);    
    % new knots along two directions (uniform)
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {[] newKnotsY};
    solid     = nrbkntins(solid,newKnots);
    uKnot     = cell2mat(solid.knots(1));
    vKnot     = cell2mat(solid.knots(2));
end

figure
hold on
nrbctrlplot(solid);
axis equal
axis off

%% convert to IGA format

[mesh] = buildIGA2DMeshForSurface (solid);

% mesh of the upper shell
mesh1  = mesh;
mesh1.globElems = mesh1.globElems + size(mesh.controlPts,1);

meshes{1}=mesh;
meshes{2}=mesh1;

noCtrPts = mesh.noPts + mesh1.noPts;

%% Bezier extraction operators

[Cxe]           = bezierExtraction2D(mesh.uKnot,mesh.vKnot,mesh.p,mesh.q);

eps = 1e-10;
fixedNodes      = find(abs(mesh.controlPts(:,1)-L) < eps);
forcedNodesL    = find(abs(mesh.controlPts(:,1)-0) < eps);

fixedNodes   = [fixedNodes; fixedNodes-1];

fixedNodes   = [fixedNodes; fixedNodes+mesh.noPts];
forcedNodesU = [forcedNodesL+mesh.noPts];

% build mesh of interface elements


lowerNodes  = 1:size(mesh.controlPts,1);
upperNodes  = size(mesh.controlPts,1)+1:size(mesh.controlPts,1)*2;

iElements   = zeros(mesh.noElemsU*mesh.noElemsV,2*(mesh.p+1)*(mesh.q+1));

for e=1:mesh.noElems
    sctr = mesh.globElems(e,:);
    iElements(e,1:(mesh.p+1)*(mesh.q+1))     = lowerNodes(sctr);
    iElements(e,(mesh.p+1)*(mesh.q+1)+1:end) = upperNodes(sctr);
end

%% write to jive mesh


fileName = '~/code/jive/bezier/delamination/dcb3d/dcbShell.mesh';

file = fopen(fileName, 'wt');

fprintf(file, '<Nodes>\n');

ip = 1;
for im=1:length(meshes)
    controlPts = meshes{im}.controlPts;
    for i=1:length(controlPts)
        fprintf(file, '  %1d %2.6f %2.6f %2.6f', ip, ...
            controlPts(i,1),controlPts(i,2),controlPts(i,3));
        fprintf(file, ';\n');
        ip = ip + 1;
    end
end

fprintf(file, '</Nodes>\n\n');

fprintf(file, '<Elements>\n');

% write shell elements

ie = 1;
for im=1:length(meshes)
    element = meshes{im}.globElems;
    for i=1:size(element,1)
        fprintf(file, '  %1d %1d', ie-1, element(i,:) ); fprintf(file, ';\n');
        ie = ie + 1;
    end
end

% write interface elements

for i=1:size(iElements,1)
    fprintf(file, '  %1d %1d', meshes{1}.noElems*2+i-1, iElements(i,:) );
    fprintf(file, ';\n');
end

fprintf(file, '</Elements>\n\n');

% write Bezier extractors
% first for shell elements
% then for interface elements

fprintf(file, '<ElementDatabase name="C">\n');

fprintf(file, ' <Column name = "irows" type = "int">\n');

ie = 1;
for im=1:length(meshes)
    C = meshes{im}.C;
    for e=1:meshes{im}.noElems
        Ce = C(:,:,e);
        [row,col] = find(Ce);
        fprintf(file, '  %1d ', ie-1); ie = ie + 1;
        for i=1:length(row)
            fprintf(file, '%1d ', row(i)-1);
        end
        fprintf(file, ';\n');
    end
end

for e=1:size(iElements,1)
    Ce = Cxe(:,:,e);
    [row,col] = find(Ce);
    fprintf(file, '  %1d ', meshes{1}.noElems*2 + e-1);
    for i=1:length(row)
        fprintf(file, '%1d ', row(i)-1);
    end
    fprintf(file, ';\n');
end

fprintf(file, ' </Column>\n');

fprintf(file, ' <Column name = "jcols" type = "int">\n');

ie = 1;
for im=1:length(meshes)
    C = meshes{im}.C;
    for e=1:meshes{im}.noElems
        Ce = C(:,:,e);
        [row,col] = find(Ce);
        
        fprintf(file, '  %d ', ie-1); ie = ie + 1;
        for i=1:length(row)
            fprintf(file, '%1d ', col(i)-1);
        end
        fprintf(file, ';\n');
    end
end

for e=1:size(iElements,1)
    Ce = Cxe(:,:,e);
    [row,col] = find(Ce);
    
    fprintf(file, '  %d ', meshes{1}.noElems*2+e-1);
    for i=1:length(row)
        fprintf(file, '%1d ', col(i)-1);
    end
    fprintf(file, ';\n');
end

fprintf(file, ' </Column>\n');

fprintf(file, ' <Column name = "values" type = "float">\n');


ie = 1;
for im=1:length(meshes)
    C = meshes{im}.C;
    for e=1:meshes{im}.noElems
        Ce = C(:,:,e);
        [row,col,val] = find(Ce);
        
        fprintf(file, '  %d ', ie-1); ie = ie + 1;
        for i=1:length(row)
            fprintf(file, '%2.1f ', val(i));
        end
        fprintf(file, ';\n');
    end
end

for e=1:size(iElements,1)
    Ce = Cxe(:,:,e);
    [row,col,val] = find(Ce);
    
    fprintf(file, '  %d ', meshes{1}.noElems*2+e-1);
    for i=1:length(row)
        fprintf(file, '%2.1f ', val(i));
    end
    fprintf(file, ';\n');
end
fprintf(file, ' </Column>\n');

% write weights

fprintf(file, ' <Column name = "weights" type = "float">\n');

ie = 1;
for im=1:length(meshes)
    weights = meshes{im}.weights;
    element = meshes{im}.locElems;
    for e=1:meshes{im}.noElems
        w = weights(element(e,:));
        fprintf(file, '  %1d ',ie-1); ie = ie + 1;
        for j=1:length(w)
            fprintf(file, '%2.4f ', w(j));
        end
        fprintf(file, ';\n');
    end
end

for e=1:size(iElements,1)
    sctr = iElements(e,:);
    w = weights(sctr(1:0.5*end));
    fprintf(file, '  %1d ',meshes{1}.noElems*2+e-1);
    for j=1:length(w)
        fprintf(file, '%2.4f ', w(j));
    end
    fprintf(file, ';\n');
end

fprintf(file, ' </Column>\n');
fprintf(file, '</ElementDatabase>\n\n');

% write element groups

fprintf(file, '<ElementGroup name="solid">\n{');

for i=1:meshes{1}.noElems*2
    fprintf(file, '  %1d', i-1);
    
end
fprintf(file, '}\n');
fprintf(file, '</ElementGroup>\n');

fprintf(file, '<ElementGroup name="interface">\n{');

for i=1+meshes{1}.noElems*2:meshes{1}.noElems*2+size(iElements,1)
    fprintf(file, '  %1d', i-1);
    
end
fprintf(file, '}\n');
fprintf(file, '</ElementGroup>\n');

% differentiate contact elements and cohesive elements

id = find(unique(uKnot)==a0/L) - 1;

surf = zeros(mesh.noElemsV,mesh.noElemsU);
count = 1;
for j=1:mesh.noElemsV
    for i=1:mesh.noElemsU
        surf(j,i) = count;
        count = count + 1;
    end
end

contacts  = surf(:,1:id);
cohesives = surf(:,id+1:end);

contacts  = contacts - 1;
cohesives = cohesives - 1;

contacts  = contacts(:);
cohesives = cohesives(:);

fprintf(file, '<ElementGroup name="contacts">\n{');

for i=1:length(contacts)
    fprintf(file, '  %1d', contacts(i) + meshes{1}.noElems*2);
    
end
fprintf(file, '}\n');
fprintf(file, '</ElementGroup>\n');

fprintf(file, '<ElementGroup name="cohesives">\n{');

for i=1:length(cohesives)
    fprintf(file, '  %1d', cohesives(i) + meshes{1}.noElems*2);
    
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

fprintf(file, '<NodeGroup name="fup">\n{');

for i=1:length(forcedNodesU)
    fprintf(file, '  %1d', forcedNodesU(i));
    
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="fdown">\n{');

for i=1:length(forcedNodesL)
    fprintf(file, '  %1d', forcedNodesL(i));
    
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fclose(file);


disp('writing jive mesh done!!!')

pause

%%
% read jem-jive result and do post-processing

e1 = 110e3;
e2 = 10e3;
e3 = 10e3;
nu12 = 0.27;
nu23 = 0.30;
nu31 = 0.27;
g12  = 5e3;
g23  = e2/(2+2*nu23);
g31  = 5e3;

dirs = [0];

for i =1:length(dirs)
    material = createOrthotropicMaterial(dirs(i),2,e1,e2,e3,nu12,nu23,nu31,g12,g23,g31);
    materials{i} = material;
end

matMap = ones(mesh.noElems,1);

vMesh=buildVisualizationMesh2D(solid);

vtuFileName = '../results/dcbShell';
resultFile = '~/code/jive/bezier/delamination/dcb3d/dcbShell.dat';
V = xml_parseany(fileread(resultFile));

noStep = size(V,2)/2;

for it=1:noStep
    uData=V{2*it-1}.Section{1};   % displacement
    dData=V{2*it}.Section{1}; % damage
    disp = str2num(uData.CONTENT);
    dam  = str2num(dData.CONTENT);
    damage = zeros(noCtrPts,1);
    damage(dam(:,1),1) = dam(:,2);
    U    = disp(:,2:end);
    for im=1:length(meshes)
        mesh = meshes{im};
        vtuFile = sprintf('../results/%s%s%d%s%d',vtuFileName,'-mesh-',im,'-',it);
        ok   = writeVTKShellsForJive(mesh,vMesh,vtuFile,U,damage,materials,matMap);
    end
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
vtuFileName = '../results/dcbShell';

fprintf(pvdFile,'<VTKFile byte_order="LittleEndian" type="Collection" version="0.1">\n');
fprintf(pvdFile,'<Collection>\n');

for it= 1:noStep
    for im=1:length(meshes)       
        vtuFile = sprintf('../results/%s%s%d%s%d%s',vtuFileName,'-mesh-',im,'-',it,'.vtu');
        fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''%d'' timestep=''%d''/>\n',...
            vtuFile,im,it);
    end
end

fprintf(pvdFile,'</Collection>\n');
fprintf(pvdFile,'</VTKFile>\n');

fclose(pvdFile);

