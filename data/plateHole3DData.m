% data file for the 3D plate hole
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

%% material

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

dirs = [0,-45,90,45];

for i =1:length(dirs)
    material = createOrthotropicMaterial(dirs(i),3,e1,e2,e3,nu12,nu23,nu31,g12,g23,g31);
    materials{i} = material;
end

%% geometry data
d = 25.4;
R = d/2;
L = 20*d;
W = 5*d;
t = 4/2; % thickness
no = 4;  % number of plies

ndim = 3;
%% geometry of one ply

x0 = R*tan(pi/8);

controlPts = zeros(4,5,2);

controlPts(1:2,1,1) = [R;0];
controlPts(1:2,2,1) = [R;x0];
controlPts(1:2,3,1) = [R*sqrt(2)/2; R*sqrt(2)/2];
controlPts(1:2,4,1) = [x0;R];
controlPts(1:2,5,1) = [0;R];

controlPts(1:2,1,2) = [W/2;0];
controlPts(1:2,2,2) = [W/2;W/4];
controlPts(1:2,3,2) = [W/2;W/2];
controlPts(1:2,4,2) = [W/4;W/2];
controlPts(1:2,5,2) = [0;W/2];


cont =1/sqrt(2);

controlPts(4,:,:)   = 1;
controlPts(4,2,1)   = cont;
controlPts(4,4,1)   = cont;

controlPts(1:2,2,1)   = controlPts(1:2,2,1) *cont;
controlPts(1:2,4,1)   = controlPts(1:2,4,1) *cont;

uKnot = [0 0 0 0.5 0.5 1 1 1];
vKnot = [0 0 1 1];

%% build NURBS object

surf1 = nrbmak(controlPts,{uKnot vKnot});

%% second patch

controlPts = zeros(4,2,2);

controlPts(1:2,1,1) = [W/2;0];
controlPts(1:2,2,1) = [L/2;0];
controlPts(1:2,1,2) = [W/2;W/2];
controlPts(1:2,2,2) = [L/2;W/2];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

%% build NURBS object for the second patch

surf2 = nrbmak(controlPts,{uKnot vKnot});

figure
hold on
nrbctrlplot(surf1);
nrbctrlplot(surf2);
axis equal
axis off

% order elevate

surf1 = nrbdegelev(surf1,[0 1]); % to bi-quadratic
surf2 = nrbdegelev(surf2,[1 1]); % to bi-quadratic

% extrusion to make solids

solid1 = nrbextrude(surf1, [0,0,t]);
solid2 = nrbextrude(surf2, [0,0,t]);

% order elevation
solid1 = nrbdegelev(solid1,[0 0 1]); % quadratic in thickness direction
solid2 = nrbdegelev(solid2,[0 0 1]); % quadratic in thickness direction

knots=[];
for ip=1:no-1
    dd = ip/no;
    knots = [knots dd dd dd];
end

solid1     = nrbkntins(solid1,{[] [] knots});
solid2     = nrbkntins(solid2,{[] [] knots});

%% h-refinement

uKnot     = cell2mat(solid1.knots(1));
vKnot     = cell2mat(solid1.knots(2));

refineCountX = 3;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
    % new knots along two directions (uniform)
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {newKnotsX newKnotsY []};
    
    % h-refinement
    solid1    = nrbkntins(solid1,newKnots);
    uKnot     = cell2mat(solid1.knots(1));
    vKnot     = cell2mat(solid1.knots(2));
end

uKnot     = cell2mat(solid2.knots(1));
vKnot     = cell2mat(solid2.knots(2));

for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
    % new knots along two directions (uniform)
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {[] newKnotsY []};
    
    % h-refinement
    solid2    = nrbkntins(solid2,newKnots);
    uKnot     = cell2mat(solid2.knots(1));
    vKnot     = cell2mat(solid2.knots(2));
end

refineCountX1 = 3;
for i=1:refineCountX1
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
    % new knots along two directions (uniform)
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {newKnotsX [] []};
    
    % h-refinement
    solid2    = nrbkntins(solid2,newKnots);
    uKnot     = cell2mat(solid2.knots(1));
    vKnot     = cell2mat(solid2.knots(2));
end

% figure
% hold on
% nrbctrlplot(solid);
% axis equal
% axis off

%% convert to IGA format

mesh1  = buildIGA3DMesh(solid1);
mesh2  = buildIGA3DMesh(solid2);
vMesh1 = buildVMesh3DInterface(solid1);
vMesh2 = buildVMesh3DInterface(solid2);
iMesh1 = buildInterfaceMesh(mesh1,no,t);
iMesh2 = buildInterfaceMesh(mesh2,no,t);

% element nodes in patch2 are numbered from mesh1.noPts;
mesh2.globElems  = mesh2.globElems + mesh1.noPts;
iMesh2.globElems = iMesh2.globElems + mesh1.noPts;

meshes{1}  = mesh1;
meshes{2}  = mesh2;
imeshes{1} = iMesh1;
imeshes{2} = iMesh2;
vmeshes{1} = vMesh1;
vmeshes{2} = vMesh2;

meshData.mesh    = meshes;
meshData.imesh   = imeshes;
meshData.vmesh   = vmeshes;
meshData.noElems = meshes{1}.noElems + meshes{2}.noElems;
meshData.noPts   = meshes{1}.noPts   + meshes{2}.noPts;

% Plot mesh and enriched nodes to check
figure
hold on
plot_mesh(vMesh1.node,vMesh1.element,'B8','b-',1.1);
plot_mesh(vMesh2.node,vMesh2.element,'B8','b-',1.1);
n5 = plot3(mesh1.controlPts(:,1),mesh1.controlPts(:,2),mesh1.controlPts(:,2),'go');
n6 = plot3(mesh2.controlPts(:,1),mesh2.controlPts(:,2),mesh2.controlPts(:,2),'ro');
set(n5,'MarkerSize',6, 'MarkerFaceColor','g','MarkerEdgeColor','k','LineWidth',0.9);
axis off
%% node set

eps = 1e-10;
fixedXNodes     = find(abs(mesh1.controlPts(:,1)-0) < eps);
fixedYNodes1    = find(abs(mesh1.controlPts(:,2)-0) < eps);
fixedYNodes2    = find(abs(mesh2.controlPts(:,2)-0) < eps);
fixedZNodes1    = find(abs(mesh1.controlPts(:,3)-0) < eps);
fixedZNodes2    = find(abs(mesh2.controlPts(:,3)-0) < eps);
forcedNodes     = find(abs(mesh2.controlPts(:,1)-L/2) < eps);

forcedNodes     = forcedNodes  + mesh1.noPts;
fixedYNodes2    = fixedYNodes2 + mesh1.noPts;
fixedZNodes2    = fixedZNodes2 + mesh1.noPts;
fixedYNodes     = [fixedYNodes1;fixedYNodes2];
fixedZNodes     = [fixedZNodes1;fixedZNodes2];

aa              = find(abs(mesh1.controlPts(:,1)-W/2) < eps);
bb              = find(abs(mesh2.controlPts(:,1)-W/2) < eps) +  mesh1.noPts;

%% element sets
%% the stacking is [45/90/-45/0]s
%A half model through the thickness was used due to the sym- metry of the layup

matMap1 = ones(mesh1.noElems,1);
matMap2 = ones(mesh2.noElems,1);

elemset1 = [];
elemset2 = [];
elemset3 = [];
elemset4 = [];
count    = 1;
for iw=1:mesh1.noElemsW
    for iv=1:mesh1.noElemsV
        for iu=1:mesh1.noElemsU
            if     iv == 1
                elemset1 = [elemset1;count];
            elseif iv == 2                 
                elemset2 = [elemset2;count];
            elseif iv == 3                 
                elemset3 = [elemset3;count];
            else                                           
                elemset4 = [elemset4;count];
            end
            count = count + 1;
        end
    end
end

matMap1(elemset1,1) = 1;
matMap1(elemset2,1) = 2;
matMap1(elemset3,1) = 3;
matMap1(elemset4,1) = 4;

s1 = [];
s2 = [];
s3 = [];
s4 = [];
count    = 1;
for iw=1:mesh2.noElemsW
    for iv=1:mesh2.noElemsV
        for iu=1:mesh2.noElemsU
            if     iv == 1
                s1 = [s1;count];
            elseif iv == 2                 
                s2 = [s2;count];
            elseif iv == 3                 
                s3 = [s3;count];
            else                                           
                s4 = [s4;count];
            end
            count = count + 1;
        end
    end
end

matMap2(s1,1) = 1;
matMap2(s2,1) = 2;
matMap2(s3,1) = 3;
matMap2(s4,1) = 4;

elemset1 = [elemset1;s1 + mesh1.noElems];
elemset2 = [elemset2;s2 + mesh1.noElems];
elemset3 = [elemset3;s3 + mesh1.noElems];
elemset4 = [elemset4;s4 + mesh1.noElems];

meshData.matMap{1}=matMap1;
meshData.matMap{2}=matMap2;


%% write to jive mesh


fileName = '~/code/jive/bezier/delamination/hole/hole.mesh';

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

% write solid elements

ie = 1;
for im=1:length(meshes)
    element = meshes{im}.globElems;
    for i=1:size(element,1)
        fprintf(file, '  %1d %1d', ie-1, element(i,:) ); 
        fprintf(file, ';\n');
        ie = ie + 1;
    end
end

% write interface elements

ie = 1;
for im=1:length(imeshes)
    element = imeshes{im}.globElems;
    for i=1:size(element,1)
        fprintf(file, '  %1d %1d', meshData.noElems+ie-1, element(i,:) );
        fprintf(file, ';\n');
        ie = ie + 1;
    end
end

fprintf(file, '</Elements>\n\n');

% write Bezier extractors
% first for solid elements
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

e = 1;
for im=1:length(meshes)
    Cxe = meshes{im}.Cxe;
    for il=1:no-1
        for ie=1:meshes{im}.noElemsU*meshes{im}.noElemsV
            Ce = Cxe(:,:,ie);
            [row,col] = find(Ce);
            fprintf(file, '  %1d ', meshData.noElems+e-1);
            for i=1:length(row)
                fprintf(file, '%1d ', row(i)-1);
            end
            fprintf(file, ';\n');
            e=e+1;
        end
    end
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

e = 1;
for im=1:length(meshes)
    Cxe = meshes{im}.Cxe;
    for il=1:no-1
        for ie=1:meshes{im}.noElemsU*meshes{im}.noElemsV
            Ce = Cxe(:,:,ie);
            [row,col] = find(Ce);
            fprintf(file, '  %1d ', meshData.noElems+e-1);
            for i=1:length(row)
                fprintf(file, '%1d ', col(i)-1);
            end
            fprintf(file, ';\n');
            e=e+1;
        end
    end
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
            fprintf(file, '%2.2f ', val(i));
        end
        fprintf(file, ';\n');
    end
end

e = 1;
for im=1:length(meshes)
    Cxe = meshes{im}.Cxe;
    for il=1:no-1
        for ie=1:meshes{im}.noElemsU*meshes{im}.noElemsV
            Ce = Cxe(:,:,ie);
            [row,col,val] = find(Ce);
            fprintf(file, '  %1d ', meshData.noElems+e-1);
            for i=1:length(row)
                fprintf(file, '%2.2f ', val(i));
            end
            fprintf(file, ';\n');
            e=e+1;
        end
    end
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

ie = 1;
for im=1:length(imeshes)
    weights = meshes{im}.weights;
    element = imeshes{im}.locElems;
    for e=1:size(element,1)
        w = weights(element(e,:));
        fprintf(file, '  %1d ',meshData.noElems+ie-1); ie = ie + 1;
        for j=1:length(w)
            fprintf(file, '%2.4f ', w(j));
        end
        fprintf(file, ';\n');
    end
end

fprintf(file, ' </Column>\n');
fprintf(file, '</ElementDatabase>\n\n');

%% write element groups

fprintf(file, '<ElementGroup name="solid">\n{');

for i=1:meshData.noElems
    fprintf(file, '  %1d', i-1);
    
end
fprintf(file, '}\n');
fprintf(file, '</ElementGroup>\n');

fprintf(file, '<ElementGroup name="solid1">\n{');

for i=1:length(elemset1)
    fprintf(file, '  %1d', elemset1(i)-1);
    
end
fprintf(file, '}\n');
fprintf(file, '</ElementGroup>\n');

fprintf(file, '<ElementGroup name="solid2">\n{');

for i=1:length(elemset2)
    fprintf(file, '  %1d', elemset2(i)-1);
    
end
fprintf(file, '}\n');
fprintf(file, '</ElementGroup>\n');

fprintf(file, '<ElementGroup name="solid3">\n{');

for i=1:length(elemset3)
    fprintf(file, '  %1d', elemset3(i)-1);
    
end
fprintf(file, '}\n');
fprintf(file, '</ElementGroup>\n');

fprintf(file, '<ElementGroup name="solid4">\n{');

for i=1:length(elemset4)
    fprintf(file, '  %1d', elemset4(i)-1);
    
end
fprintf(file, '}\n');
fprintf(file, '</ElementGroup>\n');

fprintf(file, '<ElementGroup name="interface">\n{');

for i=1+meshData.noElems : meshData.noElems+...
                           size(imeshes{1}.globElems,1)+size(imeshes{2}.globElems,1)
    fprintf(file, '  %1d', i-1);
    
end
fprintf(file, '}\n');
fprintf(file, '</ElementGroup>\n');


%% write node groups

fprintf(file, '<NodeGroup name="fixX">\n{');

for i=1:length(fixedXNodes)
    fprintf(file, '  %1d', fixedXNodes(i));
    
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="fixY">\n{');

for i=1:length(fixedYNodes)
    fprintf(file, '  %1d', fixedYNodes(i));
    
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="fixZ">\n{');

for i=1:length(fixedZNodes)
    fprintf(file, '  %1d', fixedZNodes(i));
    
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="force">\n{');

for i=1:length(forcedNodes)
    fprintf(file, '  %1d', forcedNodes(i));
    
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="aa">\n{');

for i=1:length(aa)
    fprintf(file, '  %1d', aa(i));
    
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="bb">\n{');

for i=1:length(bb)
    fprintf(file, '  %1d', bb(i));
    
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fclose(file);


disp('writing jive mesh done!!!')

pause

%%
% read jem-jive result and do post-processing


vtuFileName = '../results/plateHole3d';
resultFile = '~/code/jive/bezier/delamination/hole/hole.dat';
V = xml_parseany(fileread(resultFile));


noStep = size(V,2)/2;
% if noStep=1 then V{it}.Section does not work
for it=1:noStep        
    uData=V{2*it-1}.Section{1};   % displacement
    dData=V{2*it}.Section{1}; % damage
    disp = str2num(uData.CONTENT);
    dam  = str2num(dData.CONTENT);
    damage = zeros(meshData.noPts*ndim,1);
    damage(dam(:,1),1) = dam(:,2);
    U    = disp(:,2:end);
    for ip=1:2
        vtuFile = strcat(vtuFileName,'-mesh',num2str(ip),'-',num2str(it));
        %figure; hold on;
        ok      = plotStress3DForPatch(meshData,ip,vtuFile,U,damage,materials);
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

fprintf(pvdFile,'<VTKFile byte_order="LittleEndian" type="Collection" version="0.1">\n');
fprintf(pvdFile,'<Collection>\n');

for it= 1:noStep
    for im=1:length(meshes)       
        vtuFile = sprintf('%s%s%d%s%d%s',vtuFileName,'-mesh',im,'-',it,'.vtu');
        fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''%d'' timestep=''%d''/>\n',...
            vtuFile,im,it);
    end
end

fprintf(pvdFile,'</Collection>\n');
fprintf(pvdFile,'</VTKFile>\n');

fclose(pvdFile);


