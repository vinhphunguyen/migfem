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

dirs = [90,0,0];

for i =1:length(dirs)
    material = createOrthotropicMaterial(dirs(i),3,e1,e2,e3,nu12,nu23,nu31,g12,g23,g31);
    materials{i} = material;
end

%% geometry data
R = 8;
L = 20;
W = 20;
t = 0.5; % thickness
no = 2;  % number of plies

ndim = 3;
%% geometry of one ply

x0 = R*tan(pi/8);

controlPts = zeros(4,4,3);

controlPts(1:2,1,1) = [0;0];
controlPts(1:2,2,1) = [0;0];
controlPts(1:2,3,1) = [0;0];
controlPts(1:2,4,1) = [0;0];

controlPts(1:2,1,2) = [R;0];
controlPts(1:2,2,2) = [R;x0];
controlPts(1:2,3,2) = [x0; R];
controlPts(1:2,4,2) = [0;R];

controlPts(1:2,1,3) = [W;0];
controlPts(1:2,2,3) = [W;W];
controlPts(1:2,3,3) = [W;W];
controlPts(1:2,4,3) = [0;W];


cont =0.5*(1+1/sqrt(2));

controlPts(4,:,:)   = 1;
controlPts(4,2,2)   = cont;
controlPts(4,3,2)   = cont;

controlPts(1:2,2,2)   = controlPts(1:2,2,2) *cont;
controlPts(1:2,3,2)   = controlPts(1:2,3,2) *cont;

uKnot = [0 0 0 0.5 1 1 1];
vKnot = [0 0 0.5 1 1];

%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot});

convert2DNurbs
plotMesh (controlPts,weights,uKnot,vKnot,p,q,100,'b-','try.eps');
hold on
nrbctrlplot(solid);
axis equal
axis off

% order elevate

solid = nrbdegelev(solid,[0 1]); % to bi-quadratic


%% h-refinement

uKnot     = cell2mat(solid.knots(1));
vKnot     = cell2mat(solid.knots(2));

refineCountX = 3;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
    % new knots along two directions (uniform)
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {newKnotsX []};
    
    % h-refinement
    solid     = nrbkntins(solid,newKnots);
    uKnot     = cell2mat(solid.knots(1));
    vKnot     = cell2mat(solid.knots(2));
end

for i=1:refineCountX-1
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    id = find( uKnotVectorV == 0.5);
    subKnot = uKnotVectorV(1:id);
    % new knots along two directions (uniform)
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = subKnot(1:end-1) + 0.5*diff(subKnot);
    newKnots  = {[] newKnotsY };
    
    % h-refinement
    solid     = nrbkntins(solid,newKnots);
    uKnot     = cell2mat(solid.knots(1));
    vKnot     = cell2mat(solid.knots(2));
end

solid     = nrbkntins(solid,{[] [0.75]});

uKnot     = cell2mat(solid.knots(1));
vKnot     = cell2mat(solid.knots(2));

for i=1:refineCountX+1
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    id = find( uKnotVectorV == 0.5);
    id1 = find( uKnotVectorV == 0.75);
    subKnot = uKnotVectorV(id:id1);
    % new knots along two directions (uniform)
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = subKnot(1:end-1) + 0.5*diff(subKnot);
    newKnots  = {[] newKnotsY };
    
    % h-refinement
    solid     = nrbkntins(solid,newKnots);
    uKnot     = cell2mat(solid.knots(1));
    vKnot     = cell2mat(solid.knots(2));
end

for i=1:2
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    id = find( uKnotVectorV == 0.5);
    id1 = find( uKnotVectorV == 0.75);
    subKnot = uKnotVectorV(id1:end);
    % new knots along two directions (uniform)
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = subKnot(1:end-1) + 0.5*diff(subKnot);
    newKnots  = {[] newKnotsY };
    
    % h-refinement
    solid     = nrbkntins(solid,newKnots);
    uKnot     = cell2mat(solid.knots(1));
    vKnot     = cell2mat(solid.knots(2));
end

surf=solid;
convert2DNurbs
plotMesh (controlPts,weights,uKnot,vKnot,p,q,100,'b-','try.eps');

%extrusion to make solids
solid = nrbextrude(solid, [0,0,t]);
%order elevation
solid = nrbdegelev(solid,[0 0 1]); % quadratic in thickness direction
%knot insertion to have C^-1 at coord = 0.2 (t=0.5)
solid     = nrbkntins(solid,{[] [] [1/5 1/5 2/5 2/5 2/5]});

% figure
% hold on
% nrbctrlplot(solid);
% axis equal
% axis off

%% convert to IGA format

mesh1  = buildIGA3DMesh(solid);
vMesh1 = buildVMesh3DInterface(solid);

noElemsU = mesh1.noElemsU;
noElemsV = mesh1.noElemsV;
p        = mesh1.p;
q        = mesh1.q;
uKnot    = mesh1.uKnot;
vKnot    = mesh1.vKnot;
noPtsX   = mesh1.noPtsX;
noPtsY   = mesh1.noPtsY;
controlPts=mesh1.controlPts;


iElements   = zeros(noElemsU*noElemsV,2*(p+1)*(q+1));
iElementS   = buildIGA2DMesh (surf);

e = 1;

y0 = 0.2;
delaminationNodes  =  find(abs(controlPts(:,3) - y0 ) <1e-10);
mm                 = 0.5*length(delaminationNodes);
lowerNodes         = delaminationNodes(1:mm);
upperNodes         = delaminationNodes(mm+1:end);

for i=1:noElemsU*noElemsV
    sctr = iElementS.globElems(i,:);
    iElements(e,1:(p+1)*(q+1))     = lowerNodes(sctr);
    iElements(e,(p+1)*(q+1)+1:end) = upperNodes(sctr);
    e = e + 1;
end


iMesh.locElems  = iElements;
iMesh.globElems = iElements;


meshes{1}  = mesh1;
imeshes{1} = iMesh;
vmeshes{1} = vMesh1;


meshData.mesh    = meshes;
meshData.imesh   = imeshes;
meshData.vmesh   = vmeshes;
meshData.noElems = meshes{1}.noElems;
meshData.noPts   = meshes{1}.noPts;

%% node set

eps = 1e-10;
fixedXNodes     = find(abs(mesh1.controlPts(:,1)-0) < eps);
fixedYNodes     = find(abs(mesh1.controlPts(:,2)-0) < eps);
fixedZNodes     = find(abs(mesh1.controlPts(:,3)-0) < eps);
forcedNodes     = find(abs(mesh1.controlPts(:,1)-L) < eps);
topSurfNodes    = find(abs(mesh1.controlPts(:,3)-t) < eps);

cc = intersect(forcedNodes,find(abs(mesh1.controlPts(:,2)-W) < eps));

%% element sets
%% from bottom to top: composite, metal

matMap1 = ones(mesh1.noElems,1);


elemset1 = [];
elemset2 = [];
elemset3 = [];

count    = 1;
for iw=1:mesh1.noElemsW
    for iv=1:mesh1.noElemsV
        for iu=1:mesh1.noElemsU
            if     iw == 1
                elemset1 = [elemset1;count];
            elseif iw == 2                 
                elemset2 = [elemset2;count];
            elseif iw == 3                 
                elemset3 = [elemset3;count];
            end
            count = count + 1;
        end
    end
end

matMap1(elemset1,1) = 1;
matMap1(elemset2,1) = 2;
matMap1(elemset3,1) = 3;


meshData.matMap{1}=matMap1;



%% write to jive mesh

fileName = '~/code/jive/bezier/delamination/circle-crack/circle.mesh';

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

fprintf(file, '<ElementGroup name="interface">\n{');

for i=1+meshData.noElems : meshData.noElems+size(imeshes{1}.globElems,1)
    fprintf(file, '  %1d', i-1);
    
end
fprintf(file, '}\n');
fprintf(file, '</ElementGroup>\n');

% differentiate contact elements and cohesive elements

id = find(unique(vKnot)==0.5) - 1;

surf = zeros(noElemsV,noElemsU);
count = 1;
for j=1:noElemsV
    for i=1:noElemsU
        surf(j,i) = count;
        count = count + 1;
    end
end

contacts  = surf(1:id,:);
cohesives = surf(id+1:end,:);

contacts  = contacts - 1;
cohesives = cohesives - 1;

contacts  = contacts(:);
cohesives = cohesives(:);

fprintf(file, '<ElementGroup name="contacts">\n{');

for i=1:length(contacts)
   fprintf(file, '  %1d', contacts(i) + meshData.noElems);
  
end
fprintf(file, '}\n');
fprintf(file, '</ElementGroup>\n');

fprintf(file, '<ElementGroup name="cohesives">\n{');

for i=1:length(cohesives)
   fprintf(file, '  %1d', cohesives(i) + meshData.noElems);
  
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

fprintf(file, '<NodeGroup name="top">\n{');

for i=1:length(topSurfNodes)
    fprintf(file, '  %1d', topSurfNodes(i));
    
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="aa">\n{');

for i=1:length(cc)/2
    fprintf(file, '  %1d', topSurfNodes(2*i-1));
    
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="bb">\n{');

for i=1:length(cc)/2
    fprintf(file, '  %1d', topSurfNodes(2*i));
    
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fclose(file);


disp('writing jive mesh done!!!')

pause

%%
% read jem-jive result and do post-processing


vtuFileName = '../results/circle-crack';
resultFile = '~/code/jive/bezier/delamination/circle-crack/circle.dat';
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
    for ip=1:1
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


