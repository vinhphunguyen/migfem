addpath ../C_files/
addpath ../fem_util/
addpath ../examples/
addpath ../meshing/
addpath ../nurbs-util/
addpath ../delamination/
addpath ../structural-mechanics/
addpath ../xml_toolbox/
addpath ../post-processing/

clear all
clc

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
ndim =3;
% backtracking parameters
alpha = 0.1;
beta  = 0.7;

% initial surface

uKnot = [0 0 0 1 1 1];
vKnot = [0 0 0 1 1 1];

controlPts  = zeros(4,3,3);

controlPts(1:3,1,1) = [0;0;0];
controlPts(1:3,2,1) = [1;0.8;0];
controlPts(1:3,3,1) = [2;0;0];

controlPts(1:3,1,2) = [0;0.6;0.5];
controlPts(1:3,2,2) = [1;0.8;0.5];
controlPts(1:3,3,2) = [2;0.6;0.5];

controlPts(1:3,1,3) = [0;0;1];
controlPts(1:3,2,3) = [1;0.8;1];
controlPts(1:3,3,3) = [2;0;1];

controlPts(4,:)     = 1;

oriSurface = nrbmak(controlPts,{uKnot vKnot});

figure
hold on
nrbctrlplot(oriSurface);
axis equal
axis off

[p,q,uKnot,vKnot,noPtsX,noPtsY,weights,cPts] = convertNURBSSurface (oriSurface);

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

% generating points on offset curve
thickness = -0.05;

% computing offset points

noPts=20;

xi  = linspace(0,1,noPts);
eta = linspace(0,1,noPts);

offsetPts = zeros(noPts*noPts,3);
io = 1;
for i=1:length(eta)
    ett = eta(i);
    uspan   = findspan(noPtsY-1,q,ett,vKnot);
    for j=1:length(xi)
        xii = xi(j);
        vspan  = findspan(noPtsX-1,p,xii,uKnot);
        [R dRdxi dRdeta] = NURBS2DBasisDers([xii; ett],p,q,uKnot,vKnot,weights');
        
        pts    = cPts(1:end,:);
        x      = R    * pts;
        jacob  = [dRdxi; dRdeta] * pts; % 2x2 matrix
        
        % a1, a2 and a3 vectors (surface basis vectors)
        a1     = jacob(1,:);
        a2     = jacob(2,:);
        a3     = cross(a1,a2);
        a3     = a3/norm(a3);
        offsetPts(io,:) = x + a3*thickness;
        io = io + 1;
    end
end

figure
hold on
nrbplot(oriSurface,[12 12]);
axis equal
axis off
plot3(offsetPts(:,1),offsetPts(:,2),offsetPts(:,3),'rx');

%% gradient descent method
% initial guess curve = a line


% find indices of control points that are unknowns

id  = 1:size(cPts,1);
idf = [1 3 7 9];
idu = setdiff(id,idf);

cpoints        = cPts(:,1:3);

cpoints(1,     :) = offsetPts(1,:);
cpoints(3,     :) = offsetPts(noPts,:);
cpoints(7,     :) = offsetPts(noPts*(noPts-1)+1,:);
cpoints(end,   :) = offsetPts(end,:);

vec = cpoints(end,:) - cpoints(1,:);

% dlam = 1/(noCtrPts-1);
% for i=1:noCtrPts-2
%     cpoints(1+i,:) = cpoints(1,:)+ dlam*i*vec;
% end

cpoints0 = cpoints;


% disturbe the guess curve to compute numerically the gradient of the
% energy

samPts = computePointsOnSurface(uKnot,vKnot,p,q,noPtsX-1,noPtsY-1,cpoints,xi,eta);
energy   = getEnergy(offsetPts,samPts);

h      = 1e-8;
eps    = 0.01;
eps1   = 0.01;
error  = 10;
error1 = 10;
dgamma = 0.0001;
gamma  = 0;
iiter  = 1;
maxIter   = 40;
while error > eps
    samPts = computePointsOnSurface(uKnot,vKnot,p,q,noPtsX-1,noPtsY-1,cpoints,xi,eta);
    f        = getEnergy(offsetPts,samPts);
    
    % determine the gradients numerically
    
    for j = 1:length(idu)
        i = idu(j);
        cpoints(i,:) = cpoints(i,:) + [h 0 0];
        samPts = computePointsOnSurface(uKnot,vKnot,p,q,noPtsX-1,noPtsY-1,cpoints,xi,eta);
        energy1x = getEnergy(offsetPts,samPts);
        
        cpoints(i,:) = cpoints(i,:) - [h 0 0] + [0 h 0];
        samPts = computePointsOnSurface(uKnot,vKnot,p,q,noPtsX-1,noPtsY-1,cpoints,xi,eta);
        energy1y = getEnergy(offsetPts,samPts);
        
        cpoints(i,:) = cpoints(i,:) - [0 h 0] + [0 0 h];
        samPts = computePointsOnSurface(uKnot,vKnot,p,q,noPtsX-1,noPtsY-1,cpoints,xi,eta);
        energy1z = getEnergy(offsetPts,samPts);
        
        dir = 1/h*[energy1x-energy energy1y-energy energy1z-energy];
        grad(j,:) = dir;
        cpoints(i,:) = cpoints0(i,:);
    end
    
    s = 1;
    for k=1:10
        cpnew = cpoints0;
        cpnew(idu,:) = cpnew(idu,:) - grad*s;
        samPts = computePointsOnSurface(uKnot,vKnot,p,q,noPtsX-1,noPtsY-1,cpoints,xi,eta);
        fnew     = getEnergy(offsetPts,samPts);
        if (fnew < f + s*alpha*(-grad)'*grad)
            break;
        else
            s = s*beta;
        end
    end
    
    cpoints0(idu,:) = cpoints0(idu,:) - grad*s;
    cpoints = cpoints0;
    
    samPts   = computePointsOnSurface(uKnot,vKnot,p,q,noPtsX-1,noPtsY-1,cpoints,xi,eta);
    energy   = getEnergy(offsetPts,samPts);
    error    = energy;
    
    disp (sprintf(' %s %i %s %5.4e ', 'Iter',iiter, ':', error) );
    iiter = iiter + 1;
    
    if iiter > maxIter
        break
    end
end

%% plot result

figure
hold on
%plot3(offsetPts(:,1),offsetPts(:,2),offsetPts(:,3),'rx');
%plot3(samPts(:,1),samPts(:,2),samPts(:,3),'co');

nrbplot(oriSurface,[12 12]);
axis equal
axis off

cp        = zeros(4,noPtsX,noPtsY);
for j=1:noPtsY
  cp(1:3,:,j) = cpoints((j-1)*noPtsX+1:(j-1)*noPtsX+noPtsX,:)';
end
cp(4,:)   = weights;
offSurface = nrbmak(cp,{uKnot vKnot});
%nrbctrlplot(offSurface);
nrbplot(offSurface,[20 20]);

%% make a volume from two surfaces

volumePts = zeros(4,noPtsX,noPtsY,2);

volumePts(1:4,:,:,1) = oriSurface.coefs(:,:,:);
%volumePts(1:4,:,2,1) = oriSurface.coefs(:,:,2);

volumePts(1:4,:,:,2) = offSurface.coefs(:,:,:);
%volumePts(1:4,:,2,2) = offSurface.coefs(:,:,2);

%uKnot = uKnot;
wKnot = [0 0 1 1];

vol   = nrbmak(volumePts,{uKnot vKnot wKnot});

nrbplot(vol,[16 16 16])
axis off

%% refinement and make a mesh

%order elevation
vol = nrbdegelev(vol,[0 0 1]); % quadratic in thickness direction

uKnot     = cell2mat(vol.knots(1));
vKnot     = cell2mat(vol.knots(2));

refineCountX = 5;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
    % new knots along two directions (uniform)
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnots  = {newKnotsX newKnotsY []};
    
    % h-refinement
    vol       = nrbkntins(vol,newKnots);
    uKnot     = cell2mat(vol.knots(1));
    vKnot     = cell2mat(vol.knots(2));
end

%knot insertion to have C^-1 at coord = 0.2 (t=0.5)
vol     = nrbkntins(vol,{[] [] [1/3 1/3 2/3 2/3]});

mesh  = buildIGA3DMesh(vol);
vMesh = buildVisualizationMesh3D(vol);

noElemsU = mesh.noElemsU;
noElemsV = mesh.noElemsV;
p        = mesh.p;
q        = mesh.q;
uKnot    = mesh.uKnot;
vKnot    = mesh.vKnot;
noPtsX   = mesh.noPtsX;
noPtsY   = mesh.noPtsY;
controlPts=mesh.controlPts;


% iElements   = zeros(noElemsU*noElemsV,2*(p+1)*(q+1));
% iElementS   = buildIGA2DMesh (oriSurface);
% 
% e = 1;
% 
% y0 = 0.2;
% delaminationNodes  =  find(abs(controlPts(:,3) - y0 ) <1e-10);
% mm                 = 0.5*length(delaminationNodes);
% lowerNodes         = delaminationNodes(1:mm);
% upperNodes         = delaminationNodes(mm+1:end);
% 
% for i=1:noElemsU*noElemsV
%     sctr = iElementS.globElems(i,:);
%     iElements(e,1:(p+1)*(q+1))     = lowerNodes(sctr);
%     iElements(e,(p+1)*(q+1)+1:end) = upperNodes(sctr);
%     e = e + 1;
% end


%iMesh.locElems  = iElements;
%iMesh.globElems = iElements;


meshes{1}  = mesh;
vmeshes{1} = vMesh;


meshData.mesh    = meshes;
meshData.vmesh   = vmeshes;
meshData.noElems = meshes{1}.noElems;
meshData.noPts   = meshes{1}.noPts;

%% node set

eps = 1e-10;
fixedXNodes     = find(abs(mesh.controlPts(:,1)-0) < eps);
fixedYNodes     = find(abs(mesh.controlPts(:,1)-2) < eps);
fixedZNodes     = [fixedXNodes; fixedYNodes];

forcedNodes1     = find(abs(mesh.controlPts(:,1)-1) < eps);
forcedNodes2     = find(abs(mesh.controlPts(:,3)-0.5) < eps);

forcedNodes  = intersect(forcedNodes1,forcedNodes2);

%% element sets
%% from bottom to top: composite, metal

matMap1 = ones(mesh.noElems,1);


elemset1 = [];
elemset2 = [];
elemset3 = [];

count    = 1;
for iw=1:mesh.noElemsW
    for iv=1:mesh.noElemsV
        for iu=1:mesh.noElemsU
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

fileName = '~/code/jive/bezier/delamination/doubly-curved/doubly.mesh';

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

% ie = 1;
% for im=1:length(imeshes)
%     element = imeshes{im}.globElems;
%     for i=1:size(element,1)
%         fprintf(file, '  %1d %1d', meshData.noElems+ie-1, element(i,:) );
%         fprintf(file, ';\n');
%         ie = ie + 1;
%     end
% end

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

% e = 1;
% for im=1:length(meshes)
%     Cxe = meshes{im}.Cxe;
%     for il=1:no-1
%         for ie=1:meshes{im}.noElemsU*meshes{im}.noElemsV
%             Ce = Cxe(:,:,ie);
%             [row,col] = find(Ce);
%             fprintf(file, '  %1d ', meshData.noElems+e-1);
%             for i=1:length(row)
%                 fprintf(file, '%1d ', row(i)-1);
%             end
%             fprintf(file, ';\n');
%             e=e+1;
%         end
%     end
% end

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

% e = 1;
% for im=1:length(meshes)
%     Cxe = meshes{im}.Cxe;
%     for il=1:no-1
%         for ie=1:meshes{im}.noElemsU*meshes{im}.noElemsV
%             Ce = Cxe(:,:,ie);
%             [row,col] = find(Ce);
%             fprintf(file, '  %1d ', meshData.noElems+e-1);
%             for i=1:length(row)
%                 fprintf(file, '%1d ', col(i)-1);
%             end
%             fprintf(file, ';\n');
%             e=e+1;
%         end
%     end
% end

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

% e = 1;
% for im=1:length(meshes)
%     Cxe = meshes{im}.Cxe;
%     for il=1:no-1
%         for ie=1:meshes{im}.noElemsU*meshes{im}.noElemsV
%             Ce = Cxe(:,:,ie);
%             [row,col,val] = find(Ce);
%             fprintf(file, '  %1d ', meshData.noElems+e-1);
%             for i=1:length(row)
%                 fprintf(file, '%2.2f ', val(i));
%             end
%             fprintf(file, ';\n');
%             e=e+1;
%         end
%     end
% end

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

% ie = 1;
% for im=1:length(imeshes)
%     weights = meshes{im}.weights;
%     element = imeshes{im}.locElems;
%     for e=1:size(element,1)
%         w = weights(element(e,:));
%         fprintf(file, '  %1d ',meshData.noElems+ie-1); ie = ie + 1;
%         for j=1:length(w)
%             fprintf(file, '%2.4f ', w(j));
%         end
%         fprintf(file, ';\n');
%     end
% end

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

% fprintf(file, '<ElementGroup name="interface">\n{');
% 
% for i=1+meshData.noElems : meshData.noElems+size(imeshes{1}.globElems,1)
%     fprintf(file, '  %1d', i-1);
%     
% end
% fprintf(file, '}\n');
% fprintf(file, '</ElementGroup>\n');
% 
% % differentiate contact elements and cohesive elements
% 
% id = find(unique(vKnot)==0.5) - 1;
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
% contacts  = surf(1:id,:);
% cohesives = surf(id+1:end,:);
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
%    fprintf(file, '  %1d', contacts(i) + meshData.noElems);
%   
% end
% fprintf(file, '}\n');
% fprintf(file, '</ElementGroup>\n');
% 
% fprintf(file, '<ElementGroup name="cohesives">\n{');
% 
% for i=1:length(cohesives)
%    fprintf(file, '  %1d', cohesives(i) + meshData.noElems);
%   
% end
% fprintf(file, '}\n');
% fprintf(file, '</ElementGroup>\n');

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



fclose(file);


disp('writing jive mesh done!!!')

pause

%%
% read jem-jive result and do post-processing


vtuFileName = '../results/doubly';
resultFile = '~/code/jive/bezier/delamination/doubly-curved/doubly.dat';
V = xml_parseany(fileread(resultFile));


noStep = size(V,2)/2;
% if noStep=1 then V{it}.Section does not work
for it=1:noStep        
    uData=V{2*it-1}.Section{1};   % displacement
    dData=V{2*it}.Section{1}; % damage
    disp = str2num(uData.CONTENT);
    dam  = str2num(dData.CONTENT);
    damage = zeros(meshData.noPts*ndim,1);
    %damage(dam(:,1),1) = dam(:,2);
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








