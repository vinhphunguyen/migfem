%% Geometry data

% plate dimensions
a = 10.0; % length
b = 0.1;  % thickness
t = 1;    % width

% knots
uKnot = [0 0 1 1];
vKnot = [0 0 1 1];
wKnot = [0 0 1 1];

% control points
controlPts          = zeros(4,2,2,2);

controlPts(1:3,1,1,1) = [0;0;0];
controlPts(1:3,2,1,1) = [a;0;0];

controlPts(1:3,1,2,1) = [0;b;0];
controlPts(1:3,2,2,1) = [a;b;0];

controlPts(1:3,1,1,2) = [0;0;t];
controlPts(1:3,2,1,2) = [a;0;t];

controlPts(1:3,1,2,2) = [0;b;t];
controlPts(1:3,2,2,2) = [a;b;t];

% weights
controlPts(4,:,:)   = 1;

%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot wKnot});

%% p-refinment

solid = nrbdegelev(solid,[1 0 0]); % to cubic-linear NURBS

%% h-refinement

refineLevel = 4;
for i=1:refineLevel
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
    % new knots along two directions
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    
    newKnots  = {newKnotsX {} {}};
    solid     = nrbkntins(solid,newKnots);
    uKnot      = cell2mat(solid.knots(1));
    vKnot      = cell2mat(solid.knots(2));
end

nrbkntplot(solid)
%% 

convert3DNurbs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write jem jive input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Bezier extraction operators
generateIGA3DMesh

[C,Cxi,Cet]  = bezierExtraction3D(uKnot,vKnot,wKnot,p,q,r);
    
fixedNodes    =  find(controlPts(:,1)==0);
EPS=1e-8;
forcedNodes1  =  find(abs(controlPts(:,1)-a)<EPS);
forcedNodes2  =  find(abs(controlPts(:,2)-b)<EPS);
forcedNodes3  =  find(abs(controlPts(:,2))  <EPS);
forcedNodes   =  intersect(forcedNodes1,forcedNodes2);

%% write to jive mesh 


noElems  = size(C,3);

fileName = '~/code/jive/bezier/large-displacement/beam3D.mesh';

file = fopen(fileName, 'wt');

fprintf(file, '<Nodes>\n');

for i=1:length(controlPts)
   fprintf(file, '  %1d %2.6f %2.6f %2.6f', i, controlPts(i,1),controlPts(i,2),controlPts(i,3));
   fprintf(file, ';\n');
end

fprintf(file, '</Nodes>\n\n');

fprintf(file, '<Elements>\n');

% write solid elements

for i=1:noElems
   fprintf(file, '  %1d %1d', i-1, element(i,:) );
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


fprintf(file, ' </Column>\n');

fprintf(file, '</ElementDatabase>\n\n');

% write element groups

% write node groups

fprintf(file, '<NodeGroup name="gr1">\n{');
for i=1:length(fixedNodes)
   fprintf(file, '  %1d', fixedNodes(i));
  
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="gr2">\n{');
for i=1:length(forcedNodes)
   fprintf(file, '  %1d', forcedNodes(i));
  
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="gr3">\n{');
fprintf(file, '  %1d', forcedNodes(1)); 
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fclose(file);

disp('writing mesh file for jive is done!!!')
