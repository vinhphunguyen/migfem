% NURBS data for 3D beam
% NURBS data: control points, knots, basis orders and weights.
% Two ways to define the geometry: manually or using GeoPDEs NURBS toolbox.
% 
% Vinh Phu Nguyen, Johns Hopkins University

% control points

a = 10;
b = 4;
c = 2;

% noPtsX = 21;
% noPtsY = 7;
% noPtsZ = 5;
% 
% [controlPts,elementVV]=makeB8mesh(a,b,c,noPtsX,noPtsY,noPtsZ);
% 
% % basis order
% 
% p = 2;
% q = 2;
% r = 2;
% 
% % knot vectors
% 
% knotUTemp = linspace(0,1,noPtsX-p+1);
% knotVTemp = linspace(0,1,noPtsY-q+1);
% knotWTemp = linspace(0,1,noPtsZ-r+1);
% 
% uKnot = [0 0 knotUTemp 1 1];
% vKnot = [0 0 knotVTemp 1 1];
% wKnot = [0 0 knotWTemp 1 1];
% 
% % weights
% 
% weights = ones(1,noPtsX*noPtsY*noPtsZ)';

controlPts = zeros(4,2,2,2);

controlPts(1:3,1,1,1) = [0;0;0];
controlPts(1:3,2,1,1) = [a;0;0];

controlPts(1:3,1,2,1) = [0;b;0];
controlPts(1:3,2,2,1) = [a;b;0];

controlPts(1:3,1,1,2) = [0;0;c];
controlPts(1:3,2,1,2) = [a;0;c];

controlPts(1:3,1,2,2) = [0;b;c];
controlPts(1:3,2,2,2) = [a;b;c];

controlPts(4,:,:)     = 1; % weights

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];
zKnot = [0 0 1 1];

%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot zKnot});

figure 
hold on
nrbplot(solid,[10 10 10])

%% degree elevator
solid = nrbdegelev(solid,[1 1 1]); 

%% h-refinement

refineLevel = 0;
for i=1:refineLevel
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    
    % new knots along two directions
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    
    newKnots  = {newKnotsX newKnotsY};
    solid     = nrbkntins(solid,newKnots);
    uKnot     = cell2mat(solid.knots(1));
    vKnot     = cell2mat(solid.knots(2));
end
%% 

convert3DNurbs


%%
%% Write input file for jem-jive C++ code
%%

%% write to jive mesh 

% build connectivity ...

generateIGA3DMesh

[C,Cxi,Cet,Cze] = bezierExtraction3D(uKnot,vKnot,wKnot,p,q,r);

noElems  = size(C,3);

fileName = '~/code/jive/bezier/3D/beam.mesh';

file = fopen(fileName, 'wt');

fprintf(file, '<Nodes>\n');

for i=1:length(controlPts)
   fprintf(file, '  %1d %2.8f %2.8f %2.8f', i, ...
       controlPts(i,1),controlPts(i,2),controlPts(i,3));
   fprintf(file, ';\n');
end

fprintf(file, '</Nodes>\n\n');

fprintf(file, '<Elements>\n');

% write elements

for i=1:noElems
   fprintf(file, '  %1d %1d', i-1, element(i,:) );
   fprintf(file, ';\n');
end

fprintf(file, '</Elements>\n\n');

% write Bezier extractors 

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
        fprintf(file, '%2.6f ', val(i));
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

fprintf(file, ' </Column>\n');
fprintf(file, '</ElementDatabase>\n\n');

% write node groups

fixedNodes   = find(controlPts(:,1)==0);
forcedNodes  = find(controlPts(:,1)==a);

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

fclose(file);


          
