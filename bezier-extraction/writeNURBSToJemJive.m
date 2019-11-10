
%
% Write the Bezier extraction operators for a NURBS mesh defined
% by the given knots to jem-jive file.
%
% VP Nguyen
% Cardiff University, UK
% Feburary, 2013.


addpath ../fem_util/
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../meshing/
addpath ../nurbs-util/
addpath ../nurbs-geopdes/inst/

a = 2.0;

res = 200; % resolution for plotting NURBS

uKnot = [0 0 0 1 1 1];
vKnot = [0 0 0 1 1 1];

controlPts          = zeros(4,3,3);

controlPts(1:2,1,1) = [0;0];
controlPts(1:2,2,1) = [0;a/2];
controlPts(1:2,3,1) = [0;a];

controlPts(1:2,1,2) = [0.5*(a);0];
controlPts(1:2,2,2) = [0.5*(a); 0.5*(a)];
controlPts(1:2,3,2) = [0.5*(a);a];

controlPts(1:2,1,3) = [a;0];
controlPts(1:2,2,3) = [a;0.5*a];
controlPts(1:2,3,3) = [a;a];

controlPts(4,:,:)   = 1;

% fac                 = 1/sqrt(2);
% 
% controlPts(4,2,1) = fac;
% controlPts(4,2,2) = fac;
% controlPts(4,2,3) = fac;

% homogenous coordinates (x*w,y*w,z*w)
% 
% controlPts(1:2,2,1) = controlPts(1:2,2,1)*fac;
% controlPts(1:2,2,2) = controlPts(1:2,2,2)*fac;
% controlPts(1:2,2,3) = controlPts(1:2,2,3)*fac;

%% build NURBS object

solid = nrbmak(controlPts,{uKnot vKnot});

figure
hold on
nrbplot(solid,[40 40])

%% refinement

% refineLevel = 4;
% for i=1:refineLevel
%     uKnotVectorU = unique(uKnot);
%     uKnotVectorV = unique(vKnot);
%     
%     % new knots along two directions
%     
%     newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
%     newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
%     
%     newKnots  = {newKnotsX newKnotsY};
%     solid     = nrbkntins(solid,newKnots);
%     uKnot      = cell2mat(solid.knots(1));
%     vKnot      = cell2mat(solid.knots(2));
% end

%%

% newKnotsX = 0.5;
%     newKnotsY = 0.5;
%
%     newKnots  = {newKnotsX newKnotsY};
%     solid     = nrbkntins(solid,newKnots);

convert2DNurbs

plotMesh (controlPts,weights,uKnot,vKnot,p,q,res,'r-','try.eps');

%% Bezier extraction operators

[C,Cxi,Cet]  = bezierExtraction2D(uKnot,vKnot,p,q);

noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;


%%

% generateIGA2DMesh
% bezierPoints = [];
% 
% for e=1:noElems
%     sctr   = element(e,:);         %  element scatter vector
%     
%     Ce     = C(:,:,e);             % element Bezier extraction operator
%     we     = diag(weights(sctr,:));% element weights
%     pts    = controlPts(sctr,:);   % element nodes
%     Wb     = Ce'*weights(sctr,:);  % element Bezier weights
%     
%     bezierPts = inv(diag(Wb))*Ce'*we*pts;
%     bezierPoints = [bezierPoints;bezierPts];
% end
% 
% % remove sharing control points
% 
% bezierPoints = unique_no_sort_rows(bezierPoints);
% 
% bezierElem=[1 4 7 2 5 8 3 6 9;
%             7 16 19 8 17 20 9 18 21;
%             3 6 9 10 12 14 11 13 15;
%             9 18 21 14 22 24 15 23 25];
% 
% plot_mesh(bezierPoints,bezierElem,'Q9','b-',1);
% 
% plot(bezierPoints(:,1),bezierPoints(:,2),'rs','MarkerSize',13,'MarkerFaceColor','r');

%%


convert2DNurbs
generateIGA2DMesh

% find boundary nodes for boundary conditions

fixedXNodes1  =  find(controlPts(:,1)==0);
fixedYNodes   =  fixedXNodes1;

fixedXNodes2  =  find(controlPts(:,2)==0);

fixedXNodes   = [fixedXNodes1; fixedXNodes2];

C        = bezierExtraction2D(uKnot,vKnot,p,q);

noElems  = size(C,3);

fileName = 'curvedBeam.mesh';

file = fopen(fileName, 'wt');

fprintf(file, '<Nodes>\n');

for i=1:length(controlPts)
   fprintf(file, '  %1d %1d %2.4f', i, controlPts(i,1),controlPts(i,2));
   fprintf(file, ';\n');
end

fprintf(file, '</Nodes>\n\n');

fprintf(file, '<Elements>\n');

for i=1:size(element,1)
   fprintf(file, '  %1d %1d', i-1, element(i,:) );
   fprintf(file, ';\n');
end

fprintf(file, '</Elements>\n\n');

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

fprintf(file, ' <Column name = "weights" type = "float">\n');

for i=1:size(element,1)
    w = weights(element(i,:));
    fprintf(file, '  %1d ',i-1);
    for i=1:length(w)
        fprintf(file, '%2.4f ', w(i));
    end
    fprintf(file, ';\n');
end

fprintf(file, ' </Column>\n');

fprintf(file, '</ElementDatabase>\n\n');

fprintf(file, '<NodeGroup name="gr1">\n{');

for i=1:length(fixedXNodes1)
   fprintf(file, '  %1d', fixedXNodes1(i));
  
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');
fprintf(file, '<NodeGroup name="gr2">\n{');

for i=1:length(fixedXNodes2)
   fprintf(file, '  %1d', fixedXNodes2(i));
  
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fclose(file);



