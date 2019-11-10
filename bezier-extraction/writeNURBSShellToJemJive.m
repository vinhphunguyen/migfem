
%
% Write the Bezier extraction operators for a NURBS mesh defined
% by the given knots to jem-jive file.
%
% VP Nguyen
% Cardiff University, UK
% Feburary, 2013.

% run this after running a data file that defines the CAD.

generateIGA2DMesh

C        = bezierExtraction2D(uKnot,vKnot,p,q);

noElems  = size(C,3);

fileName = '~/code/jive/bezier/freeEndsCylinderShell.mesh';

file = fopen(fileName, 'wt');

fprintf(file, '<Nodes>\n');

for i=1:length(controlPts)
   fprintf(file, '  %1d %2.8f %2.8f %2.8f', i, ...
              controlPts(i,1),controlPts(i,2),controlPts(i,3));
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
        fprintf(file, '%2.4E ', val(i));
    end
    fprintf(file, ';\n');
end
fprintf(file, ' </Column>\n');

fprintf(file, ' <Column name = "weights" type = "float">\n');

for i=1:size(element,1)
    w = weights(element(i,:));
    fprintf(file, '  %1d ',i-1);
    for i=1:length(w)
        fprintf(file, '%2.4E ', w(i));
    end
    fprintf(file, ';\n');
end

fprintf(file, ' </Column>\n');

fprintf(file, '</ElementDatabase>\n\n');



fprintf(file, '<NodeGroup name="xfix">\n{');

for i=1:length(xConsNodes)
   fprintf(file, '  %1d', xConsNodes(i));
  
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="yfix">\n{');

for i=1:length(yConsNodes)
   fprintf(file, '  %1d', yConsNodes(i));
  
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="zfix">\n{');

for i=1:length(zConsNodes)
   fprintf(file, '  %1d', zConsNodes(i));
  
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');


fprintf(file, '<NodeGroup name="AD">\n{');

for i=1:length(nodesOnAD)
   fprintf(file, '  %1d', nodesOnAD(i));
  
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="BC">\n{');

for i=1:length(nodesOnCB)
   fprintf(file, '  %1d', nodesOnCB(i));
  
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="CD">\n{');

for i=1:length(nodesOnCD)
   fprintf(file, '  %1d', nodesOnCD(i));
  
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="next2AD">\n{');

for i=1:length(nextToADNodes)
   fprintf(file, '  %1d', nextToADNodes(i));
  
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="next2BC">\n{');

for i=1:length(nextToBCNodes)
   fprintf(file, '  %1d', nextToBCNodes(i));
  
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="next2CD">\n{');

for i=1:length(nextToCDNodes)
   fprintf(file, '  %1d', nextToCDNodes(i));
  
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');
