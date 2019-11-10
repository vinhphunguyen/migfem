function writeBezierExtractionOperator(uknot,vknot,p,q,weights,fileName)
%
% Write the Bezier extraction operators for a NURBS mesh defined
% by the given knots to jem-jive file.
%
% VP Nguyen
% Cardiff University, UK
% Feburary, 2013.

C        = bezierExtraction2D(uknot,vknot,p,q);

noElems  = size(C,3);

file = fopen(fileName, 'wt');

fprintf(file, '<ElementDatabase name="C">\n');

fprintf(file, ' <Column name = "i" type = "int">\n');

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

fprintf(file, ' <Column name = "j" type = "int">\n');

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

fprintf(file, ' <Column name = "v" type = "float">\n');
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

for i=1:length(weights)
    fprintf(file, '%2.1f ', weights(i));
end
fprintf(file, ';\n');

fprintf(file, ' </Column>\n');

fprintf(file, '</ElementDatabase>\n');

fclose(file);



