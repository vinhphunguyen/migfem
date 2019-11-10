%function ok = writeMultiPatch3DToJive(data,dirichlet,newmann,fileName)
%
% write multi-patch NURBS solid to jive *.mesh
%
% VP Nguyen, May 2013
% Cardiff University, Wales, UK


file = fopen(fileName, 'wt');
patchCount = length(data.mesh);

%% write nodes

fprintf(file, '<Nodes>\n');

i = 1;
for ip=1:patchCount
    mesh      = data.mesh{ip};
    controlPts= mesh.controlPts;
    
    for in=1:length(controlPts)
        fprintf(file, '  %1d %2.6f %2.6f %2.6f', i, controlPts(in,1),controlPts(in,2),controlPts(in,3));
        fprintf(file, ';\n');
        i = i + 1;
    end
end

fprintf(file, '</Nodes>\n\n');

fprintf(file, '<Elements>\n');

%% write solid elements

i = 1;
for ip=1:patchCount
    mesh      = data.mesh{ip};
    globElems = mesh.globElems;
    
    for ie=1:size(globElems,1)
        fprintf(file, '  %1d %1d', i-1, globElems(ie,:) );
        fprintf(file, ';\n');
        i = i + 1;
    end
end

fprintf(file, '</Elements>\n\n');

%% write Bezier extractors

fprintf(file, '<ElementDatabase name="C">\n');

fprintf(file, ' <Column name = "irows" type = "int">\n');

ie = 1;
for ip=1:patchCount
    mesh      = data.mesh{ip};
    C         = mesh.C;
    noElems   = size(C,3);
    for e=1:noElems
        Ce = C(:,:,e);
        [row,col] = find(Ce);
        fprintf(file, '  %1d ', ie-1);
        for i=1:length(row)
            fprintf(file, '%1d ', row(i)-1);
        end
        fprintf(file, ';\n');
        ie = ie + 1;
    end
end

fprintf(file, ' </Column>\n');

fprintf(file, ' <Column name = "jcols" type = "int">\n');

ie = 1;
for ip=1:patchCount
    mesh      = data.mesh{ip};
    C         = mesh.C;
    noElems   = size(C,3);
    for e=1:noElems
        Ce = C(:,:,e);
        [row,col] = find(Ce);
        
        fprintf(file, '  %d ', ie-1);
        for i=1:length(row)
            fprintf(file, '%1d ', col(i)-1);
        end
        fprintf(file, ';\n');
        ie = ie + 1;
    end
    
end

fprintf(file, ' </Column>\n');

fprintf(file, ' <Column name = "values" type = "float">\n');

ie = 1;
for ip=1:patchCount
    mesh      = data.mesh{ip};
    C         = mesh.C;
    noElems   = size(C,3);
    
    for e=1:noElems
        Ce = C(:,:,e);
        [row,col,val] = find(Ce);
        
        fprintf(file, '  %d ', ie-1);
        for i=1:length(row)
            fprintf(file, '%2.2f ', val(i));
        end
        fprintf(file, ';\n');
        ie = ie + 1;
    end
end

fprintf(file, ' </Column>\n');

% write weights

fprintf(file, ' <Column name = "weights" type = "float">\n');

ie = 1;
for ip=1:patchCount
    mesh      = data.mesh{ip};
    C         = mesh.C;
    weights   = mesh.weights;
    noElems   = size(C,3);
    elems     = mesh.locElems;
    for e=1:noElems
        w = weights(elems(e,:));
        fprintf(file, '  %1d ',ie-1);
        for j=1:length(w)
            fprintf(file, '%2.4f ', w(j));
        end
        fprintf(file, ';\n');
        ie = ie + 1;
    end
end

fprintf(file, ' </Column>\n');
fprintf(file, '</ElementDatabase>\n\n');

%% write element groups


fprintf(file, '<ElementGroup name="material0">\n{');
for i=1:length(elementSet1)
    fprintf(file, '  %1d', elementSet1(i)-1);
    
end
for i=1:length(elementSet3)
    fprintf(file, '  %1d', elementSet3(i)-1);
    
end
fprintf(file, '}\n');
fprintf(file, '</ElementGroup>\n');


fprintf(file, '<ElementGroup name="material90">\n{');
for i=1:length(elementSet2)
    fprintf(file, '  %1d', elementSet2(i)-1);
    
end
for i=1:length(elementSet4)
    fprintf(file, '  %1d', elementSet4(i)-1);
    
end
fprintf(file, '}\n');
fprintf(file, '</ElementGroup>\n');




% write node groups

xnodes = dirichlet.xnodes;
ynodes = dirichlet.ynodes;
znodes = dirichlet.znodes;

fnodes = newmann.xnodes;

fprintf(file, '<NodeGroup name="gr1">\n{');

for i=1:length(xnodes)
    fprintf(file, '  %1d', xnodes(i));
    
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fprintf(file, '<NodeGroup name="gr2">\n{');

for i=1:length(fnodes)
    fprintf(file, '  %1d', fnodes(i));
    
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

ok = 1;

