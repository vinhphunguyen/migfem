% Copyright (C) 2007 Garth N. Wells
%
% Write VTK post-processing files
%
% Modified by VP Nguyen for the IGA FEM code

function VTKPostProcess(node,elementV,ndofs,etype,vtuFile,sigma,disp)


dim      = size(node,2);
numNodes = size(node,1);
numCells = size(elementV,1);
dof      = ndofs;
x        = node;
connect  = elementV;

% Output files

outfileVTU  = strcat(vtuFile, '.vtu');
results_vtu = fopen(outfileVTU, 'wt');

if(strcmp(etype, 'Quad4') || strcmp(etype, 'Quad8') || strcmp(etype, 'Quad9'))
    numVertexesPerCell = 4;
    VTKCellCode = 9;
elseif(strcmp(etype, 'B8'))
    numVertexesPerCell = 8;
    VTKCellCode = 12;
elseif(strcmp(etype, 'Tri3') || strcmp(etype, 'Tri6'))
    numVertexesPerCell = 3;
    VTKCellCode = 5;
elseif(strcmp(etype, 'Tet4'))
    numVertexesPerCell = 4;
    VTKCellCode = 10;
else
    error('Element type not known (VTKPostProcess)')
end


dof_per_vertex = 2;


%% Write headers
fprintf(results_vtu, '<VTKFile type="UnstructuredGrid"  version="0.1"   > \n');
fprintf(results_vtu, '<UnstructuredGrid> \n');
fprintf(results_vtu, '<Piece  NumberOfPoints="  %g" NumberOfCells=" %g"> \n', numNodes, numCells);

%% Write point data
fprintf(results_vtu, '<Points> \n');
if( dof_per_vertex == 1)
    fprintf(results_vtu, '<DataArray  type="Float64"  NumberOfComponents="1"  format="ascii" > \n');
else
    fprintf(results_vtu, '<DataArray  type="Float64"  NumberOfComponents="3"  format="ascii" > \n');
end

for i=1:numNodes
    if( dim == 3)
        fprintf(results_vtu, '%f ',  x(i,1:3));
    elseif(dim == 2)
        fprintf(results_vtu, '%f ',  x(i,1:2));
        fprintf(results_vtu, '0.0 ');
    end
    fprintf(results_vtu, '\n');
end

fprintf(results_vtu, '</DataArray> \n');
fprintf(results_vtu, '</Points> \n');

%% Print cells
fprintf(results_vtu, '<Cells> \n');

%% Print cell connectivity
fprintf(results_vtu, '<DataArray  type="Int32"  Name="connectivity"  format="ascii"> \n');

for i=1:numCells
    fprintf(results_vtu, '%g ',  connect(i,1:numVertexesPerCell)-1 );
    fprintf(results_vtu, '\n');
end

fprintf(results_vtu, '</DataArray> \n');

% Print cell offsets
fprintf(results_vtu, '<DataArray  type="Int32"  Name="offsets"  format="ascii"> \n');

offset = 0;
for i=1:numCells
    offset = offset + numVertexesPerCell;
    fprintf(results_vtu, '%g ', offset);
    fprintf(results_vtu, '\n');
end

fprintf(results_vtu, '</DataArray> \n');

% Print cell types
fprintf(results_vtu, '<DataArray  type="UInt8"  Name="types"  format="ascii"> \n');

for i=1:numCells
    fprintf(results_vtu, '%g ', VTKCellCode);
    fprintf(results_vtu, '\n');
end

fprintf(results_vtu, '</DataArray> \n');
fprintf(results_vtu, '</Cells> \n');

%% Print result data

sigmaXX = sigma(:,1);
sigmaYY = sigma(:,2);
sigmaXY = sigma(:,3);

if size(sigma,2)==4
    sigmaVM = sigma(:,4); % von Mises stress
end

if( dof_per_vertex == 1)
    fprintf(results_vtu, '<PointData  Scalars="U"> \n');
    fprintf(results_vtu, '<DataArray  type="Float64"  Name="U" NumberOfComponents="1" format="ascii"> \n');
else
    fprintf(results_vtu, '<PointData  Vectors="sigma"> \n');
    
    if size(sigma,2)==4
        fprintf(results_vtu, '<DataArray  type="Float64"  Name="sigma" NumberOfComponents="4" format="ascii"> \n');
    else
        fprintf(results_vtu, '<DataArray  type="Float64"  Name="sigma" NumberOfComponents="3" format="ascii"> \n');
    end
end

for i=1:numNodes
    fprintf(results_vtu, '%f   ', sigmaXX(i) );
    fprintf(results_vtu, '%f   ', sigmaYY(i) );
    fprintf(results_vtu, '%f   ', sigmaXY(i) );
    
    if size(sigma,2)==4
        fprintf(results_vtu, '%f   ', sigmaVM(i) );
    end
    
    fprintf(results_vtu, '\n');
end

fprintf(results_vtu, '</DataArray> \n');

% print displacement field


dispX   = disp(:,1);
dispY   = disp(:,2);



fprintf(results_vtu, '<DataArray  type="Float64"  Name="U" NumberOfComponents="3" format="ascii"> \n');

if ndofs==2
    for i=1:numNodes
        fprintf(results_vtu, '%12.8f   ', dispX(i) );
        fprintf(results_vtu, '%12.8f   ', dispY(i) );
        fprintf(results_vtu, '0.0' );
        fprintf(results_vtu, '\n');
    end
else
    dispZ   = disp(:,3);
    for i=1:numNodes
        fprintf(results_vtu, '%12.8f   ', dispX(i) );
        fprintf(results_vtu, '%12.8f   ', dispY(i) );
        fprintf(results_vtu, '%12.8f   ', dispZ(i) );
        fprintf(results_vtu, '\n');
    end
end
fprintf(results_vtu, '</DataArray> \n');
fprintf(results_vtu, '</PointData> \n');

% end of VTK file

fprintf(results_vtu, '</Piece> \n');
fprintf(results_vtu, '</UnstructuredGrid> \n');
fprintf(results_vtu, '</VTKFile> \n');

fclose(results_vtu);
