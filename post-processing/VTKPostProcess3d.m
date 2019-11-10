% Copyright (C) 2007 Garth N. Wells
%
% Write VTK post-processing files
%
% Modified by VP Nguyen for the IGA FEM code

function VTKPostProcess3d(node,elementV,etype,vtuFile,sigma, disp)

dim      = size(node,2);
numNodes = size(node,1);
numCells = size(elementV,1);
dof      = 3;
x        = node;
connect  = elementV;

if size(sigma,2) == 6
    sigma = [sigma zeros(size(sigma,1),1)];
end

% Output files

outfileVTU  = strcat(vtuFile, '.vtu');
results_vtu = fopen(outfileVTU, 'wt');

if(strcmp(etype, 'B8'))
    numVertexesPerCell = 8;
    VTKCellCode = 12;
elseif(strcmp(etype, 'Tet4'))
    numVertexesPerCell = 4;
    VTKCellCode = 10;
else
    error('Element type not known (VTKPostProcess)')
end


dof_per_vertex = 3;


%% Write headers
fprintf(results_vtu, '<VTKFile type="UnstructuredGrid"  version="0.1"   > \n');
fprintf(results_vtu, '<UnstructuredGrid> \n');
fprintf(results_vtu, '<Piece  NumberOfPoints="  %g" NumberOfCells=" %g"> \n', numNodes, numCells);

%% Write point data
fprintf(results_vtu, '<Points> \n');

fprintf(results_vtu, '<DataArray  type="Float64"  NumberOfComponents="3"  format="ascii" > \n');


for i=1:numNodes
    fprintf(results_vtu, '%f ',  x(i,1:3));
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
sigmaZZ = sigma(:,3);
sigmaXY = sigma(:,4);
sigmaYZ = sigma(:,5);
sigmaZX = sigma(:,6);
sigmaVM = sigma(:,7);

dispX   = disp(:,1);
dispY   = disp(:,2);
dispZ   = disp(:,3);

fprintf(results_vtu, '<PointData  Vectors="sigma"> \n');
fprintf(results_vtu, '<DataArray  type="Float64"  Name="sigma" NumberOfComponents="7" format="ascii"> \n');


for i=1:numNodes
    fprintf(results_vtu, '%f   ', sigmaXX(i) );
    fprintf(results_vtu, '%f   ', sigmaYY(i) );
    fprintf(results_vtu, '%f   ', sigmaZZ(i) );
    fprintf(results_vtu, '%f   ', sigmaXY(i) );
    fprintf(results_vtu, '%f   ', sigmaYZ(i) );
    fprintf(results_vtu, '%f   ', sigmaZX(i) );
    fprintf(results_vtu, '%f   ', sigmaVM(i) );
    fprintf(results_vtu, '\n');
end

fprintf(results_vtu, '</DataArray> \n');

% print displacement field

fprintf(results_vtu, '<DataArray  type="Float64"  Name="U" NumberOfComponents="3" format="ascii"> \n');


for i=1:numNodes
    fprintf(results_vtu, '%f   ', dispX(i) );
    fprintf(results_vtu, '%f   ', dispY(i) );
    fprintf(results_vtu, '%f   ', dispZ(i) );
    fprintf(results_vtu, '\n');
end

fprintf(results_vtu, '</DataArray> \n');

if size(disp,2) == 4
    damage   = disp(:,4);
    
    fprintf(results_vtu, '<DataArray  type="Float64"  Name="dam" NumberOfComponents="3" format="ascii"> \n');
        
    for i=1:numNodes
        fprintf(results_vtu, '%f   ', damage(i) );                
        fprintf(results_vtu, '%f   ', 0. ); 
        fprintf(results_vtu, '%f   ', 0. ); 
        fprintf(results_vtu, '\n');
    end
    
    fprintf(results_vtu, '</DataArray> \n');    
end

fprintf(results_vtu, '</PointData> \n');


% end of VTK file

fprintf(results_vtu, '</Piece> \n');
fprintf(results_vtu, '</UnstructuredGrid> \n');
fprintf(results_vtu, '</VTKFile> \n');

fclose(results_vtu);
