function writeVTKShell(vtuFile,U,node,elementV)

global element p q index elRangeU elRangeV uKnot vKnot nElmLK...
    weights noCtrPts noDofs Ke0 controlPts W GP vIdGrd vSprGK...
    jSprRw jSprCl jSprLK noPtsX noPtsY C

noElems = size(element,1);
stress = zeros(noElems,size(elementV,2),3);
disp   = zeros(noElems,size(elementV,2),3);

Ux    = U(1:3:noDofs);
Uy    = U(2:3:noDofs);
Uz    = U(3:3:noDofs);

for e=1:size(element,1)
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);         %  element scatter vector

    nn              = length(sctr);    
    sctrB           = zeros(1,3*nn);
    sctrB(1:3:3*nn) = 3*sctr-2;
    sctrB(2:3:3*nn) = 3*sctr-1;
    sctrB(3:3:3*nn) = 3*sctr;
    
    B      = zeros(3,3*nn);
    pts    = controlPts(sctr,:);
    
    uspan = FindSpan(noPtsX-1,p,xiE(1), uKnot);
    vspan = FindSpan(noPtsY-1,q,etaE(1),vKnot);
    
    elemDisp = [Ux(sctr) Uy(sctr) Uz(sctr)];
    
    % loop over Gauss points
    
    gp = 1;
    for iv=1:2
        if (iv==2)
            xiE = sort(xiE,'descend');
        end
        for iu=1:2
            Xi  = xiE(iu);
            Eta = etaE(iv);
            [N dRdxi dRdeta] = NURBS2DBasisDersSpecial([Xi; Eta],...
                p,q,uKnot,vKnot,weights',[uspan;vspan]);
                       
            disp(e,gp,:)    = N*elemDisp;
            
            gp = gp +1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% export to VTK format to plot in Mayavi or Paraview

sigmaXX = zeros(size(node,1),1);
sigmaYY = zeros(size(node,1),1);
sigmaXY = zeros(size(node,1),1);

dispX = zeros(size(node,1),1);
dispY = zeros(size(node,1),1);
dispZ = zeros(size(node,1),1);

for e=1:size(elementV,1)
    connect = elementV(e,:);
    for in=1:4
        nid = connect(in);
        sigmaXX(nid) = stress(e,in,1);
        sigmaYY(nid) = stress(e,in,2);
        sigmaXY(nid) = stress(e,in,3);
        
        dispX(nid) = disp(e,in,1);
        dispY(nid) = disp(e,in,2);
        dispZ(nid) = disp(e,in,3);
    end
end

dim      = size(node,2);
numNodes = size(node,1);
numCells = size(elementV,1);
dof      = 3;
x        = node;
connect  = elementV;

% Output files

results_vtu = fopen(strcat(vtuFile, '.vtu'), 'wt');

numVertexesPerCell = 4;
VTKCellCode        = 9;
dof_per_vertex     = 3;


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

fprintf(results_vtu, '<PointData  Vectors="sigma"> \n');

% print displacement field

fprintf(results_vtu, '<DataArray  type="Float64"  Name="U" NumberOfComponents="3" format="ascii"> \n');
for i=1:numNodes
    fprintf(results_vtu, '%f   ', dispX(i) );
    fprintf(results_vtu, '%f   ', dispY(i) );
    fprintf(results_vtu, '%f   ', dispZ(i) );    
    fprintf(results_vtu, '\n');
end
fprintf(results_vtu, '</DataArray> \n');
fprintf(results_vtu, '</PointData> \n');

% end of VTK file

fprintf(results_vtu, '</Piece> \n');
fprintf(results_vtu, '</UnstructuredGrid> \n');
fprintf(results_vtu, '</VTKFile> \n');

fclose(results_vtu);



