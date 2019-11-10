function meshToVTK(node, sigma, disp, fname)

global uKnot vKnot wKnot

noPtsX = length(unique(uKnot))-1;
noPtsY = length(unique(vKnot))-1;
noPtsZ = length(unique(wKnot))-1;

numNodes = length(node);

sigmaXX = sigma(:,1);
sigmaYY = sigma(:,2);
sigmaZZ = sigma(:,3);
sigmaXY = sigma(:,4);
sigmaYZ = sigma(:,5);
sigmaZX = sigma(:,6);

%--------------------------------------------------------------------------
% Write the header for a structured grid:
%--------------------------------------------------------------------------
fid             = fopen([fname,'.vts'],'w','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
fprintf(fid,'<?xml version="1.0"?> \n');
fprintf(fid,'<VTKFile type="StructuredGrid" version="0.1" byte_order="BigEndian" >\n');
fprintf(fid,'  <StructuredGrid  WholeExtent="%i %i %i %i %i %i">\n', ...
    [0 noPtsX 0 noPtsY 0 noPtsZ]);
fprintf(fid,'  <Piece Extent="%i %i %i %i %i %i">\n', ...
    [0 noPtsX 0 noPtsY 0 noPtsZ]);

%--------------------------------------------------------------------------
% Add point-wise data
%--------------------------------------------------------------------------
fprintf(fid,'    <PointData Vectors="Disp"  >\n');

% STRESS -----------

fprintf(fid,'      <DataArray type="Float32" Name="Stress" NumberOfComponents="6" format="ascii">\n');

for i=1:numNodes
    fprintf(fid, '%f   ', sigmaXX(i) );
    fprintf(fid, '%f   ', sigmaYY(i) );
    fprintf(fid, '%f   ', sigmaZZ(i) );
    fprintf(fid, '%f   ', sigmaXY(i) );
    fprintf(fid, '%f   ', sigmaYZ(i) );
    fprintf(fid, '%f   ', sigmaZX(i) );
    fprintf(fid, '\n');
end

fprintf(fid,'      </DataArray>\n');
% -----------------------


% DISPLACEMENT---------------  : DISP IS A 3-component vector

fprintf(fid,'      <DataArray type="Float32" Name="Displacement" NumberOfComponents="3" format="ascii">\n');
for i=1:length(disp)
    fprintf(fid,'   %g %g %g \n',disp(i,:));
end

fprintf(fid,'      </DataArray>\n');
% -----------------------

fprintf(fid,'    </PointData>\n');
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------

fprintf(fid,'    <Celldata>\n');
fprintf(fid,'    </Celldata>\n');

%--------------------------------------------------------------------------
% Add coordinates of structured grid
%--------------------------------------------------------------------------
fprintf(fid,'    <Points>\n');

fprintf(fid,'      <DataArray type="Float32" Name="Array" NumberOfComponents="3" format="ascii">\n');
for i=1:numNodes
    fprintf(fid,' %g %g %g \n',[node(i,:)]);
end

fprintf(fid,'      </DataArray>\n');
fprintf(fid,'    </Points>\n');
%--------------------------------------------------------------------------

fprintf(fid,'  </Piece> \n');
fprintf(fid,'  </StructuredGrid> \n');
fprintf(fid,'</VTKFile>\n');
fclose(fid);
