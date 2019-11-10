% convert a 2D Tspline object created by calling
%    tspline = read_bezier_extraction (filename)
%    where filename (*,iga) is created using the Tspline plugin
%    for Rhino3d.
% to the MIGFEM data structure
%
% Vinh Phu Nguyen, February 2013
% Cardiff University, UK
% nvinhphu@gmail.com

noElems     = tspline.nel; % # of elements
controlPts  = tspline.control_points(1:2,:)';
weights     = tspline.control_points(4  ,:)';
noDofs      = tspline.ndof*tspline.rdim;
noCtrPts    = tspline.ndof;
elems       = tspline.elements;
connectivity = cellfun (@transpose, {elems.connectivity}, ...
    'UniformOutput', false);
nsh         = cellfun (@numel, connectivity);
nsh_max     = max (nsh);
degrees     = zeros(noElems,2);

% for Tsplines, the elements have different number of nodes
% hence demands the use of cell data structure

element = cell(noElems,1);
C       = cell (1,noElems);

for e=1:noElems
    el           = elems(e);
    element{e}   = el.connectivity;
    C{e}         = el.extraction;
    degrees(e,:) = el.degree;
end

%% Pre-compute Bernstein basis and derivatives for ONE Bezier element
%%  do for every Bezier element type

%  remove duplicated elements

uDegree           = unique(degrees,'rows');
noBezierElemTypes = size(uDegree,1); % # of Bezier element types in the mesh

% Gauss points
% Assumption: use the Gauss rule defined by the highest order in the
% mesh!!!

[W,Q]    = quadrature(max(max(uDegree))+1, 'GAUSS', 2 );
noGpEle  = (max(max(uDegree))+1)^2;    

for et=1:noBezierElemTypes
    orders   = uDegree(et,:); % Bernstein basis order            
    p        = orders(1);
    q        = orders(2);    
    noBasis  = (p+1)*(q+1);
    
    shapes   = zeros(noGpEle,noBasis);
    derivs   = zeros(noGpEle,noBasis,2);    
        
    for gp=1:size(W,1)
        [shapes(gp,:) derivs(gp,:,:)] = getShapeGradBernstein2D(p,q,Q(gp,1),Q(gp,2));
    end
    
    bernstein(et).basis = shapes;
    bernstein(et).ders  = derivs;
    
    quad(et).weights = W;
    quad(et).points  = Q;
end

