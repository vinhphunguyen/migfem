% function [node,elementV]=buildVisualizationMesh(controlPts,weights,...
%                         uKnot,vKnot,p,q)
%
% build a Q4 mesh from image of knots.
% This mesh is used for visualization purpose.
% Standard FE visualization function like plot_field can then be reused.
% This mesh is also useful for mesh and discontinuity interaction
% in extended IGA
% Vinh Phu Nguyen
% Delft University of Technology


%noPtsY      = length(vKnot)-q-1;
%noPtsX      = length(uKnot)-p-1;

controlPtsX = controlPts(:,1);
controlPtsY = controlPts(:,2);

% get rid of zero measure knot spans

uKnotVec = unique(uKnot);
vKnotVec = unique(vKnot);

% number of distinct knot values

noKnotsU = length(uKnotVec);
noKnotsV = length(vKnotVec);

%  using homogeneous coordinates

projcoord = nurb2proj(noPtsX*noPtsY, controlPts, weights);

dim=size(projcoord,2);

node  = zeros(noKnotsU*noKnotsV,2);
count = 1;

for vk=1:noKnotsV
    eta = vKnotVec(vk);
    for uk=1:noKnotsU
        xi = uKnotVec(uk);
        tem = SurfacePoint(noPtsX-1,p,uKnot,noPtsY-1,q,vKnot,...
                           projcoord,dim,xi,eta);
        node(count,1) = tem(1)/tem(3);
        node(count,2) = tem(2)/tem(3);
        count = count + 1;
    end
end

% build Q4 elements

nnx   = noKnotsU;
nny   = noKnotsV;
inc_u = 1;
inc_v = nnx;
node_pattern = [1 2 nnx+2 nnx+1];
elementV     = make_elem(node_pattern,nnx-1,nny-1,inc_u,inc_v);







