function jump = computeJump(N,r,e,pos,enrich_node,U)
% compute the displacement jump of a point on the crack surface
% the point belongs to element "e"
% Used in crackedMeshNURBS.m
% Vinh Phu Nguyen
% Johns Hopkins University

global controlPts element

sctr = element(e,:);
nn   = length(sctr);

jump = zeros(1,2);

for in = 1 : nn
    nodeI = sctr(in);
    posI  = pos(nodeI);
    enrI  = enrich_node(nodeI);
    
    if     (enrI == 1)     % H(x) enriched node
        jump = jump + 2*N(in)*[U(2*posI-1) U(2*posI)];
    elseif (enrI == 2) % B(x) enriched node
        jump = jump + 2*sqrt(r)*N(in)*[U(2*posI-1) U(2*posI)];
    end
end
