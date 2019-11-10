function [pos] = computePosOfEnrichedNodes(noCtrPts,enrich_node)
% We use fictitious nodes/control points to handle these
% additional dofs. At a H(x) enriched node, we add one fantom node and
% at tip enriched node, four fantom nodes are added. These fictitious nodes
% are numbered from the total number of true nodes, ie, from noCtrPts+1 ...

pos     = zeros(noCtrPts,1);
nsnode  = 0 ; % # of split nodes (H enriched)
ntnode  = 0 ; % # of tip enriched nodes (homogeneous crack)
ninode  = 0 ; % # of inclusion nodes
nitnode = 0 ; % # of bi-material tip/mat interface enriched
nItnode = 0 ; % # of bi-material tip enriched
nITnode = 0 ; % # of bi-material tip enriched with 4 homo branch funcs
              % and Moes function for weak discontinuity

for i = 1 : noCtrPts
    enrnoI = enrich_node(i);
    if     (enrnoI == 1)
        pos(i) = (noCtrPts + nsnode*1 + ntnode*4 + ninode*1 + ...
                  nitnode*13 + nItnode*12 + nITnode*5) + 1 ;
        nsnode = nsnode + 1 ;
    elseif (enrnoI == 2)
        pos(i) = (noCtrPts + nsnode*1 + ntnode*4 + ninode*1 + ...
                  nitnode*13 + nItnode*12 + nITnode*5) + 1 ;
        ntnode = ntnode + 1 ;
    elseif (enrnoI == 3)
        pos(i) = (noCtrPts + nsnode*1 + ntnode*4 + ninode*1 + ...
                  nitnode*13 + nItnode*12 + nITnode*5) + 1 ;
        ninode = ninode + 1 ;
    elseif (enrnoI == 4)
        pos(i) = (noCtrPts + nsnode*1 + ntnode*4 + ninode*1 + ...
                  nitnode*13 + nItnode*12 + nITnode*5) + 1 ;
        nitnode = nitnode + 1 ;
    elseif (enrnoI == 5)
        pos(i) = (noCtrPts + nsnode*1 + ntnode*4 + ninode*1 + ...
                  nitnode*13 + nItnode*12 + nITnode*5) + 1 ;
        nItnode = nItnode + 1 ;
    elseif (enrnoI == 6)
        pos(i) = (noCtrPts + nsnode*1 + ntnode*4 + ninode*1 + ...
                  nitnode*13 + nItnode*12 + nITnode*5) + 1 ;
        nITnode = nITnode + 1 ;
    end
end