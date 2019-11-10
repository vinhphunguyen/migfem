function K = doConstraint(K,nodesOnCD, nextToCDNodes,  nodesOnCB, nextToBCNodes,...
                            nodesOnAD, nextToADNodes, penaltyStiffness)

for i=1:length(nodesOnCD)
    sctr  = [nodesOnCD(i) nextToCDNodes(i)];
    sctrx = 3*sctr-2;
    sctry = 3*sctr-1;
    sctrz = 3*sctr;
    
    K(sctrx,sctrx) = K(sctrx,sctrx) + penaltyStiffness;
    K(sctry,sctry) = K(sctry,sctry) + penaltyStiffness;
    K(sctrz,sctrz) = K(sctrz,sctrz) + penaltyStiffness;
end

for i=1:length(nodesOnCB)
    sctr  = [nodesOnCB(i) nextToBCNodes(i)];
    sctrx = 3*sctr-2;
    sctry = 3*sctr-1;
    sctrz = 3*sctr;
    
    K(sctrx,sctrx) = K(sctrx,sctrx) + penaltyStiffness;
    K(sctry,sctry) = K(sctry,sctry) + penaltyStiffness;
    K(sctrz,sctrz) = K(sctrz,sctrz) + penaltyStiffness;
end

for i=1:length(nodesOnAD)
    sctr  = [nodesOnAD(i) nextToADNodes(i)];
    sctrx = 3*sctr-2;
    sctry = 3*sctr-1;
    sctrz = 3*sctr;
    
    K(sctrx,sctrx) = K(sctrx,sctrx) + penaltyStiffness;
    K(sctry,sctry) = K(sctry,sctry) + penaltyStiffness;
    K(sctrz,sctrz) = K(sctrz,sctrz) + penaltyStiffness;
end
