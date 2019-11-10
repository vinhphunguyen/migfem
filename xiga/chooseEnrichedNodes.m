% for one element, if max(phi)*min(phi) < 0
% and max(psi) < 0, then it is a split element
% If max(phi)*min(phi) < 0 and max(psi)*min(psi) < 0, it is
% tip element
% Vinh Phu Nguyen
% LTDS, ENISE, Saint Etienne, France

% Data structures for elements cut by crack
% Array split_elem contains the number of elements which are completely
% cut by crack. Similarly, array tip_elem stores number of tip element

enrich_node = zeros(noCtrPts,1);
crack_node  = zeros(noCtrPts,1); % which crack to which the node enriched

count1 = 0;
count2 = 0;

% loop over elements
for iel = 1 : numelem
    sctr    = elementV(iel,:);
    sctrIGA = element(iel,:);
    
    % loop over cracks
    for iCr = 1 : noCracks
        phi  = levelSets(iCr,sctr,1);
        psi  = levelSets(iCr,sctr,2);
        
        if ( max(phi)*min(phi) < 0 )
            if max(psi) < 0
                count1                 = count1 + 1 ; % ah, one split element
                split_elem(count1)     = iel;                
                tip_enr_pos = find(enrich_node(sctrIGA)==2);
                tip_enr     = sctrIGA(tip_enr_pos);
                
                if isempty(tip_enr)
                    enrich_node(sctrIGA) = 1;
                else
                    hea_enr = setdiff(sctrIGA,tip_enr);
                    enrich_node(hea_enr) = 1;
                end
                crack_node(sctrIGA) = iCr;
            elseif max(psi)*min(psi) < 0
                count2               = count2 + 1 ; % ah, one tip element
                tip_elem(count2)     = iel;
                enrich_node(sctrIGA) = 2;
                crack_node(sctrIGA)  = iCr;
            end
        end
    end
end
