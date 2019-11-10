function sctrB = assembly(e,enrich_node,pos)

% determine the scatter vector to assemble K matrix
% e: element under consideration
% element: element connectivity, seems obsolete because it is defined as a
% global variable.
% enrich_node: enrichment type for enriched nodes

global element

sctr = element(e,:);
nn   = length(sctr);

%for k = 1 : nn
%    sctrBfem(2*k-1) = 2*sctr(k)-1 ;
%    sctrBfem(2*k)   = 2*sctr(k)   ;
%end

% equivalent code with above but without loop

sctrBfem           = zeros(1,2*nn);
sctrBfem(1:2:2*nn) = 2*sctr-1;
sctrBfem(2:2:2*nn) = 2*sctr;

% determine the contribution of enriched nodes

enrnodes = enrich_node(sctr);

if ( any(enrnodes) == 0 ) % Non-enriched elements
    sctrB = sctrBfem ;
    clear sctrBfem;
else
    sn = length(find(enrnodes == 1)); % Heavise enriched
    tn = length(find(enrnodes == 2)); % homogeneous crack tip enriched
    in = length(find(enrnodes == 3)); % inclusion enriched
    cn = length(find(enrnodes == 4)); % bi-mat crack enriched
    
    sctrBxfem = zeros(1,2*(sn*1+tn*4+in*1+cn*13));
    cnt       = 0 ;
    
    % loop over nodes of the current element
    
    for k = 1 : nn
        nk     = sctr(k);
        enr_nk = enrich_node(nk);
        pnk    = pos(nk);
        pnk2   = 2 * pnk;
        
        % encounter a Heaviside enriched or inclusion enriched node
        
        if     ( enr_nk == 1) || ( enr_nk == 3)
            cnt = cnt + 1 ;
            sctrBxfem(2*cnt - 1) = pnk2 - 1;
            sctrBxfem(2*cnt    ) = pnk2    ;
            
        % encounter a homogeneous crack tip enriched node
            
        elseif ( enr_nk == 2)
            cnt = cnt + 1 ;
            sctrBxfem(2*cnt - 1) = pnk2 - 1;
            sctrBxfem(2*cnt    ) = pnk2    ;
            
            cnt = cnt + 1 ;
            sctrBxfem(2*cnt - 1) = pnk2 + 1; %2 * (pnk+1) - 1;
            sctrBxfem(2*cnt    ) = pnk2 + 2; %2 * (pnk+1)    ;
            
            cnt = cnt + 1 ;
            sctrBxfem(2*cnt - 1) = pnk2 + 3; %2 * (pnk+2) - 1;
            sctrBxfem(2*cnt    ) = pnk2 + 4; %2 * (pnk+2)    ;
            
            cnt = cnt + 1 ;
            sctrBxfem(2*cnt - 1) = pnk2 + 5; %2 * (pnk+3) - 1;
            sctrBxfem(2*cnt    ) = pnk2 + 6; %2 * (pnk+3)    ;
            
        % encounter a bi-material crack tip enriched node
            
        elseif ( enr_nk == 4)
            for i=0:12
                cnt = cnt + 1 ;
                sctrBxfem(2*cnt - 1) = pnk2 + 2*i-1;
                sctrBxfem(2*cnt    ) = pnk2 + 2*i;
            end
        elseif ( enr_nk == 5)
            for i=0:11
                cnt = cnt + 1 ;
                sctrBxfem(2*cnt - 1) = pnk2 + 2*i-1;
                sctrBxfem(2*cnt    ) = pnk2 + 2*i;
            end
        elseif ( enr_nk == 6)
            for i=0:4
                cnt = cnt + 1 ;
                sctrBxfem(2*cnt - 1) = pnk2 + 2*i-1;
                sctrBxfem(2*cnt    ) = pnk2 + 2*i;
            end
        end
    end
    sctrB = [ sctrBfem sctrBxfem ];
    
    clear sctrBfem;
    clear sctrBxfem;
end

