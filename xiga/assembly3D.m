function sctrB = assembly3D(e,element,enrich_node,pos)
% B scatter matrix for 3D XIGA problems
% Vinh Phu Nguyen
% Delft University of Technology
% The Netherlands

sctr = element(e,:);
nn   = length(sctr);

% standard part 

% for k = 1 : nn
%     sctrBfem(3*k-2) = 3*sctr(k)-2 ;
%     sctrBfem(3*k-1) = 3*sctr(k)-1 ;
%     sctrBfem(3*k)   = 3*sctr(k)   ;
% end

nn3      = 3*nn;
sctrBfem = zeros(1,nn3);

sctrBfem(1:3:nn3) = 3*sctr-2;
sctrBfem(2:3:nn3) = 3*sctr-1;
sctrBfem(3:3:nn3) = 3*sctr;

% enrich part

enrnode = enrich_node(sctr);

if ( any(enrnode) == 0 ) % Non-enriched elements
    sctrB = sctrBfem ;
else
    % find the number of H enriched nodes and
    % branch enriched nodes => size of enriched B
    
    sn        = size(find(enrnode == 1),1);
    tn        = size(find(enrnode == 2),1);
    sctrBxfem = zeros(1,3*(sn*1+tn*4));
    cnt       = 0 ;
    
    for k = 1 : nn
        nk  = sctr(k);
        pk  = pos(nk);
        ek  = enrich_node(nk);
        pk3 = 3 * pk;
        
        if     (ek == 1)   % Heaviside enriched node
            cnt = cnt + 1 ;
            sctrBxfem(3*cnt - 2) = pk3 - 2;
            sctrBxfem(3*cnt - 1) = pk3 - 1;
            sctrBxfem(3*cnt    ) = pk3    ;
        elseif (ek == 2)   % branch enriched nodes
            cnt = cnt + 1 ;
            sctrBxfem(3*cnt - 2) = pk3 - 2;
            sctrBxfem(3*cnt - 1) = pk3 - 1;
            sctrBxfem(3*cnt    ) = pk3    ;

            cnt = cnt + 1 ;
            sctrBxfem(3*cnt - 2) = pk3+1; %3 * (pk+1) - 2;
            sctrBxfem(3*cnt - 1) = pk3+2; %3 * (pk+1) - 1;
            sctrBxfem(3*cnt    ) = pk3+3; %3 * (pk+1)    ;

            cnt = cnt + 1 ;
            sctrBxfem(3*cnt - 2) = pk3+4;%3 * (pk+2) - 2;
            sctrBxfem(3*cnt - 1) = pk3+5;%3 * (pk+2) - 1;
            sctrBxfem(3*cnt    ) = pk3+6;%3 * (pk+2)    ;

            cnt = cnt + 1 ;
            sctrBxfem(3*cnt - 2) = pk3+7;%3 * (pk+3) - 2;
            sctrBxfem(3*cnt - 1) = pk3+8;%3 * (pk+3) - 1;
            sctrBxfem(3*cnt    ) = pk3+9;%3 * (pk+2)    ;
        end
    end
    sctrB = [ sctrBfem sctrBxfem ];
    clear sctrBfem;
    clear sctrBxfem;
end

% sctr
% sctrB

