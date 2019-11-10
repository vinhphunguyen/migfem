function elemDisp = element_disp(e,pos,enrich_node,U)

% From the unknowns vector U, extract the dofs
% associated with element "e"
% Then epsilon = B*U
% Used for enriched IGA
% Vinh Phu Nguyen
% Johns Hopkins University, USA, @2012

global controlPts element

sctr = element(e,:);
nn   = length(sctr);

% stdU contains true nodal displacement

stdU           = zeros(2*nn,1);
stdU(1:2:2*nn) = U(2*sctr-1);
stdU(2:2:2*nn) = U(2*sctr  );

% A contains enriched dofs
A = [];

if ( any(enrich_node(sctr)) == 0 ) % Non-enriched elements
    elemDisp = stdU ; return
else                               % having enriched DOFs
    for in = 1 : nn
        nodeI = sctr(in);
        posI  = pos(nodeI);
        enrnoI= enrich_node(nodeI);
        
        if     (enrnoI == 1)     % H(x) enriched node
            AA = [U(2*posI-1);U(2*posI)];
            A  = [A;AA];
        elseif (enrnoI == 2)     % B(x) enriched node (homo. crack)
            for i=0:3
                AA = [U(2*posI+2*i-1);U(2*posI+2*i)];
                A  = [A;AA];
            end
        elseif (enrnoI == 3)     % inclusion enriched node
            AA = [U(2*posI-1);U(2*posI)];
            A  = [A;AA];
        elseif (enrnoI == 4)     % bi-mat tip/mat interface enriched node
            for i=0:12              
                AA = [U(2*posI+2*i-1);U(2*posI+2*i)];
                A  = [A;AA];
            end
        elseif (enrnoI == 5)     % bi-mat tip enriched node
            for i=0:11             
                AA = [U(2*posI+2*i-1);U(2*posI+2*i)];
                A  = [A;AA];
            end
        elseif (enrnoI == 6)     % homo. tip enriched node/mat interface enriched node
            for i=0:4             
                AA = [U(2*posI+2*i-1);U(2*posI+2*i)];
                A  = [A;AA];
            end
        end
    end
end

% total
elemDisp = [stdU;A];
