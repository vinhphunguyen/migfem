function elemDisp = elementDisp3d(e,pos,enrich_node,U)

% From the unknowns vector u, extract the parameters
% associated with the element "e"
% Then epsilon = B*U
% Used for 3D enriched IGA or XIGA
% Vinh Phu Nguyen
% Delft University of Technology
% The Netherlands

global controlPts element

sctr = element(e,:);
nn   = length(sctr);

% stdU contains standard dofs

idx    = 0 ;
stdU   = zeros(3*nn,1);

for in = 1 : nn
    idx      = idx + 1;
    nodeI    = sctr(in) ;
    nodeI3   = 3*nodeI;
    idx3     = 3*idx;
    stdU(idx3-2) = U(nodeI3-2);
    stdU(idx3-1) = U(nodeI3-1);
    stdU(idx3)   = U(nodeI3  );
end

% A contains enriched dofs
A = [];

if ( any(enrich_node(sctr)) == 0 ) % Non-enriched elements
    elemDisp = stdU ;
else                               % having enriched DOFs
    for in = 1 : nn
        nodeI = sctr(in);
        posI  = pos(nodeI);
        posI3 = 3*posI;
        enrI  = enrich_node(nodeI);
        
        if     (enrI == 1)     % H(x) enriched node
            AA = [U(posI3-2);
                  U(posI3-1);
                  U(posI3)];
            A  = [A;AA];
        elseif (enrI == 2)     % B(x) enriched node
            AA = [U(posI3-2);
                  U(posI3-1);
                  U(posI3);
                  U(posI3+1);
                  U(posI3+2);
                  U(posI3+3);
                  U(posI3+4);
                  U(posI3+5);
                  U(posI3+6);
                  U(posI3+7);
                  U(posI3+8);
                  U(posI3+9);
                 ];
            A  = [A;AA];
        end
    end
end

% total
elemDisp = [stdU;A];

clear stdU; clear A; clear AA;
