nVoid = size(VOID,1);
nNode = size(node,1);
nCPts = size(controlPts,1);
chi   = zeros(nVoid,nNode);
CHI   = zeros(nVoid,nCPts);
for iVoid = 1:nVoid
    xc = VOID(iVoid,1);         % x-coordinate of void center
    yc = VOID(iVoid,2);         % y-coordinate of void center
    rc = VOID(iVoid,3);         % radius of void
    
    % Define the level set for the inclusion
    for iNode = 1:nNode
        xi = node(iNode,1);       % X-coordinate of current node
        yi = node(iNode,2);       % Y-coordinate of current node
        C  = sqrt((xi-xc)^2+(yi-yc)^2)-rc; % Level set value
        if abs(C) < 1e-6
            chi(iVoid,iNode) = 0;
        else
            chi(iVoid,iNode) = C;
        end
    end
    
    % Define the level set for the inclusion
    for iNode = 1:nCPts
        xi = controlPts(iNode,1);       % X-coordinate of current node
        yi = controlPts(iNode,2);       % Y-coordinate of current node
        C  = sqrt((xi-xc)^2+(yi-yc)^2)-rc; % Level set value
        if abs(C) < 1e-6
            CHI(iVoid,iNode) = 0;
        else
            CHI(iVoid,iNode) = C;
        end
    end
end


%%

enrich_node     = zeros(noCtrPts,1);
inclusion_node  = zeros(noCtrPts,1); % which inclusion to which the node enriched

count1 = 0;
count2 = 0;

% loop over elements to detect split elements
% and elements in the voids

inactiveElems = [];

for iel = 1 : length(elementV)
    sctr    = elementV(iel,:);
    sctrIGA = element(iel,:);
    for iVoid = 1:nVoid
        phi     = chi(iVoid,sctr);
        
        if    ( max(phi)*min(phi) < 0 )
            count1                = count1 + 1 ; % one split element
            splitElems(count1,:)  = [iel iVoid];
            enrich_node(sctrIGA)  = 3;
            inclusion_node(sctrIGA)  = iVoid;
        elseif max(phi) < 0
            count2                 = count2 + 1 ; % one inactive element
            inactiveElems(count2)  = iel;
        end
    end
end

eps = 1e-5;

ll    = length(splitElems);
idx   = [];

for ie=1:ll
    e    = splitElems(ie);
    sctr = elementV(e,:);
    phi  = chi(sctr);
    count = 0;
    for i=1:4
        if abs(phi(i)) < eps
            count = count + 1;
        end
    end
    tem = find(phi < 0);
    % remove ie from splitElems
    if (count==1) && length(tem) == 3
        idx = [idx; ie];
    end
end

inc_nodes = find(enrich_node == 3);

% inactiveElems   = [inactiveElems splitElems(idx)];
% splitElems(idx) = [];


% Combine separate level set functions into a single global function
CHI = min(CHI,[],1);