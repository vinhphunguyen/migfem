% compute normal and tangent level sets of cracks
% levelSets(crackId,nodeId,1): normal level set of nodeI for crackId
% levelSets(crackId,nodeId,2): tangent level set of nodeI for crackId
% Vinh Phu Nguyen
% LTDS, ENISE, Saint Etienne, France

levelSets = zeros(noCracks,numnode,2);

for iCr = 1 : noCracks
    
    xCr  = reshape(xCrack(iCr,:,:),2,2);
    xTip = xTips(iCr,:);
    
    x0   = xCr(1,1); y0 = xCr(1,2);
    x1   = xCr(2,1); y1 = xCr(2,2);
    seg  = xCr(2,:) - xCr(1,:);   % tip segment
    t    = 1/norm(seg)*seg;
    
    for i = 1 : numnode
        x                  = node(i,1);
        y                  = node(i,2);
        l                  = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) ;
        phi                = (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0);
        levelSets(iCr,i,1) = phi/l;            % normal LS
        levelSets(iCr,i,2) = ([x y]-xTip)*t';  % tangent LS
    end
end
