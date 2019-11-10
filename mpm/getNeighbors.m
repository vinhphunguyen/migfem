function neighbors = getNeighbors(elemId, numx, numy)

% find the neighbors of a given element "elemId"
% in a structured grid of numx x numy elements
% VP Nguyen
% To be used in MPM code.

row = ceil(elemId/numx);
col = rem (elemId,numx);

if col == 0
    col = numx;
end

if     row  == 1
    rows = [1 2];
elseif row == numy
    rows = [numy-1 numy];
else
    rows = [row-1 row row + 1];
end

if     col == 1
    cols = [1 2];
elseif col == numx
    cols = [numx-1 numx];
else
    cols = [col-1 col col + 1];
end

index = zeros(numy,numx);

id = 1;

for i=1:numy
    for j=1:numx
        index(i,j) = id; 
        id = id + 1;
    end
end

neighbors=index(rows,cols);

neighbors=neighbors(:);

