function mesh = buildGA1DMesh (knotVec,p)


uniqueKnots=unique(knotVec);

% --------------------------------------------------------
% ------------- Define element connectivities ------------
% --------------------------------------------------------


ne            = length(uniqueKnots)-1;   % number of elements
elRange       = zeros(ne,2);        % initialise matrices
elConn        = zeros(ne,p+1);
elKnotIndices = zeros(ne,2);

% determine our element ranges and the corresponding knot indices
element=1;
previousKnotVal=0;

for i=1:length(knotVec) 
    currentKnotVal=knotVec(i);
    if knotVec(i)~=previousKnotVal
        elRange(element,:)=[previousKnotVal currentKnotVal];
        elKnotIndices(element,:)=[i-1 i];
        element=element+1;
    end
    previousKnotVal=currentKnotVal;
end

numRepeatedKnots=0;

for e=1:ne
    indices=(elKnotIndices(e,1)-p+1):elKnotIndices(e,1);
    previousKnotVals=knotVec(indices);
    currentKnotVals=ones(1,p)*knotVec(elKnotIndices(e,1));
    if isequal(previousKnotVals,currentKnotVals) && length(nonzeros(previousKnotVals))>1;
        numRepeatedKnots=numRepeatedKnots+1;
    end
    elConn(e,:)=(elKnotIndices(e,1)-p):elKnotIndices(e,1);
end

mesh.range   = elRange;
mesh.element = elConn;



