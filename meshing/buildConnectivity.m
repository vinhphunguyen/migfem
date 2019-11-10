function [elRange,elConn] = buildConnectivity(p,knotVec,noElems)
% compute connectivity of 1D NURBS (for one direction)
% also define the element ranges i.e. [xi1,xi2]
% Adapted from the IGABEM code of Robert Simpson, Cardiff, UK
% Vinh Phu Nguyen
% Delft University of Technology, The Netherlands

elRange         = zeros(noElems,2);   
elKnotIndices   = zeros(noElems,2);
elConn          = zeros(noElems,p+1);
ne              = length(unique(knotVec))-1;   % number of elements

element         = 1;
previousKnotVal = 0;

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
