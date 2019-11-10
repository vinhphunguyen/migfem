
controlPts  = [controlPts weights'];

% -------------------------------------
% ------------ h-refinement -----------
% -------------------------------------

% the idea here is to take the original knot vector and add knots at point
% in between existing points

if refinement > 0
for c=1:refinement
    n           = length(knotVec)-1-p;
    uniqueKnots = unique(knotVec);
    newKnotsX   = uniqueKnots(1:end-1) + 0.5*diff(uniqueKnots);
    nonewkX     = size(newKnotsX,2);
    weightedPts = [controlPts(:,1).*controlPts(:,3) ...
                   controlPts(:,2).*controlPts(:,3) controlPts(:,3)];
    
    [newKnots,newControlPts] = ...
        RefineKnotVectCurve(n-1,p,knotVec,weightedPts,newKnotsX,nonewkX-1);
    
    knotVec       = newKnots;
    controlPts  = [newControlPts(:,1)./newControlPts(:,3) ...
                   newControlPts(:,2)./newControlPts(:,3) newControlPts(:,3)];
end
end

uniqueKnots=unique(knotVec);

% if refinement>1
%     oldKnotVec = knotVec;
%     newKnots   = zeros((length(uniqueKnots)-1),refinement-1);
%     
%     for i=1:size(newKnots,1)
%         distributed   = linspace(uniqueKnots(i),uniqueKnots(i+1),refinement+1);
%         newKnots(i,:) = distributed(2:end-1);
%     end
%     
%     [m n]    = size(newKnots);
%     newKnots = sort(reshape(newKnots,1,m*n));
% 
%     % and now redefine the control points
%     for knot=1:length(newKnots)                       % loop over all the new knots
%         newControlPts = zeros(length(weightedPts)+1,3); % increase the array of control points
%         knot_bar      = newKnots(knot);
%         kval          = find(knotVec>knot_bar, 1 )-1;
%         
%         for i=1:length(newControlPts)
%             if     i<= (kval - p)
%                 alpha=1;
%             elseif i>=(kval-p+1) && i<= kval
%                 alpha= ( newKnots(knot) - oldKnotVec(i) ) /...
%                        ( oldKnotVec(i+p) - oldKnotVec(i) );
%             else
%                 alpha=0;
%             end 
%             
%             newPoints=zeros(1,3);
%             
%             if i~=1
%                 newPoints=(1-alpha).*weightedPts(i-1,:);
%             end
%             
%             if i~= length(newControlPts)
%                newPoints=newPoints + alpha * weightedPts(i,:);
%             end
%             
%             newControlPts(i,:)=newPoints;
%         end
%         weightedPts = newControlPts;
%         knotVec     = sort([knotVec knot_bar]);
%         oldKnotVec  = knotVec;
%     end
% 
% controlPts=[weightedPts(:,1)./weightedPts(:,3) weightedPts(:,2)./weightedPts(:,3) weightedPts(:,3)];
% uniqueKnots=unique(knotVec);
% end


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





