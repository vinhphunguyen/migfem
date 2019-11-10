function vMesh=buildVMesh2DInterface(surface)
%
% build a Q4 mesh from image of knots.
% This mesh is used for visualization purpose.
% Target: NURBS surfaces with C^0 lines for delamination problems.
% Assumption: C^0 lines along eta direction.
%
% Vinh Phu Nguyen
% April 2013.
% Cardiff University, Wales, UK

[p,q,uKnot,vKnot,noPtsX,noPtsY,weights,controlPts] = convertNURBSSurface (surface);

% get rid of zero measure knot spans

uKnotVec = unique(uKnot);
vKnotVec = unique(vKnot);

% number of distinct knot values

noKnotsU = length(uKnotVec);
noKnotsV = length(vKnotVec);
noElemsU = noKnotsU-1; % # of elements xi dir.
noElemsV = noKnotsV-1; % # of elements eta dir.

multiplicity  = histc(vKnot,vKnotVec); % find duplicated knots
c0LineCnt     = sum(multiplicity(2:end-1)==q+1); % number of discontinuity lines

controlPts = controlPts(:,1:2);

%  using homogeneous coordinates

projcoord = nurb2proj(noPtsX*noPtsY, controlPts, weights);
dim   = size(projcoord,2);
node  = zeros(noKnotsU*(noKnotsV+c0LineCnt),2); size(node)
count = 1;
for vk=1:noKnotsV
    eta = vKnotVec(vk);
    for uk=1:noKnotsU
        xi = uKnotVec(uk);
        tem = SurfacePoint(noPtsX-1,p,uKnot,noPtsY-1,q,vKnot,projcoord,dim,xi,eta);
        node(count,1) = tem(1)/tem(3);
        node(count,2) = tem(2)/tem(3);
        count = count + 1;
    end
    % if multiplicity = order+1 => C^0 line
    if (multiplicity(vk) == q+1) && (vk ~= 1) && (vk ~= noKnotsV)
        for uk=1:noKnotsU
            xi = uKnotVec(uk);
            tem = SurfacePoint(noPtsX-1,p,uKnot,noPtsY-1,q,vKnot,projcoord,dim,xi,eta);
            node(count,1) = tem(1)/tem(3);
            node(count,2) = tem(2)/tem(3);
            count = count + 1;
        end
    end
end

% build Q4 elements

chan  = zeros(noKnotsV+c0LineCnt,noKnotsU);
count = 1;
k     = 0;
for i=1:noKnotsV
    if (multiplicity(i) == q+1) && (i ~= 1) && (i ~= noKnotsV)
        for j=1:noKnotsU
            chan(k+1,j) = count;
            count     = count + 1;
        end
        
        for j=1:noKnotsU
            chan(k+2,j) = count;
            count = count + 1;
        end
        k = k + 2;
    else
        k = k + 1;
        for j=1:noKnotsU
            chan(k,j) = count;
            count     = count + 1;
        end
        %k = i;
    end
end

chan

elConnU = zeros(noElemsU,2);
elConnV = zeros(noElemsV,2);

for i=1:noElemsU
    elConnU(i,:) = i:i+1;
end

e=1;
for i=1:noKnotsV-1
    if (multiplicity(i) == q+1) && (i ~= 1)
        e = e + 1;
    end
    elConnV(i,:) = e:e+1;
    e = e + 1;
end


noElems = noElemsU * noElemsV;
element = zeros(noElems,4);


e = 1;
for v=1:noElemsV
    vConn = elConnV(v,:);
    for u=1:noElemsU
        c = 1;
        uConn = elConnU(u,:);
        for i=1:length(vConn)
            for j=1:length(uConn)
                element(e,c) = chan(vConn(i),uConn(j));
                c = c + 1;
            end
        end
        e = e + 1;
    end
end

% swap 3rd and 4th columns to make it Q4 elements
element(:,[3 4]) = element(:,[4 3]);

vMesh.node    =  node;
vMesh.element = element;





